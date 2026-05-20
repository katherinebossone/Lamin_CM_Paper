
library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(zoo)
library(GenomicRanges)

# parameters
MIN_READS <- 20
NEW_RES   <- 20000
PSEUDO    <- 1e-20

# Region
ROI_CHR    <- "chr10"
PLOT_START <- 124000000
PLOT_END   <- 127131000

# Gene regions used
# Cdk4 chr10: 124000000 - 127131000
# Jun chr4: 94000000 - 97500000

# File locations
wt_files  <- c("P0_5_HiC/indiv/P0_5_LB1Fplus_8450_2_5kb.tsv",
               "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_4_5kb.tsv")
mut_files <- c("P0_5_HiC/indiv/P0_5_LB1Fd_8450_7_5kb.tsv",
               "P0_5_HiC/indiv/P0_5_LB1Fd_8450_8_5kb.tsv")

wt_ins_file <- "fanc_insulation/LB1Fplus_fanc.insulation_1mb.bed"

out_dir <- "P0_5_HiC/dynamic_binned"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Define functions for pipeline
# Load FANC files and call boundaries
load_insulation <- function(file) {
  ins <- fread(file, header = FALSE)
  setnames(ins, c("chr", "start", "end", "dot1", "insulation_score", "dot2"))
  ins[, insulation_score := as.numeric(insulation_score)]
  setorder(ins, chr, start)
  ins[, insulation_smooth := zoo::rollmean(insulation_score, k = 3, fill = NA),
      by = chr]
  return(ins)
}

call_boundaries <- function(ins, min_score = -0.3) {
  window_size   <- 8
  boundary_list <- list()
  for (chr_name in unique(ins$chr)) {
    chr_data <- ins[chr == chr_name]
    chr_data[, boundary := FALSE]
    n <- nrow(chr_data)
    if (n <= 2 * window_size) next
    for (i in (window_size + 1):(n - window_size)) {
      window_vals <- chr_data$insulation_smooth[(i - window_size):(i + window_size)]
      center_val  <- chr_data$insulation_smooth[i]
      if (!is.na(center_val) &&
          all(center_val < window_vals[-(window_size + 1)], na.rm = TRUE))
        chr_data$boundary[i] <- TRUE
    }
    boundary_list[[chr_name]] <- chr_data[boundary == TRUE &
                                            insulation_smooth < min_score]
  }
  boundaries <- rbindlist(boundary_list)
  GRanges(seqnames = boundaries$chr,
          ranges   = IRanges(start = boundaries$start, end = boundaries$end),
          score    = boundaries$insulation_smooth)
}

get_boundaries_in_region <- function(boundaries_gr, chrom, start_bp, end_bp) {
  chr_bounds <- boundaries_gr[seqnames(boundaries_gr) == chrom]
  positions  <- start(chr_bounds) - 1   # convert to 0-based
  positions[positions >= start_bp & positions <= end_bp]
}

# Genome-wide ICE total (cis only, all chromosomes, all replicates)

genome_wide_ice_total <- function(files) {
  total <- 0
  for (f in files) {
    hic <- fread(f, header = FALSE, sep = "\t",
                 select = c(1, 4, 8),
                 colClasses = list(character = c(1, 2)))
    setnames(hic, c("chr1", "chr2", "Count_ICE"))
    total <- total + sum(hic[chr1 == chr2, Count_ICE], na.rm = TRUE)
  }
  return(total)
}

# Load and pool region contacts
load_and_rebin <- function(file_path) {
  hic <- fread(file_path, header = FALSE, sep = "\t",
               colClasses = list(character = c(1, 4)))
  setnames(hic, c("chr1","start1","end1","chr2","start2","end2",
                  "Count_raw","Count_ICE"))
  hic <- hic[chr1 == ROI_CHR & chr2 == ROI_CHR]
  hic[, bin1 := (start1 %/% NEW_RES) * NEW_RES]
  hic[, bin2 := (start2 %/% NEW_RES) * NEW_RES]
  hic[bin1 > bin2, c("bin1","bin2") := .(bin2, bin1)]
  hic <- hic[bin1 >= PLOT_START & bin1 <= PLOT_END &
               bin2 >= PLOT_START & bin2 <= PLOT_END]
  hic <- hic[, .(Count_raw = sum(Count_raw, na.rm = TRUE),
                 Count_ICE = sum(Count_ICE, na.rm = TRUE)),
             by = .(bin1, bin2)]
  return(hic)
}

pool_replicates <- function(files) {
  combined <- rbindlist(lapply(files, load_and_rebin))
  pooled   <- combined[, .(Count_raw = sum(Count_raw, na.rm = TRUE),
                           Count_ICE = sum(Count_ICE, na.rm = TRUE)),
                       by = .(bin1, bin2)]
  return(pooled)
}

# Build matrices
build_matrices <- function(pooled) {
  bins <- sort(unique(c(pooled$bin1, pooled$bin2)))
  n    <- length(bins)
  idx  <- setNames(seq_along(bins), as.character(bins))
  raw_mat <- matrix(0, n, n, dimnames = list(bins, bins))
  ice_mat <- matrix(0, n, n, dimnames = list(bins, bins))
  xi    <- idx[as.character(pooled$bin1)]
  yi    <- idx[as.character(pooled$bin2)]
  valid <- !is.na(xi) & !is.na(yi)
  raw_mat[cbind(xi[valid], yi[valid])] <- pooled$Count_raw[valid]
  raw_mat[cbind(yi[valid], xi[valid])] <- pooled$Count_raw[valid]
  ice_mat[cbind(xi[valid], yi[valid])] <- pooled$Count_ICE[valid]
  ice_mat[cbind(yi[valid], xi[valid])] <- pooled$Count_ICE[valid]
  return(list(raw = raw_mat, ice = ice_mat, bins = bins))
}

# Align matrices to a shared bin set
align_mat <- function(mat, all_bins) {
  n   <- length(all_bins)
  out <- matrix(0, n, n, dimnames = list(all_bins, all_bins))
  idx <- match(rownames(mat), as.character(all_bins))
  valid <- !is.na(idx)
  out[idx[valid], idx[valid]] <- mat[valid, valid]
  return(out)
}

# Scale region ICE matrix by genome-wide library total
scale_to_genome_total <- function(ice_mat, genome_total, target = 1e9) {
  ice_mat * (target / genome_total)
}

# Dynamic binning such that bins have defined minimum reads + creating log2FC matrix
dynamic_bin_log2fc <- function(wt_raw, wt_ice_scaled,
                               mut_raw, mut_ice_scaled) {
  n      <- nrow(wt_raw)
  result <- matrix(NA_real_, n, n, dimnames = dimnames(wt_ice_scaled))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      
      if (wt_raw[i, j] >= MIN_READS && mut_raw[i, j] >= MIN_READS) {
        wt_val  <- wt_ice_scaled[i, j]
        mut_val <- mut_ice_scaled[i, j]
      } else {
        wt_val <- mut_val <- NA_real_
        for (radius in seq_len(n)) {
          low_i  <- max(1, i - radius); high_i <- min(n, i + radius)
          low_j  <- max(1, j - radius); high_j <- min(n, j + radius)
          # Keep upper triangle orientation
          if (j >= i) low_j  <- max(low_j,  low_i)
          else        high_j <- min(high_j, high_i)
          both_met <- sum(wt_raw[low_i:high_i,  low_j:high_j], na.rm = TRUE) >= MIN_READS &&
            sum(mut_raw[low_i:high_i, low_j:high_j], na.rm = TRUE) >= MIN_READS
          hit_edge <- low_i == 1 && high_i == n && low_j == 1 && high_j == n
          if (both_met || hit_edge) {
            wt_val  <- mean(wt_ice_scaled[low_i:high_i,  low_j:high_j], na.rm = TRUE)
            mut_val <- mean(mut_ice_scaled[low_i:high_i, low_j:high_j], na.rm = TRUE)
            break
          }
        }
      }
      
      if (!is.na(wt_val) && !is.na(mut_val))
        result[i, j] <- log2((mut_val + PSEUDO) / (wt_val + PSEUDO))
    }
  }
  return(result)
}

# Pipeline
# Calculate genome-wide library totals
wt_genome_total  <- genome_wide_ice_total(wt_files)
mut_genome_total <- genome_wide_ice_total(mut_files)

# Load and pool plot region
wt_pooled  <- pool_replicates(wt_files)
wt_mats    <- build_matrices(wt_pooled)
rm(wt_pooled); gc()

mut_pooled <- pool_replicates(mut_files)
mut_mats   <- build_matrices(mut_pooled)
rm(mut_pooled); gc()

# Align to shared bins
all_bins <- sort(unique(c(wt_mats$bins, mut_mats$bins)))

wt_raw  <- align_mat(wt_mats$raw,  all_bins)
wt_ice  <- align_mat(wt_mats$ice,  all_bins)
mut_raw <- align_mat(mut_mats$raw, all_bins)
mut_ice <- align_mat(mut_mats$ice, all_bins)
rm(wt_mats, mut_mats); gc()

# Scale by genome-wide totals 
wt_ice_scaled  <- scale_to_genome_total(wt_ice,  wt_genome_total)
mut_ice_scaled <- scale_to_genome_total(mut_ice, mut_genome_total)
rm(wt_ice, mut_ice); gc()

# Dynamic binning log2FC
log2fc_mat <- dynamic_bin_log2fc(wt_raw, wt_ice_scaled, mut_raw, mut_ice_scaled)
rm(wt_raw, mut_raw, wt_ice_scaled, mut_ice_scaled); gc()

# Load insulation and call WT boundaries
wt_ins       <- load_insulation(wt_ins_file)
wt_gr        <- call_boundaries(wt_ins, min_score = -0.3)
bound_pos_bp <- get_boundaries_in_region(wt_gr, ROI_CHR, PLOT_START, PLOT_END)

# Plot
mat_to_long <- function(mat, start_bp, end_bp) {
  bins <- as.numeric(rownames(mat))
  keep <- bins >= start_bp & bins <= end_bp
  mat  <- mat[keep, keep]
  df   <- reshape2::melt(mat)
  colnames(df) <- c("x", "y", "value")
  df$x <- as.numeric(as.character(df$x))
  df$y <- as.numeric(as.character(df$y))
  return(df)
}

fc_long <- mat_to_long(log2fc_mat, PLOT_START, PLOT_END)
plot_limit <- quantile(abs(fc_long$value), 0.99, na.rm = TRUE)
mb_breaks  <- pretty(range(c(fc_long$x, fc_long$y), na.rm = TRUE) / 1e6, n = 8) * 1e6

p <- ggplot(fc_long, aes(x, y, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(
    low      = "blue3",
    mid      = "white",
    high     = "red3",
    midpoint = 0,
    limits   = c(-plot_limit, plot_limit),
    oob      = squish,
    na.value = "grey95",
    name     = "log2FC"
  ) +
  scale_x_continuous(breaks = mb_breaks,
                     labels = function(x) paste0(x / 1e6, " Mb"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = mb_breaks,
                     labels = function(x) paste0(x / 1e6, " Mb"),
                     expand = c(0, 0)) +
  labs(
    title    = "log2(ICE LB1null / ICE LB1het)",
    subtitle = "Genome-wide ICE scaling · Dynamic binning (min 20 reads) · WT boundaries",
    x        = ROI_CHR,
    y        = ROI_CHR
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    panel.grid    = element_blank(),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    aspect.ratio  = 1
  )

if (length(bound_pos_bp) > 0) {
  p <- p +
    geom_vline(xintercept = bound_pos_bp, color = "black",
               linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = bound_pos_bp, color = "black",
               linewidth = 0.3, alpha = 0.7) +
    geom_vline(xintercept = GENE_POS, color = "cyan",
             linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = GENE_POS, color = "cyan",
               linewidth = 0.3, alpha = 0.7)
  
}
p