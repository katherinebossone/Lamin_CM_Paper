library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(zoo)
library(GenomicRanges)

# Parameters
MIN_READS   <- 20
ORIG_RES    <- 5000
NEW_RES     <- 10000
ROI_CHR     <- "chr10"
PLOT_START  <- 50.0000e6
PLOT_END    <- 57.75000e6

# Load files
wt_files  <- c("P0_5_LB1Fplus_8450_2_5kb.tsv",
               "P0_5_LB1Fplus_8450_4_5kb.tsv")
mut_files <- c("P0_5_LB1Fd_8450_7_5kb.tsv",
               "P0_5_LB1Fd_8450_8_5kb.tsv")

wt_ins_file  <- "LB1Fplus_fanc.insulation_1mb.bed"
mut_ins_file <- "LB1Fd_fanc.insulation_1mb.bed"

out_dir <- "P0_5_HiC/dynamic_binned"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Insulation score loading and boundary calling
load_insulation <- function(file) {
  ins <- fread(file, header = FALSE)
  setnames(ins, c("chr", "start", "end", "dot1", "insulation_score", "dot2"))
  ins[, insulation_score := as.numeric(insulation_score)]
  setorder(ins, chr, start)
  ins[, insulation_smooth := zoo::rollmean(insulation_score, k = 3, fill = NA), by = chr]
  return(ins)
}

call_boundaries <- function(ins, min_score = -0.3) {
  window_size <- 8
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
          all(center_val < window_vals[-(window_size + 1)], na.rm = TRUE)) {
        chr_data$boundary[i] <- TRUE
      }
    }
    chr_boundaries <- chr_data[boundary == TRUE & insulation_smooth < min_score]
    boundary_list[[chr_name]] <- chr_boundaries
  }
  boundaries <- rbindlist(boundary_list)
  GRanges(seqnames = boundaries$chr,
          ranges   = IRanges(start = boundaries$start, end = boundaries$end),
          score    = boundaries$insulation_smooth)
}

get_boundaries_in_region <- function(boundaries_gr, chrom, start_bp, end_bp) {
  chr_bounds <- boundaries_gr[seqnames(boundaries_gr) == chrom]
  positions  <- start(chr_bounds)
  positions[positions >= start_bp & positions <= end_bp]
}

# Load Hi-C data, rebin, pool replicates
load_and_rebin <- function(file_path) {
  message("  Reading: ", basename(file_path))
  hic <- fread(file_path, header = FALSE, sep = "\t",
               colClasses = list(character = c(1, 4)))
  setnames(hic, c("chr1","start1","end1","chr2","start2","end2",
                  "Count_raw","Count_ICE"))
  hic <- hic[chr1 == ROI_CHR & chr2 == ROI_CHR]
  message("    Pairs on ", ROI_CHR, ": ", nrow(hic))
  hic[, bin1 := (start1 %/% NEW_RES) * NEW_RES]
  hic[, bin2 := (start2 %/% NEW_RES) * NEW_RES]
  hic[bin1 > bin2, c("bin1","bin2") := .(bin2, bin1)]
  hic <- hic[, .(Count_raw = sum(Count_raw, na.rm = TRUE),
                 Count_ICE = sum(Count_ICE, na.rm = TRUE)),
             by = .(bin1, bin2)]
  message("    After rebinning to ", NEW_RES/1000, "kb: ", nrow(hic), " bin pairs")
  return(hic)
}

pool_replicates <- function(files) {
  reps     <- lapply(files, load_and_rebin)
  combined <- rbindlist(reps)
  pooled   <- combined[, .(Count_raw = sum(Count_raw, na.rm = TRUE),
                           Count_ICE = sum(Count_ICE, na.rm = TRUE)),
                       by = .(bin1, bin2)]
  message("  Pooled replicates: ", nrow(pooled), " unique bin pairs")
  return(pooled)
}

# Build square matrices
build_matrices <- function(pooled) {
  bins <- sort(unique(c(pooled$bin1, pooled$bin2)))
  n    <- length(bins)
  idx  <- setNames(seq_along(bins), as.character(bins))
  raw_mat <- matrix(0, n, n, dimnames = list(bins, bins))
  ice_mat <- matrix(0, n, n, dimnames = list(bins, bins))
  xi <- idx[as.character(pooled$bin1)]
  yi <- idx[as.character(pooled$bin2)]
  valid <- !is.na(xi) & !is.na(yi)
  raw_mat[cbind(xi[valid], yi[valid])] <- pooled$Count_raw[valid]
  raw_mat[cbind(yi[valid], xi[valid])] <- pooled$Count_raw[valid]
  ice_mat[cbind(xi[valid], yi[valid])] <- pooled$Count_ICE[valid]
  ice_mat[cbind(yi[valid], xi[valid])] <- pooled$Count_ICE[valid]
  message("  Matrices: ", n, " x ", n, " bins at ", NEW_RES/1000, "kb")
  return(list(raw = raw_mat, ice = ice_mat, bins = bins))
}

# Dynamic binning
dynamic_bin_single <- function(raw_mat, ice_mat) {
  n      <- nrow(raw_mat)
  result <- matrix(NA_real_, n, n, dimnames = dimnames(ice_mat))
  for (i in seq_len(n)) {
    if (i %% 50 == 0) message("    Row ", i, " / ", n)
    for (j in seq_len(n)) {
      if (raw_mat[i, j] >= MIN_READS) {
        result[i, j] <- ice_mat[i, j]
      } else {
        for (radius in seq_len(n)) {
          low_i  <- max(1, i - radius); high_i <- min(n, i + radius)
          low_j  <- max(1, j - radius); high_j <- min(n, j + radius)
          if (j >= i) low_j  <- max(low_j,  low_i)
          else        high_j <- min(high_j, high_i)
          raw_sub <- raw_mat[low_i:high_i, low_j:high_j]
          if (sum(raw_sub, na.rm = TRUE) >= MIN_READS ||
              (low_i == 1 && high_i == n && low_j == 1 && high_j == n)) {
            result[i, j] <- mean(ice_mat[low_i:high_i, low_j:high_j], na.rm = TRUE)
            break
          }
        }
      }
    }
  }
  return(result)
}

# Align matrices to shared bin set
align_mat <- function(mat, all_bins) {
  n   <- length(all_bins)
  out <- matrix(0, n, n, dimnames = list(all_bins, all_bins))
  idx <- match(rownames(mat), as.character(all_bins))
  valid <- !is.na(idx)
  out[idx[valid], idx[valid]] <- mat[valid, valid]
  return(out)
}

# Pipeline
message("Loading WT")
wt_pooled  <- pool_replicates(wt_files)
wt_mats    <- build_matrices(wt_pooled)

message("Loading Mutant")
mut_pooled <- pool_replicates(mut_files)
mut_mats   <- build_matrices(mut_pooled)

all_bins <- sort(unique(c(wt_mats$bins, mut_mats$bins)))

wt_raw  <- align_mat(wt_mats$raw,  all_bins)
wt_ice  <- align_mat(wt_mats$ice,  all_bins)
mut_raw <- align_mat(mut_mats$raw, all_bins)
mut_ice <- align_mat(mut_mats$ice, all_bins)

message("Dynamic binning: WT")
wt_dynbin <- dynamic_bin_single(wt_raw, wt_ice)

message("Dynamic binning: Mutant")
mut_dynbin <- dynamic_bin_single(mut_raw, mut_ice)

# Save full chromosome outputs
save_mat <- function(mat, fname) {
  df <- as.data.table(melt(mat))
  setnames(df, c("bin1","bin2","value"))
  fwrite(df, file.path(out_dir, fname), sep = "\t")
  message("  Saved: ", fname)
}
save_mat(wt_dynbin,  "WT_dynbin.tsv")
save_mat(mut_dynbin, "mutant_dynbin.tsv")

# Load insulation and call WT boundaries
wt_ins  <- load_insulation(wt_ins_file)
wt_gr   <- call_boundaries(wt_ins, min_score = -0.3)

# Get WT boundary positions in plot region
bound_pos_bp <- get_boundaries_in_region(wt_gr, ROI_CHR, PLOT_START, PLOT_END)

# Plot
mat_to_long <- function(mat, start_bp, end_bp) {
  df <- melt(mat)
  colnames(df) <- c("x","y","value")
  df$x <- as.numeric(as.character(df$x))
  df$y <- as.numeric(as.character(df$y))
  df <- df[df$x >= start_bp & df$x <= end_bp &
             df$y >= start_bp & df$y <= end_bp, ]
  return(df)
}

make_heatmap <- function(df, title, low, mid, high, midpoint,
                         vmin = NULL, vmax = NULL, limit_sym = FALSE,
                         boundaries = NULL) {
  lim <- if (limit_sym) {
    mx <- quantile(abs(df$value), 0.99, na.rm = TRUE); c(-mx, mx)
  } else {
    c(if (is.null(vmin)) 0 else vmin,
      if (is.null(vmax)) quantile(df$value, 0.99, na.rm = TRUE) else vmax)
  }
  mb_breaks <- pretty(range(c(df$x, df$y), na.rm = TRUE) / 1e6, n = 8) * 1e6
  
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = low, mid = mid, high = high,
                         midpoint = midpoint, limits = lim,
                         oob = squish, na.value = "grey95") +
    scale_x_continuous(breaks = mb_breaks,
                       labels = function(x) paste0(x/1e6, " Mb"),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = mb_breaks,
                       labels = function(x) paste0(x/1e6, " Mb"),
                       expand = c(0, 0)) +
    labs(title = title, x = ROI_CHR, y = ROI_CHR) +
    theme_minimal(base_size = 10) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5),
          panel.grid   = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1)
  
  # Overlay WT boundaries
  if (!is.null(boundaries) && length(boundaries) > 0) {
    p <- p +
      geom_vline(xintercept = boundaries, color = "black",
                 linewidth = 0.3, alpha = 0.7, linetype = "solid") +
      geom_hline(yintercept = boundaries, color = "black",
                 linewidth = 0.3, alpha = 0.7, linetype = "solid")
  }
  return(p)
}

wt_long  <- mat_to_long(wt_dynbin,  PLOT_START, PLOT_END)
mut_long <- mat_to_long(mut_dynbin, PLOT_START, PLOT_END)

ice_cap <- quantile(c(wt_long$value, mut_long$value), 0.95, na.rm = TRUE)

p_wt <- make_heatmap(wt_long, "WT",
                     low = "white", mid = "#FB6A4A", high = "#67000D",
                     midpoint = ice_cap / 2,
                     vmin = 0, vmax = ice_cap,
                     boundaries = bound_pos_bp)

p_mut <- make_heatmap(mut_long, "Mutant",
                      low = "white", mid = "#FB6A4A", high = "#67000D",
                      midpoint = ice_cap / 2,
                      vmin = 0, vmax = ice_cap,
                      boundaries = bound_pos_bp)
