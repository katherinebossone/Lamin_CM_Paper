library(GenomicRanges)
library(rtracklayer)

# Load data
H3K27me3 <- import("~/HiC/hilands_datasets/H3K27me3_log2fc.bw", format = "bigWig")
H3K9me3  <- import("~/HiC/hilands_datasets/WT_H3K9me3_merged.bw", format = "bigWig")
H3K4me1  <- import("~/HiC/hilands_datasets/H3K4me1_log2fc.bw", format = "bigWig")
H3K4me3  <- import("~/HiC/hilands_datasets/H3K4me3_log2fc.bw", format = "bigWig")
H3K9me2  <- import("~/cutandrun/012126/H3K9me2_log2fc.bw", format = "bigWig")
H3K27ac  <- import("~/cutandrun/012126/H3K27ac_log2fc.bw", format = "bigWig")
ATAC     <- read.table("~/HiC/hilands_datasets/ATAC_rebin_5kb.bedgraph", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
LAD      <- import("~/lamincc/LB1het_average.bw", format = "bigWig")

# Load blacklist files
blacklist1 <- import("~/low_signal_regions.bed", format = "BED")
blacklist2 <- import("~/HiC/mm10-blacklist.v2.bed", format = "BED")  # Replace with actual second file

# Combine and reduce blacklists
combined_blacklist <- GenomicRanges::reduce(c(blacklist1, blacklist2))

# Preprocessing functions
prep_log_data <- function(data, score_name, cap_min = -3, cap_max = 3) {
  data$score <- pmin(pmax(data$score, cap_min), cap_max)
  output <- data.frame(
    chrom = data$seqnames,
    start = data$start,
    end   = data$end,
    score = data$score
  )
  colnames(output)[4] <- score_name
  return(output)
}

prep_log_data2 <- function(data, score_name, cap_min = -3, cap_max = 3) {
  data$score <- log2(pmax(data$score, 1e-6))
  data$score <- pmin(pmax(data$score, cap_min), cap_max)
  output <- data.frame(
    chrom = data$seqnames,
    start = data$start,
    end   = data$end,
    score = data$score
  )
  colnames(output)[4] <- score_name
  return(output)
}
# Function to remove blacklisted regions
remove_blacklisted_regions <- function(df, blacklist) {
  gr <- GRanges(seqnames = df$chrom,
                ranges = IRanges(start = df$start, end = df$end))
  overlaps <- findOverlaps(gr, blacklist)
  keep <- setdiff(seq_along(gr), queryHits(overlaps))
  df_filtered <- df[keep, ]
  return(df_filtered)
}

# Convert to data frames
H3K9me3_df   <- as.data.frame(H3K9me3)
H3K9me2_df   <- as.data.frame(H3K9me2)
H3K27me3_df  <- as.data.frame(H3K27me3)
H3K27ac_df  <- as.data.frame(H3K27ac)
H3K4me1_df   <- as.data.frame(H3K4me1)
H3K4me3_df   <- as.data.frame(H3K4me3)
LAD_df       <- as.data.frame(LAD)
colnames(ATAC) <- c("seqnames", "start", "end", "score")

# Log-transform signals
H3K9me3_log  <- prep_log_data(H3K9me3_df, "H3K9me3_score")
H3K9me2_log  <- prep_log_data(H3K9me2_df, "H3K9me2_score")
H3K27me3_log <- prep_log_data(H3K27me3_df, "H3K27me3_score")
H3K27ac_log <- prep_log_data(H3K27ac_df, "H3K27ac_score")
H3K4me1_log  <- prep_log_data(H3K4me1_df, "H3K4me1_score")
H3K4me3_log  <- prep_log_data(H3K4me3_df, "H3K4me3_score")
ATAC_log     <- prep_log_data2(ATAC, "ATAC_score")
LAD_log      <- prep_log_data(LAD_df, "LAD_score")

# Filter to standard chromosomes
chromosomes_to_keep <- paste0("chr", c(1:19))
filter_chromosomes <- function(df) {
  df[df$chrom %in% chromosomes_to_keep, ]
}

H3K9me3_log  <- filter_chromosomes(H3K9me3_log)
H3K27me3_log <- filter_chromosomes(H3K27me3_log)
H3K4me1_log  <- filter_chromosomes(H3K4me1_log)
H3K4me3_log  <- filter_chromosomes(H3K4me3_log)
ATAC_log     <- filter_chromosomes(ATAC_log)
LAD_log      <- filter_chromosomes(LAD_log)
H3K9me2_log  <- filter_chromosomes(H3K9me2_log)
H3K27ac_log <- filter_chromosomes(H3K27ac_log)

# Trim first/last 3Mb from chromosomes
trim_chromosome_edges <- function(df, trim_size = 3e6) {
  trimmed <- do.call(rbind, lapply(split(df, df$chrom), function(chrom_df) {
    max_end <- max(chrom_df$end)
    chrom_df[chrom_df$start >= trim_size & chrom_df$end <= (max_end - trim_size), ]
  }))
  return(trimmed)
}

H3K9me3_log  <- trim_chromosome_edges(H3K9me3_log)
H3K27me3_log <- trim_chromosome_edges(H3K27me3_log)
H3K4me1_log  <- trim_chromosome_edges(H3K4me1_log)
H3K4me3_log  <- trim_chromosome_edges(H3K4me3_log)
ATAC_log     <- trim_chromosome_edges(ATAC_log)
LAD_log      <- trim_chromosome_edges(LAD_log)
H3K9me2_log  <- trim_chromosome_edges(H3K9me2_log)
H3K27ac_log <- trim_chromosome_edges(H3K27ac_log)


# Remove blacklisted regions (using combined blacklist)
H3K9me3_log  <- remove_blacklisted_regions(H3K9me3_log, combined_blacklist)
H3K27me3_log <- remove_blacklisted_regions(H3K27me3_log, combined_blacklist)
H3K4me1_log  <- remove_blacklisted_regions(H3K4me1_log, combined_blacklist)
H3K4me3_log  <- remove_blacklisted_regions(H3K4me3_log, combined_blacklist)
ATAC_log     <- remove_blacklisted_regions(ATAC_log, combined_blacklist)
LAD_log      <- remove_blacklisted_regions(LAD_log, combined_blacklist)
H3K9me2_log  <- remove_blacklisted_regions(H3K9me2_log, combined_blacklist)
H3K27ac_log <- remove_blacklisted_regions(H3K27ac_log, combined_blacklist)

# Adjust start coordinates to 0-based
H3K9me3_log$start  <- H3K9me3_log$start - 1
H3K27me3_log$start <- H3K27me3_log$start - 1
H3K4me1_log$start  <- H3K4me1_log$start - 1
H3K4me3_log$start  <- H3K4me3_log$start - 1
LAD_log$start      <- LAD_log$start - 1
H3K9me2_log$start  <- H3K9me2_log$start - 1
H3K27ac_log$start <- H3K27ac_log$start - 1

# Merge all signals by region
combined <- Reduce(function(x, y) merge(x, y, by = c("chrom", "start", "end"), all = TRUE),
                   list(H3K9me2_log,H3K9me3_log, H3K27me3_log, H3K4me1_log, H3K4me3_log, H3K27ac_log, ATAC_log, LAD_log))
combined[is.na(combined)] <- 0

# Prepare PCA input
pca_data_resampled <- combined[, c("H3K9me2_score", "H3K9me3_score", "H3K27me3_score", "H3K4me1_score", "H3K4me3_score", "H3K27ac_score", "ATAC_score", "LAD_score")]
pca_data_scaled <- scale(pca_data_resampled)

# Perform PCA
pca_result_resampled <- prcomp(pca_data_scaled, center = TRUE, scale. = TRUE)
pca_scores_df <- as.data.frame(pca_result_resampled$x[, 1:3])  # PC1–PC3

# Combine with genomic coordinates
selected_columns <- combined[, c("chrom", "start", "end")]
final_data <- cbind(selected_columns, pca_scores_df)

# Save output
write.table(final_data, file = "~/HiC/hilands_datasets/PC_file_log_cap.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
