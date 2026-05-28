#!/usr/bin/env Rscript

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(ggplot2)

# Load and adjust files
het_bed  <- "Desktop/bigwig_lamins/LB1het_LAC_2state_S2.bed"
null_bed <- "Desktop/bigwig_lamins/LB1null_LAC_2state_S2.bed"

het  <- import(het_bed,  format = "BED")
null <- import(null_bed, format = "BED")
start(het)  <- start(het)  - 1
start(null) <- start(null) - 1

if (inherits(het,  "GRangesList")) het  <- unlist(het,  use.names = FALSE)
if (inherits(null, "GRangesList")) null <- unlist(null, use.names = FALSE)

keep_chrs <- paste0("chr", 1:19)
het  <- het[seqnames(het)  %in% keep_chrs]
null <- null[seqnames(null) %in% keep_chrs]

# Merge adjacent bins into continuous LADs
het  <- GenomicRanges::reduce(het,  min.gapwidth = 1)
null <- GenomicRanges::reduce(null, min.gapwidth = 1)

# Remove first 3 Mb of each chromosome
all_chrs <- unique(c(seqlevels(het), seqlevels(null)))
exclude  <- GRanges(seqnames = all_chrs, ranges = IRanges(start = 1, end = 3e6))
het  <- het[!overlapsAny(het,  exclude)]
null <- null[!overlapsAny(null, exclude)]

# Finding overlapping LADs
overlaps <- findOverlaps(het, null)

het_shared_idx  <- unique(queryHits(overlaps))
null_shared_idx <- unique(subjectHits(overlaps))
het_unique      <- het[-het_shared_idx]

lost_sizes <- width(het_unique)

# Identifying shaved LADs
threshold  <- 0.10
ov_df      <- as.data.frame(overlaps)
het_to_null <- split(ov_df$subjectHits, ov_df$queryHits)
frac_lost  <- numeric(length(het_shared_idx))

for (k in seq_along(het_shared_idx)) {
  het_idx     <- het_shared_idx[k]
  het_range   <- het[het_idx]
  null_indices <- het_to_null[[as.character(het_idx)]]
  if (length(null_indices) == 0) { frac_lost[k] <- 1; next }
  null_ranges    <- null[null_indices]
  shared_regions <- GenomicRanges::reduce(pintersect(rep(het_range, length(null_ranges)), null_ranges))
  frac_lost[k]   <- 1 - sum(width(shared_regions)) / width(het_range)
}

shaved_idx   <- which(frac_lost > threshold)
shaved_sizes <- width(het[het_shared_idx[shaved_idx]])

# Plot
plot_df <- data.frame(
  size_Mb  = c(lost_sizes, shaved_sizes) / 1e6,
  Category = rep(c("Complete LAD loss", "LAD loss at edge"),
                 c(length(lost_sizes), length(shaved_sizes)))
)

ggplot(plot_df, aes(x = Category, y = size_Mb, fill = Category)) +
  geom_violin(aes(fill = Category), trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = c("#909090", "#4a90b8")) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = c("#909090", "#4a90b8")) +
  labs(x = NULL, y = "LAD size (log, Mb)") +
  theme_bw() +
  theme(legend.position = "none")
