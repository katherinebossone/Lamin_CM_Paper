#!/usr/bin/env Rscript

library(rtracklayer)
library(GenomicRanges)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(grid)

# Load and adjust start coordinates 
het_bed  <- "Desktop/bigwig_lamins/LB1het_LAC_2state_S2.bed"
null_bed <- "Desktop/bigwig_lamins/LB1null_LAC_2state_S2.bed"

het  <- import(het_bed,  format = "BED")
null <- import(null_bed, format = "BED")

start(het)  <- start(het)  - 1
start(null) <- start(null) - 1

if (inherits(het,  "GRangesList")) het  <- unlist(het,  use.names = FALSE)
if (inherits(null, "GRangesList")) null <- unlist(null, use.names = FALSE)

# Remove first 3 Mb of each chromosome and keep autosomes only 
keep_chrs <- paste0("chr", 1:19)

het  <- het[seqnames(het)  %in% keep_chrs]
null <- null[seqnames(null) %in% keep_chrs]

exclude <- GRanges(
  seqnames = keep_chrs,
  ranges   = IRanges(start = 1, end = 3e6)
)

het  <- het[!overlapsAny(het,  exclude)]
null <- null[!overlapsAny(null, exclude)]

## Plot venn diagram comparing unique and shared LADs ##

overlaps_bins    <- findOverlaps(het, null)
shared_bins      <- length(unique(queryHits(overlaps_bins)))
het_only_bins    <- length(het)  - shared_bins
null_only_bins   <- length(null) - shared_bins

bin_size         <- 5000
mouse_genome_bp  <- 2.65e9
total_bins       <- mouse_genome_bp / bin_size

shared_pct    <- round(100 * shared_bins    / total_bins, 2)
het_only_pct  <- round(100 * het_only_bins  / total_bins, 2)
null_only_pct <- round(100 * null_only_bins / total_bins, 2)

# plot
venn.plot <- draw.pairwise.venn(
  area1      = het_only_pct + shared_pct,
  area2      = null_only_pct + shared_pct,
  cross.area = shared_pct,
  category   = c("Het LAD genome", "Null LAD genome"),
  fill       = c("#4DBBD5FF", "#E64B35FF"),
  alpha      = 0.6,
  cex        = 2,
  cat.cex    = 1.5,
  ind        = FALSE
)

grid.newpage()
grid.draw(venn.plot)

## Compare regions of LAD loss vs shaved LADs ##
# Merge adjacent bins into continuous LADs
het_lads  <- GenomicRanges::reduce(het,  min.gapwidth = 1)
null_lads <- GenomicRanges::reduce(null, min.gapwidth = 1)

# Re-apply 3 Mb exclusion after merging (seqlevels may have expanded)
all_chrs <- unique(c(seqlevels(het_lads), seqlevels(null_lads)))
exclude2 <- GRanges(seqnames = all_chrs, ranges = IRanges(start = 1, end = 3e6))

het_lads  <- het_lads[!overlapsAny(het_lads,  exclude2)]
null_lads <- null_lads[!overlapsAny(null_lads, exclude2)]

# Overlapping LADs between het and null
overlaps_lads   <- findOverlaps(het_lads, null_lads)
het_shared_idx  <- unique(queryHits(overlaps_lads))

# Completely lost LADs (het only, no overlap with null)
het_unique  <- het_lads[-het_shared_idx]
lost_sizes  <- width(het_unique)

# Shaved LADs: shared but >10 % of het LAD lost relative to null coverage
threshold    <- 0.10
ov_df        <- as.data.frame(overlaps_lads)
het_to_null  <- split(ov_df$subjectHits, ov_df$queryHits)

frac_lost <- numeric(length(het_shared_idx))

for (k in seq_along(het_shared_idx)) {
  het_idx      <- het_shared_idx[k]
  het_range    <- het_lads[het_idx]
  null_indices <- het_to_null[[as.character(het_idx)]]
  
  if (length(null_indices) == 0) { frac_lost[k] <- 1; next }
  
  null_ranges    <- null_lads[null_indices]
  shared_regions <- GenomicRanges::reduce(
    pintersect(rep(het_range, length(null_ranges)), null_ranges)
  )
  frac_lost[k] <- 1 - sum(width(shared_regions)) / width(het_range)
}

shaved_idx   <- which(frac_lost > threshold)
shaved_sizes <- width(het_lads[het_shared_idx[shaved_idx]])

plot_df <- data.frame(
  size_Mb  = c(lost_sizes, shaved_sizes) / 1e6,
  Category = rep(
    c("Complete LAD loss", "LAD loss at edge"),
    c(length(lost_sizes), length(shaved_sizes))
  )
)

# Plot
ggplot(plot_df, aes(x = Category, y = size_Mb, fill = Category)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = c("#909090", "#4a90b8")) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = c("#909090", "#4a90b8")) +
  labs(x = NULL, y = "LAD size (log, Mb)") +
  theme_bw() +
  theme(legend.position = "none")
