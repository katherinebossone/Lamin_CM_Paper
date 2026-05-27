library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyverse)


# Load files
het_bed  <- "../LB1het_LAC_2state_S2.bed"   
null_bed <- "../LB1null_LAC_2state_S2.bed"    

cscore_het_file  <- "LB1Fplus_modified_cscore.bedgraph"  
cscore_null_file <- "LB1Fd_modified_cscore.bedgraph"     

het  <- import(het_bed, format = "BED")
null <- import(null_bed, format = "BED")

# Adjust start positions for proper comparison
start(het)  <- start(het) - 2
start(null) <- start(null) - 2

if (inherits(het, "GRangesList"))  het  <- unlist(het, use.names = FALSE)
if (inherits(null, "GRangesList")) null <- unlist(null, use.names = FALSE)

# Keep only autosomal chromosomes
keep_chrs <- paste0("chr", 1:19)
het  <- het[seqnames(het) %in% keep_chrs]
null <- null[seqnames(null) %in% keep_chrs]

# Merge adjacent bins into continuous LADs and exclude first 3Mb of each chromosome
het  <- GenomicRanges::reduce(het, min.gapwidth = 1)
null <- GenomicRanges::reduce(null, min.gapwidth = 1)

all_chrs <- unique(c(seqlevels(het), seqlevels(null)))
exclude <- GRanges(seqnames = all_chrs,
                   ranges = IRanges(start = 1, end = 3e6))
het  <- het[!overlapsAny(het, exclude)]
null <- null[!overlapsAny(null, exclude)]

# Compute unique het LAD regions
het_unique <- GenomicRanges::setdiff(het, null)

# Import cScore bedGraphs
cscore_het  <- import(cscore_het_file, format = "bedGraph")
cscore_null <- import(cscore_null_file, format = "bedGraph")


# Subset c-scores by unique het LAD coordinates
het_bins  <- subsetByOverlaps(cscore_het,  het_unique)
null_bins <- subsetByOverlaps(cscore_null, het_unique)

hits <- findOverlaps(het_bins, null_bins, type = "equal")

df_unique <- data.frame(
  het_score  = mcols(het_bins)$score[queryHits(hits)],
  null_score = mcols(null_bins)$score[subjectHits(hits)]
)

# Compute difference
eps <- 1e-6
df_unique$delta <- df_unique$null_score - df_unique$het_score

# Plot
plot_unique <- ggplot(df_unique, aes(x = "", y = delta)) +
  geom_boxplot(fill = "steelblue", alpha = 0.8, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-0.6, 0.5)) +
  theme_minimal() +
  ylab("Delta Cscore (LB1null − LB1het)") +
  xlab("") +
  ggtitle("Cscore change in unique het LAD regions")

print(plot_unique)

# Wilcoxon signed-rank test
wilcox_unique <- wilcox.test(df_unique$log2FC, mu = 0)
cat("\nWilcoxon signed-rank test:\n")
print(wilcox_unique)
