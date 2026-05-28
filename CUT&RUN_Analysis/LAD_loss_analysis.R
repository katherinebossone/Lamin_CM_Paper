library(rtracklayer)
library(GenomicRanges)
library(VennDiagram)
library(tidyverse)
library(grid)

# Load and adjust files
het_bed  <- "Desktop/bigwig_lamins/LB1het_LAC_2state_S2.bed"
null_bed <- "Desktop/bigwig_lamins/LB1null_LAC_2state_S2.bed"

het  <- import(het_bed, format="BED")
null <- import(null_bed, format="BED")

start(het)  <- start(het) -1
start(null) <- start(null) -1

if (inherits(het,"GRangesList"))  het  <- unlist(het, use.names=FALSE)
if (inherits(null,"GRangesList")) null <- unlist(null, use.names=FALSE)

# Removing first 3 Mb of chromosomes and keep only autosomes

keep_chrs <- paste0("chr",1:19)

het  <- het[seqnames(het) %in% keep_chrs]
null <- null[seqnames(null) %in% keep_chrs]

exclude <- GRanges(
  seqnames = keep_chrs,
  ranges = IRanges(start=1, end=3e6)
)

het  <- het[!overlapsAny(het, exclude)]
null <- null[!overlapsAny(null, exclude)]

# Find overlapping bins

overlaps <- findOverlaps(het, null)

shared_bins_idx_het  <- unique(queryHits(overlaps))
shared_bins_idx_null <- unique(subjectHits(overlaps))

shared_bins <- length(shared_bins_idx_het)

# Identify unique bins

het_only_bins  <- length(het)  - shared_bins
null_only_bins <- length(null) - shared_bins

# Calculate total genome bins
bin_size <- 5000

mouse_genome_bp <- 2.65e9
total_bins <- mouse_genome_bp / bin_size

# Convert to genome percentage


shared_pct    <- round(100 * shared_bins    / total_bins, 2)
het_only_pct  <- round(100 * het_only_bins  / total_bins, 2)
null_only_pct <- round(100 * null_only_bins / total_bins, 2)

# Plot
venn.plot <- draw.pairwise.venn(
  area1 = het_only_pct + shared_pct,
  area2 = null_only_pct + shared_pct,
  cross.area = shared_pct,
  category = c("Het LAD genome","Null LAD genome"),
  fill = c("#4DBBD5FF","#E64B35FF"),
  alpha = 0.6,
  cex = 2,
  cat.cex = 1.5,
  ind = FALSE
)

grid.draw(venn.plot)
