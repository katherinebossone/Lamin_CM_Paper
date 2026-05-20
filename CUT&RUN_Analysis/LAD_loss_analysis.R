library(rtracklayer)
library(GenomicRanges)
library(VennDiagram)
library(tidyverse)
library(grid)

# Load in the bed files from HMM
het_bed <- "../LB1het_LAC_2state_S2.bed"      
null_bed <- "../LB1null_LAC_2state_S2.bed"    
het <- import(het_bed, format = "BED")
null <- import(null_bed, format = "BED")
start(het) <- start(het) - 1
start(null) <- start(null) - 1

# Flatten GRange list
if (inherits(het, "GRangesList")) het <- unlist(het, use.names = FALSE)
if (inherits(null, "GRangesList")) null <- unlist(null, use.names = FALSE)
keep_chrs <- paste0("chr", 1:19)  # adjust if your seqnames do not have "chr" prefix
het <- het[seqnames(het) %in% keep_chrs]
null <- null[seqnames(null) %in% keep_chrs]

# Ensure all overlapping lads are continuous
cat("Merging adjacent bins into continuous LADs...\n")
het <- GenomicRanges::reduce(het, min.gapwidth = 1)
null <- GenomicRanges::reduce(null, min.gapwidth = 1)

# Exclude first 3Mb
cat("Removing first 3 Mb of each chromosome...\n")
all_chrs <- unique(c(seqlevels(het), seqlevels(null)))
exclude <- GRanges(seqnames = all_chrs,
                   ranges = IRanges(start = 1, end = 3e6))
het <- het[!overlapsAny(het, exclude)]
null <- null[!overlapsAny(null, exclude)]

# Find overlaps between het and null LADs
cat("Finding overlapping LADs...\n")
overlaps <- findOverlaps(het, null)

# Identify unique and shared lads
het_shared_idx <- unique(queryHits(overlaps))
null_shared_idx <- unique(subjectHits(overlaps))

het_unique <- het[-het_shared_idx]   
null_unique <- null[-null_shared_idx] 

shared_het <- het[het_shared_idx]
shared_null <- null[null_shared_idx]


## Venn diagram ##

# Number of unique/shared LADs
num_het <- length(het)
num_null <- length(null_unique)
num_shared_null <- length(null) - length(null_unique)  # nullant LADs overlapping het
num_shared_het  <- length(het) - length(het_unique)    # het LADs overlapping nullant

grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = num_het,
  area2 = num_null,
  cross.area = num_shared_het,   # Use actual shared nullant count
  category = c("het LADs", "nullant LADs"),
  fill = c("#4DBBD5FF", "#E64B35FF"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  ind = FALSE,
  print.mode = "raw"
)
grid.draw(venn.plot)


# Identify lost LADs 
lost_LADs <- het_unique
lost_sizes <- width(lost_LADs)
threshold <- 0.10

# Pre-split overlaps by het LAD index
ov_df <- as.data.frame(overlaps)
het_to_null <- split(ov_df$subjectHits, ov_df$queryHits)

frac_lost <- numeric(length(het_shared_idx))

for (k in seq_along(het_shared_idx)) {
  
  het_idx <- het_shared_idx[k]
  het_range <- het[het_idx]
  
  # nullant indices overlapping THIS wild-type LAD
  null_indices <- het_to_null[[as.character(het_idx)]]
  
  if (length(null_indices) == 0) {
    # Should not happen because these are shared LADs
    frac_lost[k] <- 1
    next
  }
  
  null_ranges <- null[null_indices]
  
  # Compute retained region
  shared_regions <- GenomicRanges::reduce(pintersect(rep(het_range, length(null_ranges)), null_ranges))
  retained_width <- sum(width(shared_regions))
  
  frac_lost[k] <- 1 - retained_width / width(het_range)
}

# Determine which LADs are shaved
shaved_idx <- which(frac_lost > threshold)
shaved_LADs <- het[het_shared_idx[shaved_idx]]
shaved_sizes <- width(shaved_LADs)

pct_shaved <- 100 * length(shaved_LADs) / length(het_shared_idx)
cat("Percentage of shaved LADs:", round(pct_shaved, 1), "%\n")


# Summarize LAD sizes
size_summary <- tibble(
  Category = c("Lost LADs", "Shaved LADs"),
  Median = c(median(lost_sizes), median(shaved_sizes)),
  Mean = c(mean(lost_sizes), mean(shaved_sizes)),
  Min = c(min(lost_sizes), min(shaved_sizes)),
  Max = c(max(lost_sizes), max(shaved_sizes)),
  Count = c(length(lost_sizes), length(shaved_sizes))
)

print(size_summary)

# Boxplot of LAD sizes
lost_sizes_Mb <- lost_sizes / 1e6
shaved_sizes_Mb <- shaved_sizes / 1e6

dev.new(width = 8, height = 6)
boxplot(lost_sizes_Mb, shaved_sizes_Mb,
        names = c("Lost LADs", "Shaved LADs"),
        ylab = "LAD size (Mb)",
        col = c("#757575", "#1f77a3"))
