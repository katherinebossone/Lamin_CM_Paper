library(rtracklayer)  
library(GenomicRanges)
library(dplyr)         
library(tidyr)
library(ggplot2)


# Read hiland file

hilands <- read.table("state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr", "start", "end", "state")

hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges   = IRanges(start = hilands$start, end = hilands$end),
  state    = hilands$state
)

cscore_het <- import("LB1Fplus_modified_cscore_negated.bedgraph")

# Aggregate signal per HiLands region
aggregate_signal <- function(cscore_gr, hilands_gr) {
  
  hits <- findOverlaps(hilands_gr, cscore_gr)
  
  df <- data.frame(
    hiland_id = queryHits(hits),
    cscore_score = mcols(cscore_gr)$score[subjectHits(hits)]
  )
  
  summary <- df %>%
    group_by(hiland_id) %>%
    summarise(mean_signal = mean(cscore_score, na.rm = TRUE)) %>%
    ungroup()
  
  # attach state info
  summary$state <- hilands_gr$state[summary$hiland_id]
  summary
}

# Aggregate for both conditions
het_summary  <- aggregate_signal(cscore_het,  hilands_gr)

# Map the wt condition
het_summary$state <- factor(
  het_summary$state,
  levels = c("P", "B", "G","Y","O","R")
)

ggplot(het_summary, aes(x = state, y = mean_signal, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c(
    "P"="#a020f0","B"="#0000ff","G"="#006400",
    "Y"="#ffd700","O"="#ffa500","R"="#ff0000"
  )) +
  theme_minimal() +
  ylab("LB1het Cscore") +
  xlab("HiLands state")
