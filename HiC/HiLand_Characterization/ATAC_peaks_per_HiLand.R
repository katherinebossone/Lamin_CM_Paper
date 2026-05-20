library(rtracklayer)   # for reading bigWig / bedGraph
library(GenomicRanges) # for genomic ranges and overlaps
library(dplyr)         # for data manipulation
library(tidyr)
library(ggplot2)

# Read files
hilands <- read.table("state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr", "start", "end", "state")

hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges   = IRanges(start = hilands$start, end = hilands$end),
  state    = hilands$state
)
peaks <- read.table("../LB1Het_ATAC_consensus_peaks.mLb.clN.bed", header = FALSE)
colnames(peaks) <- c("chr", "start", "end")

peaks_gr <- GRanges(
  seqnames = peaks$chr,
  ranges   = IRanges(start = peaks$start, end = peaks$end)
)

# Count peaks overlapping each HiLand region
hits <- findOverlaps(hilands_gr, peaks_gr)

peak_counts <- data.frame(
  hiland_id = queryHits(hits)
) %>%
  group_by(hiland_id) %>%
  summarise(n_peaks = n()) %>%
  ungroup() %>%
  mutate(state = hilands_gr$state[hiland_id])

# Summarize total peaks per state
state_peak_counts <- peak_counts %>%
  group_by(state) %>%
  summarise(total_peaks = sum(n_peaks))

state_peak_counts

# Plot
state_peak_counts$state <- factor(state_peak_counts$state,
                                  levels = c("P", "B", "G", "Y", "O", "R"))

# Calculate percentages for labels
state_peak_counts <- state_peak_counts %>%
  mutate(
    pct = total_peaks / sum(total_peaks) * 100,
    label = paste0(state, "\n", round(pct, 1), "%")
  )

ggplot(state_peak_counts, aes(x = "", y = total_peaks, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "P"="#a020f0","B"="#0000ff","G"="#006400",
    "Y"="#ffd700","O"="#ffa500","R"="#ff0000"
  )) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "ATAC peaks per HiLands state")
