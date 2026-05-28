library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

keep_chrs <- paste0("chr", 1:19)


# Load Hi-C data and scale
load_hic_data <- function(file_path) {
  hic <- fread(file_path, header = FALSE)
  setnames(hic, c("chr1","start1","end1","chr2","start2","end2","Count_raw","Count"))
  hic <- hic[chr1 == chr2 & chr1 %in% keep_chrs]
  return(hic)
}

scale_counts <- function(hic) {
  total <- sum(hic$Count, na.rm = TRUE)
  hic$Count_scaled <- hic$Count / total
  return(hic)
}

het1 <- load_hic_data("P0_5_LB1Fplus_8450_2_5kb.tsv")
het2 <- load_hic_data("P0_5_LB1Fplus_8450_4_5kb.tsv")
null1 <- load_hic_data("P0_5_LB1Fd_8450_7_5kb.tsv")
null2 <- load_hic_data("P0_5_LB1Fd_8450_8_5kb.tsv")

het1 <- scale_counts(het1)
het2 <- scale_counts(het2)
null1 <- scale_counts(null1)
null2 <- scale_counts(null2)

# Load HiLands states
hilands <- fread("state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr","start","end","state")
hilands <- hilands[chr %in% keep_chrs]
hilands$start <- as.numeric(hilands$start)
hilands$end   <- as.numeric(hilands$end)
hilands_gr <- GRanges(seqnames = hilands$chr,
                      ranges = IRanges(start = hilands$start, end = hilands$end),
                      state = hilands$state)

# Assign states to Hi-C bins
assign_states <- function(hic_df, hilands_gr) {
  gr1 <- GRanges(seqnames = hic_df$chr1, ranges = IRanges(start = hic_df$start1, end = hic_df$end1))
  gr2 <- GRanges(seqnames = hic_df$chr2, ranges = IRanges(start = hic_df$start2, end = hic_df$end2))
  overlap1 <- findOverlaps(gr1, hilands_gr)
  overlap2 <- findOverlaps(gr2, hilands_gr)
  hic_df$state1 <- NA_character_
  hic_df$state2 <- NA_character_
  hic_df$state1[queryHits(overlap1)] <- mcols(hilands_gr)$state[subjectHits(overlap1)]
  hic_df$state2[queryHits(overlap2)] <- mcols(hilands_gr)$state[subjectHits(overlap2)]
  hic_df <- hic_df[!is.na(hic_df$state1) & !is.na(hic_df$state2), ]
  return(hic_df)
}

het1 <- assign_states(het1, hilands_gr)
het2 <- assign_states(het2, hilands_gr)
null1 <- assign_states(null1, hilands_gr)
null2 <- assign_states(null2, hilands_gr)

# Summarize interactions per state pair
summarize_state_interactions <- function(hic_df) {
  hic_df %>%
    mutate(stateA = pmin(state1, state2),
           stateB = pmax(state1, state2)) %>%
    group_by(stateA, stateB) %>%
    summarise(total = sum(Count_scaled, na.rm = TRUE), .groups = "drop")
}

het1_s <- summarize_state_interactions(het1)
het2_s <- summarize_state_interactions(het2)
null1_s <- summarize_state_interactions(null1)
null2_s <- summarize_state_interactions(null2)

# Merge replicates and calculate log2FC
het <- full_join(het1_s, het2_s, by = c("stateA","stateB"), suffix = c("_1","_2")) %>%
  mutate(mean_het = rowMeans(select(., total_1, total_2), na.rm = TRUE)) %>%
  select(stateA, stateB, mean_het)

null <- full_join(null1_s, null2_s, by = c("stateA","stateB"), suffix = c("_1","_2")) %>%
  mutate(mean_null = rowMeans(select(., total_1, total_2), na.rm = TRUE)) %>%
  select(stateA, stateB, mean_null)

heatmap_df <- full_join(het, null, by = c("stateA","stateB")) %>%
  mutate(log2FC = log2(mean_null / mean_het)*100)

# Mirror matrix
heatmap_mirror <- heatmap_df
heatmap_mirror$stateA <- heatmap_df$stateB
heatmap_mirror$stateB <- heatmap_df$stateA
heatmap_full <- rbind(heatmap_df, heatmap_mirror)

# Plot
state_levels <- c("P","B","G","Y","O","R")
heatmap_full$stateA <- factor(heatmap_full$stateA, levels = state_levels)
heatmap_full$stateB <- factor(heatmap_full$stateB, levels = state_levels)

ggplot(heatmap_full, aes(x = stateA, y = stateB, fill = log2FC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(log2FC, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_fixed() +
  labs(title = "Global HiLands Interaction Changes (LB1null vs LB1het)",
       x = "State A", y = "State B", fill = "log2FC * 100") +
  theme_minimal()
