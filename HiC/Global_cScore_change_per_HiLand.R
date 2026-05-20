library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)


# Load files
hilands <- read.table("state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr", "start", "end", "state")

hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges   = IRanges(hilands$start, hilands$end),
  state    = hilands$state
)

wt  <- import("LB1Fplus_cscore.bedgraph")
mut <- import("LB1Fd_cscore.bedgraph")

# Aggregate cScore per HiLands region
get_region_score <- function(cscore_gr, hilands_gr) {
  
  hits <- findOverlaps(hilands_gr, cscore_gr)
  
  df <- data.frame(
    region_id = queryHits(hits),
    score     = mcols(cscore_gr)$score[subjectHits(hits)]
  )
  
  df %>%
    dplyr::group_by(region_id) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(state = hilands_gr$state[region_id])
}

wt_sum  <- get_region_score(wt, hilands_gr)  %>% dplyr::rename(wt_score  = mean_score)
mut_sum <- get_region_score(mut, hilands_gr) %>% dplyr::rename(mut_score = mean_score)

wt_sum$wt_score   <- -wt_sum$wt_score
mut_sum$mut_score <- -mut_sum$mut_score


# Merge + compute delta cScore
combined <- inner_join(wt_sum, mut_sum, by = c("region_id", "state")) %>%
  mutate(
    delta_cscore = mut_score - wt_score
  )

# Summary statistics per state
combined$state <- factor(
  combined$state,
  levels = c("P", "B", "G", "Y", "O", "R")
)

state_summary <- combined %>%
  group_by(state) %>%
  summarise(
    n = n(),
    mean_delta = mean(delta_cscore, na.rm = TRUE),
    median_delta = median(delta_cscore, na.rm = TRUE)
  )

print(state_summary)

# Wilcoxon test per state
wilcox_stats <- combined %>%
  group_by(state) %>%
  summarise(
    p_value = wilcox.test(delta_cscore)$p.value
  ) %>%
  mutate(q_value = p.adjust(p_value, method = "fdr"))

print(wilcox_stats)

# Plot 
ggplot(combined, aes(x = state, y = delta_cscore, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(
    "P"="#a020f0","B"="#0000ff","G"="#006400",
    "Y"="#ffd700","O"="#ffa500","R"="#ff0000"
  )) +
  theme_minimal() +
  ylim(-0.3, 0.3) +
  ylab("(mutant − WT)") +
  xlab("HiLands state")
