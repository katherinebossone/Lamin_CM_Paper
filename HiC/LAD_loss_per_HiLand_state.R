library(rtracklayer)  
library(GenomicRanges)
library(dplyr)         
library(tidyr)
library(ggplot2)


# Read hiland file

hilands <- read.table("Desktop/state_merged_sorted.bed", header = FALSE)
colnames(hilands) <- c("chr", "start", "end", "state")

hilands_gr <- GRanges(
  seqnames = hilands$chr,
  ranges   = IRanges(start = hilands$start, end = hilands$end),
  state    = hilands$state
)

lamin_het <- import("Desktop/LB1het_LAC_zscore.bedgraph")
lamin_null <- import("Desktop/LB1null_LAC_zscore.bedgraph")

# Aggregate signal per HiLands region
aggregate_signal <- function(atac_gr, hilands_gr) {
  
  hits <- findOverlaps(hilands_gr, atac_gr)
  
  df <- data.frame(
    hiland_id = queryHits(hits),
    atac_score = mcols(atac_gr)$score[subjectHits(hits)]
  )
  
  summary <- df %>%
    group_by(hiland_id) %>%
    summarise(mean_signal = mean(atac_score, na.rm = TRUE)) %>%
    ungroup()
  
  # attach state info
  summary$state <- hilands_gr$state[summary$hiland_id]
  summary
}

# Aggregate for both conditions
het_summary  <- aggregate_signal(lamin_het,  hilands_gr)
null_summary <- aggregate_signal(lamin_null, hilands_gr)

# Combine and compute log2 fold-change

het_summary2 <- het_summary %>%
  mutate(het_signal = mean_signal) %>%
  dplyr::select(hiland_id, state, het_signal)

null_summary2 <- null_summary %>%
  mutate(null_signal = mean_signal) %>%
  dplyr::select(hiland_id, state, null_signal)

combined <- het_summary2 %>%
  left_join(null_summary2, by = c("hiland_id", "state")) %>%
  mutate(
    log2FC = log2(null_signal / het_signal)
  )

# Summarize by HiLands state

state_summary <- combined %>%
  group_by(state) %>%
  summarise(
    mean_het     = mean(het_signal, na.rm = TRUE),
    mean_null    = mean(null_signal, na.rm = TRUE),
    mean_log2FC  = mean(log2FC, na.rm = TRUE)
  )


# summary statistics per state

state_summary_stats <- combined %>%
  group_by(state) %>%
  summarise(
    n            = n(),
    mean_het     = mean(het_signal, na.rm = TRUE),
    mean_null    = mean(null_signal, na.rm = TRUE),
    median_het   = median(het_signal, na.rm = TRUE),
    median_null  = median(null_signal, na.rm = TRUE),
    sd_het       = sd(het_signal, na.rm = TRUE),
    sd_null      = sd(null_signal, na.rm = TRUE),
    mean_log2FC  = mean(log2FC, na.rm = TRUE)
  )

state_summary_stats


# Wilcoxon test per state

wilcox_stats <- combined %>%
  group_by(state) %>%
  summarise(
    p_value = wilcox.test(null_signal, het_signal, paired = FALSE)$p.value
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "fdr")
  )

wilcox_stats

# plot

combined$state <- factor(
  combined$state,
  levels = c("P", "B", "G","Y","O","R")
)
combined_BP <- combined %>%
  filter(state %in% c("B", "P"))
ggplot(combined_BP, aes(x = state, y = log2FC, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c(
    "P"="#a020f0","B"="#0000ff","G"="#006400",
    "Y"="#ffd700","O"="#ffa500","R"="#ff0000"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-2.5, 2.5) +
  theme_minimal() +
  ylab("log2(LB1null / LB1het)") +
  xlab("HiLands state")
