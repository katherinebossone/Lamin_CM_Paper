library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Import files
cscore_het  <- import("Documents/LB1Fplus_cscore.bedgraph", format = "bedgraph")
cscore_null <- import("Documents/LB1Fd_cscore.bedgraph", format = "bedgraph")

# Find overlapping regions between the two conditions
hits <- findOverlaps(cscore_het, cscore_null)

# Build combined dataframe of overlapping region scores
combined <- data.frame(
  het_signal  = mcols(cscore_het)$score[queryHits(hits)],
  null_signal = mcols(cscore_null)$score[subjectHits(hits)]
) %>%
  mutate(delta = null_signal - het_signal)

# Summary statistics
summary_stats <- combined %>%
  summarise(
    n           = n(),
    mean_het    = mean(het_signal,  na.rm = TRUE),
    mean_null   = mean(null_signal, na.rm = TRUE),
    median_het  = median(het_signal,  na.rm = TRUE),
    median_null = median(null_signal, na.rm = TRUE),
    sd_het      = sd(het_signal,  na.rm = TRUE),
    sd_null     = sd(null_signal, na.rm = TRUE),
    mean_delta = mean(delta, na.rm = TRUE)
  )
summary_stats

# Wilcoxon test (global)
wilcox.test(combined$null_signal, combined$het_signal, paired = FALSE)

# Plot 
ggplot(combined, aes(x = "", y = delta)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "steelblue", width = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-0.5, 0.5) +
  theme_minimal() +
  ylab("Delta Cscore LB1null - LB1het") +
  xlab("Global")
