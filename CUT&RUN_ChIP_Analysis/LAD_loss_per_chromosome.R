library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)

# BED files locations
beds <- list(
  LAC_wt  = "LB1Fplus_LAC_2state_S2plus.bed",
  LAC_mut = "LB1Fminus_LAC_2state_S2plus.bed",
)

# mm10 chromosome sizes
mm10_chrom_sizes <- data.frame(
  chr = paste0("chr", 1:19),
  size = c(
    195471971, 182113224, 160039680, 156508116, 151834684,
    149736546, 145441459, 129401213, 124595110, 130694993,
    122082543, 120129022, 120421639, 124902244, 104043685,
    98207768, 94987271, 90702639, 61431566
  )
)

# Calculate per-chromosome LAD coverage
get_chr_coverage <- function(bed_file, chrom_sizes) {

  gr <- import(bed_file)
  gr <- gr[seqnames(gr) %in% chrom_sizes$chr]

  cov <- coverage(gr)

  cov_df <- data.frame(
    chr = chrom_sizes$chr,
    covered_bp = sapply(chrom_sizes$chr, function(chr) {
      if (chr %in% names(cov)) {
        sum(cov[[chr]] > 0)
      } else {
        0
      }
    }),
    chr_size = chrom_sizes$size
  )

  cov_df$coverage_fraction <- cov_df$covered_bp / cov_df$chr_size

  cov_df
}

## Calculate LAD coverage for all samples
coverage_df <- bind_rows(
  lapply(names(beds), function(name) {
    df <- get_chr_coverage(beds[[name]], mm10_chrom_sizes)
    df$sample <- name
    df
  })
)

# Calculate per-chromosome % LAD loss for Lamin A/C
lac_loss <- coverage_df %>%
  filter(sample %in% c("LAC_wt", "LAC_mut")) %>%
  select(chr, sample, coverage_fraction) %>%
  pivot_wider(
    names_from = sample,
    values_from = coverage_fraction
  ) %>%
  mutate(
    pct_loss = 100 * (LAC_wt - LAC_mut) / LAC_wt
  )

## Print average loss
mean_loss <- mean(lac_loss$pct_loss, na.rm = TRUE)

cat(round(mean_loss, 2))

# Plot
p <- ggplot(lac_loss, aes(x = "", y = pct_loss)) +
  geom_boxplot(fill = "magenta", alpha = 0.6, width = 0.35, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3, alpha = 0.9) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(x = NULL, y = "% LAD coverage loss per chromosome") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank())

print(p)
