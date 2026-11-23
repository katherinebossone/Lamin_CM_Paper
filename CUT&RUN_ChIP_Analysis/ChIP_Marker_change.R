library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)

# BED files locations
beds <- list(
  H3K27ac_wt  = "Desktop/E13_LB1Fplus_Nkx2_5Cre_H3K27ac_Ave_S1.bed",
  H3K27ac_mut = "Desktop/E13_LB1Fd_Nkx2_5Cre_H3K27ac_Ave_S1.bed",
  H3K9me2_wt  = "Desktop/E13_LB1Fplus_Nkx2_5Cre_H3K9me2_Ave_S2.bed",
  H3K9me2_mut = "Desktop/E13_LB1Fd_Nkx2_5Cre_H3K9me2_Ave_S2.bed"
)

mm10_chrom_sizes <- data.frame(
  chr = paste0("chr", 1:19),
  size = c(
    195471971, 182113224, 160039680, 156508116, 151834684,
    149736546, 145441459, 129401213, 124595110, 130694993,
    122082543, 120129022, 120421639, 124902244, 104043685,
    98207768, 94987271, 90702639, 61431566
  )
)

# Calculate per-chromosome marker coverage
get_chr_coverage <- function(bed_file, chrom_sizes, excl = 3e6) {
  gr <- import(bed_file)
  gr <- gr[seqnames(gr) %in% chrom_sizes$chr]
  keep <- GRanges(chrom_sizes$chr, IRanges(excl + 1, chrom_sizes$size))
  gr <- GenomicRanges::intersect(gr, keep, ignore.strand = TRUE)
  
  cov <- coverage(gr)
  
  cov_df <- data.frame(
    chr = chrom_sizes$chr,
    covered_bp = sapply(chrom_sizes$chr, function(chr) {
      if (chr %in% names(cov)) sum(cov[[chr]] > 0) else 0
    }),
    chr_size = chrom_sizes$size - excl
  )
  cov_df$coverage_fraction <- cov_df$covered_bp / cov_df$chr_size
  cov_df
}

# Calculate marker coverage for all samples
coverage_df <- bind_rows(
  lapply(names(beds), function(name) {
    df <- get_chr_coverage(beds[[name]], mm10_chrom_sizes)
    df$sample <- name
    df
  })
)

# Split sample name into isoform + genotype, so H3K27ac and H3K9me2 can be handled
coverage_df <- coverage_df %>%
  separate(sample, into = c("isoform", "genotype"), sep = "_")

# Calculate per-chromosome % marker change
change_df <- coverage_df %>%
  dplyr::select(chr, isoform, genotype, coverage_fraction) %>%
  pivot_wider(
    names_from = genotype,
    values_from = coverage_fraction
  ) %>%
  mutate(
    pct_change = 100 * (mut - wt) / wt
  )

# Print average change per isoform
mean_change <- change_df %>%
  group_by(isoform) %>%
  summarise(mean_pct_change = mean(pct_change, na.rm = TRUE))
print(mean_change)

# Wilcoxon signed-rank test
change_wide <- change_df %>%
  dplyr::select(chr, isoform, pct_change) %>%
  pivot_wider(names_from = isoform, values_from = pct_change)

wilcox_result_H3K27ac <- wilcox.test(change_wide$H3K27ac)
wilcox_result_H3K9me2 <- wilcox.test(change_wide$H3K9me2)
print(wilcox_result_H3K27ac)
print(wilcox_result_H3K9me2)

# Plot
ggplot(change_df, aes(x = isoform, y = pct_change, fill = isoform)) +
  geom_boxplot(alpha = 0.6, width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3, alpha = 0.9) +
  scale_fill_manual(values = c(H3K27ac = "magenta", H3K9me2 = "green")) +
  scale_y_continuous(limits = c(-40, 40), breaks = seq(-50, 50, 10)) +
  labs(x = NULL, y = "% LAD coverage change per chromosome") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
