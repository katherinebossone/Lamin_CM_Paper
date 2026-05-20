library(rtracklayer)
library(ggplot2)
library(ggrastr)

# Load files
cscore_file <- "LB1Fplus_cscore.bedgraph"   
lad_file    <- "LB1het_LAC_zscore.bedgraph"     

cscore_data <- import(cscore_file, format = "bedGraph")
lad_data    <- import(lad_file, format = "bedGraph")

df_cscore <- as.data.frame(cscore_data)
df_lad    <- as.data.frame(lad_data)

# Shift cscore start positions by +1
df_cscore$start <- df_cscore$start + 1
df_lad$start <- df_lad$start -1

# Merge by genomic position
merged_df <- merge(df_cscore, df_lad, by = c("seqnames", "start"), suffixes = c("_cscore", "_lad"))


ggplot(merged_df, aes(x = score_cscore, y = score_lad)) +
  stat_bin_2d(bins = 300) +
  scale_fill_gradientn(
    colors = c("darkblue", "green"),
    limits = c(0, 20),  # clamp density counts
    oob = scales::squish, # values outside limits will be squished to min/max
    name = "Density"
  ) +
  labs(title = "Scatterplot of cscore vs lad", x = "cScore", y = "Lamin-A signal (zscore)") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "right")
