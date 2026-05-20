library(dplyr)
library(ggplot2)
library(scales)

# Prior to running R script run hicPlotDistVsCounts to get the correctly normalized data
# hicPlotDistVsCounts   -m "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_2.mcool::resolutions/10000" "P0_5_HiC/indiv/P0_5_LB1Fplus_8450_4.mcool::resolutions/10000" "P0_5_HiC/indiv/P0_5_LB1Fd_8450_7.mcool::resolutions/10000" "P0_5_HiC/indiv/P0_5_LB1Fd_8450_8.mcool::resolutions/10000"  -o distance_vs_counts.pdf   --labels "LB1+/- rep1" "LB1 +/- rep2" "LB1-/- rep1" "LB1-/- rep2"   --maxdepth 100000000   --outFileData distance_data.tsv

# Read TSV output from hicPlotDistVsCounts
raw <- read.table("distance_data.tsv", header = TRUE, sep = "\t")

# Drop stray header rows and fix spacing
raw <- raw[raw$Matrix != "Matrix", ]
raw$Matrix <- trimws(raw$Matrix)

# Keep only chr1-19 (and "all" rows which are the genome-wide aggregates)
raw <- raw[raw$Chromosome %in% c(paste0("chr", 1:19), "all"), ]

# Rename columns
df_long <- raw %>%
  dplyr::rename(label = Matrix, distance = Distance, contact_freq = Contacts) %>%
  filter(distance > 0)

# Fix column types
df_long$contact_freq <- as.numeric(df_long$contact_freq)
df_long$distance <- as.numeric(df_long$distance)

# Assign groups for coloring
df_long$group <- ifelse(grepl("\\+/-", df_long$label), "LB1+/-", "LB1-/-")

# Library-size normalize within each sample
df_long <- df_long %>%
  group_by(label) %>%
  mutate(norm_by_lib_size = contact_freq / sum(contact_freq, na.rm = TRUE)) %>%
  ungroup()

# Log-bin smoothing
bins <- sapply(0:100, function(i) 10^(i / 25)) * 10000

df_smoothed <- df_long %>%
  group_by(label, group, grp = findInterval(distance, bins)) %>%
  summarise(val = mean(norm_by_lib_size, na.rm = TRUE), .groups = "drop") %>%
  mutate(Distance = sapply(grp, function(i) 10^(i / 25) * 10000))

# Plot
label_order <- c("LB1+/- rep1", "LB1 +/- rep2", "LB1-/- rep1", "LB1-/- rep2")
color_values <- setNames(
  c("#2166AC", "#6BAED6", "#D7191C", "#F4A582"))

df_smoothed$label <- factor(df_smoothed$label)

ggplot(df_smoothed, aes(x = Distance, y = val, color = label)) +
  geom_line(linewidth = 0.6) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10000, 1e8),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(base = 10, sides = "bl", linewidth = 0.3) +
  scale_color_manual(values = color_values) +
  labs(x = "Genomic distance (bp)", y = "Normalized contact frequency", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(size = 12, color = "black"),
    axis.title       = element_text(size = 13, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    legend.position  = "right"
  )
