library(GenomicRanges)
library(dplyr)
library(ggplot2)

## Promoter Density per HiLands State ##


mcols(promoters)$gene <- mcols(promoters)$gene_name  # add 'gene' column

# Find overlaps between promoters and HiLands states
hits <- findOverlaps(promoters, hilands_gr)

prom_states <- data.frame(
  gene  = mcols(promoters)$gene[queryHits(hits)],
  state = mcols(hilands_gr)$state[subjectHits(hits)],
  stringsAsFactors = FALSE
)

# Count promoters per state
prom_counts <- prom_states %>%
  group_by(state) %>%
  summarise(
    n_promoters = dplyr::n_distinct(gene),
    .groups = "drop"
  )

# Compute total genomic size of each HiLands state
state_sizes <- data.frame(
  state = mcols(hilands_gr)$state,
  width = width(hilands_gr)
) %>%
  group_by(state) %>%
  summarise(
    total_bp = sum(width),
    .groups = "drop"
  )

# Calculate promoter density (per Mb)
prom_density <- prom_counts %>%
  left_join(state_sizes, by = "state") %>%
  mutate(
    promoters_per_mb = n_promoters / (total_bp / 1e6)
  )

# Plot
prom_density$state <- factor(prom_density$state, levels = c("P","B","G","Y","O","R"))

ggplot(prom_density, aes(x = state, y = promoters_per_mb, fill = state)) +
  geom_bar(stat = "identity") +
  labs(
    x = "HiLands State",
    y = "Promoters per Mb",
    title = "Promoter Density Across HiLands States"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "P" = "#a020f0",
    "B" = "#0000ff",
    "G" = "#006400",
    "Y" = "#ffd700",
    "O" = "#ffa500",
    "R" = "#ff0000"
  ))

## Pie chart of promoter distribution across hiLands states ##

# Compute fraction of promoters per state
prom_pie <- prom_counts %>%
  mutate(
    fraction = n_promoters / sum(n_promoters),
    label = paste0(state, " (", round(fraction*100,1), "%)")
  )

# Plot
ggplot(prom_pie, aes(x = "", y = fraction, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Promoter Distribution Across HiLands States") +
  theme_void() + # remove axes
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size=3) +
  scale_fill_manual(values = c(
    "P" = "#a020f0",
    "B" = "#0000ff",
    "G" = "#006400",
    "Y" = "#ffd700",
    "O" = "#ffa500",
    "R" = "#ff0000"
  ))
