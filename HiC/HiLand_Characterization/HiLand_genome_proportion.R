library(GenomicRanges) 
library(dplyr)        
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

# Calculate genome composition per state
state_sizes <- hilands %>%
  mutate(size = end - start) %>%
  group_by(state) %>%
  summarise(total_bp = sum(size)) %>%
  mutate(
    pct = total_bp / sum(total_bp) * 100,
    label = paste0(state, "\n", round(pct, 1), "%"),
    state = factor(state, levels = c("P", "B", "G", "Y", "O", "R"))
  )

ggplot(state_sizes, aes(x = "", y = total_bp, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Genome composition by HiLands state") +
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
