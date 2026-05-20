# Load libraries
library(ggplot2)
library(scales) 


# Use output from chromatinHMM modeling
`final_data_state` <- read.delim("hilands/CM_hiland_states_6_5000.txt", header=FALSE, stringsAsFactors=FALSE)

#Use combined R object from Chromatin_marker_PCA_for_HiLands.R script
bed_combined <- data.frame(
  chrom = combined[1],
  start = combined[[2]],
  end   = combined[[3]],
  state = final_data_state[[1]]
)

df <- merge(bed_combined, combined, by = c("chrom", "start", "end"))

# Define markers and states
marker_cols <- c("H3K9me2_score","H3K9me3_score","H3K27me3_score","H3K4me1_score", 
                 "H3K4me3_score", "H3K27ac_score", "ATAC_score", "LAD_score")
states <- unique(bed_combined$state)

state_levels <- c("0", "4", "5", "2", "3", "1")
marker_levels <- c("LAD_score","H3K9me2_score","H3K9me3_score","H3K27me3_score", 
                   "H3K4me1_score", "H3K4me3_score", "H3K27ac_score","ATAC_score")


# z-score the signal per marker
for (marker in marker_cols) {
  df[[paste0(marker, "_z")]] <- scale(df[[marker]])
}

results_z <- list()

for (marker in marker_cols) {
  z_marker <- paste0(marker, "_z")
  for (s in states) {
    state_rows <- df[df$state == s, ]
    other_rows <- df[df$state != s, ]
    
    mean_in_state <- mean(state_rows[[z_marker]], na.rm = TRUE)
    mean_in_others <- mean(other_rows[[z_marker]], na.rm = TRUE)
    
    delta_z <- mean_in_state-mean_in_others # difference in standardized signal
    
    results_z[[length(results_z) + 1]] <- data.frame(
      state = s,
      marker = marker,
      delta_z = delta_z
    )
  }
}

df_z <- do.call(rbind, results_z)

df_z$state <- factor(df_z$state, levels = state_levels)
df_z$marker <- factor(df_z$marker, levels = marker_levels)

# plot
ggplot(df_z, aes(x = state, y = marker, fill = delta_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "red3",
                       midpoint = 0, , limits = c(-1, 1), oob = squish,
                       name = "Z-score") +
  theme_minimal() +
  labs(x = "ChromHMM State", y = "Marker") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
