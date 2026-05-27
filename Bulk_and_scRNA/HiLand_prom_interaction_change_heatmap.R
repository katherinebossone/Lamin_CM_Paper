library(pheatmap)
library(dplyr)
library(gridExtra)

# Use final_counts file from HiLand_promoter_interaction_change.R scripe

# Parameters
FOCUS_STATE <- "Y"      # "P", "B", "G", "Y", "O", "R"
COLOR_LIMIT <- 0.5

heatmap_colors <- colorRampPalette(c("blue3", "white", "red3"))(101)
breaks <- seq(-COLOR_LIMIT, COLOR_LIMIT, length.out = 102)

# Build matrix (ordered by log2FC)
build_matrix <- function(data, reg, focus_state) {
  mat <- data %>%
    filter(regulation == reg, flank_state == focus_state) %>%
    group_by(gene) %>%
    summarise(fold_change = mean(fold_change, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  
  rownames(mat) <- mat$gene
  mat <- as.matrix(mat[, "fold_change", drop = FALSE])
  
  # clamp values
  mat[mat >  COLOR_LIMIT] <-  COLOR_LIMIT
  mat[mat < -COLOR_LIMIT] <- -COLOR_LIMIT
  
  # highest upregulated at top
  mat <- mat[order(mat[,1], decreasing = TRUE), , drop = FALSE]
  
  return(mat)
}

# Plot 
plot_heatmap <- function(mat, title) {
  if (nrow(mat) == 0) {
    message("No data for: ", title)
    return(NULL)
  }
  
  p <- pheatmap(
    mat,
    color         = heatmap_colors,
    breaks        = breaks,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    border_color  = NA,
    main          = title,
    na_col        = "#DDDDDD"
  )
  
  return(p)
}

# Build matrices
mat_up   <- build_matrix(final_counts, "Up",   FOCUS_STATE)
mat_down <- build_matrix(final_counts, "Down", FOCUS_STATE)

# Generate heatmaps
Promup <- plot_heatmap(
  mat_up,
  title = paste0("Upregulated — state ", FOCUS_STATE, " | log2(Mut/WT)")
)

Promdown <- plot_heatmap(
  mat_down,
  title = paste0("Downregulated — state ", FOCUS_STATE, " | log2(Mut/WT)")
)
