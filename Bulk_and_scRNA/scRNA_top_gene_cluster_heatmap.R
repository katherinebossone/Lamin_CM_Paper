library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)

# Subset to non-S/G2/M ventricular compact cardiomyocytes
clusters_of_interest <- c(0, 1, 2, 3, 8)
Idents(integrated_data) <- "seurat_clusters"
subset_data <- subset(integrated_data, idents = clusters_of_interest)

# Filter to public cells
meta <- subset_data@meta.data
public_cells <- rownames(meta)[!is.na(meta$stage)]

# SCT find markers
subset_data <- PrepSCTFindMarkers(object = subset_data)
Idents(subset_data) <- "seurat_clusters"

markers <- FindAllMarkers(
  subset_data,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Top 20 marker genes per cluster
# Filter by adjusted p-value and rank by fold change
# Order by cluster so genes appear in biological order on heatmap
top_genes <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  arrange(factor(cluster, levels = c(8, 0, 2, 1, 3))) %>%
  pull(gene)

top_genes <- unique(top_genes)
top_genes <- top_genes[top_genes %in% rownames(subset_data)]

# Average expression
DefaultAssay(subset_data) <- "SCT"
avg_exp <- AverageExpression(
  subset_data,
  features = top_genes,
  group.by = "seurat_clusters",
  assays   = "SCT",
  slot     = "data"
)$SCT

# Reorder columns
col_order <- paste0("g", c(8, 0, 2, 1, 3))
avg_exp   <- avg_exp[, col_order]
colnames(avg_exp) <- c("8", "0", "2", "1", "3")

# Z-score scaling
avg_exp_scaled_uncapped <- t(scale(t(avg_exp)))
avg_exp_scaled <- pmin(pmax(avg_exp_scaled_uncapped, -2.5), 2.5)

# Preserve gene order from marker analysis
avg_exp_scaled         <- avg_exp_scaled[top_genes[top_genes %in% rownames(avg_exp_scaled)], ]
avg_exp_scaled_uncapped <- avg_exp_scaled_uncapped[rownames(avg_exp_scaled), ]

# Color scale 
col_fun <- colorRamp2(
  c(-2.5, -1.25, 0, 1.25, 2.5),
  rev(brewer.pal(5, "RdYlBu"))
)

# Assign each gene to the cluster where it has highest average z-score
gene_cluster_assignment <- apply(avg_exp_scaled_uncapped, 1, function(x) {
  colnames(avg_exp_scaled_uncapped)[which.max(x)]
})

# Order factor levels to match column order
row_split <- factor(
  gene_cluster_assignment,
  levels = c("8", "0", "2", "1", "3")
)

# Plot
ht <- Heatmap(
  avg_exp_scaled,
  name                  = "Z-score",
  col                   = col_fun,
  cluster_rows          = TRUE,
  clustering_distance_rows = function(m) dist(avg_exp_scaled_uncapped[rownames(m), ]),
  cluster_columns       = FALSE,
  row_split             = row_split,
  row_title_gp          = gpar(fontsize = 9, fontface = "bold"),
  row_gap               = unit(2, "mm"),
  show_row_dend         = TRUE,
  show_column_dend      = FALSE,
  show_row_names        = TRUE,
  show_column_names     = TRUE,
  row_names_gp          = gpar(fontsize = 6),
  column_names_gp       = gpar(fontsize = 11),
  column_names_rot      = 0,
  border                = TRUE,
  rect_gp               = gpar(col = NA),
  heatmap_legend_param  = list(
    title         = "Z-score",
    legend_height = unit(4, "cm")
  )
)

draw(ht, heatmap_legend_side = "right")
