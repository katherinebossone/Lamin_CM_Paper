library(Seurat)
library(dplyr)
library(ggplot2)

Idents(integrated_data) <- "seurat_clusters"

# Combine clusters into one identity
integrated_data$combined_cluster <- ifelse(
  integrated_data$seurat_clusters %in% c("0","1","2","3","4","6","7","8"),
  "cluster_compact",
  as.character(integrated_data$seurat_clusters)
)

# Create a combined label that merges the cluster group and genotype
integrated_data$cluster_genotype <- paste(
  integrated_data$combined_cluster,
  integrated_data$genotype,
  sep = "_"
)

# Set that as the active identity
Idents(integrated_data) <- "cluster_genotype"

# SCT find markers
integrated_data <- PrepSCTFindMarkers(integrated_data)

# Run FindMarkers comparing the same combined cluster across genotypes
markers_cluster_compact <- FindMarkers(
  object = integrated_data,
  ident.1 = "cluster_compact_LB1Fd",  # combined cluster in genotype 1
  ident.2 = "cluster_compact_LB1Fplus",     # combined cluster in genotype 2
  assay = "SCT",
  min.pct = 0.25,
  logfc.threshold = 0
)

write.table(
  markers_cluster_compact,
  file = "markers_combined_clusters_all_compact_vent_LB1Fd_vs_LB1Fplus.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

# Convert Seurat output to a plotting dataframe
df <- markers_cluster_compact
df$gene <- rownames(df)

#Plot
fdr_thresh <- 0.05
lfc_thresh <- 0.25  

df$color <- "gray"
df$color[df$p_val_adj < fdr_thresh & df$avg_log2FC >  lfc_thresh] <- "red"
df$color[df$p_val_adj < fdr_thresh & df$avg_log2FC < -lfc_thresh] <- "blue"

df$neg_log10_FDR <- -log10(df$p_val_adj)

ggplot(df, aes(x = logFC, y = neg_log10_FDR, color = color)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_identity() +
  
  geom_hline(yintercept = -log10(fdr_thresh),
             linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", linewidth = 0.4) +
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj p-value)",
    title = "Volcano Plot: cluster_compact_LB1Fd vs LB1Fplus"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
