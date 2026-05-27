library(Seurat)
library(dplyr)
library(ggplot2)

# Do not subset out any clusters prior to running.

# Normalize and scale both datasets
ShenJ_TC <- SCTransform(ShenJ_TC, verbose = FALSE)
P05Heart.combined <- SCTransform(P05Heart.combined, verbose = FALSE)
Feng_W_TC <- SCTransform(C57BL6_seurat_object, verbose = FALSE)

# Perform CCA for integration
integration_anchors <- FindIntegrationAnchors(
  object.list = list(ShenJ_TC, P05Heart.combined,C57BL6_seurat_object),
  dims = 1:12  # Adjust based on your PCA dimensions
)

# Integrate the datasets
integrated_data <- IntegrateData(anchorset = integration_anchors, dims = 1:12)

# Run PCA and UMAP
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, dims = 1:12)
integrated_data <- FindNeighbors(integrated_data, dims = 1:12)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, dims = 1:12)

# Visualize
DefaultAssay(integrated_data) <- "SCT"
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)  # Visualize different datasets

integrated_data$merged_cluster <- dplyr::recode(
  integrated_data$seurat_clusters,
  `0` = "Fibroblast", `6` = "Fibroblast", `7` = "Fibroblast", `2` = "Atria_CM", `15` = "Atria_CM",
  `9` = "Fibroblast", `14` = "Fibroblast", `18` = "Fibroblast", `1` = "Vent_CM", `13` = "Vent_CM",
  `12` = "Epicardial", `22` = "Epicardial", `5` = "Vent_CM", `19` = "Vent_CM", `8` = "Vent_CM",
  `10` = "Immune",
  `17` = "Blood",
  `4` = "Vasc_ec", `11` = "Vasc_ec", `20` = "Vasc_ec", `21` = "Vasc_ec",
  `3` = "Endocardial_ec", `16` = "Endocardial_ec"
)

integrated_data$seurat_clusters <- as.character(integrated_data$seurat_clusters)
