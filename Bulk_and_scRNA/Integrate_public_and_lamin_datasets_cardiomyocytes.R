library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

# Non-cardiomyocytes already removed prior to running script
# To focus on ventricular clusters
time_course_data <- subset(Feng_W_TC, subset = cell_type_spec == "ventricular_cm")
cells_to_remove <- WhichCells(ShenJ_TC_CM, idents = c(0,13,28))
time_course_data2 <- subset(ShenJ_TC_CM, cells = setdiff(Cells(ShenJ_TC_CM), cells_to_remove))
cells_to_remove <- WhichCells(P05_filtered, idents = c(11,5))
P05Heart_filtered_v <- subset(P05_filtered, cells = setdiff(Cells(P05_filtered), cells_to_remove))

# Normalize and scale both datasets
P05Heart_filtered_v <- SCTransform(P05Heart_filtered_v, verbose = FALSE)
time_course_data <- SCTransform(time_course_data, verbose = FALSE)
time_course_data2 <- SCTransform(time_course_data2, verbose = FALSE)

# Perform CCA for integration
integration_anchors <- FindIntegrationAnchors(
  object.list = list(P05Heart_filtered_v, time_course_data, time_course_data2),
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

# Visualize the results
DefaultAssay(integrated_data) <- "SCT"
DimPlot(subset_obj, reduction = "umap", split.by = "genotype", pt.size = 0.5, label = FALSE)  # Visualize different datasets
