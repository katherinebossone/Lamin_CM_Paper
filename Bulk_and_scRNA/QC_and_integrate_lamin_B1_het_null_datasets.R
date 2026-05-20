library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

#Read in data

dirname2 <- "~/scRNA/KB2/"
counts_matrix_filename2 = paste0(dirname2,"/filtered_feature_bc_matrix/")
counts2 <- Read10X(data.dir = counts_matrix_filename2) 

dirname3 <- "~/scRNA/KB3/"
counts_matrix_filename3 = paste0(dirname3,"/filtered_feature_bc_matrix/")
counts3 <- Read10X(data.dir = counts_matrix_filename3) 

dirname4 <- "~/scRNA/KB4/"
counts_matrix_filename4 = paste0(dirname4,"/filtered_feature_bc_matrix/")
counts4 <- Read10X(data.dir = counts_matrix_filename4) 

# check
#   counts[1:10, 1:3] 
#   dim(counts)
#   counts_per_cell <- Matrix::colSums(counts)
#   counts_per_gene <- Matrix::rowSums(counts)
#   genes_per_cell <- Matrix::colSums(counts>0)
#   hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
#   hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
#   plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
#   title('counts vs genes per cell')
#   plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

#Merge matrix
P05Heart_Null_rep1<-CreateSeuratObject(counts = counts2, min.cells = 3, min.features = 2000, project = "LB1Fd")
P05Heart_Null_rep2<-CreateSeuratObject(counts = counts3, min.cells = 3, min.features = 2000, project = "LB1Fd")
P05Heart_Het_rep1<-CreateSeuratObject(counts = counts4, min.cells = 3, min.features = 2000, project = "LB1Fplus")

P05Heart_Het_rep1$genotype <- "LB1Fplus"
P05Heart_Het_rep1$replicate <- "pRep1"

P05Heart_Null_rep1$genotype <- "LB1Fd"
P05Heart_Null_rep1$replicate <- "dRep1"

P05Heart_Null_rep2$genotype <- "LB1Fd"
P05Heart_Null_rep2$replicate <- "dRep2"


#QC to remove based on MT and read counts
P05Heart_Het_rep1[["percent.mt"]] <- PercentageFeatureSet(P05Heart_Het_rep1, pattern = "^mt")
P05Heart_Null_rep1[["percent.mt"]] <- PercentageFeatureSet(P05Heart_Null_rep1, pattern = "^mt")
P05Heart_Null_rep2[["percent.mt"]] <- PercentageFeatureSet(P05Heart_Null_rep2, pattern = "^mt")
# VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
P05Heart_Het_rep1 <- subset(P05Heart_Het_rep1, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 30)
P05Heart_Null_rep1 <- subset(P05Heart_Null_rep1, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 30)
P05Heart_Null_rep2 <- subset(P05Heart_Null_rep2, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 30)

#Combine and normalize
#Combine and normalize
P05Heart.combined <- merge(P05Heart_Het_rep1, y = c(P05Heart_Null_rep1,P05Heart_Null_rep2), add.cell.ids = c("LB1Fplus_rep1", "LB1Fd_rep1", "LB1Fd_rep2"),
                           project = "Merged_LB1Fplus_LB1Fd")

P05Heart.combined <- NormalizeData(P05Heart.combined)
P05Heart.combined <- FindVariableFeatures(P05Heart.combined)
P05Heart.combined <- ScaleData(P05Heart.combined)
P05Heart.combined <- RunPCA(P05Heart.combined)

options(future.globals.maxSize = 30 * + 1024^3)
P05Heart.combined <- SCTransform(P05Heart.combined)
P05Heart.combined <- RunPCA(P05Heart.combined, verbose = F)
P05Heart.combined <- IntegrateLayers(
  object = P05Heart.combined,
  method = CCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
P05Heart.combined <- FindNeighbors(P05Heart.combined, dims = 1:17, reduction = "integrated.dr")
P05Heart.combined <- FindClusters(P05Heart.combined, resolution = 1)
P05Heart.combined <- RunUMAP(P05Heart.combined, dims = 1:17, reduction = "integrated.dr")

FeaturePlot(P05Heart_filtered, features = c("Sln","Lmnb1","Ttn"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey","purple"), reduction = "umap",label = TRUE, pt.size = 0.5) + NoLegend()

combined <- PrepSCTFindMarkers(combined)
combined$celltype.cluster <- paste(combined$seurat_clusters, combined$genotype, sep = "_")
Idents(combined) <- "celltype.cluster"

# To find markers for all clusters run:
# all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.26)

#Remove non-CM clusters and re-cluster
cells_to_remove <- WhichCells(P05Heart.combined, idents = c(0, 2, 6, 7, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30))
P05Heart_filtered <- subset(P05Heart.combined, cells = setdiff(Cells(P05Heart.combined), cells_to_remove))

#ElbowPlot(P05Heart_filtered)

P05Heart_filtered <- FindNeighbors(P05Heart_filtered, reduction = "integrated.dr", dims = 1:17)
P05Heart_filtered <- FindClusters(P05Heart_filtered, resolution = 0.5)
P05Heart_filtered <- RunUMAP(P05Heart_filtered, dims = 1:17, reduction = "integrated.dr")

P05Heart_filtered$orig.ident <- factor(P05Heart_filtered$orig.ident,
                                       levels = c("LB1Fplus", "LB1Fd"))
DimPlot(P05Heart_filtered, reduction = "umap", group.by = c("seurat_clusters"),label = TRUE) + NoLegend()
DimPlot(P05Heart_filtered, reduction = "umap", split.by = "orig.ident",label = TRUE, pt.size = 0.5) + NoLegend()
