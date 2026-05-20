library(Seurat)
library(dplyr)
library(tidyr)
library(gprofiler2)
library("scProportionTest")
library(ggplot2)


DefaultAssay(integrated_data) <- "SCT"
prop_test <- sc_utils(integrated_data)

prop_test <- permutation_test(
  prop_test, cluster_identity = "seurat_clusters",
  sample_1 = "LB1Fplus", sample_2 = "LB1Fd",
  sample_identity = "orig.ident"
)

permutation_plot(prop_test)


## Proportion stacked bar charts ##
LB1het_null <- subset(integrated_data, subset = orig.ident == c("LB1Fplus","LB1Fd"))
E13_P1 <- subset(integrated_data, subset = stage == c("E13.5","P1"))

# Extract the cluster and genotype information
cluster_counts_LB1het_null <- table(LB1het_null$seurat_clusters, LB1het_null$orig.ident)
cluster_counts_E13_P1 <- table(E13_P1$seurat_clusters, E13_P1$stage)

cluster_counts_LB1het_null_df <- as.data.frame(cluster_counts_LB1het_null)
colnames(cluster_counts_LB1het_null_df) <- c("Cluster", "Genotype", "CellCount")

cluster_counts_E13_P1_df <- as.data.frame(cluster_counts_E13_P1)
colnames(cluster_counts_E13_P1_df) <- c("Cluster", "Stage", "CellCount")

#Calculate the proportions within each genotype
cluster_counts_LB1het_null_df <- cluster_counts_LB1het_null_df %>%
  group_by(Genotype) %>%
  mutate(Proportion = CellCount / sum(CellCount)) %>%
  ungroup()

cluster_counts_E13_P1_df <- cluster_counts_E13_P1_df %>%
  group_by(Stage) %>%
  mutate(Proportion = CellCount / sum(CellCount)) %>%
  ungroup()

cluster_counts_LB1het_null_df$Genotype <- factor(cluster_counts_LB1het_null_df$Genotype, levels = c("LB1Fplus", "LB1Fd"))
cluster_counts_E13_P1_df$Stage <- factor(cluster_counts_E13_P1_df$Stage, levels = c("E13.5", "P1"))

ggplot(cluster_counts_LB1het_null_df, aes(x = Genotype, y = Proportion, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +  # Thin black border
  geom_text(aes(label = as.factor(Cluster)), 
            position = position_stack(vjust = 0.5),  # Centers the text in each section
            color = "white", size = 3) +  # White text color
  labs(x = "Condition", y = "Value", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

ggplot(cluster_counts_E13_P1_df, aes(x = Stage, y = Proportion, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +  # Thin black border
  geom_text(aes(label = as.factor(Cluster)), 
            position = position_stack(vjust = 0.5),  # Centers the text in each section
            color = "white", size = 3) +  # White text color
  labs(x = "Condition", y = "Value", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

## cell cycle dim plot ##
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


integrated_data_cc <- CellCycleScoring(integrated_data, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)
DimPlot(integrated_data_cc)
