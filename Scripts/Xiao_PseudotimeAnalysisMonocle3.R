library(R.utils)
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(Matrix)

outputPath <- 'hepatocyte-damage-score/Data/Output/'


mergedObj <- readRDS(
  file = paste0(
  outputPath, 
  '/Results/CellFateAnalysis/SeuratObjectCellFateAnalysisMerged.rds'))

DimPlot(mergedObj, group.by = 'SenescenceDamageStatus', pt.size = 0.8, alpha = 0.7)

# not sure if this is right, or I have to explcitly provide the 'RNA' assay or the 'SCT' assay 
seuratCDS <- as.cell_data_set(mergedObj,
                              assay = 'SCT')

seuratCDS <- preprocess_cds(seuratCDS, norm_method = 'none')




####

# Reduce dimension (UMAP by default)
seuratCDS <- reduce_dimension(seuratCDS ,
                              reduction_method = "UMAP")
plot_cells(seuratCDS,color_cells_by = 'CellType', label_cell_groups = TRUE)
plot_cells(seuratCDS,color_cells_by = 'cluster', label_cell_groups = TRUE)
plot_cells(seuratCDS,color_cells_by = 'SenescenceDamageStatus', label_cell_groups = TRUE)

# Cluster the cells
seuratCDS <- cluster_cells(seuratCDS)

# Learn the trajectory graph
seuratCDS <- learn_graph(seuratCDS)
plot_cells(seuratCDS, color_cells_by = "cluster",label_groups_by_cluster=FALSE, label_leaves=FALSE)


seuratCDS <- order_cells(seuratCDS)
plot_cells(seuratCDS,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


rootcell <- rownames(mergedObj@meta.data[mergedObj@meta.data$SenescenceDamageStatus == "lowHDSroot",])
endpoint1 <- rownames(mergedObj@meta.data[mergedObj@meta.data$SenescenceDamageStatus == "HighHDS_LowSI",])
endpoint2 <- rownames(mergedObj@meta.data[mergedObj@meta.data$SenescenceDamageStatus == "HighHDS_HighSI",])


seuratCDS <- order_cells(seuratCDS, 
                          root_cells = rootcell)

partitions()
