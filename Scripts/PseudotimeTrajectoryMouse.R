# Pseudotime trajectory inference on Xiao et al. NASH mouse hepatocytes 
# 
library(Seurat)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(Matrix)
library(scater)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(viridis)

pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/'

outputPath = 'hepatocyte-damage-score/Data/Output/'

# load HS function

source('SharedFunctions.R')

# Load list HDAG

HDAG <- read.csv(
  file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')

# we need only the files from the 9 month NASH mouse
FileNames <- list.files(path = pathData)
FileNames <- FileNames[c(10,11,12)]

metaDataAnnotations <- read.table(
  paste0(pathData,
         'scitranslmed.adc9653_data_file_s1.csv'),
  sep = ';', 
  dec = ',', 
  header = TRUE)

metaDataAnnotations$cellBarcode <- 
  gsub('.*_', replacement = '', metaDataAnnotations$X)
metaDataAnnotations$condition <- 
  gsub('_.*', replacement = '', metaDataAnnotations$X)


# Read data:

counts9moNASH <- ReadMtx(mtx =  paste0(pathData, FileNames[3]), 
                       cells = paste0(pathData, FileNames[1]),
                       features = paste0(pathData, FileNames[2]),
                       feature.column = 2) 

Seurat9moNASH <- CreateSeuratObject(counts = counts9moNASH, 
                                  project = gsub('_.*', 
                                                 replacement = '', 
                                                 x = FileNames[1]),
                                  min.cells = 3)

remove(counts9moNASH)

Seurat9moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNASH , 
                                                      pattern = "^mt-")


Seurat9moNASH  <- subset(Seurat9moNASH , 
                       subset = 
                         nFeature_RNA > 500 &
                         nFeature_RNA < 8000 &
                         percent.mt < 2 )

temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH9mo',]
rownames(temp) <- temp$cellBarcode
temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNASH@meta.data)),]

Seurat9moNASH@meta.data$cellBarcode <- NA
Seurat9moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
  temp$cellBarcode
Seurat9moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
  temp$CellCluster
Seurat9moNASH@meta.data$deleteNA <- 
  rep("not na", length(Seurat9moNASH@meta.data$CellType))
Seurat9moNASH@meta.data$deleteNA[is.na(Seurat9moNASH@meta.data$CellType)] <-
  "na"

Seurat9moNASH <- subset(Seurat9moNASH,
                        subset = deleteNA == 'not na')


Seurat9moNASH <- subset(Seurat9moNASH ,
                     subset = 
                       CellType ==  "PC-Hep" |
                       CellType == "mNASH-Hep1"  |
                       CellType == "mNASH-Hep2" |
                       CellType == "PP-Hep" |
                       CellType == "Int-Hep" )

HDS <- DS_calc.func(exprMatrices = GetAssayData(Seurat9moNASH), 
                    DSignature = HDAG)

DataFrameHDS <- data.frame('HDS' = unname(HDS),
                           'CellID' = names(unlist(HDS)))

identical(DataFrameHDS$CellID, rownames(Seurat9moNASH@meta.data))
# TRUE, so:

Seurat9moNASH@meta.data$HDS <- DataFrameHDS$HDS


# Normalization & Dimension reduction & Clustering for Psuedotime 
# trajectory inference

Seurat9moNASH <- NormalizeData(Seurat9moNASH)
Seurat9moNASH <- FindVariableFeatures(Seurat9moNASH)
all.genes <- rownames(Seurat9moNASH)
Seurat9moNASH <- ScaleData(Seurat9moNASH, features = all.genes)
Seurat9moNASH <- RunPCA(Seurat9moNASH,
                        feautres =
                          VariableFeatures(Seurat9moNASH))

DimPlot(Seurat9moNASH , reduction = "pca", group.by = 'CellType')
DimPlot(Seurat9moNASH , reduction = "pca",dims = c(3,4), group.by = 'CellType')

Seurat9moNASH <- JackStraw(Seurat9moNASH, num.replicate = 100)
Seurat9moNASH <- ScoreJackStraw(Seurat9moNASH, dims = 1:15)
JackStrawPlot(Seurat9moNASH, dims = 1:15)
ggsave('XiaoMouse9moNASH_BeforePseudotime_15PCsJackStrawPlot.png',
       path = paste0(outputPath, 'Results/'))
Seurat9moNASH <- RunUMAP(Seurat9moNASH, dims = 1:14)
Seurat9moNASH  <- FindNeighbors(Seurat9moNASH , dims = 1:14)
Seurat9moNASH <- FindClusters(Seurat9moNASH, 
                              resolution = 0.4)

# checking out the cluster labels
head(Idents(Seurat9moNASH), 5)

DimPlot(Seurat9moNASH, reduction = "umap", group.by = 'CellType')

ggsave('XiaoMouse9moNASH_14dimensionsUmapColoredByCelltypeBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

DimPlot(Seurat9moNASH, reduction = "umap", group.by = 'seurat_clusters')
ggsave('XiaoMouse9moNASH_14dimensionsUmapColoredClustersBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

# users supply their own dimension reduction
# and cell clustering results to TSCAN

library(TSCAN)

scExperiment9moNASH <- as.SingleCellExperiment(Seurat9moNASH,
                                               assay = 'RNA')

colLabels(scExperiment9moNASH) <- Seurat9moNASH@meta.data$seurat_clusters

pseudo.mnn <- TSCAN::quickPseudotime( scExperiment9moNASH,
                                      use.dimred = "PCA", 
                                      dist.method = 'mnn')

plot(pseudo.mnn$mst)


minHDS <- min(Seurat9moNASH@meta.data$HDS)
maxHDS <- max(Seurat9moNASH@meta.data$HDS)


mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)

UmapPseudotime <- plotUMAP(scExperiment9moNASH, 
         colour_by=I(mnn.pseudo),
         point_size = 2,
         point_alpha = 0.7,
         point_shape = 18) + 
  scale_color_viridis(option = "G", direction = -1) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, y=UMAP_2, group=edge)) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "PT") +
  xlim(c(-8,8)) + ylim(c(-5,4))

UmapPseudotimeHDS <-plotUMAP(scExperiment9moNASH, 
                             colour_by ='HDS', 
                             text_by = "label", 
                             text_colour="black",
                             point_size = 2,
                             point_alpha = 0.7,
                             point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, 
                        y=UMAP_2, 
                        group=edge)) +    
  scale_color_viridis(option = 'B',
                      limits = c(minHDS,maxHDS),
                      direction = -1) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "HDS") +
  xlim(c(-8,8)) + ylim(c(-5,4))

UMAPcelltypes <- plotUMAP(scExperiment9moNASH, 
                          colour_by='CellType',
                          point_size = 2,
                          point_alpha = 0.7,
                          point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, y=UMAP_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'top',
         # legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16,
                                   angle = 90),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "") +
  xlim(c(-8,8)) + ylim(c(-5,4))

pannelsHorizontal <- UmapPseudotime | UmapPseudotimeHDS | UMAPcelltypes

ggsave('XiaoNASH9moMousePseudotimeHDSPannelHorizontalVersion.png',
       path = paste0(outputPath, 'Results/'))



pannelsTest <- (UmapPseudotime | (plotUMAP(scExperiment9moNASH, 
                                           colour_by='CellType',
                                           point_size = 3,
                                           point_alpha = 1,
                                           point_shape = 1) +
                                    geom_line(data = pseudo.mnn$connected$UMAP, 
                                              mapping=aes(x=UMAP_1, y=UMAP_2, 
                                                          group=edge)) +
                                    scale_color_colorblind() +
                                    theme(legend.position = 'right',
                                          legend.key.size = unit(0.1, 'cm'),
                                          legend.title = element_text(size = 20),
                                          legend.text = element_text(size = 20),
                                          axis.title.x = element_text(size = 0),
                                          axis.title.y = element_text(size = 0)) +
                                    labs(color = "") +
                                    xlim(c(-8,8)) + ylim(c(-5,4)))) / plotUMAP(scExperiment9moNASH, 
                                                                                 colour_by ='HDS', 
                                                                                 text_by = "label", 
                                                                                 text_colour="black",
                                                                                 point_size = 3,
                                                                                 point_alpha = 1,
                                                                                 point_shape = 1) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, 
                        y=UMAP_2, 
                        group=edge)) +    
  scale_color_viridis(option = 'H',
                      limits = c(minHDS,maxHDS)) +
  theme(legend.position = 'right',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "HDS") +
  xlim(c(-8,8)) + ylim(c(-5,4))

ggsave('XiaoNASH9moMousePseudotimeHDSPannelVersionI.png',
       path = paste0(outputPath, 'Results/'))

## repeat process for mouse NC 3 months

# Read data:

FileNames <- list.files(path = pathData)
FileNames <- FileNames[c(1,2,3)]

counts3moNC <- ReadMtx(mtx =  paste0(pathData, FileNames[3]), 
                         cells = paste0(pathData, FileNames[1]),
                         features = paste0(pathData, FileNames[2]),
                         feature.column = 2) 

Seurat3moNC <- CreateSeuratObject(counts = counts3moNC, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = FileNames[1]),
                                    min.cells = 3)

remove(counts3moNC)
gc()

Seurat3moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNC , 
                                                      pattern = "^mt-")


Seurat3moNC  <- subset(Seurat3moNC , 
                         subset = 
                           nFeature_RNA > 500 &
                           nFeature_RNA < 8000 &
                           percent.mt < 2 )

temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC3mo',]
rownames(temp) <- temp$cellBarcode
temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNC@meta.data)),]

Seurat3moNC@meta.data$cellBarcode <- NA
Seurat3moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- 
  temp$cellBarcode
Seurat3moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- 
  temp$CellCluster
Seurat3moNC@meta.data$deleteNA <- 
  rep("not na", length(Seurat3moNC@meta.data$CellType))
Seurat3moNC@meta.data$deleteNA[is.na(Seurat3moNC@meta.data$CellType)] <-
  "na"

Seurat3moNC <- subset(Seurat3moNC,
                        subset = deleteNA == 'not na')


Seurat3moNC <- subset(Seurat3moNC ,
                        subset = 
                          CellType ==  "PC-Hep" |
                          CellType == "mNASH-Hep1"  |
                          CellType == "mNASH-Hep2" |
                          CellType == "PP-Hep" |
                          CellType == "Int-Hep" )

HDS <- DS_calc.func(exprMatrices = GetAssayData(Seurat3moNC), 
                    DSignature = HDAG)

DataFrameHDS <- data.frame('HDS' = unname(HDS),
                           'CellID' = names(unlist(HDS)))

identical(DataFrameHDS$CellID, rownames(Seurat3moNC@meta.data))
# TRUE, so:

Seurat3moNC@meta.data$HDS <- DataFrameHDS$HDS

# Normalization & Dimension reduction & Clustering for Psuedotime 
# trajectory inference

Seurat3moNC <- NormalizeData(Seurat3moNC)
Seurat3moNC <- FindVariableFeatures(Seurat3moNC)
all.genes <- rownames(Seurat3moNC)
Seurat3moNC <- ScaleData(Seurat3moNC, features = all.genes)
Seurat3moNC <- RunPCA(Seurat3moNC,
                        feautres =
                          VariableFeatures(Seurat3moNC))

DimPlot(Seurat3moNC , reduction = "pca", group.by = 'CellType')
ggsave('XiaoMouse3moNC_BeforePseudotime_PC1and2_HepatocyteTypes.png',
       path = paste0(outputPath, 'Results/'))
DimPlot(Seurat3moNC , reduction = "pca",dims = c(3,4), group.by = 'CellType')

Seurat3moNC <- JackStraw(Seurat3moNC, num.replicate = 100)
Seurat3moNC <- ScoreJackStraw(Seurat3moNC, dims = 1:15)
JackStrawPlot(Seurat3moNC, dims = 1:15)
ggsave('XiaoMouse3moNC_BeforePseudotime_15PCsJackStrawPlot.png',
       path = paste0(outputPath, 'Results/'))
Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:14)
Seurat3moNC  <- FindNeighbors(Seurat3moNC , dims = 1:14)
Seurat3moNC <- FindClusters(Seurat3moNC, 
                              resolution = 0.5)

# checking out the cluster labels
head(Idents(Seurat3moNC), 5)

DimPlot(Seurat3moNC, reduction = "umap", group.by = 'CellType')

ggsave('XiaoMouse3moNC_14dimensionsUmapColoredByCelltypeBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

DimPlot(Seurat3moNC, reduction = "umap", group.by = 'seurat_clusters')
ggsave('XiaoMouse3moNC_14dimensionsUmapColoredClustersBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

# Pseudotime trajectory inference

scExperiment3moNC <- as.SingleCellExperiment(Seurat3moNC,
                                               assay = 'RNA')

colLabels(scExperiment3moNC) <- Seurat3moNC@meta.data$seurat_clusters

pseudo.mnn <- TSCAN::quickPseudotime( scExperiment3moNC,
                                      use.dimred = "PCA", 
                                      dist.method = 'mnn')

plot(pseudo.mnn$mst)


minHDS <- min(Seurat3moNC@meta.data$HDS)
maxHDS <- max(Seurat3moNC@meta.data$HDS)


mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)

UmapPseudotime <- plotUMAP(scExperiment3moNC, 
                           colour_by=I(mnn.pseudo),
                           point_size = 2,
                           point_alpha = 0.7,
                           point_shape = 18) + 
  scale_color_viridis(option = "B", direction = -1) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, y=UMAP_2, group=edge)) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "PT") +
  xlim(c(-8,5.5)) + ylim(c(-5.2,5.2))

UmapPseudotimeHDS <-plotUMAP(scExperiment3moNC, 
                             colour_by ='HDS', 
                             text_by = "label", 
                             text_colour="black",
                             point_size = 2,
                             point_alpha = 0.7,
                             point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, 
                        y=UMAP_2, 
                        group=edge)) +    
  scale_color_viridis(option = 'H',
                      limits = c(minHDS,maxHDS)) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
   labs(color = "HDS") +
  xlim(c(-8,5.5)) + ylim(c(-5.2,5.2))

UMAPcelltypes <- plotUMAP(scExperiment3moNC, 
                          colour_by='CellType',
                          point_size = 2,
                          point_alpha = 0.7,
                          point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=UMAP_1, y=UMAP_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'top',
        # legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16,
                                   angle = 90),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "") +
  xlim(c(-8,5.5)) + ylim(c(-5.2,5.2))

pannelsHorizontal <- UmapPseudotime | UmapPseudotimeHDS | UMAPcelltypes

ggsave('XiaoNC3moMousePseudotimeHDSPannelHorizontalVersion.png',
       path = paste0(outputPath, 'Results/'))


# Pseudotime inferred on merged da