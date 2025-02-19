# # HDS of Mouse hepatocytes from GSE189600 (Xiao et. al, 2023)
# - input is raw data, human version of HDAG list, our functions
# - subset, then normalized counts for each sample and the merged all 
# normalized counts
# - then added annotation from authors to metadata and extracted counts from 
# cells annotated as hepatocytes (four types) 
# - then I calculated the weighted HDS of those hepatocytes and made vioin
# plots of the HDS distribution per condition and hepatocyte type

# cell filtering as described in supplementary data from Xiaoo et al. 
# Before merging, each dataset was filtered with the parameters 
# nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 2 for mice 



library(Seurat)
library(rtracklayer)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(viridis)
library(glmGamPoi)
library(TSCAN)
library(scater)
library(scales)

# hepatocyte damage associated genes
HDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')

# load functions to calculate HDS
source('SharedFunctions.R')

pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/'

outputPath = 'hepatocyte-damage-score/Data/Output/'

FileNames <- list.files(path = pathData)

metaDataAnnotations <- read.table('hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/scitranslmed.adc9653_data_file_s1.csv', 
  sep = ';', dec = ',', header = TRUE)

metaDataAnnotations$cellBarcode <- 
  gsub('.*_', replacement = '', metaDataAnnotations$X)
metaDataAnnotations$condition <- 
  gsub('_.*', replacement = '', metaDataAnnotations$X)


# I. Read raw data & process & calculate all values: HDS, Pseudotime, etc
{
# read raw counts 
counts3moNC <- ReadMtx(mtx =  paste0(pathData, FileNames[3]), 
                    cells = paste0(pathData, FileNames[1]),
                    features = paste0(pathData, FileNames[2]),
                    feature.column = 2) 

Seurat3moNC <- CreateSeuratObject(counts = counts3moNC, 
                               project = gsub('_.*', 
                                              replacement = '', 
                                              x = FileNames[1]),
                               min.cells = 3)

class(Seurat3moNC)

remove(counts3moNC)

Seurat3moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNC, 
                                                      pattern = "^mt-")


Seurat3moNC  <- subset(Seurat3moNC , 
                       subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 2)

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

Seurat3moNC <- SCTransform(Seurat3moNC)


## Calculate HDS
HDSSeurat3moNC <- DS_calc.func(exprMatrices = 
                                 GetAssayData(Seurat3moNC,
                                              assay = 'SCT',
                                              layer = 'counts'),
                               DSignature = HDAG)

HDSSeurat3moNCRNA <- DS_calc.func(exprMatrices = 
                                    GetAssayData(Seurat3moNC,
                                                 assay = 'RNA',
                                                 layer = 'counts'),
                                  DSignature = HDAG)

Seurat3moNC@meta.data$condition <- '3m NC'
Seurat3moNC@meta.data$HDS <- unname(HDSSeurat3moNC)
Seurat3moNC@meta.data$HDSrawCounts <- unname(HDSSeurat3moNCRNA)


# dimension reduction & PT trajectory inference
Seurat3moNC <- RunPCA(Seurat3moNC)
Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:20)
Seurat3moNC <- FindNeighbors(Seurat3moNC, dims = 1:20)
Seurat3moNC <- FindClusters(Seurat3moNC, 
                            res = 0.3)
DimPlot(Seurat3moNC, label = TRUE)

scExp3mo <- as.SingleCellExperiment(Seurat3moNC,
                                    assay = 'SCT')
colLabels(scExp3mo) <- Seurat3moNC@meta.data$SCT_snn_res.0.3

pseudo.mnn3mo  <- TSCAN::quickPseudotime( scExp3mo,
                                          use.dimred = "PCA", 
                                          dist.method = 'mnn')
plot(pseudo.mnn3mo$mst)
mnn.pseudo3mo <- 
  averagePseudotime(pseudo.mnn3mo$ordering)

remove(temp)
gc()


# read raw counts 9moNC
counts9moNC <- ReadMtx(mtx =  paste0(pathData, FileNames[6]), 
                       cells = paste0(pathData, FileNames[4]),
                       features = paste0(pathData, FileNames[5]),
                       feature.column = 2) 

Seurat9moNC <- CreateSeuratObject(counts = counts9moNC, 
                                  project = gsub('_.*', 
                                                 replacement = '', 
                                                 x = FileNames[4]),
                                  min.cells = 3)

Seurat9moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNC, 
                                                    pattern = "^mt-")

Seurat9moNC  <- subset(Seurat9moNC , 
                       subset = 
                         nFeature_RNA > 500 &
                         nFeature_RNA < 8000 &
                         percent.mt < 2)

remove(counts9moNC)

temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC9mo',]
rownames(temp) <- temp$cellBarcode
temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNC@meta.data)),]

Seurat9moNC@meta.data$cellBarcode <- NA
Seurat9moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- 
  temp$cellBarcode
Seurat9moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- 
  temp$CellCluster
Seurat9moNC@meta.data$deleteNA <- 
  rep("not na", length(Seurat9moNC@meta.data$CellType))
Seurat9moNC@meta.data$deleteNA[is.na(Seurat9moNC@meta.data$CellType)] <-
  "na"

Seurat9moNC <- subset(Seurat9moNC,
                      subset = deleteNA == 'not na')


Seurat9moNC <- subset(Seurat9moNC ,
                      subset = 
                        CellType ==  "PC-Hep" |
                        CellType == "mNASH-Hep1"  |
                        CellType == "mNASH-Hep2" |
                        CellType == "PP-Hep" |
                        CellType == "Int-Hep" )

Seurat9moNC <- SCTransform(Seurat9moNC)

# calculate HDS:

HDSSeurat9moNC <- DS_calc.func(exprMatrices = 
                                 GetAssayData(Seurat9moNC,
                                              assay = 'SCT',
                                              layer = 'counts'),
                               DSignature = HDAG)

HDSSeurat9moNCRNA <- DS_calc.func(exprMatrices = 
                                    GetAssayData(Seurat9moNC,
                                                 assay = 'RNA',
                                                 layer = 'counts'),
                                  DSignature = HDAG)

Seurat9moNC@meta.data$condition <- '9m NC'
Seurat9moNC@meta.data$HDS <- unname(HDSSeurat9moNC)
Seurat9moNC@meta.data$HDSrawCounts <- unname(HDSSeurat9moNCRNA)



# dimension reduction & PT trajectory inference
Seurat9moNC <- RunPCA(Seurat9moNC)
Seurat9moNC <- RunUMAP(Seurat9moNC, dims = 1:20)
Seurat9moNC <- FindNeighbors(Seurat9moNC, dims = 1:20)
Seurat9moNC <- FindClusters(Seurat9moNC,
                              resolution = 0.2)
DimPlot(Seurat9moNC, label = TRUE)

scExp9mo <- as.SingleCellExperiment(Seurat9moNC,
                                    assay = 'SCT')
colLabels(scExp9mo) <- Seurat9moNC@meta.data$SCT_snn_res.0.2

pseudo.mnn9mo  <- TSCAN::quickPseudotime( scExp9mo,
                                          use.dimred = "PCA", 
                                          dist.method = 'mnn')
plot(pseudo.mnn9mo$mst)
mnn.pseudo9mo <- 
  averagePseudotime(pseudo.mnn9mo$ordering)

remove(temp)

# 3 month NASH

counts3moNASH <- ReadMtx(mtx =  paste0(pathData, FileNames[9]), 
                       cells = paste0(pathData, FileNames[7]),
                       features = paste0(pathData, FileNames[8]),
                       feature.column = 2) 

Seurat3moNASH <- CreateSeuratObject(counts = counts3moNASH, 
                                  project = gsub('_.*', 
                                                 replacement = '', 
                                                 x = FileNames[9]),
                                  min.cells = 3)

Seurat3moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNASH, 
                                                    pattern = "^mt-")

Seurat3moNASH  <- subset(Seurat3moNASH , 
                       subset = 
                         nFeature_RNA > 500 &
                         nFeature_RNA < 8000 &
                         percent.mt < 2)

remove(counts3moNASH)

temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH3mo',]
rownames(temp) <- temp$cellBarcode
temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNASH@meta.data)),]

Seurat3moNASH@meta.data$cellBarcode <- NA
Seurat3moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- 
  temp$cellBarcode
Seurat3moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- 
  temp$CellCluster
Seurat3moNASH@meta.data$deleteNA <- 
  rep("not na", length(Seurat3moNASH@meta.data$CellType))
Seurat3moNASH@meta.data$deleteNA[is.na(Seurat3moNASH@meta.data$CellType)] <-
  "na"

Seurat3moNASH <- subset(Seurat3moNASH,
                      subset = deleteNA == 'not na')


Seurat3moNASH <- subset(Seurat3moNASH ,
                      subset = 
                        CellType ==  "PC-Hep" |
                        CellType == "mNASH-Hep1"  |
                        CellType == "mNASH-Hep2" |
                        CellType == "PP-Hep" |
                        CellType == "Int-Hep" )
Seurat3moNASH <- SCTransform(Seurat3moNASH)

## Calculate HDS
HDSSeurat3moNASH <- DS_calc.func(exprMatrices = 
                                   GetAssayData(Seurat3moNASH,
                                                assay = 'SCT',
                                                layer = 'counts'),
                                 DSignature = HDAG)

HDSSeurat3moNASHRNA <- DS_calc.func(exprMatrices = 
                                      GetAssayData(Seurat3moNASH,
                                                   assay = 'RNA',
                                                   layer = 'counts'),
                                    DSignature = HDAG)

Seurat3moNASH@meta.data$condition <- '3m NASH'
Seurat3moNASH@meta.data$HDS <- unname(HDSSeurat3moNASH)
Seurat3moNASH@meta.data$HDSrawCounts <- unname(HDSSeurat3moNASHRNA)

remove(temp)
gc()

#dimension reduction & PT trajectory inference
Seurat3moNASH <- RunPCA(Seurat3moNASH)
Seurat3moNASH <- RunUMAP(Seurat3moNASH, dims = 1:14)
Seurat3moNASH <- FindNeighbors(Seurat3moNASH, dims = 1:14)
# adjusted resolution to make it more comparable to annotation and other metadata
Seurat3moNASH <- FindClusters(Seurat3moNASH, resolution = 0.1)
DimPlot(Seurat3moNASH, label = TRUE)

scExp3moNASH <- as.SingleCellExperiment(Seurat3moNASH,
                                    assay = 'SCT')
colLabels(scExp3moNASH) <- Seurat3moNASH@meta.data$SCT_snn_res.0.1

pseudo.mnn3moNASH  <- TSCAN::quickPseudotime( scExp3moNASH,
                                          use.dimred = "PCA", 
                                          dist.method = 'mnn')
plot(pseudo.mnn3moNASH$mst)
mnn.pseudo3moNASH <- 
  averagePseudotime(pseudo.mnn3moNASH$ordering)


# 9 month NASH

counts9moNASH <- ReadMtx(mtx =  paste0(pathData, FileNames[12]), 
                         cells = paste0(pathData, FileNames[10]),
                         features = paste0(pathData, FileNames[11]),
                         feature.column = 2) 

Seurat9moNASH <- CreateSeuratObject(counts = counts9moNASH, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = FileNames[12]),
                                    min.cells = 3)

Seurat9moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNASH, 
                                                    pattern = "^mt-")

Seurat9moNASH  <- subset(Seurat9moNASH , 
                         subset = 
                           nFeature_RNA > 500 &
                           nFeature_RNA < 8000 &
                           percent.mt < 2)

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

remove(counts9moNASH, temp)

Seurat9moNASH <- SCTransform(Seurat9moNASH)

## calculate HDS
HDSSeurat9moNASH <- DS_calc.func(exprMatrices = 
                                   GetAssayData(Seurat9moNASH,
                                                assay = 'SCT',
                                                layer = 'counts'),
                                 DSignature = HDAG)

HDSSeurat9moNASHRNA <- DS_calc.func(exprMatrices = 
                                      GetAssayData(Seurat9moNASH,
                                                   assay = 'RNA',
                                                   layer = 'counts'),
                                    DSignature = HDAG)


Seurat9moNASH@meta.data$condition <- '9m NASH'
Seurat9moNASH@meta.data$HDS <- unname(HDSSeurat9moNASH)
Seurat9moNASH@meta.data$HDSrawCounts <- unname(HDSSeurat9moNASHRNA)

# dimension reduction & PT trajectory inference
Seurat9moNASH <- RunPCA(Seurat9moNASH)
Seurat9moNASH <- RunUMAP(Seurat9moNASH, dims = 1:14)
Seurat9moNASH <- FindNeighbors(Seurat9moNASH, dims = 1:14)
Seurat9moNASH <- FindClusters(Seurat9moNASH, resolution = 0.3)
DimPlot(Seurat9moNASH, label = TRUE)

scExp9moNASH <- as.SingleCellExperiment(Seurat9moNASH,
                                        assay = 'SCT')
colLabels(scExp9moNASH) <- Seurat9moNASH@meta.data$SCT_snn_res.0.3

pseudo.mnn9moNASH  <- TSCAN::quickPseudotime( scExp9moNASH,
                                              use.dimred = "PCA", 
                                              dist.method = 'mnn')
plot(pseudo.mnn9moNASH$mst)
mnn.pseudo9moNASH <- 
  averagePseudotime(pseudo.mnn9moNASH$ordering)

PTresultsList <- list(PTvalues = list('3moNC' = mnn.pseudo3mo, 
                                      '9moNC'= mnn.pseudo9mo, 
                                      '3moNASH'= mnn.pseudo3moNASH, 
                                      '9moNASH' = mnn.pseudo9moNASH),
                      PTordering = list('3moNC' =pseudo.mnn3mo,
                                        '9moNC'= pseudo.mnn9mo,
                                        '3moNASH'= pseudo.mnn3moNASH,
                                        '9moNASH'= pseudo.mnn9moNASH),
                      scExpObjects = list('3moNC' = scExp3mo,
                                          '9moNC'= scExp9mo,
                                          '3moNASH'= scExp3moNASH,
                                          '9moNASH'= scExp9moNASH))

saveRDS(PTresultsList, 
        file = paste0(outputPath, 'Results/Xiao2023_PseudoTimeTSCAN.rds'))


}
#################### HERE
# ## dimension reduction 
# DefaultAssay(Seurat3moNC) <- 'SCT'
# 
# # SCT scale.data are the variable features.
# VariableFeatures(Seurat3moNC[["SCT"]]) <- 
#   rownames(Seurat3moNC[["SCT"]]@scale.data)
# Seurat3moNC <- RunPCA(Seurat3moNC,
#                       feautres = VariableFeatures(Seurat3moNC))
# Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:14)
# 
# DimPlot(Seurat3moNC , 
#         reduction = "pca", 
#         group.by = 'CellType')
# 
# Seurat3moNC <- FindNeighbors(Seurat3moNC,
#                              dims = 1:14)
# Seurat3moNC <- FindClusters(Seurat3moNC,
#                             resolution = 0.4)
# DimPlot(Seurat3moNC , 
#         reduction = "umap", 
#         group.by = 'SCT_snn_res.0.4')
# 
# scExp3mo <- as.SingleCellExperiment(Seurat3moNC,
#                                     assay = 'SCT')
# 
# colLabels(scExp3mo) <- Seurat3moNC@meta.data$SCT_snn_res.0.4
# 
# pseudo.mnn3mo  <- TSCAN::quickPseudotime( scExp3mo,
#                                       use.dimred = "PCA", 
#                                       dist.method = 'mnn')
# 
# plot(pseudo.mnn3mo$mst)
# 
# mnn.pseudo3mo <- averagePseudotime(pseudo.mnn3mo$ordering)
# 
# 
# # 9moNC
# VariableFeatures(Seurat9moNC[["SCT"]]) <- 
#   rownames(Seurat9moNC[["SCT"]]@scale.data)
# Seurat9moNC <- RunPCA(Seurat9moNC,
#                       feautres = VariableFeatures(Seurat9moNC))
# 
# DimPlot(Seurat9moNC , 
#         reduction = "umap", 
#         group.by = 'CellType')
# 
# Seurat9moNC <- RunUMAP(Seurat9moNC, dims = 1:14)
# 
# Seurat9moNC <- FindNeighbors(Seurat9moNC,
#                              dims = 1:14)
# Seurat9moNC <- FindClusters(Seurat9moNC,
#                             resolution = 0.4)
# DimPlot(Seurat9moNC, 
#         reduction = "umap", 
#         group.by = 'SCT_snn_res.0.4')
# 
# 
# scExp9mo <- as.SingleCellExperiment(Seurat9moNC,
#                                     assay = 'SCT')
# 
# colLabels(scExp9mo) <- Seurat9moNC@meta.data$SCT_snn_res.0.4
# 
# pseudo.mnn9mo  <- TSCAN::quickPseudotime( scExp9mo,
#                                           use.dimred = "PCA", 
#                                           dist.method = 'mnn')
# 
# plot(pseudo.mnn9mo$mst)
# 
# mnn.pseudo9mo <- averagePseudotime(pseudo.mnn9mo$ordering)
# 
# # 3moNASH
# VariableFeatures(Seurat3moNASH[['SCT']]) <- 
#   rownames(Seurat3moNASH[['SCT']]@scale.data)
# Seurat3moNASH <- RunPCA(Seurat3moNASH,
#                         features = VariableFeatures(Seurat3moNASH))
# 
# DimPlot(Seurat3moNASH,
#         reduction = 'umap',
#         group.by = 'CellType')
# 
# Seurat3moNASH <- RunUMAP(Seurat3moNASH,
#                          dims = 1:20)
# 
# Seurat3moNASH <- FindNeighbors(Seurat3moNASH,
#                              dims = 1:20)
# Seurat3moNASH <- FindClusters(Seurat3moNASH,
#                             resolution = 0.4)
# DimPlot(Seurat3moNASH, 
#         reduction = "umap", 
#         group.by = 'SCT_snn_res.0.4')
# 
# scExp3moNASH <- as.SingleCellExperiment(Seurat3moNASH,
#                                     assay = 'SCT')
# colLabels(scExp3moNASH) <- Seurat3moNASH@meta.data$SCT_snn_res.0.4
# pseudo.mnn3moNASH  <- TSCAN::quickPseudotime(scExp3moNASH,
#                                           use.dimred = "PCA", 
#                                           dist.method = 'mnn')
# plot(pseudo.mnn3moNASH$mst)
# 
# mnn.pseudo3moNASH <- averagePseudotime(pseudo.mnn3moNASH$ordering)
# 
# # 9moNASH
# VariableFeatures(Seurat9moNASH[['SCT']]) <- 
#   rownames(Seurat9moNASH[['SCT']]@scale.data)
# Seurat9moNASH <- RunPCA(Seurat9moNASH,
#                         features = VariableFeatures(Seurat9moNASH))
# 
# DimPlot(Seurat9moNASH,
#         reduction = 'umap',
#         group.by = 'CellType')
# 
# Seurat9moNASH <- RunUMAP(Seurat9moNASH,
#                          dims = 1:14)
# Seurat9moNASH <- FindNeighbors(Seurat9moNASH,
#                                dims = 1:14)
# Seurat9moNASH <- FindClusters(Seurat9moNASH,
#                               resolution = 0.3)
# DimPlot(Seurat9moNASH, 
#         reduction = "umap", 
#         group.by = 'SCT_snn_res.0.3')
# DimPlot(Seurat9moNASH, 
#         reduction = "umap", 
#         group.by = 'CellType')
# 
# scExp9moNASH <- as.SingleCellExperiment(Seurat9moNASH,
#                                         assay = 'SCT')
# 
# colLabels(scExp9moNASH) <- Seurat9moNASH@meta.data$SCT_snn_res.0.3
# 
# pseudo.mnn9moNASH  <- TSCAN::quickPseudotime(scExp9moNASH,
#                                              use.dimred = "PCA", 
#                                              dist.method = 'mnn')
# 
# plot(pseudo.mnn9moNASH$mst)
# 
# mnn.pseudo9moNASH <- averagePseudotime(pseudo.mnn9moNASH$ordering)

# II. Store all values in dataframe for easy plotting with ggplot
{
## bring all HDS from all samples to one dataframe
# sanity check: meta data matches HDS 
identical(names(unlist(HDSSeurat3moNC)),
          rownames(Seurat3moNC@meta.data))
identical(names(unlist(HDSSeurat9moNC)),
          rownames(Seurat9moNC@meta.data))
identical(names(unlist(HDSSeurat3moNASH)),
          rownames(Seurat3moNASH@meta.data))
identical(names(unlist(HDSSeurat9moNASH)),
          rownames(Seurat9moNASH@meta.data))
# all TRUE: good to go
DataFrameHDS <- data.frame('HDS' = unname(HDSSeurat3moNC),
                           'CellID' = names(unlist(HDSSeurat3moNC)))

DataFrameHDS$condition <- rep('NC', length(DataFrameHDS$HDS))
DataFrameHDS$age <- rep('3 months', length(DataFrameHDS$HDS))
DataFrameHDS$sample <- rep('NC 3 months', length(DataFrameHDS$HDS))
identical(DataFrameHDS$CellID, rownames(Seurat3moNC@meta.data))
DataFrameHDS$HepatocyteType <- Seurat3moNC@meta.data$CellType


DataFrameHDS <- rbind.data.frame(DataFrameHDS,
                                 data.frame('HDS' = unname(HDSSeurat9moNC),
                                            'CellID' = names(unlist(HDSSeurat9moNC)),
                                            'condition' = rep('NC', 
                                                              length(HDSSeurat9moNC)),
                                            'age' = rep('9 months', 
                                                        length(HDSSeurat9moNC)),
                                            'sample' = rep('NC 9 months', 
                                                           length(HDSSeurat9moNC)),
                                            'HepatocyteType' = Seurat9moNC@meta.data$CellType),
                                 
                                 data.frame('HDS' = unname(HDSSeurat3moNASH),
                                            'CellID' = names(unlist(HDSSeurat3moNASH)),
                                            'condition' = rep('NASH', 
                                                              length(HDSSeurat3moNASH)),
                                            'age' = rep('3 months', 
                                                        length(HDSSeurat3moNASH)),
                                            'sample' = rep('NASH 3 months', 
                                                           length(HDSSeurat3moNASH)),
                                            'HepatocyteType' = Seurat3moNASH@meta.data$CellType),
                                 
                                 data.frame('HDS' = unname(HDSSeurat9moNASH),
                                            'CellID' = names(unlist(HDSSeurat9moNASH)),
                                            'condition' = rep('NASH', 
                                                              length(HDSSeurat9moNASH)),
                                            'age' = rep('9 months', 
                                                        length(HDSSeurat9moNASH)),
                                            'sample' = rep('NASH 9 months', 
                                                           length(HDSSeurat9moNASH)),
                                            'HepatocyteType' = Seurat9moNASH@meta.data$CellType)
                                 
                                 
                                 )

Before <- DataFrameHDS

DataFrameHDS$condition <- factor(DataFrameHDS$condition,
                                 ordered = TRUE,
                                 levels = c("NC","NASH" ))

DataFrameHDS$age <- factor(DataFrameHDS$age,
                           ordered = TRUE, 
                           levels = c("3 months",
                                      "9 months"))

DataFrameHDS$sample <- factor(DataFrameHDS$sample,
                           ordered = TRUE, 
                           levels = c("NC 3 months",
                                      "NC 9 months",
                                      "NASH 3 months",
                                      "NASH 9 months"))


DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                      ordered = TRUE,
                                      levels = c("PP-Hep",
                                                 "Int-Hep",
                                                 "PC-Hep",
                                                 "mNASH-Hep1",
                                                 "mNASH-Hep2"))


remove(Before)
saveRDS(DataFrameHDS, file = paste0(outputPath, 'Xiao2023HDSsnRNAseqAllHep.rds'))

}

# III. Violin plot of HDS
####
{
PlotHDSConditions <- 
  ggplot( DataFrameHDS,
          aes( x = sample,
               y = HDS,
               color = sample
          )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind()  +
  theme( text = element_text(size=25),
         plot.title = element_text(size=25),
         legend.position = 'right',
         axis.text.x = element_blank()) +
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title ="Diet & Age"))+
  labs(y = "HDS", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")


PlotHDSHepTypes <- 
  ggplot( DataFrameHDS,
          aes( x = sample,
               y = HDS,
               color = sample
          )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Health Status")) +
  labs(y = "HDS", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~HepatocyteType)

# For manuscript:

ggplot( DataFrameHDS[DataFrameHDS$HepatocyteType != 'mNASH-Hep1' &
                         DataFrameHDS$HepatocyteType != 'mNASH-Hep2',],
          aes( x = sample,
               y = HDS,
               color = sample
          )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Health Status")) +
  labs(y = "HDS", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~HepatocyteType)

# in manuscript:
 ggsave('XiaoEtAl_HDS_ViolinPlots_OnlyZones.pdf',
       path = paste0(outputPath, 'Results/'))
 ggsave('XiaoEtAl_HDS_ViolinPlots_OnlyZones.png',
        path = paste0(outputPath, 'Results/'))
}
## IV. Plot UMAPs with CellType annotation, Pseudotime and HDS
{
# 1. Plot Pseudotime Inferece Results

UmapPT3moNC <- plotUMAP(scExp3mo,
                        colour_by=I(mnn.pseudo3mo),
                        point_size = 1,
                        point_alpha = 0.7) +
  scale_color_viridis(option = "B", 
                      direction = -1
                     # limits = c(min(I(mnn.pseudo3mo)),
                      #           max(I(mnn.pseudo3mo)))
                      ) +
  geom_line(data = pseudo.mnn3mo$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, group=edge)) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "PT") 

UmapPT9moNC <- plotUMAP(scExp9mo,
                        colour_by=I(mnn.pseudo9mo),
                        point_size = 1,
                        point_alpha = 0.7) +
  scale_color_viridis(option = "B", 
                      direction = -1
                     # limits = c(min(I(mnn.pseudo9mo)),
                    #             max(I(mnn.pseudo9mo)))
                      ) +
  geom_line(data = pseudo.mnn9mo$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, group=edge)) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "PT") 

UmapPT3moNASH <- plotUMAP(scExp3moNASH,
                        colour_by=I(mnn.pseudo3moNASH),
                        point_size = 1,
                        point_alpha = 0.7) +
  scale_color_viridis(option = "B", 
                      direction = -1
                     # limits = c(min(I(mnn.pseudo3moNASH)),
                      #           max(I(mnn.pseudo3moNASH)))
                        ) +
  geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, group=edge)) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "PT")

UmapPT9moNASH <- plotUMAP(scExp9moNASH,
                          colour_by=I(mnn.pseudo9moNASH),
                          point_size = 1,
                          point_alpha = 0.7) +
  scale_color_viridis(option = "B", 
                      direction = -1
                    #  limits = c(min(I(mnn.pseudo9moNASH)),
                    #             max(I(mnn.pseudo9moNASH)))
                      ) +
  geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, group=edge)) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "PT") 

## Plot HDS: 
# color palette will use min and max values HDS from each sample as 
# reference points --> each UMAP needs it's own HDS legend. 
# this is necessary cause the color resolution is needed to compare to pseudo
# time trajectory 

# Find global color boarders: 1%-5% top and bottom 
summary(DataFrameHDS)
length(DataFrameHDS$sample)

limitBottom1percent <- DataFrameHDS[order(DataFrameHDS$HDS),'HDS'][round((15825)*0.05)]
limitTop1percent <- DataFrameHDS[order(DataFrameHDS$HDS, decreasing = TRUE),'HDS'][round((15825)*0.01)]

UmapHDS3mo <-plotUMAP(scExp3mo,
                      colour_by ='HDS',
                      point_size = 0.5,
                      point_alpha = 1) +
  geom_line(data = pseudo.mnn3mo$connected$UMAP, 
            mapping=aes(x=umap_1, 
                        y=umap_2, 
                        group=edge)) +
  # scale_color_viridis(option = "turbo", 
  #                     direction = 1) +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1,
                        limits = c(limitBottom1percent,
                                    limitTop1percent),
                        oob = squish
                        ) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "HDS") 

UmapHDS9mo <-plotUMAP(scExp9mo,
                      colour_by ='HDS',
                      point_size = 0.5,
                      point_alpha = 1) +
  geom_line(data = pseudo.mnn9mo$connected$UMAP, 
            mapping=aes(x=umap_1, 
                        y=umap_2, 
                        group=edge)) +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1,
                        limits = c(limitBottom1percent,
                                   limitTop1percent),
                        oob = squish) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "HDS") 

UmapHDS3moNASH <-plotUMAP(scExp3moNASH,
                      colour_by ='HDS',
                      point_size = 0.5,
                      point_alpha = 1
                      ) +
  geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, 
                        y=umap_2, 
                        group=edge)) +
  # scale_color_viridis(option = "turbo", direction = 1,
  #                     limits = c(limitBottom1percent,
  #                                limitTop1percent)) +
   scale_color_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(limitBottom1percent,
                                    limitTop1percent),
                         oob = squish) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "HDS") 

UmapHDS9moNASH <-plotUMAP(scExp9moNASH,
                          colour_by ='HDS',
                          point_size = 0.5,
                          point_alpha = 1) +
  geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, 
                        y=umap_2, 
                        group=edge)) +
   # scale_color_viridis(option = "turbo", direction = 1,
   #                     limits = c(limitBottom1percent,
   #                                limitTop1percent),
   #                     oob = squish) +
   scale_color_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(limitBottom1percent,
                                    limitTop1percent),
                         oob = squish) +
  theme(legend.position = 'top',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  labs(color = "HDS") 

##### Plot Hepatocyte Zone Annotation

UMAPct3moNC <- plotUMAP(scExp3mo, 
                          colour_by='CellType',
                          point_size = 0.5,
                          point_alpha = 0.7) +
  geom_line(data = pseudo.mnn3mo$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'none',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text = element_blank()
        ) 

UMAPct9moNC <- plotUMAP(scExp9mo, 
                        colour_by='CellType',
                        point_size = 0.5,
                        point_alpha = 0.7) +
  geom_line(data = pseudo.mnn9mo$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'none',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text = element_blank()
  ) 



UMAPct3moNASH <- plotUMAP(scExp3moNASH, 
                        colour_by='CellType',
                        point_size = 0.5,
                        point_alpha = 0.7) +
  geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'none',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text = element_blank()
  ) 


UMAPct9moNASH <- plotUMAP(scExp9moNASH, 
                          colour_by='CellType',
                          point_size = 0.5,
                          point_alpha = 0.7) +
  geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'none',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text = element_blank()
  ) 

patchWorkCTs <- (UMAPct3moNC | UMAPct9moNC | UMAPct3moNASH | UMAPct9moNASH )
patchWorkHDS <- (UmapHDS3mo | UmapHDS9mo | UmapHDS3moNASH | UmapHDS9moNASH)
patchWorkPT <- (UmapPT3moNC | UmapPT9moNC | UmapPT3moNASH | UmapPT9moNASH)

allPatch <- (patchWorkCTs / patchWorkHDS / patchWorkPT ) + theme(legend.position = 'none')

allPatch
ggsave('XiaoMouseEverySampleUmapsCelltypePTHDSNewest.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 21,
       width = 35,
       units = 'cm')

# For manuscript:

Umaps3moNC <- (UMAPct3moNC | UmapHDS3mo | UmapPT3moNC)
Umaps3moNC
ggsave('XiaoMouse3mNCUmaps.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')
ggsave('XiaoMouse3mNCUmaps.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')


Umaps9moNC <- (UMAPct9moNC | UmapHDS9mo | UmapPT9moNC)
Umaps9moNC
ggsave('XiaoMouse9mNCUmaps.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')
ggsave('XiaoMouse9mNCUmaps.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')

Umaps3moNASH <- (UMAPct3moNASH | UmapHDS3moNASH | UmapPT3moNASH)
Umaps3moNASH
ggsave('XiaoMouse3mNASHUmaps.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')
ggsave('XiaoMouse3mNASHUmaps.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')


Umaps9moNASH <- (UMAPct9moNASH | UmapHDS9moNASH | UmapPT9moNASH)
Umaps9moNASH
ggsave('XiaoMouse9mNASHUmaps.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')
ggsave('XiaoMouse9mNASHUmaps.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 10,
       width = 35,
       units = 'cm')






Umap_3months <- Umaps3moNC / Umaps3moNASH

Umap_3months
ggsave('XiaoMouse3monthsSampleUmapsIndivualLegends.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 35,
       units = 'cm')

ggsave('XiaoMouse3monthsSampleUmapsIndivualLegends.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 35,
       units = 'cm')


Umap_9months <- Umaps9moNC/ Umaps9moNASH

Umap_9months 
ggsave('XiaoMouse9monthsSampleUmapsIndivualLegends.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 35,
       units = 'cm')
ggsave('XiaoMouse9monthsSampleUmapsIndivualLegends.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 35,
       units = 'cm')

# to save the celltype legend
UMAPct9moNASH + theme(legend.position = 'top',
                     legend.key.size = unit(1.5, 'cm'),
                     legend.title = element_blank(),
                     legend.text = element_text(size = 16),
                     axis.ticks = element_blank(),
                     axis.text = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6)))

ggsave('CellTypeLegend.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 30,
       units = 'cm')
ggsave('CellTypeLegend.png',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 30,
       units = 'cm')

}
# Repeat the pseudotime trajectory inference on merged samples:

## this merging step does not work after SCT transform
## merge previously normalized and filtered sample counts

mergedHepatocytes <- merge(x = Seurat3moNC, y = list(Seurat9moNC,
                                                     Seurat3moNASH,
                                                     Seurat9moNASH), 
                           add.cell.ids = c("NC 3 months", "NC 9 months", "NASH 3 months",
                                            "NASH 9 months"), 
                           project = "MouseXiaoEtAl",
                           merge.data = TRUE)

saveRDS(object = mergedHepatocytes, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllMouseHepMerged.rds')




# Dimension Reduction 
# testing this without having normalized data 
DefaultAssay(mergedHepatocytes) <- 'SCT'

# merged SCT scale.data are the variable features.

VariableFeatures(mergedHepatocytes[["SCT"]]) <- rownames(mergedHepatocytes[["SCT"]]@scale.data)

mergedHepatocytes <- RunPCA(mergedHepatocytes,
                            feautres = VariableFeatures(mergedHepatocytes))

DimPlot(mergedHepatocytes, 
        reduction = "pca", 
        group.by = 'CellType')
# PCA results 
DimPlot(mergedHepatocytes, 
        reduction = "pca", 
        group.by = 'condition')

pannelPCA1_2 <- (DimPlot(mergedHepatocytes, reduction = "pca", group.by = 'CellType') | DimPlot(mergedHepatocytes, reduction = "pca", group.by = 'condition') )

ggsave('XiaoMouseAllHepMerged_PCA.png',
       path = paste0(outputPath, 'Results/'))

DimPlot(mergedHepatocytes, 
        reduction = "pca", 
        dims = c(3,4), 
        group.by = 'CellType')

DimPlot(mergedHepatocytes, 
        reduction = "pca", 
        dims = c(4,5), 
        group.by = 'CellType')

# UMAP and clustering
mergedHepatocytes <- RunUMAP(mergedHepatocytes, 
                             dims = 1:20)

DimPlot(mergedHepatocytes, 
        reduction = "umap", 
        group.by = 'CellType')


ggsave('XiaoMouseAllHepMerged_20dimUmapColoredByCelltypeBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

mergedHepatocytes  <- FindNeighbors(mergedHepatocytes, 
                                    dims = 1:20)
mergedHepatocytes<- FindClusters(mergedHepatocytes, 
                                 resolution = 0.2)
## 0.4 generated 9 clusters 
# checking out the cluster labels
head(Idents(mergedHepatocytes), 5)

DimPlot(mergedHepatocytes, 
        reduction = "umap", 
        group.by = 'seurat_clusters')

ggsave('XiaoMouseAllHepMerged_20dimFindClusterUmapColoredClustersBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

DimPlot(mergedHepatocytes, 
        reduction = "umap", 
        group.by = 'condition')

ggsave('XiaoMouseAllHepMerged_20dimFindClusterUmapColoredSamplesBeforePseudotime.png',
       path = paste0(outputPath, 'Results/'))

# Trajectory Inference

scExperiment <- as.SingleCellExperiment(mergedHepatocytes,
                                        assay = 'SCT')

colLabels(scExperiment) <- mergedHepatocytes@meta.data$seurat_clusters

pseudo.mnn <- TSCAN::quickPseudotime( scExperiment,
                                      use.dimred = "UMAP", 
                                      dist.method = 'mnn')

plot(pseudo.mnn$mst)


minHDS <- min(mergedHepatocytes@meta.data$HDS)
maxHDS <- max(mergedHepatocytes@meta.data$HDS)


mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)

## plot results

UmapPseudotime <- plotUMAP(scExperiment, 
                           colour_by=I(mnn.pseudo),
                           point_size = 1,
                           point_alpha = 0.7,
                           point_shape = 18) + 
  scale_color_viridis(option = "B", direction = -1) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, group=edge)) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "PT") + guides(x = "none", y = "none")

UmapPseudotimeHDS <-plotUMAP(scExperiment, 
                             colour_by ='HDS', 
                             text_by = "label", 
                             text_colour="black",
                             point_size = 1,
                             point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=umap_1, 
                        y=umap_2, 
                        group=edge)) +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1,
                        limits = c(minHDS, maxHDS)) +
  theme(legend.position = 'top',
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "HDS") + guides(x = "none", y = "none")

UMAPcelltypes <- plotUMAP(scExperiment, 
                          colour_by='CellType',
                          point_size = 1,
                          point_alpha = 0.7,
                          point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'top',
        # legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16,
                                   angle = 90),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "") + guides(x = "none", y = "none")

UMAPsamples <- plotUMAP(scExperiment, 
                        colour_by='condition',
                        point_size = 1,
                        point_alpha = 0.7,
                        point_shape = 18) +
  geom_line(data = pseudo.mnn$connected$UMAP, 
            mapping=aes(x=umap_1, y=umap_2, 
                        group=edge)) +
  scale_color_colorblind() +
  theme(legend.position = 'top',
        # legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16,
                                   angle = 90),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  labs(color = "") + guides(x = "none", y = "none")

pannelsHorizontalCelltype <- UmapPseudotime | UmapPseudotimeHDS | UMAPcelltypes
pannelsHorizontalSamples <- UmapPseudotime | UmapPseudotimeHDS | UMAPsamples
allUmaps <- UmapPseudotime | UmapPseudotimeHDS | UMAPsamples | UMAPcelltypes 

pannelsHorizontalCelltype
ggsave('XiaoMousePseudotimeHDSCelltypePannelHorizontalVersion.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 30,
       units = 'cm')

pannelsHorizontalSamples
ggsave('XiaoMousePseudotimeHDSSamplePannelHorizontalVersion.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 30,
       units = 'cm')

allUmaps
ggsave('XiaoMousePseudotimeAllPannelHorizontalVersion.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 33,
       units = 'cm')


### Psuedotime violin plots

DataFramePT <- data.frame('PT' = unname(mnn.pseudo),
                          'cell_ID' = gsub(".*_", replacement = '', 
                                           x = names(mnn.pseudo)),
                          'condition' = gsub("_.*", replacement = '', 
                                             x = names(mnn.pseudo)))
mnn.pseudo


table(DataFrameHDS[,c(6,5)])
# HepatocyteType NC 3 months NC 9 months NASH 3 months NASH 9 months
# PP-Hep            1918        1051          1907          1169
# Int-Hep           2156        1196          3083          1110
# PC-Hep             195         122           343           231
# mNASH-Hep1           1           1            42           526
# mNASH-Hep2          44          29           205           496

table(DataFrameHDS[,5])
# NC 3 months   NC 9 months NASH 3 months NASH 9 months 
# 4314          2399          5580          3532 


### 

