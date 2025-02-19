### CellFatrAnalysis with other gene sets ### 

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

genefilter <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv'
  )

geneNoTInHDAG <- setdiff(genefilter$Genes, HDAG$gene_symbol)

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
  
  Seurat3moNC@meta.data$condition <- '3m NC'
  Seurat3moNC@meta.data$HDS <- unname(HDSSeurat3moNC)
  
  
  # dimension reduction
  Seurat3moNC <- RunPCA(Seurat3moNC)
  Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:20)
  DimPlot(Seurat3moNC, label = TRUE)
  
  remove(temp)
  remove(counts3moNC)
  
  
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
  Seurat9moNC@meta.data$condition <- '9m NC'
  Seurat9moNC@meta.data$HDS <- unname(HDSSeurat9moNC)
  
  # dimension reduction & PT trajectory inference
  Seurat9moNC <- RunPCA(Seurat9moNC)
  Seurat9moNC <- RunUMAP(Seurat9moNC, dims = 1:20)
  Seurat9moNC <- FindNeighbors(Seurat9moNC, dims = 1:20)
  Seurat9moNC <- FindClusters(Seurat9moNC,
                              resolution = 0.2)
  
  remove(temp)
  remove(counts9moNC)
  
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
  
  Seurat3moNASH@meta.data$condition <- '3m NASH'
  Seurat3moNASH@meta.data$HDS <- unname(HDSSeurat3moNASH)
  
  remove(temp)
  remove(counts3moNASH)
  gc()
  
  #dimension reduction
  Seurat3moNASH <- RunPCA(Seurat3moNASH)
  Seurat3moNASH <- RunUMAP(Seurat3moNASH, dims = 1:14)
  Seurat3moNASH <- FindNeighbors(Seurat3moNASH, dims = 1:14)
  Seurat3moNASH <- FindClusters(Seurat3moNASH, resolution = 0.1)
  
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
  
  Seurat9moNASH <- SCTransform(Seurat9moNASH)
  
  ## calculate HDS
  HDSSeurat9moNASH <- DS_calc.func(exprMatrices = 
                                     GetAssayData(Seurat9moNASH,
                                                  assay = 'SCT',
                                                  layer = 'counts'),
                                   DSignature = HDAG)
  
  Seurat9moNASH@meta.data$condition <- '9m NASH'
  Seurat9moNASH@meta.data$HDS <- unname(HDSSeurat9moNASH)
  
  # dimension reduction & PT trajectory inference
  Seurat9moNASH <- RunPCA(Seurat9moNASH)
  Seurat9moNASH <- RunUMAP(Seurat9moNASH, dims = 1:14)
  Seurat9moNASH <- FindNeighbors(Seurat9moNASH, dims = 1:14)
  Seurat9moNASH <- FindClusters(Seurat9moNASH, resolution = 0.3)
  
  remove(counts9moNASH, temp)
  }

# II. CellFate: Load gene sets 
{
  
  ## KRAS up-regulated: marker for cancer? 
  
  KRAS_upGenes <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_KRAS_SIGNALING_UP.v2023.2.Mm.csv',
                           header = FALSE)
  KRAS_upGenes <- KRAS_upGenes$V1
  
  
  ## 
  apoptosisFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_APOPTOSIS.v2023.2.Mm_onlyGenes.csv',
                                sep = ';',
                                header = FALSE,
                                stringsAsFactors = FALSE
  )
  apoptosisFeatures <- c(apoptosisFeatures$V1)
  
  
  # Oncogene induced cell sescence (MSibDB)
  # Mouse Gene Set: GOBP_ONCOGENE_INDUCED_CELL_SENESCENCE
  # saved at:'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/GOBP_ONCOGENE_INDUCED_CELL_SENESCENCE.v2023.2.Mm.tsv',
  
  oncoSenGenes <- c('Hmga1b',
                    'Cdkn1a',
                    'Cdkn2a',
                    'Hmga1',
                    'Hmga2',
                    'Hras',
                    'Pml',
                    'Spi1')
}

{
  DoHeatmap(
    Seurat9moNASH,
    features = oncoSenGenes,
    group.by = "CellType",
      

  )
  
  DoHeatmap(
    Seurat9moNASH,
    features = apoptosisFeatures,
    group.by = "CellType",
    
  )
}
