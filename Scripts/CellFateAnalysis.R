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

geneNotInHDAG <- setdiff(genefilter$Genes, HDAG$gene_symbol)

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


# I. Read raw data & process & calculate all values: HDS
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
  
  #Senescence markers: SenMayo
  senMayoGenes <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/SenMayoGeneSetMouse.csv',
                           sep = ';')
  mayoFeatures <- unlist(senMayoGenes$Gene.murine.)
  
  #GSEA Hallmark Apoptosis
  # apoptosisFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_APOPTOSIS.v2023.2.Mm_onlyGenes.csv',
  #                               sep = ';',
  #                               header = FALSE,
  #                               stringsAsFactors = FALSE
  # )
  # apoptosisFeatures <- c(apoptosisFeatures$V1)
  
  # Oncogene induced cell sescence (MSibDB)
  # Mouse Gene Set: GOBP_ONCOGENE_INDUCED_CELL_SENESCENCE
  # saved at:'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/GOBP_ONCOGENE_INDUCED_CELL_SENESCENCE.v2023.2.Mm.tsv',
  
  # oncoSenGenes <- c('Hmga1b',
  #                   'Cdkn1a',
  #                   'Cdkn2a',
  #                   'Hmga1',
  #                   'Hmga2',
  #                   'Hras',
  #                   'Pml',
  #                   'Spi1')
  # 
  
 # # Mouse Gene Set: HALLMARK_G2M_CHECKPOINT for proliferation? 
 #  
 #  G2CheckpointFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_G2M_CHECKPOINT.v2023.2.Mm.csv',
 #                                sep = ';',
 #                                header = FALSE,
 #                                stringsAsFactors = FALSE
 #  )
 #  G2CheckpointFeatures <- c(G2CheckpointFeatures$V1)
  # 
  # # REACTOME_NF_KB_IS_ACTIVATED_AND_SIGNALS_SURVIVAL
  # 
  # SignalSurvival <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/REACTOME_NF_KB_IS_ACTIVATED_AND_SIGNALS_SURVIVAL.csv',
  #                            sep = ';',
  #                            header = FALSE)
  # SignalSurvival <- SignalSurvival$V1
  # 
  
}

## III. Do gene sets intersect with 42 HDS genes?
{
  # Top 42 genes
  intersect(HDAG$gene_symbol[1:42], apoptosisFeatures)
  intersect(HDAG$gene_symbol[1:42], mayoFeatures)
  intersect(HDAG$gene_symbol[1:42], G2CheckpointFeatures)
  
  # All hepatocyte damage associated genes 
  intersect(HDAG$gene_symbol, apoptosisFeatures)
  intersect(HDAG$gene_symbol, mayoFeatures)
  intersect(HDAG$gene_symbol, G2CheckpointFeatures)
  intersect(apoptosisFeatures, geneNoTInHDAG)
  
  
  ## Do a Fischer-Test --> 
}

### IV. Calculate enrichment for gene sets with AUCell 
{
  cells_rankings3moNC <- AUCell_buildRankings(LayerData(Seurat3moNC, 
                                                   assay = "SCT", 
                                                   layer = "counts"), 
                                         plotStats = TRUE)
  
  SenMayocells_AUC3moNC <- AUCell_calcAUC(mayoFeatures, 
                              cells_rankings3moNC)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	44 (37% of 118)
  
  
  AUCell_exploreThresholds(SenMayocells_AUC3moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
    if(identical(colnames(SenMayocells_AUC3moNC@assays@data$AUC), 
            rownames(Seurat3moNC@meta.data)) == TRUE){
      Seurat3moNC@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC3moNC@assays@data$AUC)
    }
  
#### Apoptosis 3moNC
  
  Apoptosiscells_AUC3moNC <- AUCell_calcAUC(apoptosisFeatures, 
                                          cells_rankings3moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	26 (16% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC3moNC@assays@data$AUC), 
               rownames(Seurat3moNC@meta.data)) == TRUE){
    Seurat3moNC@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNC@assays@data$AUC)
  }
  
#### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC3moNC <- AUCell_calcAUC(G2CheckpointFeatures, 
                                            cells_rankings3moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	22 (11% of 195)
  
  AUCell_exploreThresholds(G2cells_AUC3moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC3moNC@assays@data$AUC), 
               rownames(Seurat3moNC@meta.data)) == TRUE){
    Seurat3moNC@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC3moNC@assays@data$AUC)
  }
  
 # in 3 NC Senescence signal the highest
  
  ## 9mo NC ############
  
  cells_rankings9moNC <- AUCell_buildRankings(LayerData(Seurat9moNC, 
                                                        assay = "SCT", 
                                                        layer = "counts"), 
                                              plotStats = TRUE)
  
  SenMayocells_AUC9moNC <- AUCell_calcAUC(mayoFeatures, 
                                          cells_rankings9moNC)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	41 (35% of 118)
  # 
  
  AUCell_exploreThresholds(SenMayocells_AUC9moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(SenMayocells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    Seurat9moNC@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC9moNC@assays@data$AUC)
  }
  
  # Apoptosis 
  
  Apoptosiscells_AUC9moNC <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings9moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	32 (20% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC9moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    Seurat9moNC@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC9moNC@assays@data$AUC)
  }
  
  #### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC9moNC <- AUCell_calcAUC(G2CheckpointFeatures, 
                                     cells_rankings9moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	19 (10% of 195)
  
  AUCell_exploreThresholds(G2cells_AUC9moNC, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    Seurat9moNC@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC9moNC@assays@data$AUC)
  }
  
  ## 3mo NASH ############
  
  cells_rankings3moNASH <- AUCell_buildRankings(LayerData(Seurat3moNASH, 
                                                        assay = "SCT", 
                                                        layer = "counts"), 
                                              plotStats = TRUE)
  
  SenMayocells_AUC3moNASH <- AUCell_calcAUC(mayoFeatures, 
                                          cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	33 (28% of 118)
   
  AUCell_exploreThresholds(SenMayocells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(SenMayocells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC3moNASH@assays@data$AUC)
  }
  
  # Apoptosis 
  
  Apoptosiscells_AUC3moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the data set: 
  #   geneSet: 	20 (12% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNASH@assays@data$AUC)
  }
  
  #### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC3moNASH <- AUCell_calcAUC(G2CheckpointFeatures, 
                                     cells_rankings3moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	14 (7% of 195)
  
  AUCell_exploreThresholds(G2cells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC3moNASH@assays@data$AUC)
  }
  
  ## 3mo NASH ############
  
  cells_rankings3moNASH <- AUCell_buildRankings(LayerData(Seurat3moNASH, 
                                                          assay = "SCT", 
                                                          layer = "counts"), 
                                                plotStats = TRUE)
  
  SenMayocells_AUC3moNASH <- AUCell_calcAUC(mayoFeatures, 
                                            cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	33 (28% of 118)
  
  AUCell_exploreThresholds(SenMayocells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(SenMayocells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC3moNASH@assays@data$AUC)
  }
  
  # Apoptosis 
  
  Apoptosiscells_AUC3moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                              cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the data set: 
  #   geneSet: 	20 (12% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNASH@assays@data$AUC)
  }
  
  #### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC3moNASH <- AUCell_calcAUC(G2CheckpointFeatures, 
                                       cells_rankings3moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	14 (7% of 195)
  
  AUCell_exploreThresholds(G2cells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC3moNASH@assays@data$AUC)
  }
  
  
  ## 3mo NASH ############
  
  cells_rankings3moNASH <- AUCell_buildRankings(LayerData(Seurat3moNASH, 
                                                          assay = "SCT", 
                                                          layer = "counts"), 
                                                plotStats = TRUE)
  
  SenMayocells_AUC3moNASH <- AUCell_calcAUC(mayoFeatures, 
                                            cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	33 (28% of 118)
  
  AUCell_exploreThresholds(SenMayocells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(SenMayocells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC3moNASH@assays@data$AUC)
  }
  
  # Apoptosis 
  
  Apoptosiscells_AUC3moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                              cells_rankings3moNASH)
  
  # Genes in the gene sets NOT available in the data set: 
  #   geneSet: 	20 (12% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNASH@assays@data$AUC)
  }
  
  #### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC3moNASH <- AUCell_calcAUC(G2CheckpointFeatures, 
                                       cells_rankings3moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	14 (7% of 195)
  
  AUCell_exploreThresholds(G2cells_AUC3moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC3moNASH@assays@data$AUC)
  }
  
  ## 9mo NASH ############
  
  cells_rankings9moNASH <- AUCell_buildRankings(LayerData(Seurat9moNASH, 
                                                          assay = "SCT", 
                                                          layer = "counts"), 
                                                plotStats = TRUE)
  
  SenMayocells_AUC9moNASH <- AUCell_calcAUC(mayoFeatures, 
                                            cells_rankings9moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	24 (20% of 118)
  
  AUCell_exploreThresholds(SenMayocells_AUC9moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(SenMayocells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    Seurat9moNASH@meta.data$AUCell_SenMayo <- c(SenMayocells_AUC9moNASH@assays@data$AUC)
  }
  
  # Apoptosis 
  
  Apoptosiscells_AUC9moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                              cells_rankings9moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	13 (8% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC9moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(Apoptosiscells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    Seurat9moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC9moNASH@assays@data$AUC)
  }
  
  #### Proliferation-g2checkpoint 3moNC
  
  G2cells_AUC9moNASH <- AUCell_calcAUC(G2CheckpointFeatures, 
                                       cells_rankings9moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	11 (6% of 195)
   
  AUCell_exploreThresholds(G2cells_AUC9moNASH, 
                           plotHist=TRUE, 
                           assign=TRUE)
  ## save this plot 
  
  if(identical(colnames(G2cells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    Seurat9moNASH@meta.data$AUCell_G2Checkpoint <- c(G2cells_AUC9moNASH@assays@data$AUC)
  }
  
  ## Signal Survival 
  
  SignalSurvival_AUC9moNASH <- AUCell_calcAUC(SignalSurvival,
                                              cells_rankings9moNASH)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	1 (7% of 14)
  
  AUCell_exploreThresholds(SignalSurvival_AUC9moNASH,
                           plotHist = TRUE,
                           assignCells = TRUE)
  
  ## THE WEIRDEST CURVE I HAVE SEEN!! 
  
  if(identical(colnames(SignalSurvival_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    Seurat9moNASH@meta.data$AUCell_SignalSurvival <- c(SignalSurvival_AUC9moNASH@assays@data$AUC)
  }
  
  
  

}

### V. Format all metadata into one data frame to ease plotting
{
  DataFrameAll <- data.frame(Seurat3moNC@meta.data)
  
  
  DataFrameAll <- rbind.data.frame(DataFrameAll,
                                   data.frame(Seurat9moNC@meta.data[,-c(12,13)]),
                                   data.frame(Seurat3moNASH@meta.data[,-c(12,13)]),
                                   data.frame(Seurat9moNASH@meta.data[,-c(12,13)]))
                                   
  DataFrameAll$condition <- factor(DataFrameAll$condition,
                                   ordered = TRUE,
                                   levels = c("3m NC","9m NC", 
                                              "3m NASH", "9m NASH" ))
  
  
  DataFrameAll$CellType <- factor(DataFrameAll$CellType,
                                        ordered = TRUE,
                                        levels = c("PP-Hep",
                                                   "Int-Hep",
                                                   "PC-Hep",
                                                   "mNASH-Hep1",
                                                   "mNASH-Hep2"))
  
}

### VI. Violin plot of AUCell results
{
  ggplot( DataFrameAll,
          aes( x = CellType,
               y = AUCell_SenMayo,
               color = CellType
          )) + 
    geom_violin(trim = FALSE) + 
    scale_color_colorblind() +
    theme( text = element_text(size=25),
           plot.title = element_text(size=24),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "SenMayo AUCell Score", x = '') + 
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") 
  
  ggplot( DataFrameAll,
          aes( x = condition,
               y = AUCell_SenMayo,
               color = condition
          )) + 
    geom_violin(trim = FALSE) + 
    scale_color_colorblind() +
    theme( text = element_text(size=25),
           plot.title = element_text(size=24),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "SenMayo AUCell Score", x = '') + 
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") 
   
   ggplot( DataFrameAll,
            aes( x = CellType,
                 y = AUCell_SenMayo,
                 color = CellType
            )) + 
    geom_violin(trim = FALSE) + 
    scale_color_colorblind() +
    theme( text = element_text(size=25),
           plot.title = element_text(size=24),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "SenMayo AUCell Score", x = '') + 
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") +
    facet_wrap(~condition)
   
   ggplot( DataFrameAll,
           aes( x = condition,
                y = AUCell_SenMayo,
                color = condition
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "SenMayo AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange")
   
   
   ### Apoptosis
   
   ggplot( DataFrameAll,
           aes( x = CellType,
                y = AUCell_Apoptosis,
                color = CellType
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "Apoptosis AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") 
   
   ggplot( DataFrameAll,
           aes( x = condition,
                y = AUCell_Apoptosis,
                color = condition
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "Apoptosis AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") 
   
   ggplot( DataFrameAll,
           aes( x = CellType,
                y = AUCell_Apoptosis,
                color = CellType
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "Apoptosis AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") +
     facet_wrap(~condition)
   
   ### G2-Checkpoint
   
   ggplot( DataFrameAll,
           aes( x = CellType,
                y = AUCell_G2Checkpoint,
                color = CellType
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "G2-checkp. AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") 
   
   ggplot( DataFrameAll,
           aes( x = condition,
                y = AUCell_G2Checkpoint,
                color = condition
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "G2-checkp. AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") 
   
   ggplot( DataFrameAll,
           aes( x = CellType,
                y = AUCell_G2Checkpoint,
                color = CellType
           )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_blank()) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Health Status")) +
     labs(y = "G2-checkp. AUCell Score", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") +
     facet_wrap(~condition)
  
}
## --> all of three cell fates have stronger signals in NASH 9 month sample

### VI. Correlation in all cells devided by hepatocyte zone
{
  ### Senescence & HDS
  # Scatter plots with Pearson correlation and linear regression line
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_SenMayo_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  # Scatter plots with Spearman correlation and linear regression line
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'spearman')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_SenMayo_Spearman_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  # Senescence & apoptosis
  ggplot(DataFrameAll,
         aes(x = AUCell_Apoptosis, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'pearson')
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_Apoptosis_SenMayo_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  # Senescence & G2-checkpoint
  ggplot(DataFrameAll,
         aes(x = AUCell_G2Checkpoint, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'pearson')
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_Apoptosis_G2Checkpoint_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  ### Senescence & HDS
  # Scatter plots with Pearson correlation and linear regression line
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_SenMayo_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  ### Apoptosis and HDS
  # Scatter plots with Spearman correlation and linear regression line
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_Apoptosis)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'spearman')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_Apoptosis_Spearman_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  # Scatter plots with Pearson correlation and linear regression line  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_Apoptosis)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_Apoptosis_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  # Apoptosis &v G2-Checkpoint
  # Scatter plots with Pearson correlation and linear regression line  
  ggplot(DataFrameAll,
         aes(x = AUCell_G2Checkpoint, 
             y = AUCell_Apoptosis)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_G2Checkpoint_Apoptosis_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y =  AUCell_G2Checkpoint)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_G2Checkpoint_HDS_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  
  
    
}

#### Activity scores in HDS bins - all samples
{
  hist(DataFrameAll$HDS)
  DataFrameAll <- DataFrameAll %>% mutate(HDS_bin = cut(HDS, breaks=20))
  
  ggplot(data = DataFrameAll, 
                 aes(x=HDS_bin, y=AUCell_SenMayo, fill = CellType) ) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAll_SenMayo_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y=AUCell_SenMayo, fill = CellType) ) +
    geom_boxplot(outlier.size = 0.3) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAll_SenMayo_BoxPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y=AUCell_Apoptosis, fill = CellType) ) +
    geom_boxplot(outlier.size = 0.3) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAll_Apoptosis_BoxPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y=AUCell_G2Checkpoint, fill = CellType) ) +
    geom_boxplot(outlier.size = 0.3) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAll_G2Checkpoint_BoxPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  
}


#### ALL PLOTS ONE MOUSE: 9month-NASH
{
  # HDS & cell fate markers activity score correlation plots 
  ggplot(DataFrameAll[DataFrameAll$condition == '9m NASH',],
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9mNASHScatterplots_HDS_SenMayo_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  
  ggplot(DataFrameAll[DataFrameAll$condition == '9m NASH',],
         aes(x = HDS, 
             y = AUCell_Apoptosis)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9mNASHScatterplots_HDS_Apoptosis_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  ggplot(DataFrameAll[DataFrameAll$condition == '9m NASH',],
         aes(x = HDS, 
             y = AUCell_G2Checkpoint)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9mNASHScatterplots_HDS_G2Checkpoint_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '9m NASH',], 
         aes(x=HDS_bin, 
             y=AUCell_SenMayo, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9moNASH_SenMayo_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '9m NASH',], 
         aes(x=HDS_bin, 
             y=AUCell_Apoptosis, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9moNASH_Apoptosis_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '9m NASH',], 
         aes(x=HDS_bin, 
             y=AUCell_G2Checkpoint, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9moNASH_G2Checkpoint_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
}

#### ALL PLOTS ONE MOUSE: 3month-NC
{
  # HDS & cell fate markers activity score correlation plots 
  ggplot(DataFrameAll[DataFrameAll$condition == '3m NC',],
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3mNCScatterplots_HDS_SenMayo_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  
  ggplot(DataFrameAll[DataFrameAll$condition == '3m NC',],
         aes(x = HDS, 
             y = AUCell_Apoptosis)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3mNCScatterplots_HDS_Apoptosis_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  ggplot(DataFrameAll[DataFrameAll$condition == '3m NC',],
         aes(x = HDS, 
             y = AUCell_G2Checkpoint)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3mNCScatterplots_HDS_G2Checkpoint_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '3m NC' & DataFrameAll$CellType != 'mNASH-Hep1',], 
         aes(x=HDS_bin, 
             y=AUCell_SenMayo, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3moNC_SenMayo_ViolinPlots_ByCellType.png'),
         width = 30, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll[DataFrameAll$condition == '3m NC' & DataFrameAll$CellType != 'mNASH-Hep1',], 
         aes(x=HDS_bin, 
             y=AUCell_SenMayo, 
             fill = CellType) ) +
    geom_boxplot(outlier.size = 0.1) +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3moNC_SenMayo_BoxPlots_ByCellType.png'),
         width = 30, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll[DataFrameAll$condition == '3m NC' & DataFrameAll$CellType != 'mNASH-Hep1',], 
         aes(x=HDS_bin, 
             y=AUCell_Apoptosis, 
             fill = CellType) ) +
    geom_boxplot(outlier.size = 0.1) +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3moNC_Apoptosis_BoxPlots_ByCellType.png'),
         width = 30, 
         height = 20, 
         units = 'cm')
  
  ggplot(data = DataFrameAll[DataFrameAll$condition == '3m NC' & DataFrameAll$CellType != 'mNASH-Hep1',], 
         aes(x=HDS_bin, 
             y=AUCell_G2Checkpoint, 
             fill = CellType) ) +
    geom_boxplot(outlier.size = 0.1) +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao3moNC_G2Checkpoint_BoxPlots_ByCellType.png'),
         width = 30, 
         height = 20, 
         units = 'cm')
  
  
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '9m NASH',], 
         aes(x=HDS_bin, 
             y=AUCell_Apoptosis, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9moNASH_Apoptosis_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  # separated by zone 
  ggplot(data = DataFrameAll[DataFrameAll$condition == '9m NASH',], 
         aes(x=HDS_bin, 
             y=AUCell_G2Checkpoint, 
             fill = CellType) ) +
    geom_violin() +  
    scale_fill_colorblind() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    facet_grid(~CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao9moNASH_G2Checkpoint_ViolinPlots_ByCellType.png'),
         width = 45, 
         height = 20, 
         units = 'cm')
  
  
  
}

#### 

## To do: correlate with HDS? per zone? also correlate with eachother as sanity check?

## make the HDS bins plot
