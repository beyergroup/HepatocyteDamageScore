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
  
  mergedObj <- merge(x = Seurat3moNC, 
                     y = list(Seurat9moNC, Seurat3moNASH, Seurat9moNASH), 
                     merge.data = FALSE)
  
  mergedObjExpr <- LayerData(mergedObj, 
                             assay = "SCT", 
                             layer = "counts")
  

  
  }


# I. Read in pseudo-proliferation index genes 
# from ... 
# list originally human --> converted to mouse orthologs manually using MGI

{
  proliferationIndex <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/Pseudo-Proliferation-Index-Genes_ConvertedToMouseByPaula.csv',
                                sep = ';',
                                header = TRUE,
                                stringsAsFactors = FALSE)
  
  proliferationIndex <- proliferationIndex[-c(1,2),]
  colnames(proliferationIndex) <- c('human_gene_names', 'UNIPROT_ID', 'mouse_gene' )
  
  # Do they intersect with HDS?
  
  intersect(HDAG$gene_symbol[1:42], proliferationIndex$mouse_gene)
  # no intersection
}


# Heatmap for each mouse
{
  DoHeatmap(Seurat3moNC, features = proliferationIndex$mouse_gene,
            assay = 'SCT', slot = 'counts', group.by = 'condition')
  DoHeatmap(Seurat9moNC, features = proliferationIndex$mouse_gene,
            assay = 'SCT', slot = 'counts', group.by = 'condition')
  DoHeatmap(Seurat3moNASH, features = proliferationIndex$mouse_gene,
            assay = 'SCT', slot = 'counts', group.by = 'condition')
  DoHeatmap(Seurat9moNASH, features = proliferationIndex$mouse_gene,
            assay = 'SCT', slot = 'counts', group.by = 'condition')

}
## heatmaps show there is almost no expression of those genes
# just 


# Prol. Index Aucell merged samples
{
  
  cells_rankingsAll <- AUCell_buildRankings(mergedObjExpr, 
                                            plotStats = TRUE)
  
  Prolcells_All <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                 cells_rankingsAll)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	1 (4% of 25)
  
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/AUCellExploreThresholds_ProliferationIndex_MergedXiaoSamples.pdf") )
  
  print(AUCell_exploreThresholds(Prolcells_All, 
                           plotHist = TRUE, 
                           assign = TRUE))
  
  dev.off()
  
  if(identical(colnames(Prolcells_All@assays@data$AUC), 
               rownames(mergedObj@meta.data)) == TRUE){
    print('all good')
    mergedObj@meta.data$AUCell_ProlifIndex <- c(Prolcells_All@assays@data$AUC)
  }
  
  
}

## Calculate AUCell scores  (not merged)
#  3 month NC
{
  cells_rankings3moNC <- AUCell_buildRankings(LayerData(Seurat3moNC, 
                                                        assay = "SCT", 
                                                        layer = "counts"), 
                                              plotStats = TRUE)
  cells_AUC3moNC <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                   cells_rankings3moNC)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	1 (4% of 25)
  
  
  AUCell_exploreThresholds(cells_AUC3moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  ## save this plot 
  
  if(identical(colnames(cells_AUC3moNC@assays@data$AUC), 
               rownames(Seurat3moNC@meta.data)) == TRUE){
    Seurat3moNC@meta.data$AUCell_ProliferationIndex <- c(cells_AUC3moNC@assays@data$AUC)
  print('yes')
    }
  
  
}

#  9 month NC
{
  cells_rankings9moNC <- AUCell_buildRankings(LayerData(Seurat9moNC, 
                                                        assay = "SCT", 
                                                        layer = "counts"), 
                                              plotStats = TRUE)
  
  cells_AUC9moNC <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                   cells_rankings9moNC)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	1 (4% of 25)
  # 
  
  AUCell_exploreThresholds(cells_AUC9moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  ## save this plot 
  
  if(identical(colnames(cells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    Seurat9moNC@meta.data$AUCell_ProliferationIndex <- c(cells_AUC9moNC@assays@data$AUC)
    print('yes')
  }
  
  
}

#  3 month NASH
{
  cells_rankings3moNASH <- AUCell_buildRankings(LayerData(Seurat3moNASH, 
                                                        assay = "SCT", 
                                                        layer = "counts"), 
                                              plotStats = TRUE)
  
  cells_AUC3moNASH <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                   cells_rankings3moNASH)
  

  AUCell_exploreThresholds(cells_AUC3moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  ## save this plot 
  
  if(identical(colnames(cells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    Seurat3moNASH@meta.data$AUCell_ProliferationIndex <- c(cells_AUC3moNASH@assays@data$AUC)
    print('yes')
  }
  
}

#  9 month NASH
{
  cells_rankings9moNASH <- AUCell_buildRankings(LayerData(Seurat9moNASH, 
                                                          assay = "SCT", 
                                                          layer = "counts"), 
                                                plotStats = TRUE)
  
  cells_AUC9moNASH <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                     cells_rankings9moNASH)
  
  
  AUCell_exploreThresholds(cells_AUC9moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  ## save this plot 
  
  if(identical(colnames(cells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    Seurat9moNASH@meta.data$AUCell_ProliferationIndex <- c(cells_AUC9moNASH@assays@data$AUC)
    print('yes')
  }
  
}


# Senescence 
{
  #Senescence markers: SenMayo
  senMayoGenes <- read.csv(
    file = 
      'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/SenMayoGeneSetMouse.csv',
    sep = ';')
  mayoFeatures <- unlist(senMayoGenes$Gene.murine.)
  
  intersect(HDAG$gene_symbol[1:42], mayoFeatures)

  
  SenCells_All <- AUCell_calcAUC(mayoFeatures,
                                 cells_rankingsAll)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	44 (37% of 118)
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/AUCellExploreThresholds_SenMayoIndex_MergedXiaoSamples.pdf") )
  
  # Plot AUC histogram
  AUCell_exploreThresholds(SenCells_All, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  dev.off()
  
  if(identical(colnames(SenCells_All@assays@data$AUC), 
               rownames(mergedObj@meta.data)) == TRUE){
    print('all good')
    mergedObj@meta.data$AUCell_SenMayo <- c(SenCells_All@assays@data$AUC)
  }
  

  # DoHeatmap(Seurat3moNC, features = mayoFeatures,
  #           assay = 'SCT', slot = 'counts')
  # DoHeatmap(Seurat9moNC, features = mayoFeatures,
  #           assay = 'SCT', slot = 'counts')
  # DoHeatmap(Seurat3moNASH, features = mayoFeatures,
  #           assay = 'SCT', slot = 'counts')
  # DoHeatmap(Seurat9moNASH, features = mayoFeatures,
  #           assay = 'SCT', slot = 'counts')
  # 
  Sencells_AUC3moNC <- AUCell_calcAUC(mayoFeatures, 
                                      cells_rankings3moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	44 (37% of 118)
  
  AUCell_exploreThresholds(Sencells_AUC3moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Sencells_AUC3moNC@assays@data$AUC), 
               rownames(Seurat3moNC@meta.data)) == TRUE){
    print('all good')
    Seurat3moNC@meta.data$AUCell_SenMayo <- c(Sencells_AUC3moNC@assays@data$AUC)
  }
  
  Sencells_AUC9moNC <- AUCell_calcAUC(mayoFeatures,
                                      cells_rankings9moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	41 (35% of 118)
  
  AUCell_exploreThresholds(Sencells_AUC9moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Sencells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    print('all good')
    Seurat9moNC@meta.data$AUCell_SenMayo <- c(Sencells_AUC9moNC@assays@data$AUC)
  }
  
  Sencells_AUC3moNASH <- AUCell_calcAUC(mayoFeatures, 
                                        cells_rankings3moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	33 (12% of 118)
  
  AUCell_exploreThresholds(Sencells_AUC3moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Sencells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    print('all good')
    Seurat3moNASH@meta.data$AUCell_SenMayo <- c(Sencells_AUC3moNASH@assays@data$AUC)
  }
  
  Sencells_AUC9moNASH <- AUCell_calcAUC(mayoFeatures, 
                                        cells_rankings9moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	24 (8% of 118)
  
  AUCell_exploreThresholds(Sencells_AUC9moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Sencells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    print('all good')
    Seurat9moNASH@meta.data$AUCell_SenMayo <- c(Sencells_AUC9moNASH@assays@data$AUC)
  }
  
  # Add results to dataframe
  
  DataFrameAll$AUCell_SenMayo<- c(Seurat3moNC@meta.data$AUCell_SenMayo,
                                  Seurat9moNC@meta.data$AUCell_SenMayo,
                                  Seurat3moNASH@meta.data$AUCell_SenMayo,
                                  Seurat9moNASH@meta.data$AUCell_SenMayo)
  
  saveRDS(DataFrameAll, file = paste0)
  
}

## STORE IN DATAFRAME FOR PLOTTING 
{
  ### Format Results for plotting 
  
  
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
    
    DataFrameAll <- DataFrameAll %>% mutate(HDS_bin = cut(HDS, breaks=15))
    
    View(DataFrameAll)
    
    
    ## Alternative merged Data
    mergedDataFrame <- mergedObj@meta.data
    
    mergedDataFrame$CellType <- factor(mergedDataFrame$CellType,
                                       ordered = TRUE,
                                       levels = c("PP-Hep",
                                                  "Int-Hep",
                                                  "PC-Hep",
                                                  "mNASH-Hep1",
                                                  "mNASH-Hep2"))
    
    mergedDataFrame$condition <- factor(mergedDataFrame$condition,
                                     ordered = TRUE,
                                     levels = c("3m NC","9m NC", 
                                                "3m NASH", "9m NASH" ))
    
    mergedDataFrame <- mergedDataFrame %>% mutate(HDS_bin = cut(HDS, breaks=15))
    
    mergedDataFrame <- mergedDataFrame %>% 
      mutate( deltaSmP = AUCell_SenMayo - AUCell_ProlifIndex)
    
    mergedDataFrame$deltaSmP_centered <- 
      mergedDataFrame$deltaSmP - mean(mergedDataFrame$deltaSmP)
    
    
    deltaSmP_HDS_bySample <-
      ggplot(data = mergedDataFrame,
             aes(x = HDS,
                 y = deltaSmP_centered)) + 
      geom_point(size = 0.5)   +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      facet_wrap( ~ condition) +
      stat_cor(method = 'spearman',
               cor.coef.name = 'rho')
    
    AbsDeltaSmP_HDS_bySample <-
      ggplot(data = mergedDataFrame,
             aes(x = HDS,
                 y = abs(deltaSmP_centered))) + 
      geom_point(size = 0.5)   +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      facet_wrap( ~ condition) +
      stat_cor(method = 'spearman',
               cor.coef.name = 'rho')
    
    pdf(height = 8, 
        width = 10,
        file = 
          paste0(outputPath,
                 "Results/CellFateAnalysis/",
                 "Spearman_deltaSmP_HDS_bySampleMergedXiao.pdf")
         )
    
    print(deltaSmP_HDS_bySample)
    
    dev.off()
    
    
    pdf(height = 8, 
        width = 10,
        file = 
          paste0(outputPath,
                 "Results/CellFateAnalysis/",
                 "Spearman_AbsDeltaSmP_HDS_bySampleMergedXiao.pdf")
    )
    
    print(AbsDeltaSmP_HDS_bySample)
    
    dev.off()
    
    

    

    
    ## categorize as senscence or proliferative cellfate
    mergedDataFrame$CellFate <- NA
    mergedDataFrame$CellFate[mergedDataFrame$deltaSmP_centered >= 0.02 ] <- 'S'
    mergedDataFrame$CellFate[mergedDataFrame$deltaSmP_centered <= -0.025 ] <- 'P'
    
    mergedDataFrame$S <- FALSE
    mergedDataFrame$S[mergedDataFrame$CellFate == 'S'] <- TRUE
    
    mergedDataFrame$P <- FALSE
    mergedDataFrame$P[mergedDataFrame$CellFate == 'P'] <- TRUE
    
    
    table(mergedDataFrame$CellFate)
    # P    S 
    # 1046  1474 
    
    # Plot annotated cells
    plotCategoryAssign <-
      ggplot(data = mergedDataFrame,
           aes(x = HDS,
               y = deltaSmP_centered,
               color = CellFate)) + 
      geom_point(size = 0.5)   +
      theme_bw() +
      facet_wrap( ~ condition) +
      theme(text = element_text(size = 20))
    
    pdf( height = 8, width = 12, 
         file = paste0(outputPath,
                       "Results/CellFateAnalysis/centeredDeltaSmP_HDS_CellFateMergedObjectAnnotated.pdf") )
    print(plotCategoryAssign)
    dev.off()
    
    
    
  
}


# Find genes that correlate with HDS 
{
  
  # function to calculate the correlation of genes with HDS column of 
  # provided metadata dataframe (DataFrameAll) for the S and P groups of cells
  
  corGenes <- function(exprMat, DataFrameAll){
    corGroups <- list('corSenes' = vector(),
                      'corProlif' = vector())
    toKeepS <- DataFrameAll[DataFrameAll$S,]
    which(rownames(toKeepS) == colnames(exprMat))
    toKeepExprMatrix <- exprMat[, which(colnames(exprMat) %in% rownames(toKeepS))]
    toKeepS <- toKeepS[colnames(toKeepExprMatrix),]
    
    corGroups[['corSenes']] <- apply(toKeepExprMatrix, 1,
                                     function(X){
                                       cor.test( x = X,
                                                 y = toKeepS$HDS)
                                     })
    
    toKeepP <- DataFrameAll[DataFrameAll$P,]
    which(rownames(toKeepP) == colnames(exprMat))
    toKeepExprMatrix <- exprMat[, which(colnames(exprMat) %in% rownames(toKeepP))]
    toKeepP <- toKeepP[colnames(toKeepExprMatrix),]
    
    corGroups[['corProlif']] <- apply(toKeepExprMatrix, 1,
                                      function(X){
                                        cor.test( x = X,
                                                  y = toKeepP$HDS)
                                      })
    
    return(corGroups)
    
  }
  
  corResMerged <- corGenes(mergedObjExpr, mergedDataFrame)
  
  # Reformat results
  # initiate dataframe
  nGenes <- dim(mergedObjExpr)[1]
  allCorSes <- data.frame(row.names = names(corResMerged[[1]]))
  allCorSes$pvalue_Ses <- vector(length = nGenes)
  allCorSes$ro_Ses <- vector(length = nGenes)
  allCorProlif <- data.frame(row.names = names(corResMerged[[1]]))
  allCorProlif$pvalue_Prolif <- vector(length = nGenes)
  allCorProlif$ro_Prolif <- vector(length = nGenes)
  i <- 1
  #store values in data frames
  
for (i in i:nGenes){
    allCorSes$pvalue_Ses[i] <- corResMerged[[1]][[i]]$p.value
    allCorSes$ro_Ses[i] <- corResMerged[[1]][[i]]$estimate 
    allCorProlif$pvalue_Prolif[i] <- corResMerged[[2]][[i]]$p.value
    allCorProlif$ro_Prolif[i] <- corResMerged[[2]][[i]]$estimate 
    }

  allCorSes <- na.omit(allCorSes)
  allCorProlif <- na.omit(allCorProlif)
  allCorSes$gene <- rownames(allCorSes)
  allCorProlif$gene <- rownames(allCorProlif)
  
  allCorSes$pvalue_adj_Ses <- p.adjust(allCorSes$pvalue_Ses, method = "BH")
  allCorProlif$pvalue_adj_Prolif <- p.adjust(allCorProlif$pvalue_Prolif, method = "BH")
  
  genesPassSes <- allCorSes[which(allCorSes$pvalue_adj_Ses < 0.05 & allCorSes$ro_Ses >= 0.2),]
  genesPassPro<- allCorProlif[which(allCorProlif$pvalue_adj_Prolif < 0.05 & allCorProlif$ro_Prolif >= 0.2),]
  
  ## when mergin around 700 genes are left out, why?? 
  dfForPlot <- merge(allCorSes, allCorProlif, by = 'gene')
  dfForPlot$sigGenes <- rep(FALSE, length(dfForPlot$gene))
  dfForPlot$sigGenes[which(dfForPlot$pvalue_adj_Ses < 0.05 & dfForPlot$pvalue_adj_Prolif < 0.05)] <- TRUE
  
  ggplot(data = dfForPlot, aes(x = ro_Ses, y = ro_Prolif)) + 
    geom_point(size = 0.3, aes(color = sigGenes)) + 
    scale_color_manual( values = c("FALSE" = "black", "TRUE" = "red"))

}

# delta senescence - proliferation 
{
  DataFrameAll <- DataFrameAll %>% 
    mutate( deltaSmP = AUCell_SenMayo - AUCell_ProliferationIndex)
  
  DataFrameAll$deltaSmP_centered <- DataFrameAll$deltaSmP - mean(DataFrameAll$deltaSmP)
  DataFrameAll[order(DataFrameAll$deltaSmP, decreasing = FALSE),'RankByDelta'] <- seq(1:length(DataFrameAll$orig.ident))
  
  # FIGURE 5 a
  # homoskedasticity of deltSmP across HDS 
  ggplot(data = DataFrameAll,
         aes(x = HDS,
             y = abs(deltaSmP_centered))) + 
    geom_point(size = 1)   +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    facet_wrap( ~ condition) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'rho')
  
  # shows increase in variance of delta SmP by HDS > -0.1
  # also shows only positive correlation in hepatocytes from NASH mice
  
  # assign hepatocytes to a cell fate trajectory: 
  # by visually chose threhold 
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, 
             y = deltaSmP_centered)) 
  
  # cells deltaSmP >= 0.05 will be considered on the 'senescence trajectory' 
  
  #E FIgure 5 b: categorize hepatocytes by trajectory
  DataFrameAll$CellFateTrajectory <- NA
  DataFrameAll$CellFateTrajectory[DataFrameAll$deltaSmP_centered >= 0.025 ] <- 'S'
  # cells with tdelta SmP <= 0.025 will be considred on the 'proliferation' trajectory 
  DataFrameAll$CellFateTrajectory[DataFrameAll$deltaSmP_centered <= -0.025 ] <- 'P'
  
  
   table(DataFrameAll$CellFateTrajectory)
   # P   S 
   # 986 618 
   length(DataFrameAll$CellFateTrajectory)
  # [1] 15825
  
  # Plot annotated cells
  ggplot(data = DataFrameAll,
         aes(x = HDS,
             y = deltaSmP_centered,
             color = CellFateTrajectory)) + 
    geom_point(size = 1)   +
    theme_bw() +
    facet_wrap( ~ condition) +
    theme(text = element_text(size = 20))
  
  pdf( height = 8, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/centeredDeltaSmP_HDS_CellFateTrajectoriesAnnotated.pdf") )
  print(violinLCAsamples)
  dev.off()
  
    
  # other figures 
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, 
             y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_CellRankedByDeltaSenMayoMinusProliferationByZone.png'),
         width = 30, 
         height = 15, 
         units = 'cm')
  
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, y = deltaSmP)) +
    geom_jitter(size = 0.1) +
    facet_wrap(~ condition)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_CellRankedByDeltaSenMayoMinusProliferationByCondition.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  ## add HDS dimension by using color 
  
  
  limitBottom5percent <- DataFrameAll[order(DataFrameAll$HDS),'HDS'][round((15825)*0.05)]
  limitTop5percent <- DataFrameAll[order(DataFrameAll$HDS, decreasing = TRUE),'HDS'][round((15825)*0.05)]
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, y = deltaSmP,
             colour = HDS)) +
    geom_point(size = 0.2) +
    facet_wrap(~ condition + CellType) +
    scale_color_distiller(palette = "YlOrBr",
                          direction = 1,
                          limits = c(limitBottom5percent,
                                     limitTop5percent),
                          oob = squish
    )
  # dooes not add much information so not saving it 
  
  
  
  # delta as Y axis and HDS as x axis
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByZone.png'),
         width = 30, 
         height = 15, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ condition)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByCondition.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP_centered)) +
    geom_point(size = 0.1) +
    facet_wrap(~ condition)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationCenteredByCondition.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP_centered)) +
    geom_point(size = 0.5, alpha = 0.5) 
  
  ggsave(filename = paste0(outputPath, '/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationCentered.png'),
         width = 15,
         height = 15,
         units = 'cm')
  
  # destribution for thresholds
  ggplot(DataFrameAll, aes(x=deltaSmP_centered)) + 
    geom_density() + geom_vline(aes(xintercept=median(deltaSmP_centered)),
                               color="blue", size=1) +
    geom_vline(aes(xintercept = quantile(deltaSmP_centered)[2]), color = 'red') +
    geom_vline(aes(xintercept = quantile(deltaSmP_centered)[4]), color = 'red') +
    facet_wrap(~condition)
    
      
  
  # > summary(DataFrameAll$deltaSmP)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # -0.09519  0.02957  0.03898  0.03791  0.04763  0.11430 
  
  DataFrameAll$absoluteDeltaSmPoverThres <- 'fill'
  DataFrameAll$absoluteDeltaSmPoverThres[DataFrameAll$deltaSmP >= 0.06] <- 'pass'
  DataFrameAll$absoluteDeltaSmPoverThres[DataFrameAll$deltaSmP < 0.025] <- 'pass'
  
  ggplot(data = DataFrameAll[order(DataFrameAll$deltaSmP, decreasing = FALSE),], 
         aes(x = HDS, y = deltaSmP, color = absoluteDeltaSmPoverThres)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByCellTypeCellsColoredByThreshold.png'),
         width = 20, 
         height = 25, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y= deltaSmP,
             color = absoluteDeltaSmPoverThres)) +
    geom_jitter(size = 0.3, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("HDS (bins) ") +
    stat_summary(fun = 'mean',
                 geom = "crossbar",
                 width = .5, color = "red") +
    stat_summary(fun.data = "mean_se", 
                 geom = "errorbar", 
                 width = .1, color = 'red')+
    facet_wrap(~CellType)
  
  
  
  
}

# Find genes that correlate with HDS in delSmP in each category ? 
{
  
  DataFrameAll$S <- FALSE
  DataFrameAll$S[DataFrameAll$CellFateTrajectory == 'S'] <- TRUE
  
  DataFrameAll$P <- FALSE
  DataFrameAll$P[DataFrameAll$CellFateTrajectory == 'P'] <- TRUE
  
  
  # output: vector of length(genes expessed)
  temp <- GetAssayData(Seurat3moNC,
                       assay = 'SCT',
                       layer = 'counts')
  
  corGenesHDS <- function(exprMat){
    corGroups <- list('corSenes' = vector(),
                      'corProlif' = vector())
    toKeepS <- DataFrameAll[DataFrameAll$S,]
    which(rownames(toKeepS) == colnames(exprMat))
    toKeepExprMatrix <- exprMat[, which(colnames(exprMat) %in% rownames(toKeepS))]
    toKeepS <- toKeepS[colnames(toKeepExprMatrix),]
    
   corGroups[['corSenes']] <- apply(toKeepExprMatrix, 1,
                                   function(X){
                                     cor.test( x = X,
                                               y = toKeepS$HDS)
                                   })

   toKeepP <- DataFrameAll[DataFrameAll$P,]
   which(rownames(toKeepP) == colnames(exprMat))
   toKeepExprMatrix <- exprMat[, which(colnames(exprMat) %in% rownames(toKeepP))]
   toKeepP <- toKeepP[colnames(toKeepExprMatrix),]

   corGroups[['corProlif']] <- apply(toKeepExprMatrix, 1,
                                    function(X){
                                      cor.test( x = X,
                                                y = toKeepP$HDS)
                                    })
   
   return(corGroups)
    
  }
  
  HDScorGeneExpression <- corGenes(mergedObjExpr)
  
  NC3mo_corGeneDeltaSmP <- corGenes(temp)
  NC9mo_corGeneDeltaSmP <- corGenes(GetAssayData(Seurat9moNC,
                                                 assay = 'SCT',
                                                 layer = 'counts'))
  
  NASH3mo_corGeneDeltaSmP <- corGenes(GetAssayData(Seurat3moNASH,
                                                 assay = 'SCT',
                                                 layer = 'counts'))
  
  NASH9mo_corGeneDeltaSmP <- corGenes(GetAssayData(Seurat9moNASH,
                                                   assay = 'SCT',
                                                   layer = 'counts'))
  
  allCorSes <- data.frame('corEst' = vector(),
                          'p-Value' = vector(),
                          'sample' = vector())
  
  coef(NC3mo_corGeneDeltaSmP[[1]]$Xkr4)
  
  

  
}





## don't actually neeed this
# Center AUCell Score around zero 
{
  DataFrameAll$AUCell_SenMayo_centered <- 
    (DataFrameAll$AUCell_SenMayo - mean(DataFrameAll$AUCell_SenMayo))
  
  DataFrameAll$AUCell_ProliferationIndex_centered <- 
    (DataFrameAll$AUCell_ProliferationIndex - mean(DataFrameAll$ProliferationIndex))
}

# Violin plots of SenMayo AUCell
{
  
  # check AUCell-SenMayo scores for normality, so many data points:
  qqnorm(DataFrameAll$AUCell_SenMayo,
         main = "Q-Q plot", cex.main = 0.85)
  qqline(DataFrameAll$AUCell_SenMayo, col = "blue")
  
  ggplot( DataFrameAll,
          aes( x = CellType,
               y = AUCell_SenMayo,
               color = CellType
          )) + 
    geom_boxplot(trim = FALSE,
                 outlier.size = 1.5,
                 outlier.shape = 21) + 
    scale_color_colorblind() +
    theme( text = element_text(size=18),
           plot.title = element_text(size=18),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title = "Hep. Zone")) +
    labs(y = "SenMayo\n AUCell Score",
         x = '') +
    facet_wrap(~ condition)
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ViolinPlots_AuCell_SenMayo_ColoredZones_Condition_WithMedian.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  
  ggplot( DataFrameAll,
          aes( x = condition,
               y = AUCell_SenMayo,
               colour = condition
          )) + 
    geom_boxplot(trim = FALSE,
                 outlier.size = 1.5,
                 outlier.shape = 21) + 
    scale_color_colorblind() +
    theme( text = element_text(size=18),
           plot.title = element_text(size=18),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title = "Sample")) +
    labs(y = "SenMayo\n AUCell Score",
         x = '') + 
    stat_summary(fun.data = "median_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") 
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ViolinPlots_AuCell_SenMayo_ColoredCondition_WithMedian.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
}

## make scatter plots 
{
  # Pearson Correlation
  
  # by cell type
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoAllScatterplots_HDS_ProliferationIndexGenes_Pearson_LM.png'),
         width = 20, 
         height = 15, 
         units = 'cm')
  
  # by samples + celltype 
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~condition + CellType) + 
    stat_cor()
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/XiaoBySampleByZoneScatterplots_HDS_ProliferationIndexGenes_Pearson_LM.png'),
         width = 35, 
         height = 35, 
         units = 'cm')
  
  
}

# make violin plot
{
  ggplot( DataFrameAll,
          aes( x = condition,
               y = AUCell_ProliferationIndex,
               color = condition
          )) + 
    geom_boxplot(trim = FALSE,
                outlier.size = 1.5,
                outlier.shape = 21) + 
    scale_color_colorblind() +
    theme( text = element_text(size=16),
           plot.title = element_text(size=16),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "Pseudo-Profilieration-Index-Genes\n AUCell Score", x = '') + 
    stat_summary(fun.data = "median_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") 
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ViolinPlots_AuCell_ProliferationIndexGenes_ByCondition_WithMedian.png'),
         width = 18, 
         height = 15, 
         units = 'cm')
  
  ggplot( DataFrameAll,
          aes( x = CellType,
               y = AUCell_ProliferationIndex,
               color = CellType
          )) + 
    geom_boxplot(trim = FALSE,
                 outlier.size = 1.5,
                 outlier.shape = 21) + 
    scale_color_colorblind() +
    theme( text = element_text(size=18),
           plot.title = element_text(size=18),
           legend.position = 'right',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title = "Health Status")) +
    labs(y = "Pseudo-Profilieration-Index-Genes\n AUCell Score",
         x = '') + 
    stat_summary(fun.data = "median_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") +
    facet_wrap(~ condition)
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ViolinPlots_AuCell_ProliferationIndexGenes_ColoredZones_Condition_WithMedian.png'),
         width = 20, 
         height = 20, 
         units = 'cm')
  
  
}


# Correlate Senescence and Proliferation
{
  ggplot(DataFrameAll,
         aes(x = AUCell_SenMayo, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition + CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlots_ProliferationAUCell_SenMayo_LM_PearsonCorrelation.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = AUCell_SenMayo, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType ) + 
    stat_cor()
  
  
  ggplot(DataFrameAll,
         aes(x = AUCell_SenMayo, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition ) + 
    stat_cor()
  
  ggplot(DataFrameAll,
         aes(x = AUCell_SenMayo, 
             y = AUCell_ProliferationIndex,
             color = HDS )) +
    geom_point(size = 0.5) +
    scale_color_distiller(palette = "YlOrBr",
                          direction = 1) + 
    facet_wrap(~ condition + CellType) + 
    geom_smooth() + stat_cor()
  
  ggplot(DataFrameAll,
         aes(x = AUCell_SenMayo, 
             y = AUCell_ProliferationIndex,
             color = HDS )) +
    geom_point(size = 0.5) +
    scale_color_distiller(palette = "YlOrBr",
                          direction = 1) + 
    facet_wrap(~ condition) + 
    geom_smooth() + stat_cor()

    
    
    ggplot(DataFrameAll[DataFrameAll$condition == '3m NASH' | DataFrameAll$condition == '9m NASH',],
           aes(x = AUCell_SenMayo, 
               y = AUCell_ProliferationIndex,
               color = HDS )) +
      geom_point(size = 0.8) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1) + 
      facet_wrap(~ condition + CellType) + 
      geom_smooth(color = 'red') + stat_cor()
    scale_color_colorblind()
    
    ggplot(DataFrameAll[DataFrameAll$condition == '3m NC' | DataFrameAll$condition == '9m NC',],
           aes(x = AUCell_SenMayo, 
               y = AUCell_ProliferationIndex,
               color = HDS )) +
      geom_point(size = 0.8) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1) + 
      facet_wrap(~ condition + CellType) + 
      geom_smooth(color = 'red') + stat_cor()
    scale_color_colorblind() 
    

  


}

# Correlate Senescence and Proliferation with HDS
{
  
  ## Proliferation 
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition + CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlots_Proliferation_HDS_LM_PearsonCorrelation.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlotsBySample_Proliferation_HDS_LM_PearsonCorrelation.png'),
         width = 18, 
         height = 18, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlotsByCellType_Proliferation_HDS_LM_PearsonCorrelation.png'),
         width = 18, 
         height = 20, 
         units = 'cm')
  
  
  ## looks like there is an increase in proliferation with higher damage score 
  # that is maybe hepatocyte-zone specific/dependent 
  
  #Senescence 
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition + CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlots_SenMayo_HDS_LM_PearsonCorrelation.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ condition) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlotsBySample_SenMayo_HDS_LM_PearsonCorrelation.png'),
         width = 18, 
         height = 18, 
         units = 'cm')
  
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~ CellType) + 
    stat_cor(method = 'pearson')
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_ScatterPlotsByCellType_SenMayo_HDS_LM_PearsonCorrelation.png'),
         width = 18, 
         height = 20, 
         units = 'cm')
  
  ggplot(DataFrameAll,
         aes(x = HDS, 
             y = AUCell_SenMayo,
             colour = AUCell_ProliferationIndex)) +
    scale_color_gradient(low = 'yellow', high = 'red',
                         breaks = c(0.014, 0.125)) +
    geom_point(size = 0.4) + 
    facet_wrap(~ CellType) + 
    stat_cor(method = 'pearson')
  
  
  
  ggplot(DataFrameAll,
         aes(x = AUCell_ProliferationIndex,
             color = CellType)) +
    geom_histogram(fill = 'white', binwidth = 0.005) + 
    facet_wrap(~condition)+
    geom_vline(aes(xintercept=mean(DataFrameAll$AUCell_ProliferationIndex)),
               linetype="dashed")
  
  
  
  
}

# Binned HDS 

{
  # Colored by cell types
  # ggplot(data = DataFrameAll, 
  #        aes(x=HDS_bin, y=AUCell_SenMayo, 
  #            colour = CellType, 
  #            fill = CellType) ) +
  #   geom_violin(trim = FALSE) + 
  #   scale_color_colorblind() + 
  #   scale_fill_colorblind() +
  #   theme(axis.text.x = element_text(angle = 90)) +
  #   xlab("HDS (bins) ") 
  # 
  # ggplot(data = DataFrameAll, 
  #        aes(x = HDS_bin, 
  #            y = AUCell_ProliferationIndex, 
  #            colour = CellType, fill = CellType) ) +
  #   geom_violin(trim = FALSE) + 
  #   theme(axis.text.x = element_text(angle = 90)) +
  #   xlab("HDS (bins) ")
  
  # all cells together
  
  plotA <- ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y=AUCell_SenMayo )) +
   # geom_violin() + 
    geom_jitter(size = 0.05, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    ylim(c(0, 0.125)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("HDS (bins) ") +
    stat_summary(fun = 'mean',
                 geom = "crossbar",
                 width = .5, color = "red") +
    stat_summary(fun.data = "mean_se", 
                 geom = "errorbar", 
                 width = .1, color = 'red')
  
  
  plotB <- ggplot(data = DataFrameAll, 
         aes(x = HDS_bin, y = AUCell_ProliferationIndex) ) +
    geom_jitter(size = 0.05, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    ylim(c(0, 0.125)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    stat_summary(fun = 'mean',
                 geom = "crossbar",
                 width = .5, color = "red") +
    stat_summary(fun.data = "mean_se", 
                 geom = "errorbar", 
                 width = .1, color = "red")
  AB_patch <- plotA / plotB
  
  AB_patch
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_Bins_ProliferationSenMayoNoLineWithMeanSE.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  AB_patch2 <- (plotA + geom_hline(yintercept= 0.050550, linetype="dashed", color = "blue") + annotate('text', x=2, y=0.057, label="3rd Qu.", color = 'blue')) / (plotB + geom_hline(yintercept= 0.009148, linetype="dashed", color = "blue") + annotate('text', x=2, y=0.016, label="3rd Qu.", color = 'blue'))
  
  AB_patch2
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_Bins_ProliferationSenMayoWithLineMeanSE.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  
  # bin proliferation 
  DataFrameAll <- DataFrameAll %>% mutate(Proliferation_bin = cut(AUCell_ProliferationIndex, breaks=25))
  
  ggplot(data = DataFrameAll, 
         aes(x=Proliferation_bin, y=AUCell_SenMayo )) +
    geom_violin() + 
    geom_jitter(size = 0.05, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    ylim(c(0, 0.125)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("ProliferationIndex (bins) ") +
    facet_wrap(~ CellType) 
  
  
  ggplot(data = DataFrameAll, 
         aes(x=Proliferation_bin, y=AUCell_SenMayo )) +
    geom_violin() + 
    geom_jitter(size = 0.05, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    ylim(c(0, 0.125)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("ProliferationIndex (bins) ") 
  
  
  # senescence bins 
  
  DataFrameAll <- DataFrameAll %>% mutate(SenMayo_bin = cut(AUCell_SenMayo, breaks=25))
  
  ggplot(data = DataFrameAll, 
         aes(x=SenMayo_bin, y=AUCell_ProliferationIndex )) +
    geom_jitter(size = 0.05, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    ylim(c(0, 0.125)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("SenMAyo AUC score (bins) ") 
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_SenMayo_Bins_Y_axis_Proliferation.png'),
         width = 30, 
         height = 25, 
         units = 'cm')
  
  
  

  
}

# delta senescence - proliferation 
{
  DataFrameAll <- DataFrameAll %>% mutate( deltaSmP = AUCell_SenMayo - AUCell_ProliferationIndex)
  DataFrameAll[order(DataFrameAll$deltaSmP, decreasing = FALSE),'RankByDelta'] <- seq(1:length(DataFrameAll$orig.ident))
  
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_CellRankedByDeltaSenMayoMinusProliferationByZone.png'),
         width = 30, 
         height = 15, 
         units = 'cm')
  
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, y = deltaSmP)) +
    geom_jitter(size = 0.1) +
    facet_wrap(~ condition)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_HDS_CellRankedByDeltaSenMayoMinusProliferationByCondition.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  ## add HDS dimension by using color 
  
  
  limitBottom5percent <- DataFrameAll[order(DataFrameAll$HDS),'HDS'][round((15825)*0.05)]
  limitTop5percent <- DataFrameAll[order(DataFrameAll$HDS, decreasing = TRUE),'HDS'][round((15825)*0.05)]
  
  ggplot(data = DataFrameAll, 
         aes(x = RankByDelta, y = deltaSmP,
             colour = HDS)) +
    geom_point(size = 0.2) +
    facet_wrap(~ condition + CellType) +
    scale_color_distiller(palette = "YlOrBr",
                          direction = 1,
                          limits = c(limitBottom5percent,
                                     limitTop5percent),
                          oob = squish
    )
  # dooes not add much information so not saving it 

  
  
  # delta as Y axis and HDS as x axis
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByZone.png'),
         width = 30, 
         height = 15, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x = HDS, y = deltaSmP)) +
    geom_point(size = 0.1) +
    facet_wrap(~ condition)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByCondition.png'),
         width = 15, 
         height = 15, 
         units = 'cm')
  
  # > summary(DataFrameAll$deltaSmP)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # -0.09519  0.02957  0.03898  0.03791  0.04763  0.11430 
  
  DataFrameAll$absoluteDeltaSmPoverThres <- 'fill'
  DataFrameAll$absoluteDeltaSmPoverThres[DataFrameAll$deltaSmP >= 0.06] <- 'pass'
  DataFrameAll$absoluteDeltaSmPoverThres[DataFrameAll$deltaSmP < 0.025] <- 'pass'
  
  
  
  ggplot(data = DataFrameAll[order(DataFrameAll$deltaSmP, decreasing = FALSE),], 
         aes(x = HDS, y = deltaSmP, color = absoluteDeltaSmPoverThres)) +
    geom_point(size = 0.1) +
    facet_wrap(~ CellType)
  
  ggsave(filename = paste0(outputPath,'/Results/Trajectories/Xiao_XaxisHDS_DeltaSenMayoMinusProliferationByCellTypeCellsColoredByThreshold.png'),
         width = 20, 
         height = 25, 
         units = 'cm')
  
  ggplot(data = DataFrameAll, 
         aes(x=HDS_bin, y= deltaSmP,
             color = absoluteDeltaSmPoverThres)) +
    geom_jitter(size = 0.3, 
                width = 0.3, 
                height = 0,
                alpha = 0.8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("HDS (bins) ") +
    stat_summary(fun = 'mean',
                 geom = "crossbar",
                 width = .5, color = "red") +
    stat_summary(fun.data = "mean_se", 
                 geom = "errorbar", 
                 width = .1, color = 'red')+
    facet_wrap(~CellType)
  
  
 
   
}

# preliminary apoptosis gene set
{
  #GSEA Hallmark Apoptosis
  apoptosisFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_APOPTOSIS.v2023.2.Mm_onlyGenes.csv',
                                sep = ';',
                                header = FALSE,
                                stringsAsFactors = FALSE
  )
  apoptosisFeatures <- c(apoptosisFeatures$V1)
  
  DoHeatmap(Seurat3moNC, features = apoptosisFeatures,
            assay = 'SCT', slot = 'counts')
  DoHeatmap(Seurat9moNC, features = apoptosisFeatures,
            assay = 'SCT', slot = 'counts')
  DoHeatmap(Seurat3moNASH, features = apoptosisFeatures,
            assay = 'SCT', slot = 'counts')
  DoHeatmap(Seurat9moNASH, features = apoptosisFeatures,
            assay = 'SCT', slot = 'counts')
  
  
  Apoptosiscells_AUC3moNC <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings3moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	26 (16% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Apoptosiscells_AUC3moNC@assays@data$AUC), 
               rownames(Seurat3moNC@meta.data)) == TRUE){
    print('all good')
    Seurat3moNC@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNC@assays@data$AUC)
  }
  
  Apoptosiscells_AUC9moNC <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings9moNC)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	32 (20% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC9moNC, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Apoptosiscells_AUC9moNC@assays@data$AUC), 
               rownames(Seurat9moNC@meta.data)) == TRUE){
    print('all good')
    Seurat9moNC@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC9moNC@assays@data$AUC)
  }
  
  Apoptosiscells_AUC3moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings3moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	20 (12% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC3moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Apoptosiscells_AUC3moNASH@assays@data$AUC), 
               rownames(Seurat3moNASH@meta.data)) == TRUE){
    print('all good')
    Seurat3moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC3moNASH@assays@data$AUC)
  }
  
  Apoptosiscells_AUC9moNASH <- AUCell_calcAUC(apoptosisFeatures, 
                                            cells_rankings9moNASH)
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	13 (8% of 161)
  
  AUCell_exploreThresholds(Apoptosiscells_AUC9moNASH, 
                           plotHist = TRUE, 
                           assign = TRUE)
  
  if(identical(colnames(Apoptosiscells_AUC9moNASH@assays@data$AUC), 
               rownames(Seurat9moNASH@meta.data)) == TRUE){
    print('all good')
    Seurat9moNASH@meta.data$AUCell_Apoptosis <- c(Apoptosiscells_AUC9moNASH@assays@data$AUC)
  }
  
  # Add results to dataframe
  
  DataFrameAll$AUCell_Apoptosis <- c(Seurat3moNC@meta.data$AUCell_Apoptosis,
                                     Seurat9moNC@meta.data$AUCell_Apoptosis,
                                     Seurat3moNASH@meta.data$AUCell_Apoptosis,
                                     Seurat9moNASH@meta.data$AUCell_Apoptosis)
}
  

# Correlate Apoptosis and Proliferation
{
  ggplot(DataFrameAll,
         aes(x = AUCell_Apoptosis, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) + 
    scale_color_colorblind() + geom_smooth(method = lm) +
    facet_wrap(~CellType + condition) + 
    stat_cor()
  
  ggplot(DataFrameAll,
         aes(x = AUCell_Apoptosis, 
             y = AUCell_ProliferationIndex)) +
    geom_point(size = 0.1, alpha = 0.7) +
  stat_cor()
}
