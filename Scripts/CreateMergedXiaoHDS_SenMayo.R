###################################################
######## Analysis of Cell Fate Trajectories #######
###################################################

# Load packaged, data and gene sets

{
  
  library(Seurat)
  library(rtracklayer)
  library(sctransform)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(scater)
  library(scales)
  library(tidyr)
  library(RColorBrewer)
  library(ggvenn)
  
  # load functions to calculate HDS
  source('SharedFunctions.R')
  
  # hepatocyte damage associated genes
  dataPath <- 'hepatocyte-damage-score/Data/'
  
  countPath <- paste0(dataPath,
                      'Input/scRNAseq/GSE189600/MouseData/')
  
  outputPath <- 'hepatocyte-damage-score/Data/Output/'
  
  FileNames <- list.files(path = countPath)
  
  HDAG <- 
    read.csv(
      file = paste0(dataPath,
                    'Output/HDAG.csv'))
  
  genefilter <- 
    read.csv(
      file = paste0(dataPath,
                    'Output/GenesExpressedInHepatocytes/',
                    'GenesExpressedInHepatocytes.csv'))
  
  geneNoTInHDAG <- setdiff(genefilter$Genes, HDAG$gene_symbol)
  
  
  senMayoGenes <- read.csv(
    file = 
      paste0(dataPath,
             'Input/GeneSetsCellFates/SenMayoGeneSetMouse.csv'),
    sep = ';')
  
  # load meta data for all samples 
  
  metaDataAnnotations <- 
    read.table(paste0(dataPath,
                      'Input/scRNAseq/GSE189600/MouseData/',
                      'scitranslmed.adc9653_data_file_s1.csv'), 
               sep = ';', dec = ',', header = TRUE) %>%
    mutate(cellBarcode = gsub('.*_', replacement = '', metaDataAnnotations$X),
           condition =  gsub('_.*', replacement = '', metaDataAnnotations$X))
  
  
}

# Normal Chow 3 months
{
  # Read raw read counts 
  counts3moNC <- ReadMtx(mtx =  paste0(countPath, FileNames[3]), 
                         cells = paste0(countPath, FileNames[1]),
                         features = paste0(countPath, FileNames[2]),
                         feature.column = 2) 
  
  Seurat3moNC <- CreateSeuratObject(counts = counts3moNC, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = FileNames[1]),
                                    min.cells = 3) %>%
    AddMetaData(object = .,
                metadata = PercentageFeatureSet(., 
                                                pattern = '^mt-'),
                col.name = 'percent.mt' ) %>%
    subset(., 
           subset = 
             nFeature_RNA > 500 &
             nFeature_RNA < 8000 &
             percent.mt < 2)
  
  
  # Filter cells based on external annotations
  tempAnno <- 
    metaDataAnnotations[metaDataAnnotations$condition == 'NC3mo',]
  
  Seurat3moNC@meta.data$keepCell <- ifelse(
    rownames(Seurat3moNC@meta.data) %in% tempAnno$cellBarcode, "yes", "no")
  
  # Filter cells based on external annotations
  Seurat3moNC <-  subset( Seurat3moNC, subset = keepCell == "yes")
  
  # Sanity check and update metadata
  matched_cells <- tempAnno$cellBarcode %in% rownames(Seurat3moNC@meta.data)
  
  if (identical(tempAnno$cellBarcode[matched_cells], rownames(Seurat3moNC@meta.data))) {
    cat("Sanity check passed: Metadata aligned\n")
    Seurat3moNC <- AddMetaData(Seurat3moNC, tempAnno[matched_cells, ])
  } else {
    stop("Sanity check failed: Metadata mismatch")
  }
  
  # subset cells annotated as hepatocytes 
  
  Seurat3moNC <- subset(Seurat3moNC ,
                        subset = 
                          CellCluster ==  "PC-Hep" |
                          CellCluster == "mNASH-Hep1"  |
                          CellCluster == "mNASH-Hep2" |
                          CellCluster == "PP-Hep" |
                          CellCluster == "Int-Hep" ) %>% 
    SCTransform(.) %>% AddMetaData(object = .,
                                   metadata = unname(
                                     DS_calc.func(exprMatrices =
                                                    GetAssayData(
                                                      .,
                                                      assay = 'SCT',
                                                      layer = 'counts'),
                                                  DSignature = HDAG)),
                                   col.name = 'HDS' ) 
  
  Seurat3moNC <- RenameCells(Seurat3moNC, Seurat3moNC@meta.data$condition)
  
  Seurat3moNC@meta.data$condition <- '3m NC'
  
  remove(counts3moNC, tempAnno, matched_cells)
}

# Normal 9 months 
{
  
  counts9moNC <- ReadMtx(mtx =  paste0(countPath, FileNames[6]), 
                         cells = paste0(countPath, FileNames[4]),
                         features = paste0(countPath, FileNames[5]),
                         feature.column = 2) 
  
  Seurat9moNC <- CreateSeuratObject(counts = counts9moNC, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = FileNames[4]),
                                    min.cells = 3 ) %>%
    AddMetaData(object = .,
                metadata = PercentageFeatureSet(., 
                                                pattern = '^mt-'),
                col.name = 'percent.mt' ) %>%
    subset(., 
           subset = 
             nFeature_RNA > 500 &
             nFeature_RNA < 8000 &
             percent.mt < 2)
  
  # Filter cells based on external annotations
  tempAnno <- 
    metaDataAnnotations[metaDataAnnotations$condition == 'NC9mo',]
  
  Seurat9moNC@meta.data$keepCell <- ifelse(
    rownames(Seurat9moNC@meta.data) %in% tempAnno$cellBarcode, "yes", "no")
  
  # Filter cells based on external annotations
  Seurat9moNC <-  subset( Seurat9moNC, subset = keepCell == "yes")
  
  # Sanity check and update metadata
  matched_cells <- tempAnno$cellBarcode %in% rownames(Seurat9moNC@meta.data)
  
  if (identical(tempAnno$cellBarcode[matched_cells], rownames(Seurat9moNC@meta.data))) {
    cat("Sanity check passed: Metadata aligned\n")
    Seurat9moNC <- AddMetaData(Seurat9moNC, tempAnno[matched_cells, ])
  } else {
    stop("Sanity check failed: Metadata mismatch")
  }
  
  Seurat9moNC <- subset(Seurat9moNC ,
                        subset = 
                          CellCluster ==  "PC-Hep" |
                          CellCluster == "mNASH-Hep1"  |
                          CellCluster == "mNASH-Hep2" |
                          CellCluster == "PP-Hep" |
                          CellCluster == "Int-Hep" ) %>% 
    SCTransform(.) %>% AddMetaData(object = .,
                                   metadata = unname(
                                     DS_calc.func(exprMatrices =
                                                    GetAssayData(
                                                      .,
                                                      assay = 'SCT',
                                                      layer = 'counts'),
                                                  DSignature = HDAG)),
                                   col.name = 'HDS' ) 
  
  Seurat9moNC <- RenameCells(Seurat9moNC, 
                             Seurat9moNC@meta.data$condition)

  Seurat9moNC@meta.data$condition <- '9m NC'
  
  remove(counts9moNC, tempAnno, matched_cells)
  
  gc()
  
}

# 3 months NASH mouse 
{
  
  # 3 month NASH
  
  counts3moNASH <- ReadMtx(mtx =  paste0(countPath, FileNames[9]), 
                           cells = paste0(countPath, FileNames[7]),
                           features = paste0(countPath, FileNames[8]),
                           feature.column = 2) 
  
  Seurat3moNASH <- CreateSeuratObject(counts = counts3moNASH, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = FileNames[9]),
                                      min.cells = 3) %>%
    AddMetaData(object = .,
                metadata = PercentageFeatureSet(., 
                                                pattern = '^mt-'),
                col.name = 'percent.mt' ) %>%
    subset(., 
           subset = 
             nFeature_RNA > 500 &
             nFeature_RNA < 8000 &
             percent.mt < 2)
  
  # Filter cells based on external annotations
  tempAnno <- 
    metaDataAnnotations[metaDataAnnotations$condition == 'NASH3mo',]
  
  Seurat3moNASH@meta.data$keepCell <- ifelse(
    rownames(Seurat3moNASH@meta.data) %in% tempAnno$cellBarcode, "yes", "no")
  
  # Filter cells based on external annotations
  Seurat3moNASH <-  subset( Seurat3moNASH, subset = keepCell == "yes")
  
  # Sanity check and update metadata
  matched_cells <- tempAnno$cellBarcode %in% rownames(Seurat3moNASH@meta.data)
  
  if (identical(tempAnno$cellBarcode[matched_cells], rownames(Seurat3moNASH@meta.data))) {
    cat("Sanity check passed: Metadata aligned\n")
    Seurat3moNASH <- AddMetaData(Seurat3moNASH, tempAnno[matched_cells, ])
  } else {
    stop("Sanity check failed: Metadata mismatch")
  }
  
  Seurat3moNASH <- subset(Seurat3moNASH ,
                          subset =  
                            CellCluster ==  "PC-Hep" |
                            CellCluster == "mNASH-Hep1"  |
                            CellCluster == "mNASH-Hep2" |
                            CellCluster == "PP-Hep" |
                            CellCluster == "Int-Hep" ) %>% 
    SCTransform(.) %>% AddMetaData(object = .,
                                   metadata = unname(
                                     DS_calc.func(exprMatrices =
                                                    GetAssayData(
                                                      .,
                                                      assay = 'SCT',
                                                      layer = 'counts'),
                                                  DSignature = HDAG)),
                                   col.name = 'HDS' )
  
  Seurat3moNASH <- RenameCells(Seurat3moNASH, 
                               Seurat3moNASH@meta.data$condition)


  Seurat3moNASH@meta.data$condition <- '3m NASH'
  
  remove(counts3moNASH, tempAnno, matched_cells)
  gc()
  
}

# 9 months NASH mouse 
{
  counts9moNASH <- ReadMtx(mtx =  paste0(countPath, FileNames[12]), 
                           cells = paste0(countPath, FileNames[10]),
                           features = paste0(countPath, FileNames[11]),
                           feature.column = 2) 
  
  Seurat9moNASH <- CreateSeuratObject(counts = counts9moNASH, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = FileNames[12]),
                                      min.cells = 3) %>%
    AddMetaData(object = .,
                metadata = PercentageFeatureSet(., 
                                                pattern = '^mt-'),
                col.name = 'percent.mt' ) %>%
    subset(., 
           subset = 
             nFeature_RNA > 500 &
             nFeature_RNA < 8000 &
             percent.mt < 2)
  
  # Filter cells based on external annotations
  tempAnno <- 
    metaDataAnnotations[metaDataAnnotations$condition == 'NASH9mo',]
  
  Seurat9moNASH@meta.data$keepCell <- ifelse(
    rownames(Seurat9moNASH@meta.data) %in% tempAnno$cellBarcode, "yes", "no")
  
  # Filter cells based on external annotations
  Seurat9moNASH <-  subset( Seurat9moNASH, subset = keepCell == "yes")
  
  # Sanity check and update metadata
  matched_cells <- tempAnno$cellBarcode %in% rownames(Seurat9moNASH@meta.data)
  
  if (identical(tempAnno$cellBarcode[matched_cells], rownames(Seurat9moNASH@meta.data))) {
    cat("Sanity check passed: Metadata aligned\n")
    Seurat9moNASH <- AddMetaData(Seurat9moNASH, tempAnno[matched_cells, ])
  } else {
    stop("Sanity check failed: Metadata mismatch")
  }
  
  Seurat9moNASH <- subset(Seurat9moNASH ,
                          subset = 
                            CellCluster ==  "PC-Hep" |
                            CellCluster == "mNASH-Hep1"  |
                            CellCluster == "mNASH-Hep2" |
                            CellCluster == "PP-Hep" |
                            CellCluster == "Int-Hep" ) %>% 
    SCTransform(.) %>% AddMetaData(object = .,
                                   metadata = unname(
                                     DS_calc.func(exprMatrices =
                                                    GetAssayData(
                                                      .,
                                                      assay = 'SCT',
                                                      layer = 'counts'),
                                                  DSignature = HDAG)),
                                   col.name = 'HDS' )
  
  Seurat9moNASH <- RenameCells(Seurat9moNASH, Seurat9moNASH@meta.data$condition)
  
  Seurat9moNASH@meta.data$condition <- '9m NASH'
  
  remove(counts9moNASH, tempAnno, matched_cells)
  
  } 

# Create merged object 

{
  
  DefaultAssay(Seurat3moNC) <- "RNA"
  DefaultAssay(Seurat9moNC) <- "RNA"
  DefaultAssay(Seurat3moNASH) <- "RNA"
  DefaultAssay(Seurat9moNASH) <- "RNA"
  
  mergedObj <- 
    merge(x = Seurat3moNC,
          y = list(Seurat9moNC,Seurat3moNASH,Seurat9moNASH ), 
          add.cell.ids = c("3moNC", 
                           "9moNC", 
                           "3moNASH",
                           "9moNASH"), 
          project = "xiao")
  
  mergedObj$CellCluster <- factor(mergedObj$CellCluster,
                                     ordered = TRUE,
                                     levels = c("PP-Hep",
                                                "Int-Hep",
                                                "PC-Hep",
                                                "mNASH-Hep1",
                                                "mNASH-Hep2"))
  
  mergedObj$condition <- factor(mergedObj$condition,
                                      ordered = TRUE,
                                      levels = c("3m NC","9m NC", 
                                                 "3m NASH", "9m NASH" ))
  
  
  mergedObj <- SCTransform(mergedObj)
  mergedObj <- RunPCA(mergedObj)
  mergedObj <- RunUMAP(mergedObj, 
                       dims = 1:30)
  
  # rownames(mergedObj@meta.data) <- mergedObj@meta.data$X
  DimPlot(mergedObj, reduction = "umap", 
          group.by = c('CellCluster',"condition"))
  
  (FeaturePlot(mergedObj, 'HDS' ) | DimPlot(mergedObj, reduction = "umap", 
                                            group.by = "condition"))
  
}

# Calculate SenMayo activity score 
{
  mayoFeatures <- unlist(senMayoGenes$Gene.murine.)
  
  cells_rankingsAll <- AUCell_buildRankings(mergedObjExpr, 
                                            plotStats = TRUE)
  SenCells_All <- AUCell_calcAUC(mayoFeatures,
                                 cells_rankingsAll)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	20 (17% of 118)
  
  # Plot AUC histogram
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "AUCellExploreThresholds_SenMayoIndex_MergedXiaoSamples.pdf") )
  AUCell_exploreThresholds(SenCells_All, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  
  png( height = 6, width = 12, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "AUCellExploreThresholds_SenMayoIndex_MergedXiaoSamples.png") )
  AUCell_exploreThresholds(SenCells_All, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  
  if (identical(colnames(SenCells_All@assays@data$AUC), 
                rownames(mergedObj@meta.data)) == TRUE) {
    cat("Sanity check passed: Metadata aligned\n")
    mergedObj <- AddMetaData(mergedObj,
                             metadata = c(SenCells_All@assays@data$AUC),
                             col.name = 'AUCell_SenMayo')
  } else {
    stop("Sanity check failed: Metadata mismatch")
  }
  
  
  if(identical(colnames(SenCells_All@assays@data$AUC), 
               rownames(mergedObj@meta.data)) == TRUE){
    print('all good')
    mergedObj@meta.data$AUCell_SenMayo <- c(SenCells_All@assays@data$AUC)
  }
  

}

# Save merged object with HDS and SenMayo AUCell in RDS for further steps in 
# the analysis

{
  
saveRDS(mergedObj,
        file = paste0(outputPath,
                      "MergedXiaoHDS_SenMayoHepatocytes.rds"))
  
}



