# Supplementary Figure 12 

# Panels: A, B -> Carlessi details 
{
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(patchwork)
  library(scales)
  
  # optional: load RDS with SCTransformed counts with HDS and SHGS calculated 
  # seurat_objects <- readRDS(paste0(path_to_data, "seurat_object_SCT_HDS_and_Senescence_scores.rds"))
  
  path_to_data <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/carlessi_snRNAseq_GSE200366/"
  setwd(path_to_data)
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  #SHGS 
  shgs_genes_mouse <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds")
  
  # load functions to calculate HDS
  source('SharedFunctions.R')
  
  seurat_object <- readRDS("Hep.rds")
  View(seurat_object@meta.data)
  
  DimPlot(seurat_object)
  
  seurat_object <- SCTransform(seurat_object)
  
  counts <- GetAssayData(seurat_object,
                         assay = 'SCT',
                         layer = 'counts')
  
  ## Calculate HDS
  HDS_calc <- DS_calc.func(exprMatrices = counts, DSignature = HDAG)
  
  if(identical(rownames(seurat_object@meta.data), names(HDS_calc)) == TRUE){
    seurat_object@meta.data$HDS <- unname(HDS_calc)
    print("Matching IDs!")
  }
  
  # Functions for ploting statistics
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  kruskal.test(HDS ~ Group, data = seurat_object@meta.data)
  # Kruskal-Wallis chi-squared = 8560.5, df = 2, p-value < 2.2e-16
  
  # separate zones
  # Zone_1 -> Periportal hepatocytes
  # Zone_2 -> intermediate
  # Zone_3 ->  Pericentral 
  seurat_object@meta.data <- seurat_object@meta.data %>% 
    mutate(zone = case_when(
      cell_type == "Zone_1_Hep" ~ "Periportal",
      cell_type == "Zone_2_Hep" ~ "Intermediate",
      cell_type == "Zone_3_Hep" ~ "Pericentral",
      cell_type == "daHep" ~ "damaged Hepatocytes" 
    ))
  
  seurat_object@meta.data$zone <- factor(
    seurat_object@meta.data$zone,
    levels = c("Periportal", "Intermediate","Pericentral", "damaged Hepatocytes"), 
    ordered = TRUE)
  

  # SHGS - Senescencent Hepatocytes Gene Signature 
  
  cells_rankings <- AUCell_buildRankings(counts)
  shgs_score <- AUCell::AUCell_calcAUC(shgs_genes_mouse, cells_rankings)
  
  
  if(identical(colnames(shgs_score@assays@data$AUC), rownames(seurat_object@meta.data)) == TRUE) {
    print("Matching IDs!")
    seurat_object@meta.data$aucell_shgs <- c(shgs_score@assays@data$AUC)
  }
  
  # PANEL A: Boxplots 
  
  ggplot(seurat_object@meta.data, aes(y = aucell_shgs, x = Group, color = Group)) + 
    geom_boxplot(lwd = 1.5) + scale_color_colorblind() + 
    stat_summary(geom = "linerange", fun.data = median_IQR, size = 1.5, show.legend = FALSE, position = position_dodge(0.95)) +
    stat_summary(fun.y = "median", geom = "crossbar", show.legend = FALSE, size = 1.5, width = 0.2,position = position_dodge(0.95)) +
    theme_base() +
    theme( legend.position = "bottom",
           text = element_text(size = 20),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("SHGS activity")
  
  # Wilcoxon test
  
  pairwise.wilcox.test(seurat_object@meta.data$aucell_shgs,
                       g = seurat_object@meta.data$Group, alternative = "greater",
                       p.adjust.method = "none")
  
  pairwise.wilcox.test(seurat_object@meta.data$HDS,
                       g = seurat_object@meta.data$Group, alternative = "greater", 
                       p.adjust.method = "none")
 
   
  # PANEL B: Scatter plots HDS vs. SHGS activity per group
  ggplot(data = seurat_object@meta.data,
         aes(x = HDS, y = aucell_shgs)) + facet_grid(~Group) +
    geom_point(size =1,alpha = 0.8, color = "grey") + 
    stat_smooth(span = 0.20, fill = "orange", color = "black", lwd = 1.5) + 
    theme_classic() + theme(axis.text = element_text(size = 12),text = element_text(size=12)) +
    ylab("SHGS activity")
 
}

# Panels: C, D,E -> Xiao mouse
{
  library(Seurat)
  library(readxl)
  library(rstatix)
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
  #library(scater)
  library(scales)
  
  pathData <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/Xiao_snRNAseq_GSE189600/"
  setwd(pathData)
  
  path_to_output <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/output/"
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  #SHGS 
  shgs_genes_mouse <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds")
  
  # load functions to calculate HDS
  source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')
  
  FileNames <- list.files(path = pathData)
  
  metaDataAnnotations <- read.table("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/annotation_files/Xiao_snRNAseq_GSE189600_scitranslmed.adc9653_data_file_s1.csv", 
                                    sep = ';', dec = ',', header = TRUE)
  
  metaDataAnnotations$cellBarcode <- 
    gsub('.*_', replacement = '', metaDataAnnotations$X)
  metaDataAnnotations$condition <- 
    gsub('_.*', replacement = '', metaDataAnnotations$X)
  
  
  # I. Read raw data & process & calculate all values: HDS, Pseudotime, etc
  
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
  
  
  # Create merged object:
  DefaultAssay(Seurat3moNC) <- "RNA"
  DefaultAssay(Seurat9moNC) <- "RNA"
  DefaultAssay(Seurat3moNASH) <- "RNA"
  DefaultAssay(Seurat9moNASH) <- "RNA"
  
  mergedObj <- merge(x = Seurat3moNC,
                     y = list(Seurat9moNC,Seurat3moNASH,Seurat9moNASH ), 
                     add.cell.ids = c("3moNC", 
                                      "9moNC", 
                                      "3moNASH",
                                      "9moNASH"), 
                     project = "xiao")
  
  mergedObj$CellType <- factor(mergedObj$CellType,
                               ordered = TRUE,
                               levels = c("PP-Hep",
                                          "Int-Hep",
                                          "PC-Hep",
                                          "mNASH-Hep1",
                                          "mNASH-Hep2"))
  
  
  mergedObj <- SCTransform(mergedObj)
  mergedObj <- RunPCA(mergedObj)
  mergedObj <- RunUMAP(mergedObj, 
                       dims = 1:30)
  
  # rownames(mergedObj@meta.data) <- mergedObj@meta.data$X
  DimPlot(mergedObj, reduction = "umap", 
          group.by = c('CellCluster',"condition"))
  
  (FeaturePlot(mergedObj, 'HDS' ) | DimPlot(mergedObj, reduction = "umap", group.by = "condition"))
  
  mergedObjExpr <- LayerData(mergedObj, assay = "SCT", layer = "counts")
  
  cells_rankingsAll <- AUCell_buildRankings(mergedObjExpr, plotStats = TRUE)
  shgs_aucell <- AUCell_calcAUC(shgs_genes_mouse, cells_rankingsAll)#
  #Genes in the gene sets NOT available in the dataset: 
  #	geneSet: 	4 (4% of 100)
  
  # Plot AUC histogram
  
  pdf( height = 6, width = 12,
       file = paste0(path_to_output,
                     "AUCellExploreThresholds_SHGS_MergedXiaoSamples.pdf") )
  AUCell_exploreThresholds(shgs_aucell,plotHist = TRUE,
                           assign = TRUE)
  dev.off()
  
  if (identical(colnames(shgs_aucell@assays@data$AUC),rownames(mergedObj@meta.data)) == TRUE) {
    print("Sanity check passed: Metadata aligned")
    mergedObj <- AddMetaData(mergedObj,
                             metadata = c(shgs_aucell@assays@data$AUC),col.name = 'AUCell_SHGS')
  } else {
    stop("Sanity check failed: Metadata mismatch")}
  
  mergedDataFrame <- mergedObj@meta.data
  
  
  mergedDataFrame$CellType <- factor(mergedDataFrame$CellType,
                                     ordered = TRUE,
                                     levels = c("PP-Hep",
                                                "Int-Hep",
                                                "PC-Hep",
                                                "mNASH-Hep1",
                                                "mNASH-Hep2"))
  
  mergedDataFrame$condition<- factor(mergedDataFrame$condition,
                                     ordered = TRUE,
                                     levels = c("3m NC","9m NC", 
                                                "3m NASH", "9m NASH" ))
  
  ggplot(mergedDataFrame , aes(x = AUCell_SHGS, color = condition)) + 
    geom_density(lwd = 1.5, trim = TRUE) + scale_color_colorblind() + theme_bw() + 
    theme( text = element_text(size = 18), legend.position = "top")  
  
  
  ### Sliding window approach 
  
  all_samples_scatterplot_glm <- ggplot(data = mergedDataFrame,aes(x = HDS, y = AUCell_SHGS)) + 
    geom_point(size =1,alpha = 0.8, color = "grey") + 
    stat_smooth(fill = "orange", color = "black", lwd = 1.5) + theme_classic() + 
    theme(axis.text = element_text(size = 18),text = element_text(size = 18)) +
    ylab("SHGS Activity") + xlim(c(min(mergedDataFrame$HDS),max(mergedDataFrame$HDS)))
  
  # Visualize with window approach 
  
  # Sliding window visualization approach  
  dataSW <- mergedDataFrame
  dataSW <- dataSW[order(dataSW$HDS),] 
  
  # Initialize variables
  window_size <- 2000
  medians <- numeric()  # Store medians
  q1s <- numeric()
  q3s <- numeric() # Store IQRs
  positions <- numeric() # Store the middle position of the window
  
  # Perform sliding window calculations
  for (i in 1:(nrow(dataSW) - window_size + 1)) {
    # Get the current window
    window <- dataSW[i:(i + window_size - 1), ]
    
    # Calculate statistics for SenMayo
    medians <- c(medians, median(window$AUCell_SHGS))
    q1s <- c(q1s, quantile(window$AUCell_SHGS, 0.25))
    q3s <- c(q3s, quantile(window$AUCell_SHGS, 0.75))
    
    # Store the middle position of the window for plotting
    positions <- c(positions, mean(window$HDS))
  }
  
  # Combine results into a data frame
  resultsSWsmoothMedian <- data.frame(
    Position = positions,
    Median = medians,
    Q1 = q1s,
    Q3 = q3s
  )
  
  resultsSWsmoothMedian$IQR <- (abs(resultsSWsmoothMedian$Q3)- 
                                  abs(resultsSWsmoothMedian$Q1))

  
  IQR_MeanHDS_all_samples <-  ggplot(resultsSWsmoothMedian, aes(x = Position, y = IQR)) + 
    geom_line(size = 1.5) + xlab("Mean HDS (2000 cells sliding window)") + 
    theme_classic() + theme(axis.text = element_text(size = 18),
                            text = element_text(size = 18)) + xlim(c(min(dataSW$HDS),max(dataSW$HDS)))
  
  plot_all_samples_glm <-cowplot::plot_grid(plotlist = list(all_samples_scatterplot_glm ,
                                                            IQR_MeanHDS_all_samples), rel_heights = c(0.75,0.25), nrow = 2) 
  
  
  plot_all_samples_glm
  
  # Function for ploting statistics
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  ggplot(mergedDataFrame, aes(y = AUCell_SHGS, x = condition, color = condition)) + 
    geom_boxplot(lwd = 1.5) + scale_color_colorblind() + 
    stat_summary(geom = "linerange", fun.data = median_IQR, size = 1.5, show.legend = FALSE, position = position_dodge(0.95)) +
    stat_summary(fun.y = "median", geom = "crossbar", show.legend = FALSE, size = 1.5, width = 0.2,position = position_dodge(0.95)) +
    theme_base() +
    theme( legend.position = "bottom",
           text = element_text(size = 20),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("SHGS activity")
  
  
  
  ## PLOTS 
  
  # 1. UMAPS: conditions, HDS, SHGS, and cellfate assignment
  
  # floor and ceiling values for applying limits to the color scale
  # Squish values outside the limit to the nearest boundary
  floor_valueHDS <- quantile(mergedObj@meta.data$HDS, 0.01)  
  ceiling_valueHDS <- quantile(mergedObj@meta.data$HDS, 0.99)
  
  floor_valueSHGS <- quantile(mergedObj@meta.data$AUCell_SHGS, 0.01)  
  ceiling_valueSHGS <- quantile(mergedObj@meta.data$AUCell_SHGS, 0.99)
  
  mergedObj@meta.data$condition <- factor(mergedObj@meta.data$condition,
                                          ordered = TRUE,
                                          levels = c("3m NC","9m NC", 
                                                     "3m NASH", "9m NASH" ))
  
  set.seed(42)
  umapHDSmerged <- FeaturePlot(mergedObj, features = "HDS",
                               cells = sample(rownames(mergedObj@meta.data),
                                              replace = FALSE), pt.size = 1) + 
    theme(legend.position = "top", legend.key.size = unit(15, "mm"),
          legend.key.height = unit(5,"mm"),
          text = element_text(size=20),
          axis.title=element_blank(),axis.text=element_blank(),
          axis.ticks=element_blank()) +
    scale_colour_distiller(palette = "YlOrBr", 
                           direction = 1,
                           limits = c(floor_valueHDS, 
                                      ceiling_valueHDS),  
                           oob = scales::squish) + ggtitle("")
  
  umapSHGSmerged <- FeaturePlot(mergedObj,
                                features = "AUCell_SHGS",
                                cells = sample(rownames(mergedObj@meta.data),
                                               replace = FALSE),
                                pt.size = 1) + 
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           limits = c(floor_valueSHGS, 
                                      ceiling_valueSHGS),  
                           oob = scales::squish) +
    theme(legend.position = "top", 
          legend.key.height=unit(5,"mm"),
          legend.key.width = unit(15, "mm"),
          text = element_text(size=20),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) + ggtitle("")
  
  
  densityPlotsHDS <- ggplot(data = mergedObj@meta.data,
                            aes(x = HDS, color = condition))+
    geom_density(lwd = 1, trim = TRUE)+ theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("HDS")
  
  densityPlotsSHGS <- ggplot(data = mergedObj@meta.data,
                             aes(x = AUCell_SHGS, color = condition))+
    geom_density(lwd = 1, trim = TRUE) + theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("SHGS")
  
  grid_umaps <- cowplot::plot_grid(plotlist = list(umapHDSmerged, 
                                                   umapSHGSmerged,
                                                   densityPlotsHDS,
                                                   densityPlotsSHGS), 
                                   rel_heights = c(1,0.35), nrow = 2) 
  
  set.seed(42)
  umap_conditions <- DimPlot(mergedObj, group.by = 'condition', 
                             pt.size = 1, cells = sample(rownames(mergedObj@meta.data))) + 
    scale_color_colorblind() + 
    theme(legend.position = "bottom", axis.text.y = element_blank(),
          axis.text.x = element_blank(),axis.text = element_blank(), 
          axis.title = element_blank(), axis.ticks = element_blank(), 
          plot.title = element_blank(), legend.key.height=unit(5,"mm"),
          text = element_text(size=12)) 
  
  grid_umaps

  
  
  umap_conditions

  
  umap_cellfate <- DimPlot(merged_data, 
                           group.by = 'cellfate', pt.size = 1, cells = 
                             sample(rownames(merged_data@meta.data), replace = FALSE)) + 
    scale_color_manual(values = c('undamaged' =  '#CCBB44',
                                  'transition' = '#228833',
                                  'damaged-not-senescent' = '#4477AA',
                                  'damaged-unresolved' = '#EE6677',
                                  'damaged-senescent' = '#AA3377'))+
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          legend.key.height=unit(5,"mm"),
          legend.text = element_text(size = 6))
  
  umap_cellfate

  
  
  
}

# Panels F and G: Xiao human 
{
  library(Seurat)
  library(stringr)
  library(GSEABase)
  library(dplyr)
  library(patchwork)
  library(AUCell)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(biomaRt)
  library(rtracklayer)
  library(tidyverse)
  library(gridExtra)
  library(readxl)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  
  path_to_data <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/human/Xiao_snRNAseq_GSE189600/"
  
  setwd(path_to_data)
  
  
  meta_data_annotations <- read.table(paste0(path_to_data,
                                             'humanSnRNAseqSeuratAnnotations.csv'), 
                                      sep = ';', dec = ',', header = TRUE) %>% 
    filter(., CellCluster %in% c("hPC-Hep","hPP-Hep","hNASH-Hep","hInt-Hep"))
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"Normal1_", "")
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"Normal2_", "")
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"Normal3_", "")
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"NASH1_", "")
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"NASH2_", "")
  meta_data_annotations$X <- str_replace(meta_data_annotations$X ,"NASH3_", "")
  
  original_identities_samples <- unique(meta_data_annotations$orig.ident)
  
  # Hepatocyte Damage Associated Genes and function for calculating HDS 
  humanHDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesManuallyModified",
    sep = '\t')
  # dummy values, because rank should as in list
  humanHDAG$mean_rank <- 1:42
  
  source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')
  
  sample_directories <- list.files(path_to_data) 
  sample_directories <- sample_directories[-c(1,2)]
  
  
  sample_list <- list()
  sample_i <- 1
  for(sample_i in 1:length(sample_directories)){
    
    counts <- ReadMtx(
      mtx = paste0(path_to_data,sample_directories[sample_i],"/","matrix.mtx.gz"),
      cells = paste0(path_to_data,sample_directories[sample_i],"/","barcodes.tsv.gz"),
      features = paste0(path_to_data,sample_directories[sample_i],"/","features.tsv.gz"),
      feature.column = 2)
    
    temp_seurat <- CreateSeuratObject(counts = counts, 
                                      project = sample_directories[sample_i],
                                      min.cells = 3)
    remove(counts)
    gc()
    
    # keeps id's of hepatocytes only 
    ids_to_keep <- meta_data_annotations$X[meta_data_annotations$orig.ident == original_identities_samples[[sample_i]]]
    temp_seurat <- subset(temp_seurat, cells = ids_to_keep)
    temp_seurat@meta.data$percent_mt <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")
    sample_list[[sample_i]] <-  subset(temp_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent_mt < 10)
    
    remove(temp_seurat)
    gc()
    sample_list[[sample_i]] <- SCTransform(sample_list[[sample_i]])
  }
  
  # calculate HDS 
  
  for(sample_i in 1:length(sample_list)){
    
    HDS_temp  <- DS_calc.func(
      exprMatrices = GetAssayData(sample_list[[sample_i]],
                                  assay = 'SCT',
                                  layer = 'counts'),
      DSignature = humanHDAG,
      geneIDname = 'HumanGeneID')
    
    if(identical(names(HDS_temp), 
                 rownames(sample_list[[sample_i]]@meta.data) ) == TRUE){
      sample_list[[sample_i]]@meta.data$HDS <- c(unlist(HDS_temp))
    }
    
  }
  
  gc()
  
  # Load SHGS:
  
  shgs_genes <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025.rds")
  shgs_genes <- shgs_genes$SHGS

  # Create a merged Seurat object 
  merged_data_xiao <- merge(
    x = sample_list[[1]],
    sample_list[2:6])
  merged_data_xiao <- SCTransform(merged_data_xiao)
  merged_data_xiao <- RunPCA(merged_data_xiao)
  merged_data_xiao <- RunUMAP(merged_data_xiao, dims = 1:30)
  
  # Compute geneset enrichment with AUCell 
  cells_rankings <- AUCell_buildRankings(LayerData(merged_data_xiao, 
                                                   assay = "SCT", layer = "counts"), plotStats = TRUE)
  shgs_aucell <- AUCell_calcAUC(shgs_genes, cells_rankings)
  
  # Add AUCell results to metadata 
  if (identical(colnames(shgs_aucell@assays@data$AUC),rownames(merged_data_xiao@meta.data)) == TRUE) {
    print("Sanity check passed: Metadata aligned")
    merged_data_xiao <- AddMetaData(merged_data_xiao,
                                    metadata = c(shgs_aucell@assays@data$AUC),col.name = 'AUCell_SHGS')
  } else {
    stop("Sanity check failed: Metadata mismatch")}
  
  
  # Prepare data to plot 
  data_to_plot <- merged_data_xiao@meta.data
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., patient_status= case_when(
      (orig.ident == "snRNA-seq_human_healthy_rep1"  | orig.ident == "snRNA-seq_human_healthy_rep2"| orig.ident == "snRNA-seq_human_healthy_rep3") ~ "healthy",
      (orig.ident == "snRNA-seq_human_NASH_rep1"  | orig.ident == "snRNA-seq_human_NASH_rep2"| orig.ident == "snRNA-seq_human_NASH_rep3") ~ "NASH"))
  
  data_to_plot$patient_status <- factor(data_to_plot$patient_status, levels = c("healthy", "NASH"))
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., patient = case_when(
      orig.ident == "snRNA-seq_human_healthy_rep1" ~ "p1, f:47",
      orig.ident == "snRNA-seq_human_healthy_rep2" ~ "p2, m:50",
      orig.ident == "snRNA-seq_human_healthy_rep3" ~ "p3, m:50",
      orig.ident == "snRNA-seq_human_NASH_rep1" ~ "p1, m:49",
      orig.ident == "snRNA-seq_human_NASH_rep2" ~ "p2, f:64", 
      orig.ident == "snRNA-seq_human_NASH_rep3" ~ "p3, m:31"))
  
  data_to_plot$patient <- factor(data_to_plot$patient, levels = c("p1, f:47", "p2, m:50", "p3, m:50","p3, m:31", "p1, m:49", "p2, f:64"))
  
  
  # Functions for ploting statistics
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  # Copy paste results to plot with violin plots
  x <- data_to_plot$AUCell_SHGS[data_to_plot$patient_status == "NASH"]
  y <- data_to_plot$AUCell_SHGS[data_to_plot$patient_status == "healthy"]
  wilcox.test(x = x, y = y, data = data_to_plot, p.adjust.method = "BH", alternative = "greater")
  
  #data:  x and y
  #W = 251939390, p-value < 2.2e-16
  # alternative hypothesis: true location shift is greater than 0
  
  ggplot(data_to_plot, aes(y = AUCell_SHGS, 
                           x = patient_status, colour = patient_status)) +
    geom_violin(trim = TRUE, lwd = 1.5) + 
    scale_color_colorblind() +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR, 
                 size = 1.5 ,
                 show.legend = FALSE, 
                 position = position_dodge(0.95)) +
    stat_summary(
      fun.y = "median",
      size = 1,
      geom = "crossbar",
      show.legend = FALSE,
      width = 0.2,
      position = position_dodge(0.95)) +
    theme_classic() +
    theme( text = element_text(size = 18) ,
           axis.title.x = element_blank()) +
    annotate("text",label= "One-Sided Wilcoxon-Test,\n B.H. adjusted p-value < 2.2e-16 \n n = 3", y = 0.1, x = 1)
  
  # SHGS vs HDS 
  
  xiao_scatterplot_glm <- ggplot(data = merged_data_xiao@meta.data, aes(x = HDS, y = AUCell_SHGS)) + 
    geom_point(size =1,alpha = 0.8, color = "grey") + 
    stat_smooth(fill = "orange", color = "black", lwd = 1.5) + theme_classic() + 
    theme(axis.text = element_text(size = 18),text = element_text(size = 18)) 
  
  # Sliding window visualization approach  
  dataSW <- merged_data_xiao@meta.data
  dataSW <- dataSW[order(dataSW$HDS),] 
  
  # Initialize variables
  window_size <- 2000
  medians <- numeric()  # Store medians
  q1s <- numeric()
  q3s <- numeric() # Store IQRs
  positions <- numeric() # Store the middle position of the window
  
  # Perform sliding window calculations
  for (i in 1:(nrow(dataSW) - window_size + 1)) {
    # Get the current window
    window <- dataSW[i:(i + window_size - 1), ]
    
    # Calculate statistics for SenMayo
    medians <- c(medians, median(window$AUCell_SHGS))
    q1s <- c(q1s, quantile(window$AUCell_SHGS, 0.25))
    q3s <- c(q3s, quantile(window$AUCell_SHGS, 0.75))
    
    # Store the middle position of the window for plotting
    positions <- c(positions, mean(window$HDS))
  }
  
  # Combine results into a data frame
  resultsSWsmoothMedian <- data.frame(
    Position = positions,
    Median = medians,
    Q1 = q1s,
    Q3 = q3s
  )
  
  resultsSWsmoothMedian$IQR <- (abs(resultsSWsmoothMedian$Q3)- 
                                  abs(resultsSWsmoothMedian$Q1))

  
  
  IQR_MeanHDS_Xiao <-  ggplot(resultsSWsmoothMedian, aes(x = Position, y = IQR)) + 
    geom_line(size = 1.5) + xlab("Mean HDS (2000 cells sliding window)") + 
    theme_classic() + theme(axis.text = element_text(size = 18),
                            text = element_text(size = 18)) + xlim(c(min(dataSW$HDS),max(dataSW$HDS)))
  
  plot_IQR_glm_Xiao <-cowplot::plot_grid(plotlist = list(xiao_scatterplot_glm,
                                                         IQR_MeanHDS_Xiao), rel_heights = c(0.75,0.25), nrow = 2)
  
  
  plot_IQR_glm_Xiao
  
}

# Panels H and I: ORA differentially expressed genes damaged-senescent hepatocytes
{
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(patchwork)
  library(scales)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(Orthology.eg.db)
  
  # Run this after running chunk for PANELS A and B
  # Compute thresholds for cell-fate assignment 
  
  quantiles_hds <- quantile(seurat_object@meta.data$HDS, probs = c(0.25, 0.75))
  quantiles_shgs <- quantile(seurat_object@meta.data$aucell_shgs[seurat_object@meta.data$Group == "TAA"], probs = c(0.25, 0.75))

  
  # Add or overwrite the 'group' column based on the conditions
  seurat_object@meta.data <-  seurat_object@meta.data %>%
    mutate(cellfate = case_when(
      HDS <= quantiles_hds[1]  ~ "undamaged",
      HDS > quantiles_hds[1] & HDS <= quantiles_hds[2] ~ "transition",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[2]  ~ "damaged_senescent",
      HDS > quantiles_hds[2] & aucell_shgs <= quantiles_shgs[1] ~ "damaged_not_senescent",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[1] & aucell_shgs <= quantiles_shgs[2] ~ "damaged_unresolved",
      TRUE ~ NA_character_  # optional: assign NA if no condition is met
    ))
  
  # Control number of cells assigned per cell fate
  table( seurat_object@meta.data$cellfate)
  # Control if any cells have not been assigned a cell fate 
  sum(is.na( seurat_object@meta.data$cellfate))
  
  # Cell fates as ordered factors 
  seurat_object@meta.data$cellfate <- factor( seurat_object@meta.data$cellfate, 
                                              ordered = TRUE,
                                              levels = c('undamaged', 
                                                         'transition',
                                                         'damaged_not_senescent',
                                                         'damaged_unresolved',
                                                         'damaged_senescent'),
                                              labels = c('undamaged',
                                                         'transition',
                                                         'damaged-not-senescent',
                                                         'damaged-unresolved',
                                                         'damaged-senescent'))
  
  
  # Summarize the data: count the number of cells for each fate in each sample
  summarized_cellfates <-  seurat_object@meta.data %>%
    group_by(Group, cellfate) %>%
    summarize(Cell_Count = n(), .groups = "drop") %>%
    group_by(Group) %>% 
    mutate(proportion = Cell_Count / sum(Cell_Count))
  
  
  # Prepare a table-like data frame for annotations
  label_data <- summarized_cellfates %>% 
    group_by(Group) %>%
    mutate(n = sum(Cell_Count)) %>%
    dplyr::select(Group, cellfate, Cell_Count,n) %>%
    pivot_wider(names_from = cellfate, values_from = Cell_Count, values_fill = 0)
  
  # Convert the data frame to a long format for plotting as a table
  label_data_long <- label_data[,-2] %>%
    pivot_longer(
      cols = -Group,
      names_to = "cellfate",
      values_to = "Cell_Count"
    ) %>% 
    mutate(cellfate = factor(cellfate,
                             levels = c('undamaged',
                                        'transition',
                                        'damaged-not-senescent',
                                        'damaged-unresolved',
                                        'damaged-senescent')))
  
  
  seurat_object <- PrepSCTFindMarkers( seurat_object)
  Idents( seurat_object) <- "cellfate"
  
  # p-values obtained from this analysis should be interpreted with caution, 
  # because these tests treat each cell as an independent replicate 
  # and ignore inherent correlations between cells originating from the same sample
  # there can be many false positives!
  
  
  senescence_damage_markers <- FindMarkers( seurat_object, ident.1 =  'damaged-senescent',
                                            ident.2 = 'damaged-not-senescent')
  
  
  
  
  positive_markers <- filter(senescence_damage_markers,
                             p_val_adj < 0.05 & avg_log2FC > 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  
  positive_markers <- positive_markers[order(positive_markers$avg_log2FC, decreasing = TRUE),]
  positive_markers$gene <- rownames(positive_markers)
  
  intersect(positive_markers$gene, shgs_genes_mouse)
  # 30 genes overlap:
  # [1] "Fosl1"     "Rab31"     "Abhd2"     "Osbpl3"    "Sema3e"    "Iqgap1"    "Cd9"       "Pawr"      "Enc1"      "Tnfrsf12a"
  #[11] "Gls"       "Igdcc4"    "Sytl5"     "Bmp8b"     "Srxn1"     "Krt18"     "Tgfbr2"    "Itgav"     "Tpm1"      "Slc48a1"  
  #[21] "Tanc1"     "Synj2"     "Ppfibp1"   "Cadm1"     "Sntb2"     "Pam"       "Epha2"     "Ifngr1"    "Palmd"     "Angptl4"  
  
  intersect(positive_markers$gene, HDAG$gene_symbol[1:42])
  #[1] "Iqgap1" "Srxn1"  "Rtn4"   "Nrg1"   "Ano6
  
  # subset significant negative marker genes
  negative_markers <- filter(senescence_damage_markers, p_val_adj < 0.05 & avg_log2FC < 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  
  negative_markers <- negative_markers[order(negative_markers$avg_log2FC, decreasing = FALSE),]
  negative_markers$gene <- rownames(negative_markers)
  
  intersect(negative_markers$gene, shgs_genes_mouse)
  # none
  intersect(negative_markers$gene, HDAG$gene_symbol[1:42])
  # [1] "Ces3b"     "Selenbp2"  "Gnmt"      "Nudt7"     "Slc22a30"  "Serpina1c" "Serpinc1" 
  
  egoNegBP <- enrichGO(gene = negative_markers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       universe = rownames(merged_data))
  
  egoNegALL <- enrichGO(gene = negative_markers$gene,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
  
  egoNegBP_simplified <- simplify(
    egoNegBP,
    cutoff = 0.7,             # Similarity threshold (0.7 is a common default)
    by = "p.adjust",          # Retain the most significant term
    select_fun = min,         # Select term with the smallest p-value within each cluster
    measure = "Wang"          # Semantic similarity measure ("Wang" is common for GO terms)
  )
  # plot BP GO terms
  mutate(egoNegBP, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 40, label_format = 60, font.size = 15) 
  
  # plot simplified version of enriched BP GO 
  mutate(egoNegBP_simplified, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 30, label_format = 60, font.size = 15) 
 
  write.csv(egoNegBP@result[egoNegBP@result$p.adjust < 0.01,], 
            file = paste0(path_to_output, "carlessi_ds_vs_dns_negative_markers_GO_BP_significant_padjust_0_01.csv") )
  




  egoPosBP <- enrichGO(gene = positive_markers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       universe = rownames(merged_data))
  
  
  egoPosBP_simplified <- simplify(
    egoPosBP,
    cutoff = 0.7,             # Similarity threshold (0.7 is a common default)
    by = "p.adjust",          # Retain the most significant term
    select_fun = min,         # Select term with the smallest p-value within each cluster
    measure = "Wang"          # Semantic similarity measure ("Wang" is common for GO terms)
  )
  
  # plot BP GO terms
  mutate(egoPosBP, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 40, label_format = 60, font.size = 15) 
  
  # plot simplified version of enriched BP GO 
  mutate(egoPosBP_simplified , qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 30, label_format = 60, font.size = 15) 
 

  
  write.csv(egoPosBP@result[egoPosBP@result$p.adjust < 0.01,], 
            file = paste0(path_to_output, "carlessi_ds_vs_dns_positive_markers_GO_BP_significant_padjust_0_01.csv") )
  


  
}

# Panel J and K: Wunderlich mouse
{
  # Applying HDS to Thomas scRNAseq data
  library(Seurat)
  library(AUCell)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  
  path_to_data <- "Data/mouse/wunderlich_NIK_NASH_CDA/"
  data <- readRDS(paste0(path_to_data,"nik_samples_annotated.rds"))
  data <- subset(data, subset = class_annotation == "Hepatocyte")
  data <- subset(data, subset = (Diet == "Nash" | Diet == "Cdaa") & Genetic == "NIK-FL")
  metadata_hepatocytes <- data@meta.data
  data <- SCTransform(data)
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  #SHGS 
  shgs_genes_mouse <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds")
  
  # load functions to calculate HDS
  source('SharedFunctions.R')
  
  
  HDS <- DS_calc.func( exprMatrices = GetAssayData(data,assay = 'SCT', layer = 'counts'), 
                       DSignature = HDAG )
  
  if(identical(rownames(data@meta.data), names(HDS)) == TRUE){
    data@meta.data$HDS <- unname(HDS)
    print("Matching IDs!")}
  
  
  # SHGS - Senescencent Hepatocytes Gene Signature 
  cells_rankings <- AUCell_buildRankings(GetAssayData(data,assay = 'SCT', layer = 'counts'))
  shgs_score <- AUCell::AUCell_calcAUC(shgs_genes_mouse, cells_rankings)
  
  png( height = 6, width = 12, units = 'in',res = 300,
       file = paste0(path_to_data,"shgs_aucell_histogram.png"))
  AUCell_exploreThresholds(shgs_score, plotHist = TRUE, assign = TRUE)
  dev.off()
  dev.off()
  
  if(identical(colnames(shgs_score@assays@data$AUC), rownames(data@meta.data)) == TRUE) {
    print("Matching IDs!")
    data@meta.data$aucell_shgs <- c(shgs_score@assays@data$AUC)}
  
  data@meta.data$Diet <- factor(data@meta.data$Diet,
                                ordered = TRUE,
                                levels = c("Nash", "Cdaa"))
  
  # Ploting 
  
  # Function for ploting statistics
  ggplot(data@meta.data, aes(y = HDS, x = Diet, color = Diet)) + 
    geom_boxplot(lwd = 1.5) + 
    scale_color_manual(values = c("Nash" = "#56B4E9", "Cdaa" = "#E69F00")) +
    theme_base() +
    theme( legend.position = "bottom",
           text = element_text(size = 20),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("HDS")
  
  ggplot(data@meta.data, aes(y = aucell_shgs, x = Diet, color = Diet)) + 
    geom_boxplot(lwd = 1.5) +
    scale_color_manual(values = c("Nash" = "#56B4E9", "Cdaa" = "#E69F00")) +
    theme_base() +
    theme( legend.position = "bottom",
           text = element_text(size = 20),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("SHGS activity")
  
}

