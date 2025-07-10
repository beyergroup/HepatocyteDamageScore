
# I. Subset Liver Cell Atlas expression counts to find genes expressed 
# in hepatocytes: 
# - subset all cells annotated as hepatocytes
# - keep only nucSeq data
# - input is raw data from two experiments
# - output: merged filtered, subset and transformed data set 

library(Seurat)
library(sctransform)
library(ggplot2)
library(tidyverse)

# Input: single cell RNA couts from StSt Liver Cell Atlas experiment and from
# NAFLD Liver Cell Atlas experiment

# Standard Diet: 

pathData = 
   '/hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/rawData_mouseStSt/countTable_mouseStSt/'
pathAnnotation = 
  '/hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/annot_mouseStStAll.csv'

outputPath = '/hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/'

annotation.mouseStSt <- read.csv(
  file = pathAnnotation)

# read raw counts 
counts <- ReadMtx(mtx =  paste0(pathData,'matrix.mtx.gz'), 
                  cells = paste0(pathData,'barcodes.tsv.gz'),
                  features = paste0(pathData,'features.tsv.gz'),
                  feature.column = 1) 

# create Seurat object 
lca.mouse.hep.nucSEq.stst <- CreateSeuratObject(counts = counts, 
                                                project = "LCAStdDiet",
                                                min.cells = 3)

remove(counts)

# save indices to keep only single nuclei hepatocyte reads
cellsToKeep <- annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes' & 
                                      annotation.mouseStSt$typeSample == 'nucSeq', 
                                    'cell'] 
# subset seurat object
lca.mouse.hep.nucSEq.stst <- subset(lca.mouse.hep.nucSEq.stst, 
                                    cells = cellsToKeep)
print(cellsToKeep)

remove(cellsToKeep)

# add selected metadata to seurat object
annotation.mouseStSt <- 
  annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes' &
                         annotation.mouseStSt$typeSample == 'nucSeq', ]
annotation.mouseStSt <- 
  annotation.mouseStSt[is.element(annotation.mouseStSt$cell,
                                  names(lca.mouse.hep.nucSEq.stst$orig.ident)),]

lca.mouse.hep.nucSEq.stst@meta.data <- 
  lca.mouse.hep.nucSEq.stst@meta.data[annotation.mouseStSt$cell,]

if (identical(rownames(lca.mouse.hep.nucSEq.stst@meta.data),
              annotation.mouseStSt$cell) == TRUE){
  lca.mouse.hep.nucSEq.stst@meta.data <- 
    cbind(lca.mouse.hep.nucSEq.stst@meta.data,
          annotation.mouseStSt)
}

# repeat steps for NAFLD counts
pathData = 
  '/data/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/'
pathAnnotation = 
  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseNafldAll.csv'

annotation.mouseNAFLD <- read.csv(
  file = pathAnnotation)

# read raw counts 
counts <- ReadMtx(mtx =  paste0(pathData,'matrix.mtx.gz'), 
                  cells = paste0(pathData,'barcodes.tsv.gz'),
                  features = paste0(pathData,'features.tsv.gz'),
                  feature.column = 1) 

# create seurat object 
lca.mouse.hep.nucSEq.nafld <- CreateSeuratObject(counts = counts, 
                                                project = "LCANafld",
                                                min.cells = 3)

remove(counts)

# save indices to keep only single nuclei hepatocyte reads
cellsToKeep <- annotation.mouseNAFLD[annotation.mouseNAFLD$annot == 'Hepatocytes' &
                                       annotation.mouseNAFLD$typeSample == 'nucSeq', 
                                    'cell'] 
print(cellsToKeep)
# subset seurat object
lca.mouse.hep.nucSEq.nafld <- subset(lca.mouse.hep.nucSEq.nafld, 
                                    cells = cellsToKeep)
remove(cellsToKeep)

# add selected metadata to seurat object
annotation.mouseNAFLD <- annotation.mouseNAFLD[annotation.mouseNAFLD$annot == 'Hepatocytes' &
                                               annotation.mouseNAFLD$typeSample == 'nucSeq', ]
annotation.mouseNAFLD <- annotation.mouseNAFLD[is.element(annotation.mouseNAFLD$cell,
                                                        names(lca.mouse.hep.nucSEq.nafld$orig.ident)),]
lca.mouse.hep.nucSEq.nafld@meta.data <- lca.mouse.hep.nucSEq.nafld@meta.data[annotation.mouseNAFLD$cell,]

if (identical(rownames(lca.mouse.hep.nucSEq.nafld@meta.data),
              annotation.mouseNAFLD$cell) == TRUE){
  lca.mouse.hep.nucSEq.nafld@meta.data <- cbind(lca.mouse.hep.nucSEq.nafld@meta.data,
                                               annotation.mouseNAFLD)
}

# merge Sdt. Diet and NAFLD counts 

merged.mouse.lca.hep.nucSeq <-
  merge(x = lca.mouse.hep.nucSEq.stst,
        y = lca.mouse.hep.nucSEq.nafld, 
        add.cell.ids = c("StSt", "Nafld"), 
        project = "LiverCellAtlas")

# short QC and filtering


merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- 
  PercentageFeatureSet(merged.mouse.lca.hep.nucSeq,
                       pattern = "^mt-")

DimsBeforeFiltering <- dim(merged.mouse.lca.hep.nucSeq)

print(DimsBeforeFiltering)

# Plots pre Filtering
prefilterFeauturePlot <- VlnPlot(merged.mouse.lca.hep.nucSeq, 
                                 features = c("nFeature_RNA", 
                                              "nCount_RNA",
                                              "percent.mt" ), 
                                 ncol = 3 ,pt.size = FALSE) 

# Filtering
merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = nFeature_RNA > 500 &
                                        nFeature_RNA < 6000 &
                                        percent.mt <= 2.5 & 
                                        nCount_RNA < 30000 )

# post-filtering plots
prefilterFeauturePlot <- VlnPlot(merged.mouse.lca.hep.nucSeq, 
                                 features = c("nFeature_RNA", 
                                              "nCount_RNA",
                                              "percent.mt" ), 
                                 ncol = 3,pt.size = FALSE) 

# SCT transform counts

merged.mouse.lca.hep.nucSeq <- SCTransform(merged.mouse.lca.hep.nucSeq,
                                           verbose = FALSE)


# Save filtered filtered Seurat object
# Seurat object
saveRDS(merged.mouse.lca.hep.nucSeq, 
        paste0(outputPath,'countTable_mouseMerged_Hepatocytes_nucSeq.rds'))


ggsave(
  'prefilterFeauturePlot.png',
  plot = prefilterFeauturePlot,
  path = outputPath
)

ggsave(
  'prefilterFeauturePlot.png',
  plot = prefilterFeauturePlot,
  path = outputPath
)


