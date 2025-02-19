# HDS of Human hepatocytes from GSE189600 (Xiao et. al, 2023)
# - input is raw data, human version of HDAG list, our functions
# - subset, then normalized counts for each sample and the merged all 
# normalized counts
# - then added annotation from authors to metadata and extracted counts from 
# cells annotated as hepatocytes (four types) 
# - then I calculated the weighted HDS of those hepatocytes and made vioin
# plots of the HDS distribution per condition and hepatocyte type


library(Seurat)
library(rtracklayer)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)


pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/HumanData/'

outputPath = 'hepatocyte-damage-score/Data/Output/Results/'

samplesNames <- c('GSM6808755_Y21822A1', 'GSM6808756_Y31122A5',
                  'GSM6808757_Y31122A8', 'GSM6808758_Y21822A2',
                  'GSM6808759_Y31122A6', 'GSM6808760_Y31122A7')

# read raw counts 
countsA1 <- ReadMtx(mtx =  paste0(pathData, samplesNames[1],'_matrix.mtx'), 
                  cells = paste0(pathData, samplesNames[1],'_barcodes.tsv'),
                  features = paste0(pathData, samplesNames[1],'_features.tsv'),
                  feature.column = 2) 

A1Seurat <- CreateSeuratObject(counts = countsA1, 
                               project = samplesNames[1],
                               min.cells = 3)

A1Seurat[["percent.mt"]] <- PercentageFeatureSet(A1Seurat, 
                                                 pattern = "^MT-")


A1Seurat <- subset(A1Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A1Seurat <- SCTransform(A1Seurat)


remove(countsA1)
gc()

###

countsA2 <- ReadMtx(mtx =  paste0(pathData, samplesNames[2],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[2],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[2],'_features.tsv'),
                    feature.column = 2) 

A2Seurat <- CreateSeuratObject(counts = countsA2, 
                               project = samplesNames[2],
                               min.cells = 3)

A2Seurat[["percent.mt"]] <- PercentageFeatureSet(A2Seurat, 
                                                 pattern = "^MT-")
A2Seurat <- subset(A2Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A2Seurat <- SCTransform(A2Seurat)

remove(countsA2)
gc()


## 

countsA3 <- ReadMtx(mtx =  paste0(pathData, samplesNames[3],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[3],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[3],'_features.tsv'),
                    feature.column = 2) 

A3Seurat <- CreateSeuratObject(counts = countsA3, 
                               project = samplesNames[3],
                               min.cells = 3)

A3Seurat[["percent.mt"]] <- PercentageFeatureSet(A3Seurat, 
                                                 pattern = "^MT-")

A3Seurat <- subset(A3Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A3Seurat <- SCTransform(A3Seurat)

remove(countsA3)
gc()

##

countsA4 <- ReadMtx(mtx =  paste0(pathData, samplesNames[4],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[4],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[4],'_features.tsv'),
                    feature.column = 2) 

A4Seurat <- CreateSeuratObject(counts = countsA4, 
                               project = samplesNames[4],
                               min.cells = 3)

A4Seurat[["percent.mt"]] <- PercentageFeatureSet(A4Seurat, 
                                                 pattern = "^MT-")

A4Seurat <- subset(A4Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A4Seurat <- SCTransform(A4Seurat)

remove(countsA4)
gc()
##
countsA5 <- ReadMtx(mtx =  paste0(pathData, samplesNames[5],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[5],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[5],'_features.tsv'),
                    feature.column = 2) 

A5Seurat <- CreateSeuratObject(counts = countsA5, 
                               project = samplesNames[5],
                               min.cells = 3)

A5Seurat[["percent.mt"]] <- PercentageFeatureSet(A5Seurat, 
                                                 pattern = "^MT-")

A5Seurat <- subset(A5Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A5Seurat <- SCTransform(A5Seurat)


remove(countsA5)
gc()

countsA6 <- ReadMtx(mtx =  paste0(pathData, samplesNames[6],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[6],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[6],'_features.tsv'),
                    feature.column = 2) 

A6Seurat <- CreateSeuratObject(counts = countsA6, 
                               project = samplesNames[6],
                               min.cells = 3)

A6Seurat[["percent.mt"]] <- PercentageFeatureSet(A6Seurat, 
                                                 pattern = "^MT-")

A6Seurat <- subset(A6Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)


remove(countsA6)
gc()

# save before because memory exhausted
allSamplesBut6 <- list(A1Seurat, A2Seurat, A3Seurat,
                   A4Seurat, A5Seurat)

saveRDS(object = allSamplesBut6, 
        file = 'hepatocyte-damage-score/Data/Output/tempFile.rds')
remove(allSamplesBut6)
remove(A1Seurat, A2Seurat, A3Seurat, A4Seurat, A5Seurat)

# now memory free and it worked 
A6Seurat <- SCTransform(A6Seurat)

allSamplesBut6 <- readRDS(file = 'hepatocyte-damage-score/Data/Output/tempFile.rds')


# use meta data to subset hepatocytes 

metaDataAnnotations <- read.table(
  paste0(pathData,
         'scitranslmed.adc9653_data_files_s1_to_s6/humanSnRNAseqSeuratAnnotations.csv'), 
  sep = ';', dec = ',', header = TRUE)

# A1
allSamplesBut6[[1]]@meta.data$cell_id <- rownames(allSamplesBut6[[1]]@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "Normal.A1" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[1]]@meta.data)),]

allSamplesBut6[[1]]@meta.data$cellBarcode <- NA
allSamplesBut6[[1]]@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[1]]@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
allSamplesBut6[[1]]@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[1]]@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
allSamplesBut6[[1]]@meta.data$deleteNA <- 
  rep("not na", length(allSamplesBut6[[1]]@meta.data$CellType))
allSamplesBut6[[1]]@meta.data$deleteNA[is.na(allSamplesBut6[[1]]@meta.data$CellType)] <-
  "na"

allSamplesBut6[[1]]@meta.data <- subset(allSamplesBut6[[1]]@meta.data,
                                  subset = deleteNA == 'not na')

# A2
allSamplesBut6[[2]]@meta.data$cell_id <- rownames(allSamplesBut6[[2]]@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "Normal.A5" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[2]]@meta.data)),]

allSamplesBut6[[2]]@meta.data$cellBarcode <- NA
allSamplesBut6[[2]]@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[2]]@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
allSamplesBut6[[2]]@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[2]]@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
allSamplesBut6[[2]]@meta.data$deleteNA <- 
  rep("not na", length(allSamplesBut6[[2]]@meta.data$CellType))
allSamplesBut6[[2]]@meta.data$deleteNA[is.na(allSamplesBut6[[2]]@meta.data$CellType)] <-
  "na"

allSamplesBut6[[2]]@meta.data <- subset(allSamplesBut6[[2]]@meta.data,
                                        subset = deleteNA == 'not na')

# A3
allSamplesBut6[[3]]@meta.data$cell_id <- rownames(allSamplesBut6[[3]]@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "Normal.A8" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[3]]@meta.data)),]

allSamplesBut6[[3]]@meta.data$cellBarcode <- NA
allSamplesBut6[[3]]@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[3]]@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
allSamplesBut6[[3]]@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[3]]@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
allSamplesBut6[[3]]@meta.data$deleteNA <- 
  rep("not na", length(allSamplesBut6[[3]]@meta.data$CellType))
allSamplesBut6[[3]]@meta.data$deleteNA[is.na(allSamplesBut6[[3]]@meta.data$CellType)] <-
  "na"

allSamplesBut6[[3]]@meta.data <- subset(allSamplesBut6[[3]]@meta.data,
                                        subset = deleteNA == 'not na')

# A4

allSamplesBut6[[4]]@meta.data$cell_id <- rownames(allSamplesBut6[[4]]@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "NASH.A2" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[4]]@meta.data)),]

allSamplesBut6[[4]]@meta.data$cellBarcode <- NA
allSamplesBut6[[4]]@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[4]]@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
allSamplesBut6[[4]]@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[4]]@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
allSamplesBut6[[4]]@meta.data$deleteNA <- 
  rep("not na", length(allSamplesBut6[[4]]@meta.data$CellType))
allSamplesBut6[[4]]@meta.data$deleteNA[is.na(allSamplesBut6[[4]]@meta.data$CellType)] <-
  "na"

allSamplesBut6[[4]]@meta.data <- subset(allSamplesBut6[[4]]@meta.data,
                                        subset = deleteNA == 'not na')

#A5
allSamplesBut6[[5]]@meta.data$cell_id <- rownames(allSamplesBut6[[5]]@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "NASH.A6" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[5]]@meta.data)),]

allSamplesBut6[[5]]@meta.data$cellBarcode <- NA
allSamplesBut6[[5]]@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[5]]@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
allSamplesBut6[[5]]@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(allSamplesBut6[[5]]@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
allSamplesBut6[[5]]@meta.data$deleteNA <- 
  rep("not na", length(allSamplesBut6[[5]]@meta.data$CellType))
allSamplesBut6[[5]]@meta.data$deleteNA[is.na(allSamplesBut6[[5]]@meta.data$CellType)] <-
  "na"

allSamplesBut6[[5]]@meta.data <- subset(allSamplesBut6[[5]]@meta.data,
                                        subset = deleteNA == 'not na')

#A6
A6Seurat@meta.data$cell_id <- rownames(A6Seurat@meta.data)

temp <- metaDataAnnotations[metaDataAnnotations$orig.ident == "NASH.A7" ,]
rownames(temp) <- gsub(".*_", '', x = temp$X)
temp <- temp[intersect(gsub(".*_", '', x = temp$X), rownames(A6Seurat@meta.data)),]

A6Seurat@meta.data$cellBarcode <- NA
A6Seurat@meta.data$cellBarcode[match(gsub(".*_", '', x = temp$X), rownames(A6Seurat@meta.data))] <- 
  gsub(".*_", '', x = temp$X)
A6Seurat@meta.data$CellType[match(gsub(".*_", '', x = temp$X), rownames(A6Seurat@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
A6Seurat@meta.data$deleteNA <- 
  rep("not na", length(A6Seurat@meta.data$CellType))
A6Seurat@meta.data$deleteNA[is.na(A6Seurat@meta.data$CellType)] <-
  "na"

A6Seurat@meta.data <- subset(A6Seurat@meta.data,
                             subset = deleteNA == 'not na')

# subset hepatocytes

for(i in 1:5){
  allSamplesBut6[[i]] <- subset(x = allSamplesBut6[[i]] , 
                                subset = CellType ==  "hPC-Hep" |
                                  CellType == "hNASH-Hep"  |
                                  CellType == "hPP-Hep" |
                                  CellType == "hInt-Hep" )
}

A6Seurat <- subset(x = A6Seurat ,
                   subset = CellType ==  "hPC-Hep" |
                     CellType == "hNASH-Hep"  |
                     CellType == "hPP-Hep" |
                     CellType == "hInt-Hep" )






#If you want to merge the normalized data matrices as well as the raw count matrices
# pass merge.data = TRUE

mergedHepatocytes <- merge(allSamplesBut6[[1]], y = list(allSamplesBut6[[2]],
                                                     allSamplesBut6[[3]],
                                                     allSamplesBut6[[4]],
                                                     allSamplesBut6[[5]],
                                                     A6Seurat), 
                       add.cell.ids = c("Normal1", "Normal2", "Normal3",
                                        "NASH1", "NASH2", "NASH3"), 
                       project = "HumanXiaoEtAl",
                       merge.data = TRUE)



remove(allSamplesBut6, A6Seurat)

saveRDS(object = mergedHepatocytes, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')

mergedHepatocytes <- readRDS(file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')

# calculate HDS with HDAG list mapped to human genes

humanHDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesOnlyHumanGenesManual.csv',
    sep = ';')

# load functions to calculate HDS
source('SharedFunctions.R')

# create Expression Matrix

ExprMatrix <- GetAssayData(mergedHepatocytes)
remove(mergedHepatocytes)

# dummy values, because rank should as in list
humanHDAG$mean_rank <- 1:42

l <- dim(ExprMatrix)[2]

HDS <- DS_calc.func(ExprMatrix[,1:floor((l/2))], DSignature = humanHDAG,
                    geneIDname = 'HumanGeneID'
                    
                    )
gc()

HDS <- c(HDS, DS_calc.func(ExprMatrix[,(floor(l/2)+1):l], DSignature = humanHDAG,
                                geneIDname = 'HumanGeneID'
                                ) 
              )

length(HDS) == l 

# Prepare HDS results for plotting

mergedHepatocytes <- 
  readRDS(file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')


DataFrameHDS <- data.frame('HDS' = unname(HDS),
                           'CellID' = names(unlist(HDS)))

identical(DataFrameHDS$CellID, rownames(mergedHepatocytes@meta.data))

mergedHepatocytes@meta.data$sample <- 
  gsub('_.*', replacement = '', rownames(mergedHepatocytes@meta.data))

DataFrameHDS$sample <- mergedHepatocytes@meta.data$sample

DataFrameHDS$condition <-gsub('[1-3]_.*', 
                              replacement = '', 
                              rownames(mergedHepatocytes@meta.data))


DataFrameHDS$HepatocyteType <- mergedHepatocytes@meta.data$CellType


DataFrameHDS$condition <- factor(DataFrameHDS$condition,
                           ordered = TRUE,
                           levels = c('Normal', 'NASH'))

DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                 ordered = TRUE,
                                 levels = c('hPP-Hep', 'hInt-Hep',
                                            'hPC-Hep','hNASH-Hep'))



# Save HDS and cell meta data

saveRDS(DataFrameHDS, file = paste0(outputPath, 'Xiao2023_Human_HDS_Dataframe.rds'))
test <- readRDS(paste0(outputPath, 'Xiao2023_Human_HDS_Dataframe.rds'))
# Ploting 


PlotHDSConditions <- 
  ggplot( DataFrameHDS,
          aes( x = condition,
               y = HDS,
               color = condition
               )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=25),
         legend.position = 'right',
         axis.text.x = element_text(angle = 0, 
                                    vjust = 0.5, hjust=0.5
                                    )) +
   guides(colour = guide_legend(override.aes = aes(label = ""),
                                title="Diet"))+
   labs(y = "HDS", x = '') + 
   stat_summary(fun.data = "mean_sdl",
                fun.args = list(mult = 1),
                geom = "pointrange") +
  stat_compare_means(paired = FALSE, 
                     method = 'wilcox.test',
                     label =  "p.signif",
                     # label.y = 0.12,
                     label.x = 1.25,
                     size = 10,
                     symnum.args = list(
                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                       symbols = c("****", "***", "**", "*", "ns")))

PlotHDSConditions
ggsave(file = 'HumanGSE189600WeightedHDSPerCondition.png',
       path = outputPath)

   PlotHDSHepTypes <- 
     ggplot( DataFrameHDS,
             aes( x = condition,
                  y = HDS,
                  color = condition
             )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=24),
            legend.position = 'right',
            axis.text.x = element_text(angle = 45, 
                                       vjust = 0.5, hjust=0.5
            )) + 
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                title="Health Status")) +
     labs(y = "HDS", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange")  + 
     stat_compare_means(ref.group = 'Normal',
                        paired = FALSE, 
                        method = 'wilcox.test',
                        label =  "p.signif",
                        # label.y = 0.12,
                        label.x = 1.25,
                        size = 10,
                        symnum.args = list(
                          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                          symbols = c("****", "***", "**", "*", "ns"))) +
     facet_wrap(~HepatocyteType)
   
   PlotHDSHepTypes
   ggsave(file = 'HumanGSE189600WeightedHDSPerConditionAndCelltype.png',
          path = outputPath,
          width = 20,
          height = 25,
          units = 'cm')
   
   table(DataFrameHDS[,c(4,5)])
   
   # condition hPP-Hep hInt-Hep hPC-Hep hNASH-Hep
   # Normal    4934     3588    4888      3787
   # NASH      3193     3442    5821      5566

   kruskal.test(HDS ~ sample, data = DataFrameHDS )
   
   # Kruskal-Wallis rank sum test
   # 
   # data:  HDS by sample
   # Kruskal-Wallis chi-squared = 5200.5, df = 5, p-value < 2.2e-16
   
   
   PlotHDSSamples <- 
     ggplot( DataFrameHDS,
             aes( x = sample,
                  y = HDS,
                  color = condition
             )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=25),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=0.5
            )) +
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Diet"))+
     labs(y = "HDS", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange")  +
     stat_kruskal_test(significance = list(cutpoints = 
                                             c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("****", "***", "**", "*", "ns")
                              ),
                       label =  "p.signif",
                       size = 10
                       
                       )

   PlotHDSSamples
   ggsave(file = 'HumanGSE189600WeightedHDSPerSample.png',
          path = outputPath) 
   
   PlotHDSEvery <- 
     ggplot( DataFrameHDS,
             aes( x = HepatocyteType,
                  y = HDS,
                  color = condition
             )) + 
     geom_violin(trim = FALSE) + 
     scale_color_colorblind() +
     theme( text = element_text(size=25),
            plot.title = element_text(size=25),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=0.5
            )) +
     guides(colour = guide_legend(override.aes = aes(label = ""),
                                  title="Diet"))+
     labs(y = "HDS", x = '') + 
     stat_summary(fun.data = "mean_sdl",
                  fun.args = list(mult = 1),
                  geom = "pointrange") + 
     facet_wrap(~sample) 
   
   PlotHDSEvery
   ggsave(file = 'HumanGSE189600WeightedHDSPerCellTypeConditionSample.png',
          path = outputPath,
          width = 25,
          height = 30,
          units = 'cm') 
##### so far: human data has been SCTransformed (before subseting hepatocytes), then merged, meta data 
# from authors has been added matching the cell IDs. 

## no integration was done from my side, but their cell type labeling was done 
# after integration with Harmony

# still to do: calculate HDS on integrated samples to see if there is a difference???  

## not used anymore
## QC A1 Seurat
QCViolins<- VlnPlot(A1Seurat, features = c("nFeature_RNA", 
                               "nCount_RNA", 
                               "percent.mt"), 
        ncol = 3 ) & theme(axis.title.x = element_blank(),
                           axis.text.x  = element_blank())

QC1 <- FeatureScatter(A1Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme(legend.position = "none")
QC2 <- FeatureScatter(A1Seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
  theme(legend.position = "none")
QC3 <- FeatureScatter(A1Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position = "none")

QCPlots <- QCViolins / (QC1 | QC2 | QC3)
QCPlots
ggsave(file = 'QualityControlValuesBeforeFilteringHuman_Y21822A1.png',
       path = paste0(outputPath, 'HumanQC_GSE189600/'))

