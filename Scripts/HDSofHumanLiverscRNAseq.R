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

A6Seurat <- SCTransform(A6Seurat)


remove(countsA6)
gc()

allSamples <- list(A1Seurat, A2Seurat, A3Seurat,
                   A4Seurat, A5Seurat, A6Seurat)

#If you want to merge the normalized data matrices as well as the raw count matrices
# pass merge.data = TRUE

mergedSamples <- merge(allSamples[[1]], y = c(allSamples[[2]],
                                              allSamples[[3]],
                                              allSamples[[4]],
                                              allSamples[[5]],
                                              allSamples[[6]]), 
                       add.cell.ids = c("Normal1", "Normal2", "Normal3",
                                        "NASH1", "NASH2", "NASH3"), 
                       project = "HumanXiaoEtAl",
                       merge.data = TRUE)

remove(allSamples)
remove(A1Seurat, A2Seurat, A3Seurat, A4Seurat, A5Seurat, A6Seurat)

metaDataAnnotations <- read.table(
  paste0(pathData,
         'scitranslmed.adc9653_data_files_s1_to_s6/humanSnRNAseqSeuratAnnotations.csv'), 
  sep = ';', dec = ',', header = TRUE)

mergedSamples@meta.data$cell_id <- rownames(mergedSamples@meta.data)


temp <- metaDataAnnotations
rownames(temp) <- temp$X
temp <- temp[intersect(temp$X, rownames(mergedSamples@meta.data)),]

mergedSamples@meta.data$cellBarcode <- NA
mergedSamples@meta.data$cellBarcode[match(temp$X, rownames(mergedSamples@meta.data))] <- 
  temp$X
mergedSamples@meta.data$CellType[match(temp$X, rownames(mergedSamples@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
mergedSamples@meta.data$deleteNA <- 
  rep("not na", length(mergedSamples@meta.data$CellType))
mergedSamples@meta.data$deleteNA[is.na(mergedSamples@meta.data$CellType)] <-
  "na"

mergedSamples@meta.data <- subset(mergedSamples@meta.data,
                        subset = deleteNA == 'not na')

saveRDS(object = mergedSamples, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanSamplesMerged.rds')

mergedHepatocytes <- subset(x = mergedSamples, 
                            subset = CellType ==  "hPC-Hep" |
                         CellType == "hNASH-Hep"  |
                         CellType == "hPP-Hep" |
                         CellType == "hInt-Hep" )

saveRDS(object = mergedHepatocytes, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')
# Very suspicious: after subsetting the subset has as many cells as genes 

# calculate HDS with HDAG list mapped to human genes
mergedHepatocytes <- readRDS( file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')

humanHDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesOnlyHumanGenesManual.csv',
    sep = ';')

# load functions to calculate HDS
source('SharedFunctions.R')

# create Expression Matrix

ExprMatrix <- GetAssayData(mergedHepatocytes)

# genesMatrix <- rownames(ExprMatrix)
#
# matchingEnsemblGeneNames <- match(genesMatrix, gtfSelected$gene_id)
# 
# matchingGeneNames <- gtfSelected[matchingEnsemblGeneNames,]
# 
# genesNotMapped <- genesMatrix[is.na(matchingGeneNames$hgnc_gene)]
# genesNotMapped <- match(genesNotMapped, rownames(ExprMatrix))
# 
# # remove counts from 80 genes that dont match 
# 
# ExprMatrix <- ExprMatrix[-genesNotMapped,]
# 
# # remove NAs from mapped gene list
# matchingGeneNames <- matchingGeneNames
# matchingGeneNames <- matchingGeneNames[!(is.na(matchingGeneNames$hgnc_gene)),]
# 
# 
# if(identical(rownames(ExprMatrix), matchingGeneNames$gene_id)){
#   rownames(ExprMatrix) <- matchingGeneNames$hgnc_gene
#   print('hooray!')
# }

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

DataFrameHDS <- data.frame('HDS' = unname(HDS),
                           'CellID' = names(unlist(HDS)))

identical(DataFrameHDS$CellID, rownames(mergedHepatocytes@meta.data))

mergedHepatocytes@meta.data$sample <- 
  gsub('_.*', replacement = '', mergedHepatocytes@meta.data$cellBarcode)

DataFrameHDS$sample <- mergedHepatocytes@meta.data$sample

DataFrameHDS$condition <-gsub('[1-3]_.*', 
                              replacement = '', 
                              mergedHepatocytes@meta.data$cellBarcode)


DataFrameHDS$HepatocyteType <- mergedHepatocytes@meta.data$CellType


DataFrameHDS$condition <- factor(DataFrameHDS$condition,
                           ordered = TRUE,
                           levels = c('Normal', 'NASH'))

DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                 ordered = TRUE,
                                 levels = c('hPP-Hep', 'hInt-Hep',
                                            'hPC-Hep','hNASH-Hep'))




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
          path = outputPath)
   
   table(DataFrameHDS[,c(3,4)])
   # HepatocyteType
   # condition hInt-Hep hNASH-Hep hPC-Hep hPP-Hep
   # Normal     3588      3787    4888    4934
   # NASH       3442      5566    5821    3193
   #

   kruskal.test(HDS ~ sample, data = DataFrameHDS )
   
   # Kruskal-Wallis rank sum test
   # 
   # data:  HDS by sample
   # Kruskal-Wallis chi-squared = 6279, df = 5, p-value < 2.2e-16
   
   
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
          path = outputPath) 
##### so far: human data has been normalized, then merged, meta data 
# from authors has been added matching the cell IDs. 

## no integration was done from my side, but their cell type labeling was done 
# after integration with Harmony

# still to do: calculate HDS on integrated samples to see if there is a difference 

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

