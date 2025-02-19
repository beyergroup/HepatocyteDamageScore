
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
library(ProliferativeIndex)
library(AUCell)


pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/'

outputPath = 'hepatocyte-damage-score/Data/Output/'

mergedHepatocytes<- readRDS(file =
                              'Repositories/cell-damage-score/hepatocyte-damage-score/Data/Output/GSE189600SeuratAllMouseHepMerged.rds')


VariableFeatures(mergedHepatocytes[["SCT"]]) <- rownames(mergedHepatocytes[["SCT"]]@scale.data)

mergedHepatocytes <- RunPCA(mergedHepatocytes,
                            feautres = VariableFeatures(mergedHepatocytes))

DimPlot(mergedHepatocytes, 
        reduction = "pca", 
        group.by = 'CellType')

mergedHepatocytes <- RunUMAP(mergedHepatocytes, 
                             dims = 1:20)

mergedHepatocytes  <- FindNeighbors(mergedHepatocytes, 
                                    dims = 1:20)
mergedHepatocytes<- FindClusters(mergedHepatocytes, 
                                 resolution = 0.2)

# senescence
# p21 and p53
VlnPlot(mergedHepatocytes, 
        features = c("Cdkn1a"),
        group.by = 'CellType',
        split.by = 'condition')


VlnPlot(mergedHepatocytes, 
        features = c("Trp53"),
        group.by = 'CellType',
        split.by = 'condition')

FeaturePlot(mergedHepatocytes, 
            features = c('Cdkn1a','Trp53'))
# try p16 (Cdkn2a) 

FeaturePlot(mergedHepatocytes, 
            features = c('Glipr1','Clec12a', 'Phlda3'))


# HCC

VlnPlot(mergedHepatocytes, 
        features = c("Afp"),
        group.by = 'CellType',
        split.by = 'condition')

VlnPlot(mergedHepatocytes, 
        features = c("Tgfb1"),
        group.by = 'CellType',
        split.by = 'condition')

FeaturePlot(mergedHepatocytes, 
            features = c('Afp','Tgfb1'))

# necroptosis
FeaturePlot(mergedHepatocytes, 
            features = c('Ripk1','Ripk3', 'Mlkl'))

# apoptosis
FeaturePlot(mergedHepatocytes, 
            features = c('Tnfsf10','Tnfrsf10b', 'Fadd', 'Casp8'))

#ferroptosis


# Cell death 

minHDS <- min(mergedHepatocytes@meta.data$HDS)
maxHDS <- max(mergedHepatocytes@meta.data$HDS)

UmapHDS <- FeaturePlot(mergedHepatocytes,
                             features ='HDS') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1,
                        limits = c(minHDS, maxHDS)) 


## SenMayo gene set for senescence

senMayoGenes <- read.csv(file = 'hepatocyte-damage-score/Data/Input/SenMayoGeneSetMouse.csv',
                         sep = ';')

mayoFeatures <- list(senMayoGenes$Gene.murine.)

# which SenMayo genes are expressed? it's only 39 from 125
intersectGenes <- 
  intersect(rownames(mergedHepatocytes[["SCT"]]@scale.data), 
            unlist(mayoFeatures))


mergedHepatocytes <- AddModuleScore(
  object = mergedHepatocytes,
  features = mayoFeatures,
  ctrl = 50,
  name = 'Mayo_Features'
)
View(mergedHepatocytes@meta.data)

UmapMayoModuleScore <- FeaturePlot(mergedHepatocytes,
                       features ='Mayo_Features1') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

ggplot(mergedHepatocytes@meta.data, aes(x=HDS, y=Mayo_Features1)) +
  geom_point()


## AUCell approach

cells_rankings <- AUCell_buildRankings(LayerData(mergedHepatocytes, 
                                                 assay = "SCT", 
                                                 layer = "counts"), 
                                       plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(unlist(mayoFeatures), 
                            cells_rankings)

AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, 
                         assign=TRUE,
                         smallestPopPercent = 0.01,
                         thrP = 0.10)

identical(colnames(cells_AUC@assays@data$AUC), 
          rownames(mergedHepatocytes@meta.data))

mergedHepatocytes@meta.data$AUCell_SenMayo <- c(cells_AUC@assays@data$AUC)


# try on one sample

cells_rankings9moNASH <- AUCell_buildRankings(LayerData(subset(mergedHepatocytes, 
                                                        subset = condition == '9m NASH'), 
                                                 assay = "SCT", 
                                                 layer = "counts"), 
                                       plotStats=TRUE)

cells_AUC9moNASH <- AUCell_calcAUC(unlist(mayoFeatures), 
                            cells_rankings9moNASH)

AUCell_exploreThresholds(cells_AUC9moNASH, plotHist=TRUE, assign=TRUE,
                         smallestPopPercent = 0.01)

ggplot(mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '3m NC',], aes(x=HDS, y=Mayo_Features1)) +
  geom_point()


cor.test(mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '3m NC',]$HDS,
         mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '3m NC',]$Mayo_Features1)


##
FeaturePlot(mergedHepatocytes,features = 'Mayo_Features1')
FeaturePlot(mergedHepatocytes, features = 'AUCell_SenMayo')

VlnPlot(mergedHepatocytes, features = "Mayo_Features1", split.by = "CellType",
        group.by = 'condition')
VlnPlot(mergedHepatocytes, features = "AUCell_SenMayo", split.by = "CellType",
        group.by = 'condition')
VlnPlot(mergedHepatocytes, features = "HDS", split.by = "CellType",
        group.by = 'condition')

DoHeatmap(mergedHepatocytes,features = unlist(mayoFeatures), group.by = "CellType")
DoHeatmap(subset(mergedHepatocytes, 
                 subset = condition == '9m NASH'),
          features = unlist(mayoFeatures), group.by = "CellType")
DoHeatmap(subset(mergedHepatocytes, 
                 subset = condition == '3m NASH'),
          features = unlist(mayoFeatures), group.by = "CellType")


#AUCell correlation 
ggplot(mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '9m NASH',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point()

cor.test(mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '9m NASH',]$HDS,
         mergedHepatocytes@meta.data[mergedHepatocytes@meta.data$condition == '9m NASH',]$AUCell_SenMayo)

## DataFrame

DataFrameHDS <- data.frame('HDS' = mergedHepatocytes@meta.data$HDS,
                           'CellID' = rownames(mergedHepatocytes@meta.data),
                           'sample' = mergedHepatocytes@meta.data$condition,
                           'HepatocyteType' = mergedHepatocytes@meta.data$CellType,
                           'AUCell_SenMayo' = mergedHepatocytes@meta.data$AUCell_SenMayo,
                           'MayoFeatures1' = mergedHepatocytes@meta.data$Mayo_Features1)


DataFrameHDS$sample <- factor(DataFrameHDS$sample,
                              ordered = TRUE, 
                              levels = c("3m NC",
                                         "9m NC",
                                         "3m NASH",
                                         "9m NASH"))


DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                      ordered = TRUE,
                                      levels = c("PP-Hep",
                                                 "Int-Hep",
                                                 "PC-Hep",
                                                 "mNASH-Hep1",
                                                 "mNASH-Hep2"))

#SenMayo AUCell Result

PlotAUCellSenMayoHepTypes <- 
  ggplot( DataFrameHDS,
          aes( x = sample,
               y = AUCell_SenMayo,
               color = sample
          )) + 
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Sample")) +
  labs(y = "AUCellSenMayo", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~HepatocyteType)

PlotAUCellSenMayoHepTypes
ggsave('XiaoMouseSenMayoAUCellFacetHepZones.pdf',
       path = paste0(outputPath, 'Results/Trajectories'),
       height = 16,
       width = 25,
       units = 'cm')


PlotAUCellSenMayoSamples <- 
  ggplot( DataFrameHDS,
          aes( x = HepatocyteType,
               y = AUCell_SenMayo,
               color = HepatocyteType
          )) + 
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Hepatocyte Zone")) +
  labs(y = "AUCellSenMayo", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~sample)

PlotAUCellSenMayoSamples
ggsave('XiaoMouseSenMayoAUCellFacetSamples.pdf',
       path = paste0(outputPath, 'Results/Trajectories'),
       height = 16,
       width = 20,
       units = 'cm')


PlotSenMayoFeaturesHepTypes <- 
  ggplot( DataFrameHDS,
          aes( x = sample,
               y = MayoFeatures1,
               color = sample
          )) + 
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Sample")) +
  labs(y = "SenMayoModule", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~HepatocyteType)


PlotSenMayoFeaturesSamples <- 
  ggplot( DataFrameHDS,
          aes( x = HepatocyteType,
               y = MayoFeatures1,
               color = HepatocyteType
          )) + 
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Hepatocyte Zone")) +
  labs(y = "SenMayoModule", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange")  +
  facet_wrap(~sample)

## Correlation plots 

corPlot1 <- ggplot(DataFrameHDS[DataFrameHDS$sample == '9m NASH' & 
                      DataFrameHDS$HepatocyteType == 'PP-Hep',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point(size = 0.3) +
  geom_smooth(method=lm)


corPlot2 <- ggplot(DataFrameHDS[DataFrameHDS$sample == '9m NASH' & 
                      DataFrameHDS$HepatocyteType == 'Int-Hep',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point(size = 0.3) +
  geom_smooth(method=lm)


corPlot3 <- ggplot(DataFrameHDS[DataFrameHDS$sample == '9m NASH' & 
                      DataFrameHDS$HepatocyteType == 'PC-Hep',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point(size = 0.3) +
  geom_smooth(method=lm)


corPlot4 <- ggplot(DataFrameHDS[DataFrameHDS$sample == '9m NASH' & 
                      DataFrameHDS$HepatocyteType == 'mNASH-Hep1',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point(size = 0.3) +
  geom_smooth(method=lm)

corPlot5 <- ggplot(DataFrameHDS[DataFrameHDS$sample == '9m NASH' & 
                      DataFrameHDS$HepatocyteType == 'mNASH-Hep2',], 
       aes(x=HDS, y=AUCell_SenMayo)) +
  geom_point(size = 0.3) +
  geom_smooth(method=lm)

patchCor <- (corPlot1 | corPlot2 | corPlot3) / (corPlot4|corPlot5)

# HALLMARK APOPTOSIS

apoptosisFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/HALLMARK_APOPTOSIS.v2023.2.Mm_onlyGenes.csv',
                         sep = ';',
                         header = FALSE,
                         stringsAsFactors = FALSE
                        )
apoptosisFeatures <- c(apoptosisFeatures$V1)

intersectGenesApoptosis <- 
  intersect(rownames(mergedHepatocytes[["SCT"]]@scale.data), 
            unlist(apoptosisFeatures))
length(intersectGenesApoptosis)
# [1] 72

# AUCell both genes sets 
geneSets <- list('SenMayo' = unlist(mayoFeatures),
                 'ApoptosisHallmark' = apoptosisFeatures)

cells_AUC_both <- AUCell_calcAUC(geneSets, cells_rankings)
# Genes in the gene sets NOT available in the dataset: 
#   SenMayo: 	44 (37% of 118)
# ApoptosisHallmark: 	26 (16% of 161)
par(mfrow=c(1,2)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC_both, 
                                             plotHist=TRUE, 
                                             assign=TRUE,
                                             smallestPopPercent = 0.01,
                                             thrP = 0.10)

identical(colnames(cells_AUC_both@assays@data$AUC), 
          rownames(mergedHepatocytes@meta.data))

mergedHepatocytes@meta.data$AUCell_ApoptosisHallmark <- c(cells_AUC_both@assays@data$AUC)


