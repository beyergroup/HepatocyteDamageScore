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
library(AUCell)

# hepatocyte damage associated genes
HDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')

# load functions to calculate HDS
source('SharedFunctions.R')

pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE200366/'

outputPath = 'hepatocyte-damage-score/Data/Output/'

## SenMayo gene set for senescence

senMayoGenes <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/SenMayoGeneSetMouse.csv',
                         sep = ';')

mayoFeatures <- list(senMayoGenes$Gene.murine.)

# HALLMARK APOPTOSIS

apoptosisFeatures <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/HALLMARK_APOPTOSIS.v2023.2.Mm_onlyGenes.csv',
                              sep = ';',
                              header = FALSE,
                              stringsAsFactors = FALSE
)
apoptosisFeatures <- c(apoptosisFeatures$V1)

# daHep Markers (pre-cancerous cells)
daHep <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/Carlessi_daHepMakers.csv')
daHepUp <- daHep$X[daHep$avg_log2FC > 0 & abs(daHep$avg_log2FC) > 0.5 ]
daHepDown <- daHep$X[daHep$avg_log2FC < 0 & abs(daHep$avg_log2FC) > 0.5]
#

AllHep<- readRDS('hepatocyte-damage-score/Data/Input/scRNAseq/Hep.rds')
View(test@meta.data)

AllHep@meta.data$HDS <- DS_calc.func(exprMatrices = 
                                 GetAssayData(AllHep,
                                              assay = 'RNA',
                                              layer = 'counts'),
                               DSignature = HDAG)
# calc AUCell SenMayo
cells_rankings <- AUCell_buildRankings(LayerData(AllHep, 
                                                 assay = "RNA", 
                                                 layer = "counts"), 
                                       plotStats = TRUE)

cells_AUC <- AUCell_calcAUC(unlist(mayoFeatures), 
                            cells_rankings)

AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, 
                         assign=TRUE)

identical(colnames(cells_AUC@assays@data$AUC), 
          rownames(AllHep@meta.data))

AllHep@meta.data$AUCell_SenMayo <- c(cells_AUC@assays@data$AUC)

# calc AUCell Apoptosis Hallmarks
cells_rankingsApo <- AUCell_buildRankings(LayerData(AllHep, 
                                                 assay = "RNA", 
                                                 layer = "counts"), 
                                       plotStats = TRUE)

cells_AUCApo <- AUCell_calcAUC(apoptosisFeatures, 
                            cells_rankingsApo)

AUCell_exploreThresholds(cells_AUCApo, plotHist=TRUE, 
                         assign=TRUE)

identical(colnames(cells_AUCApo@assays@data$AUC), 
          rownames(AllHep@meta.data))

AllHep@meta.data$AUCell_Apoptosis <- c(cells_AUCApo@assays@data$AUC)


#daHep


cells_rankingsDaHep <- AUCell_buildRankings(LayerData(AllHep, 
                                                    assay = "RNA", 
                                                    layer = "counts"), 
                                          plotStats = TRUE,
                                          nCores = 1,
                                          verbose = progStat)

cells_AUCdaHepUp <- AUCell_calcAUC( daHepUp, cells_rankingsDaHep , 
                             aucMaxRank = ceiling( 0.5 * nrow(cells_rankingsDaHep))
)
cells_AUCdaHepUp <- getAUC( cells_AUCdaHepUp ) 

cells_AUCdaHepDown <- AUCell_calcAUC( daHepDown, cells_rankingsDaHep , 
                                    aucMaxRank = ceiling( 0.5 * nrow(cells_rankingsDaHep))
)
cells_AUCdaHepDown <- getAUC( cells_AUCdaHepDown ) 

identical(names(cells_AUCdaHepDown[1,]), names(cells_AUCdaHepUp[1,]))
cells_AUC_Sub <- cells_AUCdaHepUp[1,] - cells_AUCdaHepDown[1,]

par(mfrow=c(1,3))
hist(cells_AUCdaHepUp[1,])
hist(cells_AUCdaHepDown[1,])
hist(cells_AUC_Sub)
## saved

identical(names(cells_AUC_Sub), 
          rownames(AllHep@meta.data))

AllHep@meta.data$AUCell_daHepSub<- unname(cells_AUC_Sub)

CarUMAP_AUCelldaHepSub <- FeaturePlot(AllHep, features = 'AUCell_daHepSub') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

## Testing other MsigDB Gene Sets ## 
### Oncogene induced cell sescence (MSibDB)
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

intersect(oncoSenGenes, rownames(AllHep))
# 1/8 gene not expressed

cells_rankings <- AUCell_buildRankings(LayerData(AllHep,
                                                 assay = "RNA",
                                                 layer = "counts"),
                                       plotStats = TRUE,
                                       nCores = 1,
                                       verbose = progStat)

cells_AUCtemp <- AUCell_calcAUC( oncoSenGenes, 
                                 cells_rankings,
                                 aucMaxRank = 
                                   ceiling( 0.5 * nrow(cells_rankings))
)

cells_AUCtemp <- getAUC(cells_AUCtemp) 

# AUCell_exploreThresholds(cells_AUCtemp, plotHist=TRUE, 
#                          assign=TRUE)

if(identical(names(cells_AUCtemp[1,]), 
          rownames(AllHep@meta.data)) == TRUE){
  AllHep@meta.data$AUCell_OncogeneInducedSenescence <- unname(cells_AUCtemp[1,])  
}

# MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE
increasedHCCgenes <- read.csv('hepatocyte-damage-score/Data/Input/GeneSetsCellFates/MP_INCREASED_HEPATOCELLULAR_CARCINOMA_INCIDENCE.v2023.2.Mm.csv',
                              header = FALSE)
increasedHCCgenes <- increasedHCCgenes[,1]

cells_AUCtemp <- AUCell_calcAUC( increasedHCCgenes, 
                                 cells_rankings,
                                 aucMaxRank = 
                                   ceiling( 0.5 * nrow(cells_rankings))
)

# Genes in the gene sets NOT available in the dataset: 
#   geneSet: 	1 (2% of 55)
cells_AUCtemp <- getAUC(cells_AUCtemp) 

if(identical(names(cells_AUCtemp[1,]), 
             rownames(AllHep@meta.data)) == TRUE){
  AllHep@meta.data$AUCell_increasedHCCincidence <- unname(cells_AUCtemp[1,])  
}

# MP_INCREASED_LIVER_TUMOR_INCIDENCE
increasedHCCgenes <- read.csv('hepatocyte-damage-score/Data/Input/GeneSetsCellFates/MP_INCREASED_LIVER_TUMOR_INCIDENCE.v2023.2.Mm.csv',
                              header = FALSE)
increasedHCCgenes <- increasedHCCgenes[,1]

cells_AUCtemp <- AUCell_calcAUC( increasedHCCgenes, 
                                 cells_rankings,
                                 aucMaxRank = 
                                   ceiling( 0.5 * nrow(cells_rankings))
)

# Genes in the gene sets NOT available in the dataset: 
#   geneSet: 	1 (2% of 40)

cells_AUCtemp <- getAUC(cells_AUCtemp) 

if(identical(names(cells_AUCtemp[1,]), 
             rownames(AllHep@meta.data)) == TRUE){
  AllHep@meta.data$AUCell_increasedLiverTumorIncidence <- unname(cells_AUCtemp[1,])  
}

### MsigDB: Human Family of Oncogenes converde to Mouse Gene Symbols 

hOncoGenes <- read.csv('hepatocyte-damage-score/Data/Input/GeneSetsCellFates/MsigDBOncogenesFamilyHumanConvertedToMouse.csv',
                       header = FALSE)
hOncoGenes <- hOncoGenes[,1]

cells_AUCtemp <- AUCell_calcAUC( hOncoGenes, 
                                 cells_rankings,
                                 aucMaxRank = 
                                   ceiling( 0.5 * nrow(cells_rankings))
)

# 14/314 genes not in expression matrix (4%)

AUCell_exploreThresholds(cells_AUCtemp, plotHist=TRUE, assign=TRUE)
# curve doesn't look like any particular group of cells is enriched for the gene set

cells_AUCtemp <- getAUC(cells_AUCtemp) 

if(identical(names(cells_AUCtemp[1,]), 
             rownames(AllHep@meta.data)) == TRUE){
  AllHep@meta.data$AUCell_humanOncoGenes <- 
    unname(cells_AUCtemp[1,])  
}


###


## store in Dataframe

CarlessiDataFrame <- AllHep@meta.data

# Oncogene Induced Senscence correlation with HDS?

ggplot(CarlessiDataFrame,
       aes(x = AUCell_OncogeneInducedSenescence, 
           y = HDS,
           colour = cell_type)) +
  geom_point(size = 0.3, alpha = 0.5) + scale_color_colorblind()

ggplot(CarlessiDataFrame,
       aes(x = AUCell_OncogeneInducedSenescence, 
           y = HDS,
           colour = Group)) +
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_colorblind() + geom_smooth(method = lm)

ggplot(CarlessiDataFrame,
       aes(x = AUCell_OncogeneInducedSenescence, 
           y = HDS)) +
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_colorblind() + geom_smooth(method = lm)

cor.test(x = CarlessiDataFrame$HDS,
         y = CarlessiDataFrame$AUCell_OncogeneInducedSenescence)


cor.test(x = CarlessiDataFrame$HDS,
         y = CarlessiDataFrame$AUCell_increasedHCCincidence)
ggplot(CarlessiDataFrame,
       aes(x = AUCell_increasedHCCincidence, 
           y = HDS)) +
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_colorblind() + geom_smooth(method = lm)

cor.test(x = CarlessiDataFrame$HDS,
         y = CarlessiDataFrame$AUCell_increasedLiverTumorIncidence)
ggplot(CarlessiDataFrame,
       aes(x = AUCell_increasedLiverTumorIncidence, 
           y = HDS)) +
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_colorblind() + geom_smooth(method = lm)


cor.test(x = CarlessiDataFrame$HDS,
         y = CarlessiDataFrame$AUCell_humanOncoGenes)
ggplot(CarlessiDataFrame,
       aes(x = AUCell_humanOncoGenes, 
           y = HDS,
           )) +
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_colorblind() + geom_smooth(method = lm)


# Plot UMAPs
CarUMAP_HDS <- FeaturePlot(AllHep, features = 'HDS') + 
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

CarUMAP_Celltype <-DimPlot(AllHep) + scale_color_colorblind()

CarUMAP_AUCellSenMayo <- FeaturePlot(AllHep, features = 'AUCell_SenMayo') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

CarUMAP_AUCellSApoptosis <- FeaturePlot(AllHep, features = 'AUCell_Apoptosis') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

CarUMAP_AUCellOncoSen <- FeaturePlot(AllHep, features = 'AUCell_OncogeneInducedSenescence') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

CarUMAP_AUCellincreasedHCC <- FeaturePlot(AllHep, features = 'AUCell_increasedHCCincidence') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

CarUMAP_AUCellincreaseLiverTumor <- FeaturePlot(AllHep, features = 'AUCell_increasedLiverTumorIncidence') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)
CarUMAP_AUCellHumanOncoGenes <- FeaturePlot(
  AllHep, features = 'AUCell_humanOncoGenes') +
  scale_color_distiller(palette = "YlOrBr",
                        direction = 1)

toSave <- (CarUMAP_Celltype | CarUMAP_HDS)

CarlessiUmaps <- (DimPlot(AllHep,group.by = 'SampleID')|DimPlot(AllHep,group.by = 'Group') | CarUMAP_Celltype) / (CarUMAP_HDS |CarUMAP_AUCellSenMayo | CarUMAP_AUCellSApoptosis)
CarlessiUmaps
ggsave('CarlessiHDSRNAcountsSenMayoApoptosisSampleIDGroup.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 35,
       units = 'cm')

CarUMAP_HDS | CarUMAP_AUCelldaHepSub
ggsave('CarlessiHDSRNAcountsAucelldaHepSubUMAPs.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 25,
       units = 'cm')


View(AllHep@meta.data)


CarlessiHDSViolinCellTypes <- 
  ggplot( CarlessiDataFrame,
          aes( x = SampleID,
               y = HDS,
               color = Group
          )) +
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Group")) +
  labs(y = "HDS", x = 'Samples') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") + 
  facet_wrap(~cell_type)

CarlessiHDSViolinCellTypes
ggsave('CarlessiHDSofRNAcountsFacetHepZones.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 25,
       units = 'cm')

CarlessiHDSViolinSamples <- 
  ggplot( CarlessiDataFrame,
          aes( x = cell_type,
               y = HDS,
               color = cell_type
          )) +
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Cell type")) +
  labs(y = "HDS", x = 'Samples') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") + 
  facet_wrap(~SampleID)

CarlessiHDSViolinSamples
ggsave('CarlessiHDSofRNAcountsFacetSamples.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 25,
       units = 'cm')

CarlessiHDSViolinGroup <- 
  ggplot( CarlessiDataFrame,
          aes( x = cell_type,
               y = HDS,
               color = cell_type
          )) +
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Cell type")) +
  labs(y = "HDS", x = 'Samples') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") + 
  facet_wrap(~Group)

CarlessiHDSViolinGroup
ggsave('CarlessiHDSofRNAcountsFacetSamples.pdf',
       path = paste0(outputPath, 'Results/'),
       height = 16,
       width = 25,
       units = 'cm')

CarlessiSenMayoViolinCellTypes <- 
  ggplot( CarlessiDataFrame,
          aes( x = SampleID,
               y = AUCell_SenMayo,
               color = Group
          )) +
  geom_violin(trim = FALSE) +
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=24),
         legend.position = 'right',
         axis.text.x = element_blank()) + 
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Group")) +
  labs(y = "AUCell_SenMayo", x = 'Samples') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") + 
  facet_wrap(~cell_type)

CarlessiSenMayoViolinCellTypes

### CORRELATION HDS & SenMayo per Sample
allSamples <- unique(CarlessiDataFrame$SampleID)

allCorPlotsSamples <- lapply(1:9, function(i){
  return( ggplot(CarlessiDataFrame[CarlessiDataFrame$SampleID == allSamples[i],],
                 aes(x = AUCell_SenMayo, 
                     y = HDS,
                     colour = cell_type)) +
            geom_point(size = 0.3, alpha = 0.5) + scale_color_colorblind()
  )
})

allCorPlotsSamplesWithLM <- lapply(1:9, function(i){
 return( ggplot(CarlessiDataFrame[CarlessiDataFrame$SampleID == allSamples[i],],
         aes(x = AUCell_SenMayo, 
             y = HDS,
             colour = cell_type)) +
    geom_point(size = 0.3) + scale_color_colorblind() +
      geom_smooth(method=lm)
  )
})

allCorPlotsSamplesApop <- lapply(1:9, function(i){
  return( ggplot(CarlessiDataFrame[CarlessiDataFrame$SampleID == allSamples[i],],
                 aes(x = AUCell_Apoptosis, 
                     y = HDS,
                     colour = cell_type)) +
            geom_point(size = 0.3) + scale_color_colorblind()
  )
})

allCorPlotsSamplesApopWithLM <- lapply(1:9, function(i){
  return( ggplot(CarlessiDataFrame[CarlessiDataFrame$SampleID == allSamples[i],],
                 aes(x = AUCell_Apoptosis, 
                     y = HDS,
                     colour = cell_type)) +
            geom_point(size = 0.3) + scale_color_colorblind() +
            geom_smooth(method=lm)
  )
})

HDS_SenMayo_Samples <- (allCorPlotsSamples[[1]] | allCorPlotsSamples[[2]] | allCorPlotsSamples[[3]]) / (allCorPlotsSamples[[4]] | allCorPlotsSamples[[5]] | allCorPlotsSamples[[6]]) /  (allCorPlotsSamples[[7]] | allCorPlotsSamples[[8]] | allCorPlotsSamples[[9]]) 
HDS_Apoptosis_Samples <- (allCorPlotsSamplesApop[[1]] | allCorPlotsSamplesApop[[2]] | allCorPlotsSamplesApop[[3]]) / (allCorPlotsSamplesApop[[4]] | allCorPlotsSamplesApop[[5]] | allCorPlotsSamplesApop[[6]]) /  (allCorPlotsSamplesApop[[7]] | allCorPlotsSamplesApop[[8]] | allCorPlotsSamplesApop[[9]]) 

HDS_SenMayo_Samples 
ggsave('CarlessisHDS_SenMayo_SampleScatterPlots.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 20,
       width = 30,
       units = 'cm')

###

hist(CarlessiDataFrame$HDS)
CarlessiDataFrame <- CarlessiDataFrame %>% mutate(HDS_bin = cut(HDS, breaks=25))

tempA <- ggplot(data = CarlessiDataFrame, aes(x=HDS_bin, y=AUCell_daHepSub) ) +
  geom_boxplot(fill="#69b3a2") +
  xlab("HDS (bins) ")

tempA
ggsave('CarlessisHDS25BinsBoxplotAUCelldaHepSub.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 35,
       units = 'cm')

tempB<- ggplot(data = CarlessiDataFrame, aes(x=HDS_bin, y=AUCell_SenMayo) ) +
  geom_boxplot(fill="#69b3a2") +
  xlab("HDS (bins) ")


tempB
ggsave('CarlessisHDS25BinsBoxplotAUCellSenMayo.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 35,
       units = 'cm')

tempC <-ggplot(data = CarlessiDataFrame, aes(x=HDS_bin, y=AUCell_Apoptosis) ) +
  geom_boxplot(fill="#69b3a2") +
  xlab("HDS (bins) ")

tempC
ggsave('CarlessisHDS25BinsBoxplotAUCellApoptosis.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 15,
       width = 35,
       units = 'cm')

allTemp <- (tempA) / tempB / tempC
allTemp
ggsave('CarlessisHDS25BinsBoxplotAllAucellScores.pdf',
       path = paste0(outputPath, 'Results/Trajectories/'),
       height = 30,
       width = 35,
       units = 'cm')

###

tempA <- ggplot(data = CarlessiDataFrame, aes(x=HDS_bin, y=AUCell_daHepSub) ) +
  geom_boxplot(fill="#69b3a2") +
  xlab("HDS (bins) ")

