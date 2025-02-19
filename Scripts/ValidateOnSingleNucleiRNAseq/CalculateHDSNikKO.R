# Applying HDS to Thomas scRNAseq data
library(Seurat)
library(AUCell)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(viridis)

# Load Wunderlich Group data from NIK-KO experiment 
# - snRNA-seq
# - pre-processed, filtered and normalized, annotated, etc. by them

data <- readRDS('~/Desktop/HepatocyteDamageScore/scRNAseq/Thomas_Wunderlich_snRNAseq/nik_samples_annotated.rds')
data <- subset(data, 
               subset= class_annotation == "Hepatocyte")
gc()

# old PCA values: PCA had been calculated before subsetting hepatocytes
DimPlot(data, reduction = 'pca',group.by = 'Diet', 
        shuffle = TRUE,pt.size = 0.3,
        dims = c(5,6), split.by ='Genetic') 

head(Stdev(data, reduction = "pca"), 10)

# Run PCA on the hepatocyte subset
# don't know if this is mathematically correct! (ask Andreas!)

# Finding most variable genes is necessary before PCA (top 2000 default)
data <- FindVariableFeatures(data)
VariableFeaturePlot(data)
plotVariableFeatures <- 
  LabelPoints(plot = VariableFeaturePlot(data), 
              points = head(VariableFeatures(data), 10), 
              repel = TRUE)

# before PCA data has to be rescaled, the problem is that before subsetting the 
# cells by celltype == 'hepatocyte', the counts had been already scaled
# so I don't know if scaling them ONE MORE TIME is problematic
# ALSO: do I then calculate the HDS on these scaled values or on the "normal" ones

data <- ScaleData(data, 
                  features = rownames(data))
gc()

# perform PCA:
data <- RunPCA(data, 
               features = VariableFeatures(object = data))


DimPlot(data, reduction = 'umap',
        
        group.by = 'Diet', 
        split.by = 'Genetic',
        shuffle = TRUE, 
        pt.size = 0.1) + scale_color_colorblind(alpha(0.3))


metadata_hepatocytes <- data@meta.data

# PC Eigenvalues:
head(Stdev(data, reduction = "pca"), 10)
# 8.681295 8.247859 5.642825 4.562431 3.870043 3.480020 3.406927 2.799780 2.569180 2.479066
# Conclusion: each PC explain only small fractions of the variance, 
# meaning there is no clear strong trend separating the cells. 

# Plot PCA:
# PC1 & PC2
DimPlot(data, reduction = 'pca', 
        group.by = 'SampleName', 
        shuffle = TRUE,pt.size = 0.3,
        dims = c(1,2),
        split.by = 'Diet') 

DimPlot(data, reduction = 'pca', 
        group.by = 'SampleName', 
        shuffle = TRUE, pt.size = 0.3,
        dims = c(3,4),
        split.by = 'Diet') 

DimPlot(data, reduction = 'pca', 
        group.by = 'SampleName', 
        shuffle = TRUE, pt.size = 0.3,
        dims = c(5,6),
        split.by = 'Diet')

# Nice separation of diets and genotypes
PC1_6 <- DimPlot(data, reduction = 'pca', 
                 group.by = 'SampleName', 
                 shuffle = TRUE,pt.size = 0.3,
                 dims = c(1,6)) + scale_color_colorblind() + ggtitle('Conditions')



# load functions
source('SharedFunctions.R')


# Load list of Hepatocyte Damage Associated Genes
HDAG <- 
  read.csv(file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')

DS_nash <- DS_calc.func( 
  exprMatrices = GetAssayData( object = subset(data, subset = Diet == 'Nash'),
                               slot = 'counts'),
  ceilThrsh = 0.05 ,
  DSignature = HDAG ,
  ntop = 42, useOrder ="mean_rank" )

gc()

DS_cdaa <- DS_calc.func( 
  exprMatrices = GetAssayData(object = subset(data, subset = Diet == 'Cdaa'),
                              slot = 'counts'),
  ceilThrsh = 0.05 ,
  DSignature= HDAG ,
  ntop = 42, useOrder ="mean_rank" )
gc()

DS_hcc <- DS_calc.func(
  exprMatrices = GetAssayData(object = subset(data, subset = Diet == 'Hcc'),
                                     slot = 'counts'),
  ceilThrsh = 0.05 ,
  DSignature= HDAG ,
  ntop = 42, useOrder ="mean_rank" )
gc()


#Format data for plotting

DS_all <- data.frame('Cell_ID' = names(DS_nash),
                     'DS' = unname(DS_nash),
                     'Diet' = rep('NASH', length(names(DS_nash))))

DS_all <- rbind(DS_all, data.frame('Cell_ID' = names(DS_cdaa),
                                   'DS' = unname(DS_cdaa),
                                   'Diet' = rep("CDAA", length(names(DS_cdaa)))))

DS_all <- rbind(DS_all, data.frame('Cell_ID' = names(DS_hcc),
                                   'DS' = unname(DS_hcc),
                                   'Diet' = rep("HCC", length(names(DS_hcc)))))

DS_all$Genotype <- unlist(lapply(DS_all$Cell_ID, function(ii){
  return(metadata_hepatocytes[match(ii, metadata_hepatocytes$Cell_ID), 'Genetic']) 
}))

DS_all$Genotype <- factor(DS_all$Genotype,
                      ordered = TRUE,
                      levels = c("NIK-FL", "NIK-KO"))

if(identical(data@meta.data$Cell_ID, DS_all$Cell_ID) == TRUE){
  data@meta.data$HDS <- DS_all$DS
}else(print('Cell IDs do not match!'))

# Plot PCA

PC1_6_HDS <- FeaturePlot(data, 'HDS',
            dims = c(1, 6), reduction = "pca",
            repel = TRUE,
            pt.size = 0.3,
            order = FALSE,
            keep.scale = 'all') + 
  scale_color_viridis(option = "B", direction = -1) + 
  ggtitle('Hepatocyte Damage Score')

panelPCAs <- (PC1_6 | PC1_6_HDS) + 
  plot_annotation(tag_levels = 'A')


ggsave(filename = 'NIKKO_PC1_PC6_weightedHDS.png', 
       path = 'hepatocyte-damage-score/Data/Output/Results/'
)

# Violin plots

plotNASH <- ggplot( DS_all[DS_all$Diet == 'NASH',], 
        aes( x = Genotype  , y = DS,
             color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("NASH") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") +
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange"
  ) 

plotCDAA <- ggplot( DS_all[DS_all$Diet == 'CDAA',], 
                    aes( x = Genotype  , y = DS,
                         color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("CDAA") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") +
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

plotDEN <- ggplot( DS_all[DS_all$Diet == 'HCC',], 
                    aes( x = Genotype  , y = DS,
                         color = Genotype)) + 
  geom_violin(trim = FALSE) + #theme_bw() + 
  scale_color_colorblind() +
  ggtitle("DEN") +
  theme( text = element_text(size=12),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(size=12),
         legend.position = 'right') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS") + 
  coord_cartesian(ylim = c(-0.5, 0.5)) + 
  scale_y_continuous(breaks=c(-0.5,-0.25, 0, 0.25, 0.5)) +
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

panelWunderlich <- (plotNASH | plotCDAA | plotDEN) + 
  plot_annotation(tag_levels = 'A')

ggsave(filename = 'ViolinPlotsHDSNIKallHepsWeightedFunction.png', 
       path = 'hepatocyte-damage-score/Data/Output/Results/'
       )

## summary statistics

summary(DS_all)

# NASH
summary(DS_all[DS_all$Diet == 'NASH' &
                 DS_all$Genotype == 'NIK-FL', ])

summary(DS_all[DS_all$Diet == 'NASH' &
                 DS_all$Genotype == 'NIK-KO', ])

# CDAA 
summary(DS_all[DS_all$Diet == 'CDAA' &
                 DS_all$Genotype == 'NIK-FL', ])
summary(DS_all[DS_all$Diet == 'CDAA' &
                 DS_all$Genotype == 'NIK-KO', ])

# DEN
summary(DS_all[DS_all$Diet == 'HCC' &
                 DS_all$Genotype == 'NIK-FL', ])
summary(DS_all[DS_all$Diet == 'HCC' &
                 DS_all$Genotype == 'NIK-KO', ])

