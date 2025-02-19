library(Seurat)
library(AUCell)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(viridis)

#### load functions written by us
{
  source('SharedFunctions.R')
}

outputPath <- 'hepatocyte-damage-score/Data/Output/Results/'


# Hepatocyte Universal Damage Signature
HDAG <- read.csv(file = 
                   'hepatocyte-damage-score/Data/Output/HDAG.csv')

# load snRNAseq counts 

#load Seurat object
mergedLCAnucSeq <- 
  readRDS('hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/countTable_mouseMerged_Hepatocytes_nucSeq.rds')

# sample information needed for subsetting properly 
sampleInfo <- 
  read.csv('hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/GSE192740_sampleInfo_modified.csv', 
           sep = ";")

# calculate percentage of mitochondrial RNA reads

mergedLCAnucSeq[["percent.mt"]] <- 
  PercentageFeatureSet(mergedLCAnucSeq,
                       pattern = "^mt-")



QCPlotmergedsnRNAseqLCABeforeFiltering <- 
  VlnPlot(mergedLCAnucSeq,
          features = c("nFeature_RNA", "nCount_RNA","percent.mt" ),
          ncol = 3, pt.size = F) 


# filter reads

mergedLCAnucSeq <- subset(mergedLCAnucSeq, 
                      subset = nFeature_RNA > 200 &
                        nFeature_RNA < 6000 &
                        percent.mt <= 5 & 
                        nCount_RNA < 30000 &
                        typeSample == 'nucSeq')

# load expression matrix

annotation <- mergedLCAnucSeq@meta.data 

# 
ExprMatrix <- GetAssayData(mergedLCAnucSeq,
                           layer = 'counts')
#remove(mergedLCAnucSeq)
gc()

# identical(colnames(ExprMatrix), rownames(annotation))
# TRUE

# cut calculation in two steps because of vector memory exhausted

HDS <- DS_calc.func(ExprMatrix[,1:floor(ncol(ExprMatrix)/2)], 
                    DSignature = HDAG,
                    useOrder = 'mean_rank',
                    ntop = 42)

HDS <- append(HDS, 
              DS_calc.func(ExprMatrix[,floor(ncol(ExprMatrix)/2 + 1): ncol(ExprMatrix)], 
                           DSignature = HDAG,
                           useOrder = 'mean_rank',
                           ntop = 42))

identical(colnames(ExprMatrix), names(HDS))
#[1] TRUE

HDSdataframe <- data.frame('Cell_ID' = names(HDS),
                           'DS' = unname(HDS))

HDSdataframe$sample <- annotation$sample

## These are the labels for the cohorts not the diet, 
# so NAFLD cohort includes also healthy mice

HDSdataframe$cohort <- factor(annotation$orig.ident)

# renames cohort names 
levels(HDSdataframe$cohort)[levels(HDSdataframe$cohort) == "liver_cell_atlas_Nafld"] <- 'NAFLD Cohort'
levels(HDSdataframe$cohort)[levels(HDSdataframe$cohort) == "liver_cell_atlas_stst"] <- 'Standard Diet Cohort'

# NAFLD cohort has 24 and 36 weeks old mice fed either St. Diet or Western Diet

HDSdataframe$diet <- NA 

HDSdataframe[HDSdataframe$sample == "CS197" | 
               HDSdataframe$sample ==  "CS196" |
               HDSdataframe$sample == "CS193" |
               HDSdataframe$sample == "CS192" |
               HDSdataframe$sample == "CS191" |
               HDSdataframe$sample == "CS185"  |
               HDSdataframe$sample == "CS183" , 'diet'] <- 'Western Diet'

HDSdataframe$diet[is.na(HDSdataframe$diet)] <- 'Standard Diet'

HDSdataframe$diet <- factor(HDSdataframe$diet,
                                  ordered = TRUE,
                                  levels = c("Standard Diet", 
                                             "Western Diet"))

# add age

HDSdataframe$age <- NA

HDSdataframe[HDSdataframe$sample == "CS197" | 
               HDSdataframe$sample ==  "CS196" |
               HDSdataframe$sample == "CS195" |
               HDSdataframe$sample == "CS194" , 'age'] <- '36 weeks'

HDSdataframe[HDSdataframe$sample == "CS193" | 
               HDSdataframe$sample ==  "CS192" |
               HDSdataframe$sample == "CS191" |
               HDSdataframe$sample == "CS185" |
               HDSdataframe$sample == "CS183" |
               HDSdataframe$sample == "CS190" |
               HDSdataframe$sample == "CS189" |
               HDSdataframe$sample == "CS188" |
               HDSdataframe$sample == "CS184" |
               HDSdataframe$sample == "CS182", 'age'] <- '24 weeks'

HDSdataframe$age <- factor(HDSdataframe$age,
                           ordered = TRUE,
                           levels = c('24 weeks', '36 weeks'))

saveRDS(HDSdataframe, 
        file = 
          'hepatocyte-damage-score/Data/Output/weithedHDSLiverCellAtlasAllHep.rds')

## Plot results of Nafld data set 
plotNafldSampleFacetAges <- ggplot( HDSdataframe[HDSdataframe$cohort == "NAFLD Cohort",], 
                                    aes( x = diet, 
                                         y = DS,
                                         color = diet
                                    )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=12),
         plot.title = element_text(size=12),
         legend.position = 'right',
         axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5
         )) +
  guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="Diet"))+
  labs(y = "HDS", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") +
  facet_wrap(~age)

## plot distribution of data set healthy mouse (ABU samples)
plotBatchSt <- ggplot( HDSdataframe[HDSdataframe$cohort == "Standard Diet Cohort",], 
                       aes( x = diet, 
                            y = DS,
                            color = diet
                       )) + 
  geom_violin(trim = FALSE) + 
  scale_color_colorblind() +
  theme( text = element_text(size=12),
         plot.title = element_text(size=12),
         legend.position = 'none') +
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  labs(y = "HDS", x = '') + 
  stat_summary(fun.data = "mean_sdl",
               fun.args = list(mult = 1),
               geom = "pointrange") 

# save plots

panelViolinPlotsHDSAll <- (plotBatchSt | plotNafldSampleFacetAges) + 
  plot_annotation(tag_levels = 'A')

ggsave('weightedHDSViolinPlotsLiverCellAtlasBothCohorts.png', path = outputPath)


# counts of cells for reporting

table(HDSdataframe$cohort,HDSdataframe$diet, HDSdataframe$age)
# , ,  = 24 weeks
# 
# 
# Standard Diet Western Diet
# NAFLD Cohort                  1976         3987
# Standard Diet Cohort             0            0
# 
# , ,  = 36 weeks
# 
# 
# Standard Diet Western Diet
# NAFLD Cohort                  9351         4290
# Standard Diet Cohort             0            0
# 
# 

table(HDSdataframe$cohort,HDSdataframe$diet)

# Standard Diet Western Diet
# NAFLD Cohort                 11327         8277
# Standard Diet Cohort         10345            0

###########################################################
####### Part two:  Calculate values for UMAP only 24 weeks NAFLD cohort
############################################################

if(identical(HDSdataframe$Cell_ID, rownames(mergedLCAnucSeq@meta.data))){
  mergedLCAnucSeq@meta.data$age <- HDSdataframe$age
  mergedLCAnucSeq@meta.data$diet <- HDSdataframe$diet
  print(TRUE)
}
 

NAFLDcohort <- subset(mergedLCAnucSeq,
                      idents= "liver_cell_atlas_Nafld")

 
unique(NAFLDcohort@meta.data$sample)

# keep only 24 week old mice

NAFLDcohort <- subset(NAFLDcohort,
                      subset = 
                        age == '24 weeks')

remove(mergedLCAnucSeq)
gc()

# Find Variable Features
NAFLDcohort <- FindVariableFeatures(NAFLDcohort)
gc()


plotVariableFeatures <- 
  LabelPoints(plot = VariableFeaturePlot(NAFLDcohort), 
              points = head(VariableFeatures(NAFLDcohort), 10),
              repel = TRUE,
              ynudge = 0
              )
gc()

NAFLDcohort <- ScaleData(NAFLDcohort, 
                             features = rownames(NAFLDcohort))

NAFLDcohort <- RunPCA(NAFLDcohort, 
                      features = VariableFeatures(object = NAFLDcohort))

# explore Eigenvalues
head(Stdev(NAFLDcohort, reduction = "pca"), 8)


DimPlot(NAFLDcohort, reduction = 'pca', 
        group.by = 'diet', 
        shuffle = TRUE,pt.size = 0.5,
        dims = c(1,2)) 
DimPlot(NAFLDcohort, reduction = 'pca', 
        group.by = 'diet', 
        shuffle = TRUE,pt.size = 0.5,
        dims = c(3,4)) 

NAFLDcohort <- RunUMAP(NAFLDcohort, dims = c(1,2,3,4))

UMAP <- DimPlot(NAFLDcohort, reduction = 'umap', 
        group.by = 'diet', 
        shuffle = TRUE, pt.size = 0.5,
        repel = TRUE) + scale_color_colorblind() + ggtitle('Diet') +
  theme(legend.position = 'bottom',
       legend.key.size = unit(1, 'cm'),
       legend.title = element_text(size = 18),
       axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
       axis.title.y=element_blank(),
       axis.text.y=element_blank(),
       axis.ticks.y=element_blank())

# Calculate HDS 

NAFLDcohortHDS <- DS_calc.func(GetAssayData(object = NAFLDcohort, 
                                            slot = 'counts'), 
                               DSignature = HDAG,
                               useOrder = 'mean_rank',
                               ntop = 42)

# store score in Seurat object meta data for plotting 

if(identical(rownames(NAFLDcohort@meta.data), names(NAFLDcohortHDS)) == TRUE){
  NAFLDcohort@meta.data$HDS <- unname(NAFLDcohortHDS)
}else(print('Cell IDs do not match!'))

remove(NAFLDcohortHDS)
gc()

minHDS <- min(NAFLDcohort@meta.data$HDS,)
maxHDS <- max(NAFLDcohort@meta.data$HDS)


UMAP_HDS <- FeaturePlot(NAFLDcohort, 'HDS',
                        reduction = "umap",
                        pt.size = 0.5) + 
  scale_colour_distiller(palette = "YlOrBr",
                       direction = 1,
                       limits = c(minHDS, maxHDS)) + 
  ggtitle('HDS')+
  theme(legend.position = 'bottom',
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

panelUMAPs <- (UMAP | UMAP_HDS) + 
  plot_annotation(tag_levels = 'A')

panelUMAPs
ggsave(filename = 'UMAPLiverCellAtlas24WeeksNAFLDcohortHepatocytesWeightedHDSYlOrBr.png', 
       path = 'hepatocyte-damage-score/Data/Output/Results/'
)

densityPlotNoAxis <- 
  ggplot(NAFLDcohort@meta.data, 
         aes(x=HDS, color=diet)) +
  geom_density(show.legend = FALSE, 
               linewidth = 1.5) + 
  scale_color_colorblind() + theme_void()

densityPlotWithAxis <- 
  ggplot(NAFLDcohort@meta.data, 
         aes(x=HDS, color=diet)) +
  geom_density(show.legend = FALSE, 
               linewidth = 1.5) + 
  scale_color_colorblind() + theme_classic()

densityPlotNoAxis
ggsave(filename = 'DensityPlotNoAxisLiverCellAtlas24WeeksNAFLDcohortHepatocytesWeightedHDS.png', 
       path = 'hepatocyte-damage-score/Data/Output/Results/'
)

densityPlotWithAxis
ggsave(filename = 'DensityPlotWithAxisLiverCellAtlas24WeeksNAFLDcohortHepatocytesWeightedHDS.png', 
       path = 'hepatocyte-damage-score/Data/Output/Results/'
)





