
# Figure 1 

# Data for Panel E 
#load functions our functions and libraries
{
  source('SharedFunctions.R')
  library(Seurat)
  library(AUCell)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(cowplot)
  library(reshape2)
  library(plyr)
  library(RColorBrewer)
  library(dplyr)
  library(scales)
  
  
  # functions for ploting statistics
  median_IQR <- function(x) {
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  stat_box_data <- function(y) {
    return( 
      data.frame(
        y = 0.1,  #may need to modify this depending on your data
        label = paste('count =', length(y), '\n',
                      'median =', round(median(y), digits = 2), '\n')
      )
    )
  }
  
  
}

# Load and process Liver Cell Atlas data and calculate HDS 
{
  outputPath <- 'hepatocyte-damage-score/Data/Output/Results/'
  
  # Hepatocyte Universal Damage Signature
  HDAG <- read.csv(file = 
                     'hepatocyte-damage-score/Data/Output/HDAG.csv')
  
  #load  merged Seurat object, output from: merge_snRNAseq_LiverCellAtlas_murine.R
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
  
  # filter cells by quality control metrics
  mergedLCAnucSeq <- subset(mergedLCAnucSeq, 
                            subset = nFeature_RNA > 200 &
                              nFeature_RNA < 6000 &
                              percent.mt <= 5 & 
                              nCount_RNA < 30000 &
                              typeSample == 'nucSeq')
  
  annotation <- mergedLCAnucSeq@meta.data 
  
  # create expression matrix for calculating HDS
  ExprMatrix <- GetAssayData(mergedLCAnucSeq,
                             layer = 'counts')
  gc()
  # sanity check
  identical(colnames(ExprMatrix), rownames(annotation))
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
  
  # sanity check
  identical(colnames(ExprMatrix), names(HDS))
  #[1] TRUE
  
  HDSdataframe <- data.frame('Cell_ID' = names(HDS),
                             'HDS' = unname(HDS))
  
  HDSdataframe$sample <- annotation$sample
  
  ## These are the labels for the cohorts not the diet, 
  # so NAFLD cohort includes also healthy mice
  
  HDSdataframe$cohort <- factor(annotation$orig.ident)
  
  # renames cohort names 
  levels(HDSdataframe$cohort)[levels(HDSdataframe$cohort) == "liver_cell_atlas_Nafld"] <- 'NAFLD Cohort'
  levels(HDSdataframe$cohort)[levels(HDSdataframe$cohort) == "liver_cell_atlas_stst"] <- 'Standard Diet Cohort'
  
  # NAFLD cohort has 24 and 36 weeks old mice fed either St. Diet or Western Diet
  HDSdataframe$diet <- NA 
  # adding sample diet manually
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
  
 # add data to seurat object meta dats for plotting 
  if(identical(rownames(mergedLCAnucSeq@meta.data), HDSdataframe$Cell_ID)){
    mergedLCAnucSeq@meta.data$diet <- HDSdataframe$diet
    mergedLCAnucSeq@meta.data$age <- HDSdataframe$age
    mergedLCAnucSeq@meta.data$cohort <- HDSdataframe$cohort
    mergedLCAnucSeq@meta.data$HDS <- HDSdataframe$HDS
  }
  
  # Find Variable Features
  NAFLDcohort24weeks <- subset(mergedLCAnucSeq, subset = cohort == "NAFLD Cohort" & age == '24 weeks') 
  NAFLDcohort24weeks <-  FindVariableFeatures(NAFLDcohort24weeks)
  NAFLDcohort24weeks  <- ScaleData(NAFLDcohort24weeks, 
                           features = rownames(NAFLDcohort24weeks ))
  gc()
  NAFLDcohort24weeks <- RunPCA(NAFLDcohort24weeks, 
                        features = VariableFeatures(object = NAFLDcohort24weeks))
  NAFLDcohort24weeks <- RunUMAP(NAFLDcohort24weeks, dims = c(1,2,3,4))

}

# Plot Panel E
{
  colorPaletteHDS <- brewer.pal(n = 9, name = 'YlOrBr')
  minHDS <- min(NAFLDcohort24weeks@meta.data$HDS)
  maxHDS <- max(NAFLDcohort24weeks@meta.data$HDS)
  
  # NAFLD Cohort, 24-week-old mice
  gg1 <- FeaturePlot(NAFLDcohort24weeks,
                     features = "HDS",
                     cells = sample(rownames(NAFLDcohort24weeks@meta.data)),
                     pt.size = 1.5
                     ) + 
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           limits = c(minHDS, maxHDS)) +
    theme(legend.position = "bottom", 
          legend.key.height=unit(5,"mm"),
          text = element_text(size=20),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) 

  
  gg2 <- DimPlot( NAFLDcohort24weeks , 
                  group.by = "diet",
                  pt.size = 1.5 , 
                  cells = sample(rownames(NAFLDcohort24weeks@meta.data))) +
    scale_colour_colorblind()+
    theme(legend.position = "bottom",
          legend.key.height = unit(15,"mm"),
          text = element_text(size = 20),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  
  ## densito plot for PDS
  gg3 <- ggplot(data = NAFLDcohort24weeks@meta.data,
                aes(x = HDS, color = diet))+
    geom_density(lwd = 3)+ theme_bw()+
    scale_colour_colorblind()+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=20))
  
  ggpnl <-cowplot::plot_grid(plotlist = list(gg2, gg1, gg3), 
                             rel_heights = c(1,0.25), nrow = 2)
  
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "FiguresManuscript/LCA_NAFLDcohort24weeks.snRNAseq_HDS.umap.NOdim.pdf"))
  ggpnl
  dev.off()
  
  png(height = 6, width = 8,
      units = 'in',
      res = 1200,
      file = paste0(outputPath,
                    "FiguresManuscript/PNGs/LCA_NAFLDcohort24weeks.snRNAseq_HDS.umap.NOdim.png"))
  ggpnl
  dev.off() 
}

# Data Panel F: (part of this boxplot has been put in Supplementary Figure 2 C)
{ 
  # Output of cross validation script
  FirstGroupHDSListDataFrames <- 
    read.csv(file = paste0(outputPath,
                           'LeaveOneOutCrossValidationFirstGroupWeightedHDSForPlotting.csv'),
             sep = ';',
             dec = '.', stringsAsFactors = TRUE)
  FirstGroupHDSListDataFrames <- 
    FirstGroupHDSListDataFrames[FirstGroupHDSListDataFrames$conditionGroupPlotting != 'not included',]
  
  cvPlot1 <- ggplot( FirstGroupHDSListDataFrames, 
                   aes( x = conditionGroupPlotting ,
                        y = score , 
                        color = conditionGroupPlotting)) +
    scale_color_colorblind() + 
    geom_boxplot(lwd = 1.5, outlier.shape = NA) + 
    geom_jitter(size = 5,  width = 0.2)+
    theme_bw() + theme( text = element_text(size = 24) , 
                        axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.title=element_blank()) +
    stat_compare_means( size = 6 ,label = 'p.format') +
    facet_grid(cols = vars(model), scales = "free") 
  
  pdf( height = 6, width = 14, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_cv.models_firstGroup.pdf") )
  print( cvPlot1 )
  dev.off()
  
  png( height = 600, width = 1400, file = paste0(outputPath,
                                             "FiguresManuscript/PNGs/HDS_cv.models_firstGroup.png"),
       units = 'px')
  print( cvPlot1 )
  dev.off()
  
  # sample number for figure description
  table(FirstGroupHDSListDataFrames$model,FirstGroupHDSListDataFrames$conditionGroupPlotting )
  
}

#Preparing Data for Panel G 
{

  # Load data 
  spatialObjectMergedFiltered <- readRDS(
    file =
      'hepatocyte-damage-score/Data/Output/LiverCellAtlasMouseSpatialMergedFilteredSCT.rds'
    )
  
  ####################
  # Apply HDS
  ####################
  
  # Load HUDS (Hepatocyte Universal Damage Score)
  geneListHDAG <- read.csv(
    'hepatocyte-damage-score/Data/Output/HDAG.csv')
  
  # extract expression matrix 
  exprMatrix <- GetAssayData(spatialObjectMergedFiltered,
                             assay = 'LiverCellAtlasSpatial')
  
  HDSmergedfiltered <- DS_calc.func( exprMatrices = exprMatrix, 
                                     ceilThrsh = 0.05 , 
                                     DSignature = geneListHDAG, 
                                     ntop = 42 , 
                                     useOrder = 'mean_rank')
  
  HDSmergedfiltered <- data.frame('cellID' = names(HDSmergedfiltered),
                                  'HDS' = HDSmergedfiltered )
  
  # add HDS to SueratObjectMetadata for plotting 
  if (identical(row.names(spatialObjectMergedFiltered@meta.data), 
            HDSmergedfiltered$cellID) == TRUE) {
    spatialObjectMergedFiltered@meta.data$HDS <- HDSmergedfiltered$HDS
    }

  ############
  # Plot HDS #
  ############
  
  # format and center data for box-plots
  
  HDSmergedfiltered$orig.ident <- 
    spatialObjectMergedFiltered@meta.data$orig.ident
  # 
  # HDSmergedfiltered$centeredHDS <- scale(HDSmergedfiltered$HDS, 
  #                                        scale = FALSE)
  HDSmergedfiltered$zonationGroup <- 
    spatialObjectMergedFiltered@meta.data$zonationGroup
  
  HDSmergedfiltered$diet <- NA
  HDSmergedfiltered[HDSmergedfiltered$orig.ident == 'StSpatialSeuratM1S1', 'diet'] <- 
    'Standard Diet' 
  HDSmergedfiltered[HDSmergedfiltered$orig.ident != 'StSpatialSeuratM1S1', 'diet'] <- 
    'Western Diet'
  HDSmergedfiltered$diet <- factor(HDSmergedfiltered$diet,
                                   ordered = TRUE, 
                                   levels = c("Standard Diet", "Western Diet"))
  
  if(identical(HDSmergedfiltered$cellID, 
               rownames(spatialObjectMergedFiltered@meta.data))){
    print('yes')
    spatialObjectMergedFiltered@meta.data$diet <- HDSmergedfiltered$diet
  }
  
  
  

}

#Panel G
{
  maxHDS <- -0.09
  minHDS <- -0.37
  
  p3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS",
                           images = 'WDMouse002Sample001',
                           pt.size.factor = 1,
                           crop = FALSE,
                           image.alpha = 1) +
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(minHDS,maxHDS),
                         name = 'HDS', 
                         oob = squish) +
    theme(legend.key.size = unit(2, 'cm'),
          text = element_text(size = 27))
  
  p3Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                          images = 'WDMouse002Sample001',
                          alpha = 0, 
                          crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
 
  
  pdf( height = 15, width = 30, 
       file = paste0(outputPath,
                     "FiguresManuscript/HDS_LCAspatialRNAseq_WDM2S1.pdf") )
  print((wrap_plots(p3 | p3Tissue)))
  dev.off()

  
  
}



