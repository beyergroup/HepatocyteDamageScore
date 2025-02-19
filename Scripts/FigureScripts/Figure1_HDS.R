
# Figure 1 

# Data for Figure 1 E load functions our functions and libraries
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

  #load  merged Seurat object 
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

# Plot Figure 1 E
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

# Data Figure 1 F
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

# Figure 1 F
{
  SecondGroupHDS <- 
    read.csv(file = paste0(outputPath,
                           'LeaveOneOutCrossValidationSecondGroupWeightedHDSForPlotting.csv'),
             sep = ';',
             dec = '.', 
             stringsAsFactors = TRUE)
  
  SecondGroupHDS <- 
    SecondGroupHDS[SecondGroupHDS$conditionGroupPlotting != 'not included',]
  
  cvPlot2 <- ggplot(SecondGroupHDS, 
                  aes(x = conditionGroupPlotting , 
                      y = score , 
                      color = conditionGroupPlotting)) +
    scale_color_colorblind() + 
    geom_boxplot(lwd=1.5, outlier.shape = NA) + 
    geom_jitter(size = 5, width = 0.2)+
    theme_bw() + theme( text = element_text(size = 24) , 
                        axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.title=element_blank()) +
    stat_compare_means( size = 6 ,label = 'p.format') +
    facet_grid(cols = vars(model), scales = "free") 
  
  pdf( height = 6, width = 12, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_cv.models_SecondGroup.pdf") )
  print( cvPlot2 )
  dev.off()
}

# Preparing Data for Figures 1 G and Supplementary Figures

{

  # READ
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

# Figure G

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

# Supplementary Figure 1 D (rest of Liver Cell Atlas samples)
{
  HDS_LCAsamples <- na.omit(mergedLCAnucSeq@meta.data)
  
 violinLCAsamples <-ggplot(data = HDS_LCAsamples,
                           aes(x = age,
                               y = HDS,
                               color = diet)) +
   geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    ylab("HDS") + 
   ylim(c(-0.32, 0.15)) +
   scale_color_colorblind() +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR,
                 size = 1.5 ,
                 show.legend = FALSE,
                 position = position_dodge(0.95)) +
    stat_summary( fun.y = "median",  
                  size = 1, 
                  geom = "crossbar",
                  show.legend = FALSE,
                  width = 0.2,
                  position = position_dodge(0.95)) +
    theme_classic() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, 
                                      hjust = 1),
           legend.title=element_blank()) +
   stat_summary(
     fun.data = stat_box_data, 
     geom = "text",
     position = position_dodge(1),
     size = 5
     )
 
 pdf( height = 8, width = 12, 
      file = paste0(outputPath,
                    "FiguresManuscript/HDS_LCA_allsamples.pdf") )
 print(violinLCAsamples)
 dev.off()
   
 
}



# Supplementary Figure 1 E (cross-validation)

{
  # plot acetominophen on it's own because it's a time series: GSE111828 
  
  AcetaminophenGroupHDS <- 
    read.csv(file = paste0(outputPath,
                           'LeaveOneOutCrossValidationAcetaminophen.csv'),
             sep = ',',
             dec = '.', 
             stringsAsFactors = TRUE)[,-1]
  
  colorsAcetaminophen <-c('0hr' = 'black', '12hr' = '#E69F00', '24hr' = '#E69F00', 
                          '36hr' = '#E69F00', '48hr' = '#E69F00', '72hr' = '#E69F00')
  
  
  cvAcetaminophen <- ggplot(AcetaminophenGroupHDS, 
                            aes(x = condition,
                                y = score,
                                color = condition)) + 
    geom_boxplot(lwd = 1.5, outlier.shape = NA) +
    geom_jitter(size = 5, width = 0.1) + 
    theme_bw() + 
    scale_color_manual(values = colorsAcetaminophen) +
    theme( text = element_text(size = 24) , 
           axis.title.x=element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           legend.title=element_blank()) + 
    stat_compare_means(method = "kruskal.test", 
                       size = 6, 
                       label = 'p.format') + # Add global p-value
    labs(y = "HDS") 
  
  pdf( height = 6, width = 12, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_cv.modelAcetaminophen.pdf") )
  print(cvAcetaminophen)
  dev.off()
  
  table(AcetaminophenGroupHDS$model, AcetaminophenGroupHDS$condition)
  # n = 4 
  
  # Tabula Muris 
  
  AgingGroupHDS <- 
    read.csv(file = paste0(outputPath,
                           'LeaveOneOutCrossValidationAging.csv'),
             sep = ',',
             dec = '.', 
             stringsAsFactors = TRUE)[,-1]
  
  cvAgingVersion1 <- ggplot(AgingGroupHDS,
                            aes(x = factor(condition),
                                y = score)) + 
    geom_boxplot(lwd = 1.5, outlier.shape = NA) +
    geom_jitter(size = 5, width = 0.1) + 
    theme_bw() +
    theme( text = element_text(size = 24) , 
           axis.title.x=element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           legend.title=element_blank()) + 
    stat_compare_means(method = "kruskal.test", 
                       size = 8, 
                       label = 'p.format') + # Add global p-value
    labs(y = "HDS")
  
  pdf( height = 6, width = 12, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_cv.modelAgingKruskalBoxplots.pdf") )
  print(cvAgingVersion1)
  dev.off()
  
  
  cvAgingVersion2 <- ggplot(AgingGroupHDS, 
                            aes(x = as.numeric(condition),
                                y = score)) + 
    geom_jitter(size = 5, width = 0.1) + 
    theme_bw() +
    theme( text = element_text(size = 24) , 
           axis.title.x=element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           legend.title=element_blank()) +
    labs(x = "age in months" , y = "HDS") + 
    stat_smooth(method = 'lm', 
                linetype = 'dashed', 
                color= '#E69F00', 
                lwd = 1.5) + 
    stat_cor(method = 'spearman', 
             size = 8, 
             cor.coef.name = "rho") + 
    scale_x_continuous(breaks = c(1, seq(3,27,3)))
  
  pdf( height = 6, width = 12, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_cv.modelAgingSpearmanCorLM.pdf") )
  print(cvAgingVersion2)
  dev.off()
  
}

# Supplementary Figure 1 F (size & robustness to permutation)
{
  
  # Signature Size Test:
  DSsize <- readRDS('hepatocyte-damage-score/Data/Output/testSizesmouseMergedHepatocytesnucSeqQCfitleredSCT.rds')
  
  size_vec <- c(5, 10, 15, 20, 30 , 42, 50, 100, 200, 300, 500, nrow(HDAG))
  
  DSsizeListDataFrames <- 
    lapply( seq(DSsize) , function(ii){
      datt <- cbind.data.frame( score = DSsize[[ii]] ,
                                setSize = as.factor(size_vec[ii]), 
                                condition = as.factor(
                                  NAFLDcohort24weeks@meta.data$diet[
                                    match(names(DSsize[[ii]]),
                                          rownames(NAFLDcohort24weeks@meta.data))
                                    ]
                                  )
      )
      
      return(datt)
    })
  
  # y axis limits includes most extreme values 
  minMax <- c(min(DSsize[[11]]), max(DSsize[[12]]))
  
  DSsizeListViolin <- lapply( 
    seq(DSsizeListDataFrames) , function(ii){
      ppg <- ggplot( DSsizeListDataFrames[[ii]] , 
                     aes(x = condition, 
                         y = score,
                         color = condition)) + 
        geom_violin(trim = TRUE,
                    lwd = 1.5) + 
        theme_classic()+
        scale_colour_colorblind() +
         theme(
           legend.position = "none",
           legend.title = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           text = element_text(size=20)
        )  +
        ggtitle(paste(DSsizeListDataFrames[[ii]]$setSize[1], "Genes")) +
        guides(colour = guide_legend(override.aes = aes(label = ""))) +
        labs(y = "HDS") + 
        stat_summary(geom = "linerange",
                     fun.data = median_IQR,
                     size = 1.5 ,
                     show.legend = FALSE) +
        stat_summary( fun.y = "median",  
                      size = 1, 
                      geom = "crossbar",
                      show.legend = FALSE,
                      width = 0.2) +
        ylim(minMax)
      return( ppg )
    }
  )
  
  
  SizePannelFigure <- (DSsizeListViolin[[1]] | DSsizeListViolin[[2]] | DSsizeListViolin[[3]] | DSsizeListViolin[[4]] | DSsizeListViolin[[5]] | DSsizeListViolin[[6]] ) / ( DSsizeListViolin[[7]] | DSsizeListViolin[[8]] | DSsizeListViolin[[9]] | DSsizeListViolin[[10]] | DSsizeListViolin[[11]]  | DSsizeListViolin[[12]] + theme(legend.position = 'bottom' ) ) 
  
  pdf(height = 15,
      width = 20, 
       file = paste0(outputPath,
                     "FiguresManuscript/HDS_SizeTest_LCANAFLDcohort24weeks.pdf") 
      )
  print(SizePannelFigure)
  dev.off()
  
  pdf(height = 15,
      width = 20, 
      file = paste0(outputPath,
                    "FiguresManuscript/HDS_SizeTest_LCANAFLDcohort24weeks.pdf") 
  )
  print(SizePannelFigure)
  dev.off()
  
 # permutation plots
  permuResults <- readRDS(
    file = paste0(outputPath,
                  'RandomizationTest/RandomizationWeightedHDSLCAmergedNuq24WeeksNAFLDCohort.rds'))
  
  # Format Results for better plotting
  toPlot <- Reduce(rbind, lapply(permuResults ,function(X) {
    X$value <- scale(X$value, scale = FALSE)
    return(X)
  }))
  
  toPlot$HDS <- toPlot$value
  toPlot <- toPlot[,c(-2, -3)]
  
  # add cohort and diet annotation 
  toPlot$diet <- as.factor(
    NAFLDcohort24weeks@meta.data$diet[match(toPlot$cell,
                          rownames(NAFLDcohort24weeks@meta.data))])
  
  toPlot$cohort <- as.factor(
    NAFLDcohort24weeks@meta.data$cohort[match(toPlot$cell,
                                            rownames(NAFLDcohort24weeks@meta.data))])
  # plotting 
  permViolin <- ggplot(toPlot, 
                aes(x = random.percent, 
                    y = HDS, 
                    color = diet)) +
    geom_violin(trim = TRUE, lwd = 1.5) +
    theme_bw() +
    scale_colour_colorblind() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size=20)
    ) +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR,
                 size = 1.5 ,
                 show.legend = FALSE,
                 position = position_dodge(0.95)) +
    stat_summary( fun.y = "median",  
                  size = 1, 
                  geom = "crossbar",
                  show.legend = FALSE,
                  width = 0.3,
                  position = position_dodge(0.95)) +
    labs(x = "percent of studies randomized [%]")
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "FiguresManuscript/permutationTestHDS_LCAsnRNAseq_24weekNAFLDcohort.pdf") )
  print(permViolin)
  dev.off()
  

}

# Supplementary Figure 1 G (HDS spatial violin plots and rest of tissue sections)

{
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS", 
                           images = 'StMouse001Sample001', 
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
  
  p1Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                          images = 'StMouse001Sample001', 
                          alpha = 0, crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
  
  p2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS", 
                           images = 'WDMouse001Sample001', 
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
  
  p2Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                          images = 'WDMouse001Sample001',
                          alpha = 0, 
                          crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
  
  
  p4 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS",
                           images = 'WDMouse002Sample002',
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
  
  p4Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                          images = 'WDMouse002Sample002',
                          alpha = 0, 
                          crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
  
  
  layerHDSPlotSup <- wrap_plots(p1 | p2 | p4)
  layerTissuePlotSup<- wrap_plots(p1Tissue | p2Tissue | p4Tissue)
  
  pdf( height = 30, width = 30, 
       file = paste0(outputPath,
                     "FiguresManuscript/HDS_LCAspatialRNAseq_stM1_WDM1S1_WDM2S2.pdf") )
  print(layerHDSPlotSup/layerTissuePlotSup)
  dev.off()
  
  # In case they ask me for boxplots instead of violin plots 
  # spatialLCA_Boxplot <- ggplot(data = HDSmergedfiltered,
  #                         aes(x = orig.ident, y = centeredHDS,
  #                             color = diet)) +
  #   geom_boxplot(lwd = 1.5, outlier.size = NA) + 
  #   theme_bw() + 
  #   theme( text = element_text(size = 24) ,
  #          axis.title.x = element_blank(),
  #          axis.text.x = element_text(angle = 45, hjust = 1),
  #          legend.title=element_blank()) +
  #   scale_color_colorblind() +
  #   stat_compare_means(method = "kruskal.test",
  #                      ref.group = 'StSpatialSeuratM1S1',
  #                      paired = FALSE,
  #                      label =  "p.format",
  #                      size = 10)
  
  spatialLCA_ViolinPlots <-ggplot(data = HDSmergedfiltered,
                                  aes(x = orig.ident, 
                                      y = HDS,
                                      color = diet)) +
    geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    xlab("samples") + 
    ylab("HDS") + 
    stat_summary(geom = "linerange",
                 fun.data = median_IQR,
                 size = 1.5 ,
                 show.legend = FALSE) +
    stat_summary( fun.y = "median",  
                  size = 1, 
                  geom = "crossbar",
                  show.legend = FALSE,
                  width = 0.2) +
    scale_x_discrete(labels = c("StD \n M1S1", "WD \n M1S1",
                                "WD \n M2S1", "WD \n M2S2")) +
    theme_bw() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           legend.title=element_blank()) + 
    scale_color_colorblind() + 
    stat_compare_means(method = 'wilcox.test',
                       ref.group = 'StSpatialSeuratM1S1',
                       paired = FALSE,
                       label =  "p.signif",
                       label.y = -0.03,
                       size = 10,
                       symnum.args = list(
                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                         symbols = c("****", "***", "**", "*", "ns"))) 
  
  pdf( height = 8, width = 10, file = paste0(outputPath,
                                             "FiguresManuscript/HDS_spatialRNAseqLCA_ViolinWilcoxTestIQRMedian.pdf") )
  print(spatialLCA_ViolinPlots)
  dev.off()
}


