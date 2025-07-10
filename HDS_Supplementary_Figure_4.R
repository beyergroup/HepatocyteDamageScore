# Supplementary Figure 4 

# Prepare libraries and data for spatial plots 
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

  # Load data 
  spatialObjectMergedFiltered <- readRDS(
    file = 'hepatocyte-damage-score/Data/Output/LiverCellAtlasMouseSpatialMergedFilteredSCT.rds')
  
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
  
  
  HDSmergedfiltered$orig.ident <- 
    spatialObjectMergedFiltered@meta.data$orig.ident

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

# Panel A:
# HDS violin plots calculated on murine spatial RNA-seq data from Liver Cell Atlas   
{    
    maxHDS <- -0.09
    minHDS <- -0.37
    
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
      scale_color_colorblind() 

        print(spatialLCA_ViolinPlots)
    
 
  }
  
# Panel B: Distribution of HDS on tissue, samples not shown in Figure 1 
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
    

    print(layerHDSPlotSup/layerTissuePlotSup)
  

}

# Panel C: Calculate HDS on proteomics data Williams et al. 
{
  library(stringr)
  library(AUCell)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  
  path_to_data <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/williams_proteomics/"
  data <- read.csv(paste0(path_to_data,"aTableS2_LIVER_DATA_Small.csv"), sep = ",")
  
  metadata <- data[1:52,]
  sum(str_detect(data$InjectionOrder, "Prot_*"))
  data_prot <- data[str_detect(data$InjectionOrder, "Prot_*"),]
  
  #> dim(data_prot)
  #[1] 3968  630
  rownames(data_prot) <- data_prot$InjectionOrder.1
  data_prot <- data_prot[-c(1:3)]
  data_prot_col <- str_detect(colnames(data_prot), "Injection*")
  sum(data_prot_col)
  # 315 
  # select the colums with protein quantification data 
  data_prot <- data_prot[,data_prot_col]
  metadata_prot <- metadata[,data_prot_col]
  first_row_protein <- data_prot[1,]
  data_prot <- data_prot[-1,]
  proteins <- rownames(data_prot)
  data_prot <- data.frame(lapply(data_prot,as.numeric))
  rownames(data_prot) <- proteins
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  # load functions to calculate HDS
  source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')
  
  
  HDS <- DS_calc.func(exprMatrices = data_prot, 
                      DSignature = HDAG )
  
  rownames_metadata <- colnames(metadata_prot)
  metadata_prot <- t(metadata_prot) 
  colnames(metadata_prot) <- metadata_prot[1,]
  metadata_prot <- as.data.frame(metadata_prot)
  rownames(metadata_prot) <- rownames_metadata
  metadata_prot <- metadata_prot[-(1:3),]
  View(metadata)
  
  length(unique(metadata_prot$OmicsEarTag))
  
  unique(metadata_prot$Diet)
  unique(metadata_prot$Age)
  
  metadata_prot$HDS <- as.numeric(unlist(HDS))
  metadata_prot$Diet <- factor(metadata_prot$Diet, levels = c("CD", "HF"), ordered = TRUE)
  
  
  wilcox.test(x =  metadata_prot$HDS[metadata_prot$Diet == "HF"], y = metadata_prot$HDS[metadata_prot$Diet == "CD"],  alternative = "greater")
  
  ggplot(metadata_prot, aes(x = Diet, y = HDS, color = Diet)) + geom_boxplot(lwd = 1.5) + scale_color_colorblind() + theme_classic() + 
    theme( text = element_text(size = 18), axis.title.x = element_blank()) 
  
}