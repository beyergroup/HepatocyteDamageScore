### Applying HDS on Spatial Transcriptomics ###
library(Seurat)
library(hdf5r)
options(Seurat.object.assay.version = "v5")
library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(gridExtra)
library(AUCell)
library(GSEABase)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(viridis)
library(colorBlindness)

# functions
source(file='SharedFunctions.R')
# options(future.globals.maxSize = 1e+09)

#############################
# 1. Liver Cell Atlas 
# Visium 36 Weeks Mouse Livers
#############################

{
  ###########
  # Load data 
  ###########
  
  toData <- 'hepatocyte-damage-score/Data/Input/SpatialTranscriptomicData/LiverCellAtlasVisium/'
  dataPath <- paste0(toData,'GSE192742_36WeeksMouseLiver/')
  plotsPath <- 'hepatocyte-damage-score/Data/Output/Results/'
    
  
  # GEO data info: 
  #vstrain: C57Bl/6
  # platform: 10x Visium
  # condition: SD, number of added abs: 0, number of spots: 1279
  # Raw sequencing data for each sample was converted 
  # to matrices of expression counts using the Space Ranger software provided by 10X Genomics (version 1.0).
  # Genome_build: mm10 or hg19
  
  # Sample Info
  # StMouse001Sample001 -> shortfilename: JBO6
  # WDMouse001Sample001 -> shortfilename: JBO9
  # WDMouse002Sample001 -> shortfilename: JBO10
  # WDMouse002Sample002 -> shortfilename: JBO12
  
  # there is another batch of 4 standard diet mouse livers -> GSE192741
  # this are al
  

  annotationNAFLD <- read.csv(file = paste0(dataPath,"annot_mouseNafldVisium.csv"))

 
  spatialSamples <- list(
    Load10X_Spatial(data.dir = paste0(dataPath,'StSt/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "StMouse001Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "StSt/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse001Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse001Sample001", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse001Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample001/'),
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "LiverCellAtlasSpatial",
                    slice = "WDMouse002Sample001",
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample001/"))),
    Load10X_Spatial(data.dir = paste0(dataPath,'WD/Mouse002Sample002/'),
                    filename = "filtered_feature_bc_matrix.h5", 
                    assay = "LiverCellAtlasSpatial", 
                    slice = "WDMouse002Sample002", 
                    image = Read10X_Image(image.dir = paste0(
                      dataPath, "WD/Mouse002Sample002/")))
    )
  
  sampleNames <- c("StSpatialSeuratM1S1",
                   "WDSpatialSeuratM1S1",
                   "WDSpatialSeuratM2S1", 
                   "WDSpatialSeuratM2S2")
  
  names(spatialSamples) <- sampleNames

  sampleCoordinates <- list(
    spatialSamples[[1]]@images$StMouse001Sample001@coordinates,
    spatialSamples[[2]]@images$WDMouse001Sample001@coordinates,
    spatialSamples[[3]]@images$WDMouse002Sample001@coordinates,
    spatialSamples[[4]]@images$WDMouse002Sample002@coordinates
    )
  
  ### Functions from Liver Cell Atlas Git Hub
  ##### Function drawVlnPlot
  
  drawVlnPlot <- function(toPlot, colsToColor){
    toPlot$staticNr=1
    
    toPlot<-toPlot[order(toPlot[,colsToColor[1]]),]
    p_nGene <- ggplot(toPlot, aes(staticNr, nFeature_LiverCellAtlasSpatial)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    toPlot<-toPlot[order(toPlot[,colsToColor[2]]),]
    p_nUMI <- ggplot(toPlot, aes(staticNr, nCount_LiverCellAtlasSpatial)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    toPlot<-toPlot[order(toPlot[,colsToColor[3]]),]
    p_mito <- ggplot(toPlot, aes(staticNr, percent.mt)) + 
      geom_violin(fill="gray80") + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic() +
      theme(legend.position = "none", axis.title.x = element_blank())
    
    return(grid.arrange(p_nGene,p_nUMI,p_mito,ncol=3))
  }
  
  
  ##### Function drawVlnPlot_split
  drawVlnPlot_split <- function(toPlot, colsToColor){
    toPlot$staticNr=1
    
    columnNr<-which(colnames(toPlot)==colsToColor[1])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_nGene <- ggplot(toPlot, aes(orig.ident, nFeature_LiverCellAtlasSpatial, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
      geom_violin() + 
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    columnNr<-which(colnames(toPlot)==colsToColor[2])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_nUMI <- ggplot(toPlot, aes(orig.ident, nCount_LiverCellAtlasSpatial, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
      geom_violin() + 
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    columnNr<-which(colnames(toPlot)==colsToColor[3])
    toPlot<-toPlot[order(toPlot[,columnNr]),]
    p_mito <- ggplot(toPlot, aes(orig.ident, percent.mt, fill=orig.ident)) + 
      geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
      geom_violin() +
      scale_color_manual(values=c("#00BFC4", "#F8766D")) +
      theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
            axis.text.y = element_text(size = 8), axis.title.x = element_blank())
    
    thePlot<-grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
    return(thePlot)
  }
  
  ########################################  
  #### Check for outlier sports (from liver cell atlas git hub)
  ########################################
  
  ##### Sample: StSpatialSeurateM1S1 #####
  
  # p1 <- SpatialDimPlot(spatialSamples[[1]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[1]])+
  #   geom_point(size=0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol=2)
  # 
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #        file='/PATH/TO/results/image_sample1.png', 
  #        width=10, height = 5)
  # 
  # Select outlier spots
  badSpots <- sampleCoordinates[[1]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagecol < 5000 | imagerow < 4500) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[1]]$badSpot <- FALSE
  sampleCoordinates[[1]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[1]]) + 
  #   geom_point(size=0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  
  # Remove bad spots
  spatialSamples[[1]] <- spatialSamples[[1]][,setdiff(colnames(spatialSamples[[1]]),badSpots)]
  dim(spatialSamples[[1]])
  # 31053  1274
  
  # p3 <- SpatialDimPlot(spatialSamples[[1]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  
  
  ##### Sample: WDSpatialSeuratM1S1 #####

  # p1 <- SpatialDimPlot(spatialSamples[[2]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[2]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  
  # grid.arrange(p1, p2, ncol=2)
  
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # No outlier spots
  
  #### Sample: WDSpatialSeuratM2S1 ####
  # 
  # p1 <- SpatialDimPlot(spatialSamples[[3]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[3]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol = 2)
  
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # Select outlier spots
  badSpots <- sampleCoordinates[[3]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagerow > 16300) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[3]]$badSpot <- FALSE
  sampleCoordinates[[3]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[3]]) + 
  #   geom_point(size=0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  
  # Remove bad spots
  spatialSamples[[3]] <- spatialSamples[[3]][,setdiff(colnames(spatialSamples[[3]]),badSpots)]
  dim(spatialSamples[[3]])
  
  # 31053  2074
  
  #p3 <- SpatialDimPlot(spatialSamples[[3]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  

  
  #### Sample: WDSpatialSeuratM2S2 ####
  
  # p1 <- SpatialDimPlot(spatialSamples[[4]])
  # p2 <- ggplot(aes(x=imagecol, y=imagerow), 
  #              data = sampleCoordinates[[4]])+
  #   geom_point(size = 0.5) +
  #   theme_classic() +
  #   scale_y_reverse()
  # 
  # grid.arrange(p1, p2, ncol = 2)
  # 
  # ggsave(grid.arrange(p1, p2, ncol = 2), 
  #       file='/PATH/TO/results/image_sample1.png', 
  #       width=10, height = 5)
  
  # Select outlier spots
  badSpots <- sampleCoordinates[[4]] %>% mutate('spot' = rownames(.)) %>% 
    dplyr::filter(imagecol < 5000 | imagerow > 16000) %>%
    pull(spot)
  
  # Visualize bad spots
  sampleCoordinates[[4]]$badSpot <- FALSE
  sampleCoordinates[[4]][badSpots,'badSpot'] <- TRUE
  # p2 <- ggplot(aes(x = imagecol, y = imagerow, color = badSpot), 
  #              data = sampleCoordinates[[4]]) + 
  #   geom_point(size = 0.5) +
  #   theme_classic() + 
  #   scale_y_reverse()
  # 
  # grid.arrange( p1, p2, ncol = 2)
  # 
   # Remove bad spots
  spatialSamples[[4]] <- 
    spatialSamples[[4]][, setdiff(colnames(spatialSamples[[4]]),badSpots)]
  
  dim(spatialSamples[[4]])
  # 31053  2252
  
  #p3 <- SpatialDimPlot(spatialSamples[[4]])
  # ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
  #       width=15, height = 5)
  
  
  ########################################
  # Merge 
  ########################################
  
  spatialSamples[[1]]@meta.data$orig.ident <- sampleNames[1]
  spatialSamples[[2]]@meta.data$orig.ident <- sampleNames[2]
  spatialSamples[[3]]@meta.data$orig.ident <- sampleNames[3]
  spatialSamples[[4]]@meta.data$orig.ident <- sampleNames[4]
  


  spatialObjectMerged <- merge(spatialSamples[[1]], y = c(spatialSamples[[2]], 
                                                      spatialSamples[[3]],
                                                      spatialSamples[[4]]), 
                               add.cell.ids = sampleNames, merge.data = TRUE)
  
  # checking merge
  dim(spatialObjectMerged)
  head(spatialObjectMerged@meta.data)
  tail(spatialObjectMerged@meta.data)
  table(spatialObjectMerged@meta.data$orig.ident)
  
  # StSpatialSeuratM1S1 WDSpatialSeuratM1S1 WDSpatialSeuratM2S1 WDSpatialSeuratM2S2 
  # 1274                2348                2074                2252 
  
  
  #########################################
  # Filtering
  #########################################
  
  spatialObjectMerged[["percent.mt"]] <- 
    PercentageFeatureSet(spatialObjectMerged, pattern = "^mt-")
  spatialObjectMerged
  
  # Select spots to be removed: kept same thresholds as Liver Cell Atlas git hub
  
  toPlot <- spatialObjectMerged@meta.data
  toPlot$badCell1 <- FALSE
  toPlot[ toPlot$nFeature_LiverCellAtlasSpatial < 1000, 'badCell1' ] <- TRUE
  toPlot$badCell2 <- FALSE
  toPlot[toPlot$nCount_LiverCellAtlasSpatial < 5000, 'badCell2'] <- TRUE
  toPlot[toPlot$nCount_LiverCellAtlasSpatial> 60000, 'badCell2'] <- TRUE
  toPlot$badCell3 <- FALSE
  toPlot[toPlot$percent.mt > 15, 'badCell3'] <- TRUE
  
  # Here I use functions from Liver Cell Atlas GitHub
  
  p1 <- drawVlnPlot(toPlot,colsToColor = c('badCell1','badCell2','badCell3'))
  print(p1)
  ggsave(p1,file='SpotsOverThresholdSCT.png')
  

  # Do filtering
  spatialObjectMergedFiltered <-subset(spatialObjectMerged, 
                         nFeature_LiverCellAtlasSpatial > 1000 & 
                           nCount_LiverCellAtlasSpatial > 5000 & 
                           nCount_LiverCellAtlasSpatial < 60000 & 
                           percent.mt < 15)
  
  dim(spatialObjectMergedFiltered)
  # From   31053  7947 to 31053  7829
  
  ########################################
  # Normalization with SCT 
  ########################################

  ## the reason for the normalization pre filtering is that if we normalize
  ## post merging we normalize all sections together 
  
  spatialObjectMergedFiltered <- SCTransform(spatialObjectMergedFiltered , 
                                             verbose = TRUE, 
                                             assay = "LiverCellAtlasSpatial")
  gc()
  remove(spatialObjectMerged)
  
  # this would be the alternative: to filter and transform each sample individually
  # I will compare the UMAPs of both methods to check if there is an effect of 
  # normalizing the data when merged or not-merged
  
  # do filtering and SCT normalization with separate Seurat objects in case
  # we need them 
  
  # i = 1
  # while(i<=4){
  #   spatialSamples[[i]][["percent.mt"]] <- 
  #     PercentageFeatureSet(spatialSamples[[i]], pattern = "^mt-")
  #   spatialSamples[[i]] <- subset(spatialSamples[[i]], 
  #                                 nFeature_LiverCellAtlasSpatial > 1000 & 
  #                                   nCount_LiverCellAtlasSpatial > 5000 & 
  #                                   nCount_LiverCellAtlasSpatial < 60000 & 
  #                                   percent.mt < 15)
  #   spatialSamples[[i]] <- SCTransform(spatialSamples[[i]], 
  #                                      verbose = TRUE, 
  #                                      assay = "LiverCellAtlasSpatial")
  #   
  #   i <- i + 1
  #   
  # }

  
  ###########################
  # Add annotation info to meta data (from GEO) of merged and filtered object
  ###########################
  
  # the annotation file from GEO includes a zonation parameter, a zonation group 
  # and cluster identity for each cell 
  # here I noticed that there are less spots annotated than there are spots in 
  # spots in the data set, so probably they got removed
  
  # reformat meta data
  spatialObjectMergedFiltered@meta.data$barcode <- 
    sub( ".*_", "", rownames(spatialObjectMergedFiltered@meta.data))
  
  # barcodes are repeated between samples! i.e. not unique, so need sample
  # names added to barcode to differentiate between spots
  
  spatialObjectMergedFiltered@meta.data$sample_barcode <- 
    rownames(spatialObjectMergedFiltered@meta.data)
  
  #do this for each samples
  length(which(annotationNAFLD$sample == 'JBO6'))
  # [1] 1121
  length(which(spatialObjectMergedFiltered$orig.ident == 'StSpatialSeuratM1S1'))
  # [1] 1244
  
  ## there more spots in Seurat object than in the annotation file
  
  annotationNAFLD$sample_barcode <- NA
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO6'),'sample_barcode'] <-
    paste0('StSpatialSeuratM1S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO6'),'spot'])
  
  length(which(annotationNAFLD$sample == 'JBO9'))
  # 2002
  length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM1S1'))
  # 2340
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO9'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM1S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO9'),'spot'])
  
  length(which(annotationNAFLD$sample == 'JBO10'))
  # 1780
  length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM2S1'))
  # 2021
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO10'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM2S1_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO10'),'spot'])
  
  length(which(annotationNAFLD$sample == 'JBO12'))
  # 1768
  length(which(spatialObjectMergedFiltered$orig.ident == 'WDSpatialSeuratM2S2'))
  # 2224
  
  annotationNAFLD[which(annotationNAFLD$sample == 'JBO12'),'sample_barcode'] <-
    paste0('WDSpatialSeuratM2S2_', annotationNAFLD[
      which(annotationNAFLD$sample == 'JBO12'),'spot'])
  
  annotationNAFLD$sample_barcode <- gsub(pattern = '*_1|_2|_3|_4', 
                                         x = annotationNAFLD$sample_barcode, 
                                         replacement = '')
  
  # add annotation file data to Seurat object meta data
  # have to map spot barcodes from one file to the other 
  
  temp <- annotationNAFLD
  rownames(temp) <- temp$sample_barcode
  # because they dont have the same size
  temp <- temp[intersect(annotationNAFLD$sample_barcode, 
                         spatialObjectMergedFiltered@meta.data$sample_barcode),]
  
  spatialObjectMergedFiltered@meta.data$BarcodeFromAnnotation <- NA
  spatialObjectMergedFiltered@meta.data$BarcodeFromAnnotation[match(
    temp$sample_barcode, rownames(spatialObjectMergedFiltered@meta.data))] <- 
    temp$sample_barcode
  
  # spatialObjectMergedFiltered@meta.data$annot_cluster[match(
  #   temp$sample_barcode, rownames(spatialObjectMergedFiltered@meta.data))] <- 
  #   annotationNAFLD$cluster
  
  spatialObjectMergedFiltered@meta.data$deleteNA <- 
    rep("not na", length(spatialObjectMergedFiltered@meta.data$BarcodeFromAnnotation))
  spatialObjectMergedFiltered@meta.data$deleteNA[is.na(
    spatialObjectMergedFiltered@meta.data$BarcodeFromAnnotation)] <-
    "na"
  
  spatialObjectMergedFiltered <-subset(spatialObjectMergedFiltered, 
                                        deleteNA == 'not na')
  
  #  sum(!(spatialObjectMergedFiltered$sample_barcode %in% annotationNAFLD$sample_barcode))
  # says all barcodes from seurat object in annotation file
  
  # match indices to add info to meta data
  spatialObjectMergedFiltered@meta.data$zonationGroup <- NA
  spatialObjectMergedFiltered@meta.data$zonationGroup[match(
    temp$sample_barcode, rownames(spatialObjectMergedFiltered@meta.data))] <- 
    temp$zonationGroup
  # check for sanity manually: all good
  
  # match indices to add info to meta data
  spatialObjectMergedFiltered@meta.data$zonation <- NA
  spatialObjectMergedFiltered@meta.data$zonation[match(
    temp$sample_barcode, rownames(spatialObjectMergedFiltered@meta.data))] <- 
    temp$zonation
  
  # check for sanity manually: all good
  
  spatialObjectMergedFiltered@meta.data$annot_cluster <- NA
  spatialObjectMergedFiltered@meta.data$annot_cluster[match(
    temp$sample_barcode, rownames(spatialObjectMergedFiltered@meta.data))] <- 
    temp$cluster
  
  saveRDS(spatialObjectMergedFiltered,
          file = 'hepatocyte-damage-score/Data/Output/LiverCellAtlasMouseSpatialMergedFilteredSCT.rds')
}

{ 
  # When replotting, we can start from here: 
  # merged, filtered, SCT transformed samples in one Seurat object
  # remember: AUCell does not used transformed counts 
  # start here:
  spatialObjectMergedFiltered <- readRDS(
    file = 'hepatocyte-damage-score/Data/Output/LiverCellAtlasMouseSpatialMergedFilteredSCT.rds')
  
  

  ###########################
  # UMAP
  ###########################
  # 
  # # UMAP of merged and then transformed "spatialObjectMergedFiltered"
  # 
  # # alternative 1: merged and SCT transformed
  # 
  # # this step is done in the github
  # ### !!! not sure about this step 
  # VariableFeatures(spatialObjectMergedFiltered ) <- c(VariableFeatures(spatialSamples[[1]]), 
  #                                       VariableFeatures(spatialSamples[[2]]),
  #                                       VariableFeatures(spatialSamples[[3]]),
  #                                       VariableFeatures(spatialSamples[[4]]))
  # 
  # spatialObjectMergedFiltered  <- RunPCA(spatialObjectMergedFiltered, 
  #                                        verbose = FALSE)
  # 
  # spatialObjectMergedFiltered  <- FindNeighbors(spatialObjectMergedFiltered , 
  #                                               dims = 1:30)
  # 
  # spatialObjectMergedFiltered  <- FindClusters(spatialObjectMergedFiltered , 
  #                                              resolution = 0.6)
  # 
  # spatialObjectMergedFiltered  <- RunUMAP(spatialObjectMergedFiltered , dims = 1:30)
  # 
  # 
  # # calculates percetage based on total mt in all samples merged (right?)
  # spatialObjectMergedFiltered [["percent.mt_merged"]] <- 
  #   PercentageFeatureSet(spatialObjectMergedFiltered, pattern = "^mt-")
  # 
  # ## RBC stands for 'red blood cells', hemoglobin (different ones) expressing cells
  # # I would have to look at this in paper 
  # spatialObjectMergedFiltered [["percent.rbc_merged"]] <- 
  #   PercentageFeatureSet(spatialObjectMergedFiltered, pattern = "Hb[ab]-")
  # 
  # plotMito <- FeaturePlot(spatialObjectMergedFiltered, c("percent.mt"))
  # plotMitoMerged <- FeaturePlot(spatialObjectMergedFiltered, c("percent.mt_merged"))
  # # plots look similar 
  # 
  # plotRBC <- FeaturePlot(spatialObjectMergedFiltered, c("percent.rbc_merged"))
  # 
  ########################################
  # Plot on UMAP (all of this just like in github from Liver Cell Atlas)
  ########################################
  
  # ### Plot clusters
  # umapPlot <- DimPlot(spatialObjectMergedFiltered, 
  #                     reduction = "umap",
  #                     order = c(4,6),
  #                     label = T)
  # ggsave(filename = 'UMAPFilteredMergedSCT.png',
  #       path = plotsPath)
  # 
  # umapPlotSplit <- DimPlot(spatialObjectMergedFiltered, 
  #                          reduction = "umap",
  #                          order = c(4,6),
  #                          group.by = 'orig.ident')
  # 
  # pcaPlot <- DimPlot(spatialObjectMergedFiltered, 
  #         reduction = "pca",
  #         pt.size = 0.1, 
  #         group.by = 'orig.ident',
  #         shuffle = TRUE) + scale_color_colorblind() +
  #   theme(legend.justification=c(-0.1,1), legend.position=c(0,0.95)) + 
  #   labs(title = NULL)
  #   
  # umapPlot
  # ggsave(filename = 'PCAFilteredMergedSCT.png',
  #        path = plotsPath)
  # 
  # 
  # # plotting clusters predominant in WD samples last to see if they are also 
  # # present in left cluster -> it seems like there are not clsuters 4 or 5 in left clsuter
  # 
  # DimPlot(spatialObjectMergedFiltered, 
  #         reduction = "umap",
  #         order = c(4,6),
  #         pt.size = 0.7
  #         )
  # 
  # PCAPlot(spatialObjectMergedFiltered,dims = c(3,4))
  # 
  # grid.arrange(umapPlot, 
  #              umapPlotSplit, 
  #              nrow=2)
  # 
  # ggsave(filename = 'UMAPandSampleIdentitiesFilteredMergedSCT.png',
  #        path = plotsPath)
  # 
  # ### Plot mitochondrial genes percentage on UMAP
  # # when explaining plots, dont forget to mention cut-offs
  # plotMito <- FeaturePlot(spatialObjectMergedFiltered,'percent.mt_merged',
  #                       min.cutoff = 'q2',  
  #                       max.cutoff = 'q98')
  # 
  # grid.arrange(umapPlot, umapPlotSplit, plotMito, ncol=3)
  # 
  # ### Plot RBC
  # # when explaining plots, dont forget to mention cut-offs
  # plotRBC <- FeaturePlot(spatialObjectMergedFiltered, 
  #                        'percent.rbc_merged', 
  #                        min.cutoff = 'q2', 
  #                        max.cutoff = 'q98')
  # 
  # grid.arrange(umapPlot, umapPlotSplit, plotMito, plotRBC, ncol=4)
  # 
  # ggsave(grid.arrange(umapPlot, umapPlotSplit, plotMito, 
  #                     plotRBC, ncol=2, nrow=2),
  #        file=paste0(dataPath, '/spatialMergedFilteredUMAPall.png'), 
  #        width = 24, height = 6)
  # 
  # 
  # ### Plot some genes
  # FeaturePlot(spatialObjectMergedFiltered, c("Glul", "Cyp2f2"), 
  #             min.cutoff = 'q2', max.cutoff = 'q98')
  # 
  # 
  # ### Plot clusters
  # SpatialDimPlot(spatialObjectMergedFiltered,
  #                idents)
  # 
  # 
  # #grid.arrange(SpatialDimPlot(spatialObjectMergedFiltered), ncol = 2,
  # #             nrow = 2)
  # SpatialDimPlot(spatialObjectMergedFiltered,label = TRUE, label.size = 3)
  # 
  # 
  # #Plot certain clusters (this are numbers on git hub, probably other numbers in here)
  # # clustersOI<-c(6,9,5,8,3,11)
  # # 
  # # SpatialDimPlot(spatialObjectMergedFiltered, 
  # #                cells.highlight = CellsByIdentities(
  # #                  object = spatialObjectMergedFiltered, 
  # #                  idents = clustersOI), 
  # #                facet.highlight = TRUE, 
  # #                images = 'StMouse001Sample001', 
  # #                ncol=6)
  # # 
  # #  
  # # SpatialDimPlot(spatialObjectMergedFiltered, 
  # #                cells.highlight = CellsByIdentities(
  # #                  object = spatialObjectMergedFiltered, 
  # #                  idents = clustersOI), 
  # #                facet.highlight = TRUE, 
  # #                images = 'WDMouse001Sample001',
  # #                ncol=6)
  # # 
  # # SpatialDimPlot(spatialObjectMergedFiltered, 
  # #                cells.highlight = CellsByIdentities(
  # #                  object = spatialObjectMergedFiltered, 
  # #                  idents = clustersOI), 
  # #                facet.highlight = TRUE, 
  # #                images = 'WDMouse002Sample001', 
  # #                ncol=6)
  # # 
  # # ### Plot mito
  # # SpatialFeaturePlot(spatialObjectMergedFiltered, 
  # #                    features = "percent.mt") 
  # 
  # 
  # ### Plot genes (why theses genes?)
  # SpatialFeaturePlot(spatialObjectMergedFiltered, 
  #                    features = c("Glul", "Cyp2f2"))
  # 
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
  identical(row.names(spatialObjectMergedFiltered@meta.data), 
            HDSmergedfiltered$cellID)
  # TRUE
  
  spatialObjectMergedFiltered@meta.data$HDS <- HDSmergedfiltered$HDS
  
  ########################################
  # Plot HDS
  ########################################

  # format and center data for box-plots

  HDSmergedfiltered$orig.ident <- 
    spatialObjectMergedFiltered@meta.data$orig.ident
  
  HDSmergedfiltered$centeredHDS <- scale(HDSmergedfiltered$HDS, 
                                            scale = FALSE)
  HDSmergedfiltered$zonationGroup <- 
    spatialObjectMergedFiltered@meta.data$zonationGroup
  
  HDSmergedfiltered$zonation <- 
    spatialObjectMergedFiltered@meta.data$zonation
  
  HDSmergedfiltered$diet <- NA
  HDSmergedfiltered[HDSmergedfiltered$orig.ident == 'StSpatialSeuratM1S1', 'diet'] <- 
    'Standard Diet' 
  HDSmergedfiltered[HDSmergedfiltered$orig.ident != 'StSpatialSeuratM1S1', 'diet'] <- 
    'Western Diet'
  HDSmergedfiltered$diet <- factor(HDSmergedfiltered$diet,
                                      ordered = TRUE, 
                                      levels = c("Standard Diet", "Western Diet"))

  
  # Not centered spot HDS distribution per sample 
  
  comparisons <- c('StSpatialSeuratM1S1')
  HDSViolinPlot <- ggplot(data = HDSmergedfiltered,
                          aes(x = orig.ident, y = HDS,
                              color = diet)) +
    geom_violin(trim = TRUE) + xlab("samples") + ylab("HDS") +
    scale_x_discrete(labels = c("StD \n M1S1", "WD \n M1S1",
                                "WD \n M2S1", "WD \n M2S2")) +
    theme(axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18),
          axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          legend.position = 'none')  + 
    scale_color_colorblind() +
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") +
    stat_compare_means(method = 'wilcox.test',
                       ref.group = 'StSpatialSeuratM1S1',
                       paired = FALSE,
                       label =  "p.signif",
                       label.y = -0.04,
                       size = 10,
                       symnum.args = list(
                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                         symbols = c("****", "***", "**", "*", "ns"))) +
    ylim(c(-0.45, 0))
  
    HDSViolinPlot
    
  ggsave(HDSViolinPlot, 
         path = plotsPath,
         file= 'SpatialLiverCellAtlasHDSViolinPlotWilcoxonTestStRefGroup.png')
 
 # DONT THINKS I'LL DO THIS: Plot of spatial distribution of HDS centered around zero

  # add diet column to meta deta
  
  if(identical(HDSmergedfiltered$cellID, 
               rownames(spatialObjectMergedFiltered@meta.data))){
    print('yes')
    spatialObjectMergedFiltered@meta.data$diet <- HDSmergedfiltered$diet
  }
  

  # Plot HDS on spots 

  # plots for thesis:
  
  maxHDS <- max(spatialObjectMergedFiltered@meta.data$HDS)
  minHDS <- -0.3
  spatialObjectMergedFiltered@meta.data$HDSforPlots <- spatialObjectMergedFiltered@meta.data$HDS
  spatialObjectMergedFiltered@meta.data$HDSforPlots[spatialObjectMergedFiltered@meta.data$HDSforPlots <= -0.3] <- -0.3
  # for plotting save all values 
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots", 
                           images = 'StMouse001Sample001', 
                           pt.size.factor = 2) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3,maxHDS)) +
    # scale_fill_gradient(low = 'white', high = 'black',
    #                     aesthetics = "fill",
    #                     limits = c(minHDS,maxHDS),
    #                     name = "HDS") +
    # scale_fill_gradientn( colours =  viridis(n = 100,
    #                                          direction = -1),
    #                       aesthetics = "fill",
    #                       limits = c(minHDS,maxHDS),
    #                       name = "HDS") +
    theme(legend.position = 'none',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
    ) +
    ggtitle("SD M1S1")
  
  
  p2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots",
                           images = 'WDMouse001Sample001',
                           pt.size.factor = 1.8) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3, maxHDS)) +
    theme(legend.position = 'none',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
          
    ) +
    ggtitle("WD M1S1")
  
  p3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots",
                           images = 'WDMouse002Sample001',
                           pt.size.factor = 2) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3, maxHDS)) +
    theme(legend.position = 'none',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
    ) +
    ggtitle("WD M2S1")
  
  p4 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots",
                           images = 'WDMouse002Sample002',
                           pt.size.factor = 2) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3,maxHDS)) +
    theme(legend.position = 'right',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
    ) +
    ggtitle("WD M2S2")
  
  plotJustImages <- wrap_plots(SpatialPlot(spatialObjectMergedFiltered, 
                                           alpha = 0)  & 
                                 theme(legend.position = "none") & ggtitle(''))
  
  
  # combine HDS on tissue, tissue without spots and box-plot of HDS distribution
  
  layer1plot <- wrap_plots((p1 | p2 | p3 | p4))
  panelplot <-  layer1plot / plotJustImages
  panelplot

  ggsave('SpatialHDSLiverCellAtlasYlOrBr-03threshold.png',
         path = plotsPath)
  
  # test colorblind friendliness 
  cvdPlot(layer1plot & theme(plot.title = element_text(size = 5),
                             legend.title = element_text(size = 5),
                             legend.key.size = unit(0.5, 'cm')))
  
  ggsave('SpatialHDSLiverCellAtlasYlOrBr-03thresholdColorBlindTest.png',
         path = plotsPath)
  
  # all plots without tissue
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots", 
                           images = 'StMouse001Sample001', 
                           pt.size.factor = 1,
                           crop = FALSE,
                           image.alpha = 0) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3,maxHDS),
                         name = 'HDS') +
    theme(legend.position = 'none',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
    ) +
    ggtitle("SD M1S1") 
    
  p1Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                          images = 'StMouse001Sample001', 
                          alpha = 0, crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
  
  p2 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDSforPlots", 
                           images = 'WDMouse001Sample001', 
                           pt.size.factor = 1,
                           crop = FALSE,
                           image.alpha = 0) + 
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1,
                         limits = c(-0.3,maxHDS),
                         name = 'HDS') +
    theme(legend.position = 'none',
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 18)
    ) + 
    ggtitle("WD M1S1") 
 
   p2Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                           images = 'WDMouse001Sample001',
                           alpha = 0, 
                           crop = FALSE)  & 
    theme(legend.position = "none") & ggtitle('')
   
   p3 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                            features = "HDSforPlots",
                            images = 'WDMouse002Sample001',
                            pt.size.factor = 1,
                            crop = FALSE,
                            image.alpha = 0) + 
     scale_fill_distiller(palette = "YlOrBr",
                          direction = 1,
                          limits = c(-0.3, maxHDS)) +
     theme(legend.position = 'none',
           legend.key.size = unit(1, 'cm'),
           legend.title = element_text(size = 18),
           plot.title = element_text(size = 18)
     ) +
     ggtitle("WD M2S1")
   
   p3Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                           images = 'WDMouse002Sample001',
                           alpha = 0, 
                           crop = FALSE)  & 
     theme(legend.position = "none") & ggtitle('')
   
   p4 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                            features = "HDSforPlots",
                            images = 'WDMouse002Sample002',
                            pt.size.factor = 1,
                            crop = FALSE,
                            image.alpha = 0) + 
     scale_fill_distiller(palette = "YlOrBr",
                          direction = 1,
                          limits = c(-0.3,maxHDS)) +
     theme(legend.position = 'right',
           legend.key.size = unit(1, 'cm'),
           legend.title = element_text(size = 18),
           plot.title = element_text(size = 18)
     ) +
     ggtitle("WD M2S2")
   
   p4Tissue <- SpatialPlot(spatialObjectMergedFiltered, 
                           images = 'WDMouse002Sample002',
                           alpha = 0, 
                           crop = FALSE)  & 
     theme(legend.position = "none") & ggtitle('')
   
   layer1plot <- wrap_plots((p1 | p2 | p3 | p4))
   plotJustImages <- wrap_plots((p1Tissue | p2Tissue | p3Tissue | p4Tissue ))
   panelplot <-  layer1plot / plotJustImages
   panelplot
   
   ggsave('SpatialHDSLiverCellAtlasYlOrBr-03thresholdNoCropping.svg',
          device = 'svg',
          dpi = 'print',
          path = plotsPath)
  
  
  # plot zonation values from the annotation file (from the authors of liver cell atlas)
  
  zP1 <- SpatialFeaturePlot(spatialObjectMergedFiltered,
                            features = 'zonation',
                            images = 'StMouse001Sample001',
                            pt.size.factor = 2)
  zP2 <- SpatialFeaturePlot(spatialObjectMergedFiltered,
                            features = 'zonation',
                            images = 'WDMouse001Sample001',
                            pt.size.factor = 2)
  zP3 <- SpatialFeaturePlot(spatialObjectMergedFiltered,
                            features = 'zonation',
                            images = 'WDMouse002Sample001',
                            pt.size.factor = 2)
  zP4 <- SpatialFeaturePlot(spatialObjectMergedFiltered,
                            features = 'zonation',
                            images = 'WDMouse002Sample002',
                            pt.size.factor = 2)
  
  zonationPlot <- zP1 | zP2 | zP3 | zP4
  zonationPlot
  ggsave( 'LCAspatialPlotZonationAnnotation.png', 
          path = plotsPath)
  
  # only Std. Diet and WD M1 S1
  
  p1 <- SpatialFeaturePlot(spatialObjectMergedFiltered, 
                           features = "HDS", 
                           images = 'StMouse001Sample001', 
                           pt.size.factor = 2) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(minHDS,maxHDS),
                          name = "HDS") +
    theme(legend.position = 'none',
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          plot.title = element_text(size = 20)) +
    ggtitle("Standard Diet\nMouse 1 Sample 1")
  
  p2temp <-  SpatialFeaturePlot(spatialObjectMergedFiltered, 
                                      features = "HDS",
                                      images = 'WDMouse001Sample001',
                                      pt.size.factor = 1.8) + 
    scale_fill_gradientn( colours =  turbo(n = 100),
                          aesthetics = "fill", 
                          limits = c(minHDS,maxHDS),
                          name = "HDS") +
    ggtitle("Western Diet\nMouse 1 Sample 1") +  
    
    theme(legend.position = 'right',
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          plot.title = element_text(size = 20)
    ) 
  
  
  tissuePlot1 <- SpatialPlot(spatialObjectMergedFiltered, 
                             images = 'StMouse001Sample001',
                             alpha = 0) + 
    theme(legend.position = "none") + 
    ggtitle('')
  
  tissuePlot2 <- SpatialPlot(spatialObjectMergedFiltered, 
                                    images = 'WDMouse001Sample001',
                                    alpha = 0) +
    theme(legend.position = "none") + 
    ggtitle('')
  
  panel2Sections <- wrap_plots((p1 | p2temp)) / wrap_plots((tissuePlot1 | tissuePlot2))
  
  panel2Sections
  ggsave('HDSTissueFilteredMergedLiverCellAtlasOnlyTwoSections.png',
         path = plotsPath)
  
  
  ## plot annotation clusters 
  
  pp1ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'StMouse001Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  pp2ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse001Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom',  legend.title = element_text('cluster ID'))
  pp3ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse002Sample001',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  pp4ClusterIdents <- SpatialDimPlot(spatialObjectMergedFiltered, 
                                     images = 'WDMouse002Sample002',
                                     pt.size.factor = 1.8) + 
    theme(legend.position = 'bottom', legend.title = element_text('cluster ID'))
  
  
  ClusterIdsTissue <- wrap_plots((pp1ClusterIdents | pp2ClusterIdents | pp3ClusterIdents | pp4ClusterIdents))
  
  ggsave('ClusterIdsPlusMitoPercNnTissue.png', path = plotsPath)

  
  
  # plotClusters <- SpatialPlot(spatialObjectMergedFiltered,
  #                            features = 'seurat_clusters',
  #                               alpha )  & 
  #   theme(legend.position = "none") & ggtitle('') 
  # 
  # plot
  # 
  # combine HDS on tissue, tissue without spots and box-plot of HDS distribution
  # 
  # layer1plot <- wrap_plots((p1 | p2 | p3 | p4))
  # panelplot <- plotJustImages / layer1plot / plotMito 
  # panelplot + plot_annotation(
  #     title = 'HDS Spatial Distribution',
  #     subtitle = 'Centered score was calculated after merging, filtering and transforming (SCT)\nLiver Cell Atlas spatial transcriptomic data'
  #   )  
  # ggsave('CenteredHDSonTissuePlusTissueSCTFilteredMerged.png',
  #        path = plotsPath)
  # 
  # # Plot for thesis: 
  # 
  # layer1plot <- wrap_plots(((p1 | p2) / (p3 | p4)))
  
      
  # idea for interpretation, only clusters 0, 3, 4, 7, and 11 are found in 
  # healthy tissue. we can say we associate the other clusters with disease 
  # and look at how the score is distributed in those clusters not present in 
  # healthy tissue 
  

  # how can I plot the zones on the tissue ??
  #!!!!!!
  # violin plots of HDS per zone group!! per slice and also grouped 
  
  # # First Western Diet Samples All Together
  # ZonationHDSViolinPlotWesternDiet <- 
  #   ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident != "StSpatialSeuratM1S1",], 
  #          aes(x = zonationGroup, y = HDS)) +
  #   geom_violin() + xlab("Zones") + ylab("centered HDS") +
  #   theme(axis.text.x = element_text(size=12),
  #         axis.title = element_text(size = 12),
  #         title = element_text(size = 15)
  #         )  + 
  #   geom_boxplot(width=0.1, outlier.size = 0) +
  #   ggtitle('Distribution of Centered HDS per Liver Zone in Western Diet Samples')
  # 
  # 
  ZonationHDSViolinPerSample <- 
    ggplot(data = HDSmergedfilteredSCT, 
           aes(x = zonationGroup, y = HDS,
               color = diet)) +
    geom_violin(trim = FALSE) + xlab("Hepatocyte zonation") + ylab("HDS") +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          title = element_text(size = 15))  + 
    scale_color_colorblind() +
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange") + 
    facet_wrap(~orig.ident)
  
  
  ZonationHDSViolinPerSample
  ggsave(ZonationHDSViolinPerSample, 
         path = plotsPath,
         file= 'ZonationHDSViolinPerSample.png')
  
  # # HDS Violin Plot St Mouse 1 Sample 1
  # ZonationHDSViolinPlotM1S1StandardDiet <- 
  #   ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "StSpatialSeuratM1S1",], 
  #          aes(x = zonationGroup, y = HDS)) + 
  #   geom_violin() + xlab("Zones") + ylab("centered HDS") +
  #   theme(axis.text.x = element_text(size=12),
  #         axis.title = element_text(size = 12),
  #         title = element_text(size = 12)
  #   )  + 
  #   geom_boxplot(width=0.1, outlier.size = 0) +
  #   ggtitle('Distribution of HDS per Liver Zone\nin Standard Diet Mouse 1 Sample 1')
  # 
  # ZonationHDSViolinPlotM1S1StandardDiet
  # 
  # 
  # # HDS Violin Plot WD Mouse 1 Sample 1
  # ZonationHDSViolinPlotM1S1WesternDiet <- 
  #   ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM1S1",], 
  #          aes(x = zonationGroup, y = centeredHDS)) +
  #   geom_violin() + xlab("Zones") + ylab("centered HDS") +
  #   theme(axis.text.x = element_text(size=12),
  #         axis.title = element_text(size = 12),
  #         title = element_text(size = 12)
  #   )  + 
  #   geom_boxplot(width=0.1, outlier.size = 0) +
  #   ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 1 Sample 1')
  # 
  # # HDS Violin Plot WD Mouse 2 Sample 1
  # ZonationHDSViolinPlotM2S1WesternDiet <- 
  #   ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM2S1",],
  #          aes(x = zonationGroup, y = centeredHDS)) +
  #   geom_violin() + xlab("Zones") + ylab("centered HDS") +
  #   theme(axis.text.x = element_text(size=12),
  #         axis.title = element_text(size = 12),
  #         title = element_text(size = 12)
  #   )  + 
  #   geom_boxplot(width=0.1, outlier.size = 0) +
  #   ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 2 Sample 1')
  # 
  # # HDS Violin Plot WD Mouse 2 Sample 2
  # ZonationHDSViolinPlotM2S2WesternDiet <- 
  #   ggplot(data = HDSmergedfilteredSCT[HDSmergedfilteredSCT$orig.ident == "WDSpatialSeuratM2S2",], 
  #          aes(x = zonationGroup, y = centeredHDS)) +
  #   geom_violin() + xlab("Zones") + ylab("centered HDS") +
  #   theme(axis.text.x = element_text(size=12),
  #         axis.title = element_text(size = 12),
  #         title = element_text(size = 12)
  #   )  + 
  #   geom_boxplot(width=0.1, outlier.size = 0) +
  #   ggtitle('Distribution of Centered HDS per Liver Zone\nin Western Diet Mouse 2 Sample 2')
  # 
  # HDSZonesAllViolins <- (ZonationHDSViolinPlotM1S1StandardDiet|ZonationHDSViolinPlotM1S1WesternDiet)/
  #   (ZonationHDSViolinPlotM2S1WesternDiet|ZonationHDSViolinPlotM2S2WesternDiet)
  # 
  # ggsave(HDSZonesAllViolins, 
  #        path = plotsPath,
  #        file= 'HDSZonesAllViolinPlots.png')
  # 

  
}
