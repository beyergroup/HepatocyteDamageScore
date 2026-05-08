# Supplementary Figure 1 

# Panel C and D: gene signature size test, Liver Cell Atlas - NAFLD cohort
{
  
  # ============================================#
  ############# Test Signature Sizes ############ 
  # ============================================#
  
  # Input: 
  # - load functions in SharedFunctions.R
  # - Liver Cell Atlas nucSeq filtered merged SCT transformed RNA Expression Matrix
  # - corresponding annotation
  # - List of Hepatocyte Damage Associated Genes 
  
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
  
  HDSdataframe <- mergedLCAnucSeq@meta.data 
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
  if(identical(rownames(mergedLCAnucSeq@meta.data), rownames(HDSdataframe))){
    mergedLCAnucSeq@meta.data$diet <- HDSdataframe$diet
    mergedLCAnucSeq@meta.data$age <- HDSdataframe$age
    mergedLCAnucSeq@meta.data$cohort <- HDSdataframe$cohort
    mergedLCAnucSeq@meta.data$HDS <- HDSdataframe$HDS
  }
  
  # Find Variable Features
  NAFLDcohort24weeks <- subset(mergedLCAnucSeq, subset = cohort == "NAFLD Cohort" & age == '24 weeks') 
  NAFLDcohort24weeks  <- ScaleData(NAFLDcohort24weeks, 
                                   features = rownames(NAFLDcohort24weeks ))
  gc()
  exprMatrix <- GetAssayData(NAFLDcohort24weeks,
                             layer = 'counts')

  # Test sizes for violin plots 
  
  size_vec <- c(5, 10, 15, 20, 30 , 42, 50, 100, 200)
  
  DSsize <- lapply( seq(size_vec) , 
                    
                    function(ii)
                    {
                      print(size_vec[ii])
                      ntop <- size_vec[ii]
                      DS_calc.func( exprMatrices = exprMatrix,
                                    DSignature = HDAG[1:ntop,]
                      )
                    } )
  
  DSsizeListDataFrames <- 
    lapply( seq(size_vec) , function(ii){
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
  
  # Make violin plots
  
  # y axis limits includes most extreme values 
  minMax <- c(min(DSsize[[8]]), max(DSsize[[9]]))
  
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
  
  print(DSsizeListViolin)
  

 
  ############
  # Panel D:
  ############
  # Test if significance in difference of means between conditions increases with
  # larger amount of top genes in gene set 
  
    size_vec <- c(2,12,22,32,42,52,62,72,82,92,102,152,202)
    
     # iterate over set sizes and calculate the disease score
      
      DSsize <- lapply( seq(size_vec) , 
                        
                        function(ii)
                        {
                          print(size_vec[ii])
                          ntop <- size_vec[ii]
                          DS_calc.func( exprMatrices = exprMatrix,
                                        DSignature = HDAG[1:ntop,]
                          )
                        } )
      
      
  
  DSsizeListDataFrames <- 
    lapply( seq(size_vec) , function(ii){
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
  
  
  all_hds_sizes <- cbind(DSsizeListDataFrames[[1]]$score, DSsizeListDataFrames[[2]]$score)
  
  for(i in 3:length(DSsizeListDataFrames)){
    all_hds_sizes <- cbind(all_hds_sizes,
                           DSsizeListDataFrames[[i]]$score)
  }
  
  all_hds_sizes <- cbind(as.data.frame(all_hds_sizes), "sample" = NAFLDcohort24weeks@meta.data$sample[match(rownames(DSsizeListDataFrames[[1]]),
                                                                                                            rownames(NAFLDcohort24weeks@meta.data))])
  
  all_hds_sizes <- cbind(as.data.frame(all_hds_sizes), "diet" =  NAFLDcohort24weeks@meta.data$diet[match(rownames(DSsizeListDataFrames[[1]]),
                                                                                                         rownames(NAFLDcohort24weeks@meta.data))])
  
  
  
  t_test_per_size  <- c(t.test( V1 ~ diet, all_hds_sizes[,c("V1","diet")] )$statistic,
                        t.test( V2 ~ diet, all_hds_sizes[,c("V2","diet")] )$statistic,
                        t.test( V3 ~ diet, all_hds_sizes[,c("V3","diet")] )$statistic,
                        t.test( V4 ~ diet, all_hds_sizes[,c("V4","diet")] )$statistic,
                        t.test( V5 ~ diet, all_hds_sizes[,c("V5","diet")] )$statistic,
                        t.test( V6 ~ diet, all_hds_sizes[,c("V6","diet")] )$statistic,
                        t.test( V7 ~ diet, all_hds_sizes[,c("V7","diet")] )$statistic,
                        t.test( V8 ~ diet, all_hds_sizes[,c("V8","diet")] )$statistic,
                        t.test( V9 ~ diet, all_hds_sizes[,c("V9","diet")] )$statistic,
                        t.test( V10 ~ diet, all_hds_sizes[,c("V10","diet")] )$statistic,
                        t.test( V11 ~ diet, all_hds_sizes[,c("V11","diet")] )$statistic,
                        t.test( V12 ~ diet, all_hds_sizes[,c("V12","diet")] )$statistic,
                        t.test( V13 ~ diet, all_hds_sizes[,c("V13","diet")] )$statistic)
  
  size_test_df <- DataFrame("signature_size" = size_vec,
                                "t_statistic" = t_test_per_size)
  
  t_statistic_per_size <- ggplot(size_test_df, aes( x=signature_size, y= abs(t_statistic ), group=1)) +
    geom_line()+
    geom_point( alpha=0.6) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18),
          axis.text = element_text(size = 18), 
          text = element_text(size = 18)) +
    scale_x_continuous(
      breaks = size_test_df$signature_size,
      labels = size_test_df$signature_size
    ) + ylab("T-statistic") + xlab("Gene signature size")
  
  t_statistic_per_size

}

# Panel F: permutation test - Liver Cell Atlas
{
  # Need output from: Randomization_test.R
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
  
  print(permViolin)
  
}

# Panel H: cell type specificity Test - Xiao et al
# (data in workstation)
{
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

path_to_data <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/Xiao_snRNAseq_GSE189600/"
setwd(path_to_data)

# hepatocyte damage associated genes
HDAG <- read.csv(
  file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
  sep = ',')

# load functions to calculate HDS
source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')

FileNames <- list.files(path = path_to_data)

metaDataAnnotations <- read.table("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/annotation_files/Xiao_snRNAseq_GSE189600_scitranslmed.adc9653_data_file_s1.csv", 
                                  sep = ';', dec = ',', header = TRUE)

metaDataAnnotations$cellBarcode <- 
  gsub('.*_', replacement = '', metaDataAnnotations$X)
metaDataAnnotations$condition <- 
  gsub('_.*', replacement = '', metaDataAnnotations$X)

options(future.globals.maxSize = 8000 * 1024^2)
# I. Read raw data & process & calculate all values: HDS, Pseudotime, etc
file_names <- list(FileNames[1:3], FileNames[7:9],FileNames[4:6], FileNames[10:12])
sample_i <- 1
seurat_list <- list()
HDS_list <- list()
conditions_metadata <- unique(metaDataAnnotations$condition )
for(sample_i in 1:length(conditions_metadata)){
  # read raw counts 
  
  
  temp <- ReadMtx(mtx =  paste0(path_to_data, file_names[[sample_i]][3]), 
                  cells = paste0(path_to_data, file_names[[sample_i]][1]),
                  features = paste0(path_to_data, file_names[[sample_i]][2]),
                  feature.column = 2) 
  
  seurat_list[[sample_i]] <- CreateSeuratObject(counts = temp, 
                                                project = conditions_metadata[sample_i],
                                                min.cells = 3)
  
  class(seurat_list[[sample_i]] )
  
  remove(temp)
  
  seurat_list[[sample_i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[sample_i]] , 
                                                                  pattern = "^mt-")
  
  
  seurat_list[[sample_i]]  <- subset(seurat_list[[sample_i]]  , 
                                     subset = 
                                       nFeature_RNA > 500 &
                                       nFeature_RNA < 8000 &
                                       percent.mt < 2)
  
  temp_metadata <- metaDataAnnotations[metaDataAnnotations$condition == conditions_metadata[sample_i] ,]
  
  rownames(temp_metadata) <- temp_metadata$cellBarcode
  
  temp_metadata <- temp_metadata[intersect(temp_metadata$cellBarcode, rownames(seurat_list[[sample_i]]@meta.data)),]
  
  seurat_list[[sample_i]]@meta.data$cellBarcode <- NA
  seurat_list[[sample_i]]@meta.data$CellType<- NA
  
  seurat_list[[sample_i]]@meta.data$cellBarcode[match(temp_metadata$cellBarcode, rownames(seurat_list[[sample_i]]@meta.data))] <- 
    temp_metadata$cellBarcode
  
  seurat_list[[sample_i]]@meta.data$CellType[match(temp_metadata$cellBarcode, rownames(seurat_list[[sample_i]]@meta.data))] <- 
    temp_metadata$CellCluster
  
  seurat_list[[sample_i]]@meta.data$deleteNA <- 
    rep("not na", length(seurat_list[[sample_i]]@meta.data$CellType))
  seurat_list[[sample_i]]@meta.data$deleteNA[is.na(seurat_list[[sample_i]]@meta.data$CellType)] <-
    "na"
  
  seurat_list[[sample_i]] <- subset(seurat_list[[sample_i]],
                                    subset = deleteNA == 'not na')
  
  seurat_list[[sample_i]] <- subset(seurat_list[[sample_i]],
                                    subset = CellType != "Meso.Cell")
  
  
  seurat_list[[sample_i]] <- SCTransform(seurat_list[[sample_i]])
  
  seurat_list[[sample_i]]$sample <- conditions_metadata[sample_i] 
  
  ## calculate HDS
  HDS_list[[sample_i]] <- DS_calc.func(exprMatrices = 
                                         GetAssayData(seurat_list[[sample_i]],
                                                      assay = 'SCT',
                                                      layer = 'counts'),
                                       DSignature = HDAG)
  
  gc()
  
  if(identical(
    names(HDS_list[[sample_i]]), 
    rownames((seurat_list[[sample_i]]@meta.data))
  )){
    seurat_list[[sample_i]]@meta.data$HDS <- unname(HDS_list[[sample_i]])
  }
  
  
  
  
}

data_to_plot <- rbind(
  seurat_list[[1]]@meta.data,
  seurat_list[[2]]@meta.data,
  seurat_list[[3]]@meta.data,
  seurat_list[[4]]@meta.data)

# summarize hepatocyte zones into one category -> hepatocytes 
data_to_plot$celltype_sum <- data_to_plot$CellType
data_to_plot <- data_to_plot %>%
  mutate(., celltype_sum = case_when(
    (celltype_sum == "PP-Hep" | celltype_sum == "Int-Hep" | celltype_sum == "PC-Hep" | celltype_sum == "mNASH-Hep1" | celltype_sum == "mNASH-Hep2") ~ "Hepatocyte",
    (celltype_sum == "EC1" | celltype_sum == "EC2") ~ "Endothelial cell", 
    (celltype_sum == "T cell" |celltype_sum == "B cell") ~ "Lymphocyte",
    celltype_sum == "Cholangiocyte"  ~ "Cholangiocyte",
    (celltype_sum == "Mac1" | celltype_sum == "Mac2" | celltype_sum == "Mac3" ) ~ "Macrophage",
    (celltype_sum == "Stellate1" | celltype_sum == "Stellate2") ~ "Stellate cell",
    
    
    
  ) )

data_to_plot$sample <- factor(data_to_plot$sample, levels = c("NC3mo", "NC9mo", "NASH3mo", "NASH9mo"))

data_to_plot$celltype_sum <- factor(data_to_plot$celltype_sum, levels = unique(data_to_plot$celltype_sum))

mean_per_group <- data_to_plot %>%
  group_by(sample, celltype_sum) %>%
  summarize(mean = mean(HDS))

# horizontal and devided by cell type
ggplot(data_to_plot, aes(x = HDS, color = sample)) + geom_density(lwd = 1.5, trim = TRUE) + facet_grid(cols = vars(celltype_sum)) + 
  scale_color_colorblind() + theme_bw() + theme( text = element_text(size = 18), legend.position = "top") +
  geom_vline(data =  mean_per_group, aes(xintercept = mean, color = sample), lwd = 1.5, linetype = "dashed") 

}



