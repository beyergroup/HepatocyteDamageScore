# Load libraries and functions 

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
        y = -0.25,  #may need to modify this depending on your data
        label = paste('n =', length(y), '\n',
                      'M =', round(median(y), digits = 2), '\n')
      )
    )
  }
  
}



# FIgure 2 HDS

# Figure 2 D: Hepatocyte zones and HDS
{
  outputPath <- "hepatocyte-damage-score/Data/Output/"
  dfHDS <- readRDS(file = paste0(outputPath, 
                                 'Xiao2023HDSsnRNAseqAllHep.rds'))
  
  dfHDS$sample <- ordered(dfHDS$sample, 
                          levels = c("NC 3 months", 
                                     "NASH 3 months", 
                                     "NC 9 months", 
                                     "NASH 9 months"))
  
  my_comparisons <- list( c("NC 3 months", "NASH 3 months"), 
                          c("NC 9 months", "NASH 9 months"), 
                          c("NC 3 months", "NC 9 months"),
                          c("NASH 3 months","NASH 9 months" ))
  
  
  HDSbyHepType <- 
    ggplot(dfHDS[dfHDS$HepatocyteType != 'mNASH-Hep1' & dfHDS$HepatocyteType != 'mNASH-Hep2', ],
           aes( x = sample,
                 y = HDS,
                 color = condition
            )) + 
    geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    ylab("HDS") +
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
    theme_bw() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, 
                                      hjust = 1),
           legend.title=element_blank()) + 
    stat_compare_means(comparisons = my_comparisons) +
    facet_wrap(~HepatocyteType) +
    stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    position = position_dodge(1),
    size = 3
  )
  
  pdf( height = 8, width = 12, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/HDSXiao2023HDSbyHepZone.pdf") )
  print(HDSbyHepType)
  dev.off()
}


# Figure 2 F: Pseudotime vs. HDS
{
  # load TSCAN pseudotime inference results 
  # to plot against results shown in Figure 2 D
  
  
  dataPT <- readRDS(file = paste0(outputPath,
                                  'Results/Xiao2023_PseudoTimeTSCAN.rds'))
  
  # Plot pseudotime results on UMAP of hepatocytes of each sample
  listPTumaps <- sapply(names(dataPT[[1]]),
                        function(x){
                          plotUMAP(dataPT[['scExpObjects']][[x]],
                                            colour_by=I(dataPT[['PTvalues']][[x]]),
                                            point_size = 1,
                                            point_alpha = 1) +
                            scale_color_viridis(option = "B", 
                                                direction = -1) +
                            geom_line(data = dataPT[['PTordering']][[x]]$connected$UMAP, 
                                      mapping = aes(x = umap_1, y = umap_2, group = edge)) +
                            theme(legend.position = 'top',
                                  axis.ticks = element_blank(),
                                  axis.text = element_blank(),
                                  axis.title = element_text(size = 20),
                                  legend.key.size = unit(1.5, 'cm'),
                                  legend.title = element_text(size = 20),
                                  legend.text = element_text(size = 20),) +
                            labs(color = "PT") 
                        }, 
                        simplify = FALSE)
  names(listPTumaps) <- names(dataPT[[1]])
  
  
  ## Plot HDS: 
  # color palette will use min and max values HDS from each sample as 
  # reference points --> each UMAP needs it's own HDS legend. 
  # this is necessary cause the color resolution is needed to compare to pseudo
  # time trajectory 
  
  # Find global color boarders: 1%-5% top and bottom 

  
  limitBottom<- DataFrameHDS[order(DataFrameHDS$HDS),'HDS'][round((dim(DataFrameHDS)[1])*0.025)]
  limitTop <- DataFrameHDS[order(DataFrameHDS$HDS, decreasing = TRUE),'HDS'][round((dim(DataFrameHDS)[1])*0.025)]
  
  listHDSumaps <- sapply(names(dataPT[[1]]),
                         function(x){
                           # limitTop <- max(dataPT[['scExpObjects']][[x]]$HDS) * 0.95
                           # limitBottom <- min(dataPT[['scExpObjects']][[x]]$HDS) *1.05
                           
                          plotUMAP(dataPT[['scExpObjects']][[x]],
                                   colour_by ='HDS',
                                   point_size = 1,
                                   point_alpha = 1) +
                            geom_line(data = dataPT[['PTordering']][[x]]$connected$UMAP, 
                                      mapping = aes(x = umap_1, 
                                                  y = umap_2, 
                                                  group = edge)) +
                            scale_color_distiller(palette = "YlOrBr",
                                                  direction = 1,
                                                  limits = c(min(dataPT[['scExpObjects']][[x]]$HDS),
                                                             max(dataPT[['scExpObjects']][[x]]$HDS))
                                                  # limits = c(limitBottom,
                                                  #            limitTop),
                                                  #oob = squish
                                                  ) +
                            theme(legend.position = 'top',
                                  axis.ticks = element_blank(),
                                  axis.text = element_blank(),
                                  axis.title = element_text(size = 20),
                                  legend.key.size = unit(1.5, 'cm'),
                                  legend.title = element_text(size = 20),
                                  legend.text = element_text(size = 20)) +
                            labs(color = "HDS") 
                        }, 
                        simplify = FALSE)
  names(listHDSumaps) <- names(dataPT[[1]])
  
  # Plot Hepatocyte Zones
  
  listZoneumaps <- sapply(names(dataPT[[1]]), 
                          function(x){
                            plotUMAP(dataPT[['scExpObjects']][[x]], 
                                     colour_by='CellType',
                                     point_size = 1,
                                     point_alpha = 1) +
                              geom_line(data = dataPT[['PTordering']][[x]]$connected$UMAP, 
                                        mapping = aes(x = umap_1,
                                                      y = umap_2,
                                                      group = edge)) +
                              scale_color_colorblind() +
                              theme(legend.position = 'top',
                                    axis.ticks = element_blank(),
                                    axis.text = element_blank(),
                                    axis.title = element_text(size = 20),
                                    legend.key.size = unit(1.5, 'cm'),
                                    legend.title = element_blank(),
                                    legend.text = element_text(size = 8))
                            }, 
                          simplify = FALSE)
  

  umaps3moNASH <- (listHDSumaps[[3]] | listPTumaps[[3]] | listZoneumaps[[3]])
  
  pdf( height = 5, width = 15, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/HDSXiao2023_UMAPsPT_HDS_Zone_Mouse3moNASH.pdf") )
  print(umaps3moNASH)
  dev.off()

  
}

# Figure 2 G: Tranfer HDS to human data
{
  # read output of script 'HDSofHumanLiverscRNAseqForLowMemory.R"
  humanXiaoHDS <- readRDS(file = paste0(outputPath,
                                        'Results/Xiao2023_Human_HDS_Dataframe.rds'))
  
  violinXiaoHuman <- ggplot(data = humanXiaoHDS,
                            aes(x = factor(sample,
                                           level = c('Normal1', 'Normal2', 'Normal3',
                                  'NASH1', 'NASH2', 'NASH3')), 
                                y = HDS,
                                color = condition)) +
    geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    ylab("HDS") + 
    #ylim(c(-0.32, 0.15)) +
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
    theme_bw() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, 
                                      hjust = 1),
           legend.title=element_blank())  +
    stat_compare_means(label.x = 1.5)

  pdf( height = 6, width = 8, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/Xiao2023_HDS_human_snRNAseqHep.pdf") )
  print(violinXiaoHuman)
  dev.off()
  
}

# Supplementary Figure 2 G: the other pseudotimes
{
  umaps3moNC <- (listHDSumaps[[1]] | listPTumaps[[1]] | listZoneumaps[[1]])
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/HDSXiao2023_UMAPsPT_HDS_Zone_Mouse3moNC.pdf") )
  print(umaps3moNC)
  dev.off()
  
  umaps9moNC <- (listHDSumaps[[2]] | listPTumaps[[2]] | listZoneumaps[[2]])
  
  pdf( height = 5, width = 15, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/HDSXiao2023_UMAPsPT_HDS_Zone_Mouse9moNC.pdf") )
  print(umaps9moNC)
  dev.off()
  
  umaps9moNASH <- (listHDSumaps[[4]] | listPTumaps[[4]] | listZoneumaps[[4]])
  
  pdf( height = 5, width = 15, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/HDSXiao2023_UMAPsPT_HDS_Zone_Mouse9moNASH.pdf") )
  print(umaps9moNASH)
  dev.off()
  
}

# Supplementary Figure H
{
  
  violinXiaoHumanByHep <- ggplot(data = humanXiaoHDS,
                            aes(x = factor(sample,
                                           level = c('Normal1', 'Normal2', 'Normal3',
                                                     'NASH1', 'NASH2', 'NASH3')), 
                                y = HDS,
                                color = condition)) +
    geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    ylab("HDS") + 
    #ylim(c(-0.32, 0.15)) +
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
    theme_bw() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, 
                                      hjust = 1),
           legend.title=element_blank())  +
    stat_compare_means(label.x = 1.5) +
    facet_wrap(~HepatocyteType) 
  
  pdf( height = 10, width = 10, 
       file = paste0(outputPath,
                     "Results/FiguresManuscript/Xiao2023_HDS_human_snRNAseqHepByZONE.pdf") )
  print(violinXiaoHumanByHep)
  dev.off()
  
}
