# Supplementary Figure 2 -- HDS 

# Load libraries, functions and path
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

outputPath <- 'hepatocyte-damage-score/Data/Output/Results/'

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

# Panel C : (part of this boxplot has been put in Figure 1 F)
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


# Panel C: second part 
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


# Panel D & E : acetaminophen and aging data 

{
  # plot acetaminophen on its own because it's a time series: GSE111828 
  # plot aging data also on its own to show increase with age
  
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
  
 # pdf( height = 6, width = 12, file = paste0(outputPath,
  #                                           "FiguresManuscript/HDS_cv.modelAcetaminophen.pdf") )
  print(cvAcetaminophen)
#  dev.off()
  
  table(AcetaminophenGroupHDS$model, AcetaminophenGroupHDS$condition)
  
  AgingGroupHDS <- 
    read.csv(file = paste0(outputPath,
                           'LeaveOneOutCrossValidationAging.csv'),
             sep = ',',
             dec = '.', 
             stringsAsFactors = TRUE)[,-1]
  
  cv_Aging <- ggplot(AgingGroupHDS, 
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
  
 # pdf( height = 6, width = 12, file = paste0(outputPath,
  #                                           "FiguresManuscript/HDS_cv.modelAgingSpearmanCorLM.pdf") )
  print(cv_Aging )
  #dev.off()
  

  
}


