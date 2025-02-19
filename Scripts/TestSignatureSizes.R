# ============================================#
############# Test Signature Sizes ############ 
# ============================================#

# Input: 
# - load functions in SharedFunctions.R
# - Liver Cell Atlas nucSeq filtered merged SCT transformed RNA Expression Matrix
# - corresponding annotation
# - List of Hepatocyte Damage Associated Genes 
# Output: 
# - RDS with test results
# - Plot

library(AUCell)
library(BiocGenerics)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(plyr)
library(biomaRt)
library(Seurat)
library(ggpubr)
library(cowplot)
library(magrittr)
library(tidyverse)
library(dplyr)

#### load functions written by us
# generate damage signature
{
  source("SharedFunctions.R")

}

#### Load Data
{
pathToData <- 'hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/'

data <- readRDS(paste0(pathToData,
                       'ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv'))

annotation <- 
  read.csv(paste0(pathToData,
                  'ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv'), 
           row.names = 1) 

#### Load list of Hepatocyte Damage Associated Genes
HDAG <- read.csv(file = 
                   'hepatocyte-damage-score/Data/Output/HDAG.csv')

#### Load list of Genes Expressed in Hepatocytes 
geneFilter <- 
  read.csv(file ='hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')

}

#### Subset data set for testing ##
{
  # keep only NAFLD Cohort samples 24 weeks old

  annotation$diet <- NA
  annotation[annotation$Sample == "CS197" | 
               annotation$Sample ==  "CS196" |
               annotation$Sample == "CS193" |
               annotation$Sample == "CS192" |
               annotation$Sample == "CS191" |
               annotation$Sample == "CS185"  |
               annotation$Sample == "CS183" , 'diet'] <- 'Western Diet'
  
  annotation$diet[is.na(annotation$diet)] <- 'Standard Diet'
  
  toKeep <- annotation$cell_id[annotation$Sample == "CS193" | 
                                 annotation$Sample ==  "CS192" |
                                 annotation$Sample == "CS191" |
                                 annotation$Sample == "CS185" |
                                 annotation$Sample == "CS183" |
                                 annotation$Sample == "CS190" |
                                 annotation$Sample == "CS189" |
                                 annotation$Sample == "CS188" |
                                 annotation$Sample == "CS184" |
                                 annotation$Sample == "CS182"] 
    
  
   exprMatrix <- data[, toKeep]
 
}

#### Test effect of a set size i.e. number of top genes to calculate HDS 
{
  size_vec <- c(5, 10, 15, 20, 30 , 42, 50, 100, 200, 300, 500, nrow(HDAG))
  
  {
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
    

  }
 
  
  # save results in RDS file containing list of names vectors 
  
   saveRDS(DSsize,
           'hepatocyte-damage-score/Data/Output/testSizesmouseMergedHepatocytesnucSeqQCfitleredSCT.rds')
  


  
####  Plotting  
  # about the plotting functions:
  # The function mean_sdl is used for adding mean and standard deviation. 
  # It computes the mean plus or minus a constant times the standard 
  # deviation. In the R code above, the constant is
  # specified using the argument mult (mult = 1). 
  # By default mult = 2. The mean +/- SD can be added as a crossbar or
  # a pointrange.
{
  
  # Prepare data for plotting
  
  DSsizeListDataFrames <- 
    lapply( seq(DSsize) , function(ii){
      datt <- cbind.data.frame( score = DSsize[[ii]] ,
                                setSize = as.factor(size_vec[ii]), 
                                condition = as.factor(
                                  annotation$diet[match(names(DSsize[[ii]]),
                                                        annotation$cell_id)])
      )
      datt$condition <- relevel( datt$condition , ref = "Standard Diet")
      
      return(datt)
    })
  
  # plot
  
  DSsizeListViolin <- lapply( 
    seq(DSsizeListDataFrames) , function(ii){
      ppg <- ggplot( DSsizeListDataFrames[[ii]] , 
                     aes( x=condition , y=score,
                          color=condition)) + 
        geom_violin(trim = FALSE) + #theme_bw() + 
        scale_color_colorblind() +
        ggtitle(paste(DSsizeListDataFrames[[ii]]$setSize[1], "Genes")) +
        theme( text = element_text(size=16),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.x = element_blank(),
               plot.title = element_text(size=16),
               legend.position = 'none') +
        guides(colour = guide_legend(override.aes = aes(label = "")))+
        labs(y = "HDS") + 
        stat_summary(fun.data = "mean_sdl",
                     fun.args = list(mult = 1),
                     geom = "pointrange",
                     width =0.2,
                     size = 0.1
                     )  +
        ylim(c(-0.22,0.06 ))
      return( ppg )
    }
  )
}

SizePannelFigure <- (DSsizeListViolin[[1]] | DSsizeListViolin[[2]] | DSsizeListViolin[[3]] | DSsizeListViolin[[4]]) / (DSsizeListViolin[[5]] | DSsizeListViolin[[6]] | DSsizeListViolin[[7]] | DSsizeListViolin[[8]]) / (DSsizeListViolin[[9]] | DSsizeListViolin[[10]] | DSsizeListViolin[[11]] | DSsizeListViolin[[12]] + theme(legend.position = 'bottom'))

SizePannelFigure
ggsave(filename = 'SizeTestLiverCellAtlasPanelFigureOnly24WeeksNAFLDCohortWeightedDS.png', 
       units = 'px',
       path = 'hepatocyte-damage-score/Data/Output/Results/',
       width = 2600, height = 2800)
}
