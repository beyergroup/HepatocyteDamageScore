#### ==================#####
#### Randomization Test #### 
#### ==================#####

#### Load packages, data and functions ####
{
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
  library(tidyr)
  library(reshape2)
  
  
  # load functions written by us
  source("SharedFunctions.R")
  
  # paths to DEG tables of bulk studies
  
  pathToData <- 'hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/'
  path <- 'hepatocyte-damage-score/Data/Input/DEGTables/'
  studies <- c('GSE99010_CCl4_deg.csv',
            'GSE99010_WesternDiet_deg.csv',
            'GSE97234_AH_deg.csv',
            'GSE97234_ASH_deg.csv',
            'GSE83240_DEN_deg.csv',
            'GSE153580_STZ_deg.csv',
            'GSE148849_fastfooddiet_deg.csv',
            'GSE138419_AMLN_deg.csv',
            'GSE137449_CDAHFD_diet_deg.csv',
            'GSE135050_hfcfdiet_deg.csv',
            'GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv',
            'GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
            'GSE119953_HFD_deg.csv',
            'GSE119953_DDC_deg.csv',
            'GSE119441_HFD_deg.csv',
            'GSE119441_PFOA_deg.csv',
            'GSE114261_STAM_deg.csv',
            'GSE111828_Acetaminophen_deg.csv')

  
  data <- readRDS(paste0(pathToData,
                         'ExprMatrixmouseMergedHepatocytesnucSeqQCfitleredSCT.csv'))
  
  annotation <- 
    read.csv(paste0(pathToData,
                    'ExprMatrixAnnotationMouseMergedHepatocytesnucSeqQCfitleredSCT.csv'), 
             row.names = 1) 
  
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
  

  # loads list of all DEG tables
  DElist <- lapply( seq(studies), function(ii){
    read.csv(paste0(path,studies[[ii]]))}
    )
  names(DElist) <- basename(studies)
  
  #### Load list of Hepatocyte Damage Associated Genes
  HDAG <- 
    read.csv(file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')
  
  #### Load list of Genes Expressed in Hepatocytes 
  geneFilter <- 
    read.csv(file ='hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')
  
  geneFilter <- geneFilter$Genes
  
  # Load Liver Cell Atlas Hepatocytes expression matrix for testing #  
 
  # sampling 2000 random nuclei per condition 
  # set.seed(2)
  # use only cells 24 weeks NAFLD cohort

    #random sampling
    # c(sample(annotation$cell_id[annotation$condition == 'StSt'],
    #          size = 1000, replace = FALSE),
    #   sample(annotation$cell_id[annotation$condition == 'Nafld'], 
    #          size = 1000, replace = FALSE)
    # )



#### Randomization Step ####   
# Randomize a percentage of p-values and LFCs in DEGs and select genes 
# for HUDS, then calculate HDS and compare conditions.
perturbationTestFunction <- function(pert.percent = seq(10, 100, 10),
                                     ntop = 42, exprDat, DEdat, k = 50, 
                                     hepGenes){
  
   
  perturb_test <-lapply( seq(pert.percent) , function(ii) {

    
    print(paste("percent perturbed", pert.percent[ii] ))
    print(paste("these many studies:", length(DEdat)*( pert.percent[ii]/100)))
          # subset data 50 times to avoid data set bias (k times)

          HDS_pert.lvl <- Reduce( cbind.data.frame , lapply( 1:k, function(kk)
          {
            
            print(paste("randomisation round", kk ))
            
            datt <- DEdat
            # choose dataset to randomize
            whichR <- sample( 1:length(datt), length(datt)*( pert.percent[ii]/100) )
            
            # randomize p-value and lfc of chosen datasets
            print(whichR)
            datt[whichR]  <- lapply( datt[whichR]  , function(XX){
              XX$log2FoldChange <- sample( XX$log2FoldChange)
              XX$pvalue <- sample( XX$pvalue)
              return(XX) })
            
            HDS <- damage_signature.func( datt, 
                                          geneFilter = hepGenes)
            
            # calculate disease score for one dataset, sc Liver Cell Atlas, 2000 cells per condition
            XXX <- DS_calc.func( exprMatrices = exprDat, 
                                 DSignature = HDS, ntop = ntop)
            
            
            
            return(XXX)
          } )) 
          
          # prepare for plotting with ggplot
          # XX <- reshape::melt.list( data=HDS_pert.lvl )
          # XX$random.percent <- as.factor(pert.percent[ii])
          # return(XX)
          HDS_pert.lvl$cell <- rownames(HDS_pert.lvl)
          XX <- reshape::melt( data=HDS_pert.lvl , id.vars = "cell")
          XX$random.percent <- as.factor(pert.percent[ii])
          return(XX)
  })
    
  
    # returns a list of dataframes  
   return(perturb_test) 
    
}

testResult <- perturbationTestFunction(pert.percent = seq(0,100, 10),
                         ntop = 42, exprDat = exprMatrix, 
                         DEdat = DElist, k = 50, hepGenes = geneFilter)


saveRDS(testResult, 
        file = 'hepatocyte-damage-score/Data/Output/Results/RandomizationTest/RandomizationWeightedHDSLCAmergedNuq24WeeksNAFLDCohort.rds')


# Format Results for better plotting
toPlot <- Reduce( rbind, lapply( testResult , function(X) {
  X$value <- scale(X$value, scale = FALSE)
  return(X)
}))

# toPlot$condition <- gsub(pattern = "_.*", replacement = "", toPlot$cell)
# toPlot$cell_id <- gsub(pattern = ".*_", replacement = "", toPlot$cell)
 toPlot$HDS <- toPlot$value
 toPlot <- toPlot[,c(-2, -3)]

#idea 
toPlot$condition <- as.factor(annotation$diet[match(toPlot$cell,
                        annotation$cell_id)])

toPlot$condition <- factor(toPlot$condition,
                           ordered = TRUE,
                           levels = c("Standard Diet", "Western Diet"))

# plotting 
ggp <- ggplot( toPlot , aes( x=random.percent , y=HDS , color=condition)) +
   geom_boxplot(outlier.size = 0.1) +
  theme_bw()  +
  scale_color_colorblind() + 
  theme( text = element_text(size=18)) +
  labs(x = "percent of studies randomized [%]" , 
       y = "centered HDS"
       )

ggdens <- ggplot( toPlot , aes( x=random.percent , y=HDS , color=condition)) +
  geom_violin(trim = FALSE, size = 0.5) + 
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1),
               geom = "pointrange",
               position = position_dodge(0.95),
               size = 0.2
              ) + 
  scale_color_colorblind() + 
  theme( text = element_text(size=18)) +
  labs(x = "percent of studies randomized [%]" ,
       y = "centered HDS")


ggp 
ggsave(filename = 'RandomizationBoxplotsWeightedHDS24WeeksNAFLDCohort.png', 
       units = 'cm',
       path = 'hepatocyte-damage-score/Data/Output/Results/RandomizationTest/',
       width = 20, height = 25)


ggdens 
ggsave(filename = 'RandomizationViolinPlotWeightedHDS24WeeksNAFLDCohort.png', 
       units = 'cm',
       path = 'hepatocyte-damage-score/Data/Output/Results/RandomizationTest/',
       width = 20, height = 25)

}

# alternative code randomization
{
  perturbationTestFunction <- function(pert.percent = seq(10, 100, 10), 
                                       ntop = 42, exprDat, DEdat, k = 50, 
                                       hepGenes) {
    # Loop over percentages
    perturb_test <- lapply(pert.percent, function(percent) {
      message("Percent perturbed: ", percent)
      num_datasets <- length(DEdat) * (percent / 100)
      message("Number of datasets to randomize: ", num_datasets)
      
      # Perform k randomization steps
      HDS_pert_lvl <- lapply(1:k, function(step) {
        message("Randomization round: ", step)
        datt <- DEdat
        whichR <- sample(1:length(datt), num_datasets)
        datt[whichR] <- lapply(datt[whichR], function(dataset) {
          dataset$log2FoldChange <- sample(dataset$log2FoldChange)
          dataset$pvalue <- sample(dataset$pvalue)
          return(dataset)
        })
        
        # Calculate damage signature and disease score
        HDS <- damage_signature.func(datt, geneFilter = hepGenes)
        DS_calc.func(exprMatrices = exprDat, DSignature = HDS, ntop = ntop)
      })
      
      # Combine the randomization results into one data frame
      HDS_combined <- do.call(cbind, HDS_pert_lvl)
      HDS_combined$cell <- rownames(HDS_combined)
      print(head(HDS_combined))
      reshaped <- reshape2::melt(HDS_combined, 
                                 id.vars = "cell", 
                                 variable.name = "randomization", 
                                 value.name = "score")
      reshaped$random_percent <- factor(percent)
      
      print(head(reshaped))
      
      # Calculate the mean score for each cell
      averaged_scores <- reshaped %>%
        dplyr::group_by(cell, random_percent) %>%
        dplyr::summarise(mean_score = mean(score, 
                                           na.rm = TRUE))
      
      print(head(averaged_scores))
      
      return(averaged_scores)
    })
    
    # Combine results into a single dataframe
    combined_results <- do.call(rbind, perturb_test)
    return(combined_results)
  }
  
  perturbationTestFunction <- function(pert.percent = seq(10, 100, 10), 
                                       ntop = 42, exprDat, DEdat, k = 50, 
                                       hepGenes) {
    library(tidyr)
    library(dplyr)
    
    # Loop over percentages
    perturb_test <- lapply(pert.percent, function(percent) {
      message("Percent perturbed: ", percent)
      num_datasets <- length(DEdat) * (percent / 100)
      message("Number of datasets to randomize: ", num_datasets)
      
      # Perform k randomization steps
      HDS_pert_lvl <- lapply(1:k, function(step) {
        message("Randomization round: ", step)
        datt <- DEdat
        whichR <- sample(1:length(datt), num_datasets)
        datt[whichR] <- lapply(datt[whichR], function(dataset) {
          dataset$log2FoldChange <- sample(dataset$log2FoldChange)
          dataset$pvalue <- sample(dataset$pvalue)
          return(dataset)
        })
        
        # Calculate damage signature and disease score
        HDS <- damage_signature.func(datt, geneFilter = hepGenes)
        DS_calc.func(exprMatrices = exprDat, DSignature = HDS, ntop = ntop)
      })
      
      # Combine the randomization results into one data frame
      HDS_combined <- do.call(cbind, HDS_pert_lvl)
      
      # Add cell IDs as a column
      HDS_combined <- data.frame(HDS_combined, cell = rownames(HDS_combined))
      print(head(HDS_combined))
      
      # Reshape data using tidyr::pivot_longer
      reshaped <- HDS_combined %>%
        pivot_longer(
          cols = -cell,  # All columns except 'cell'
          names_to = "randomization", 
          values_to = "score"
        ) %>%
        mutate(random_percent = factor(percent))
      print(head(reshaped))
      

            
      # Calculate the mean score for each cell
      averaged_scores <- reshaped %>%
        group_by(cell, random_percent) 
    #    summarise(mean_score = mean(score, na.rm = TRUE))
      
      return(averaged_scores)
    })
    print(averaged_scores)
    
    # Combine results into a single dataframe
    combined_results <- bind_rows(perturb_test)
    return(combined_results)
  }
  
  
  
  testResultModified <- perturbationTestFunction(pert.percent = seq(0,100, 10),
                                         ntop = 42, exprDat = exprMatrix, 
                                         DEdat = DElist, k = 2, 
                                         hepGenes = geneFilter)
  
}