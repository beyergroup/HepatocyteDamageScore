###################################################
######## Analysis of Cell Fate Trajectories #######
###################################################

# I. Load packages, data and gene sets for analysis

{
  
  
  library(Seurat)
  library(rtracklayer)
  library(sctransform)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(ggthemes)
  library(reshape2)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(glmGamPoi)
  library(TSCAN)
  library(scater)
  library(scales)
  library(tidyr)
  library(RColorBrewer)
  library(VennDiagram)
  library(org.Mm.eg.db)
  library(Orthology.eg.db)

# hepatocyte damage associated genes
HDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')

genefilter <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv'
  )


senMayoGenes <- read.csv(
  file = 
    paste0(dataPath,
           'Input/GeneSetsCellFates/SenMayoGeneSetMouse.csv'),
  sep = ';')

geneNoTInHDAG <- setdiff(genefilter$Genes, HDAG$gene_symbol)

# load functions to calculate HDS
source('SharedFunctions.R')

pathData <- 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/'

outputPath <- 'hepatocyte-damage-score/Data/Output/'

FileNames <- list.files(path = pathData)

metaDataAnnotations <- read.table(
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/scitranslmed.adc9653_data_file_s1.csv', 
  sep = ';', dec = ',', header = TRUE)

metaDataAnnotations$cellBarcode <- 
  gsub('.*_', replacement = '', metaDataAnnotations$X)
metaDataAnnotations$condition <- 
  gsub('_.*', replacement = '', metaDataAnnotations$X)

# read preprocessed merged seurat object
# SCTransformed, HDS and AUCell SenMayo calculated 

mergedObj <- readRDS(
  file = paste0(outputPath,
                "MergedXiaoHDS_SenMayoHepatocytes.rds"))

}



# Intersection 
{
    venn.plot <- venn.diagram(
    x = list(Set1 = HDAG[1:42,"gene_symbol"], Set2 = senMayoGenes$Gene.murine.),
    category.names = c("SenMayo", "HDS 42 genes"),
    filename = NULL,  # Set to NULL to plot in R instead of saving to a file
    output = TRUE,
    main = "Intersection of gene sets"
  )
  dev.off()
  grid.draw(venn.plot)
}



# preprare data frame with all meta data 

{
  
  mergedDataFrame <- mergedObj@meta.data
  
  
  mergedDataFrame$CellType <- factor(mergedDataFrame$CellCluster,
                                     ordered = TRUE,
                                     levels = c("PP-Hep",
                                                "Int-Hep",
                                                "PC-Hep",
                                                "mNASH-Hep1",
                                                "mNASH-Hep2"))
  
  mergedDataFrame$condition <- factor(mergedDataFrame$condition,
                                      ordered = TRUE,
                                      levels = c("3m NC","9m NC", 
                                                 "3m NASH", "9m NASH" ))
  
  mergedDataFrame <- mergedDataFrame %>% 
    mutate(HDS_bin = cut(HDS, breaks=15))
  
  
  
  
  
}

# Plottin HDS vs. SenMayo 

{
  
  # 1. Bin cells by HDS
  

  # 2.a. Plot 'SenMayo' AUCellScore per bin
  
  HDSbinsSenMayo <- ggplot(data = mergedDataFrame, 
                  aes(x=HDS_bin, y=AUCell_SenMayo )) +
    geom_jitter(size = 0.08, 
                width = 0.3, 
                height = 0,
                alpha = 1) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    stat_summary(fun = 'mean',
                 geom = "crossbar",
                 width = .5, color = "red") +
    stat_summary(fun.data = "mean_se", 
                 geom = "errorbar", 
                 width = .1, color = 'red') +
    geom_hline(yintercept= median(mergedDataFrame$AUCell_SenMayo), 
               linetype="dashed", 
               color = "blue") + 
    ylab('SenMayo Activity Score')
  
  pdf(height = 6, 
      width = 8,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/",
               "Bins_HDS_SenMayo_Mean_SE.pdf"))
  print(HDSbinsSenMayo)
  dev.off()
  
  png(height = 6, 
      width = 8,
      units = 'in',
      res = 300,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/PNGs/",
               "Bins_HDS_SenMayo_Mean_SE.png"))
  print(HDSbinsSenMayo)
  dev.off()

  
  # 2.b. SenMayo AUCellScore x HDS per sample
  # with Spearmann correlation coefficient 
  
  SenMayo_HDS_bySample <-
    ggplot(data = mergedDataFrame,
           aes(x = HDS,
               y = AUCell_SenMayo)) + 
    geom_point(size = 0.6,
               alpha = 0.8)   +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    facet_wrap( ~ condition) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'rho') +
    ylab('SenMayo Activity Score')
  
  
  pdf(height = 6, 
      width = 6,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/",
               "Spearman_SenMayo_HDS_bySampleMergedXiao.pdf"))
  print(SenMayo_HDS_bySample)
  dev.off()
  
  png(height = 6, 
      width = 6,
      units = 'in',
      res = 300,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/PNGs/",
               "Spearman_SenMayo_HDS_bySampleMergedXiao.png"))
  print(SenMayo_HDS_bySample)
  dev.off()
  
  # 2.c. SenMayo AUCellScore x HDS per zone (all samples mixed)
  # zones as annotated by data set authors 
  # with Spearmann correlation coefficient 
  
  SenMayo_HDS_byZone <-
    ggplot(data = mergedDataFrame,
           aes(x = HDS,
               y = AUCell_SenMayo)) + 
    geom_point(size = 0.6,
               alpha = 0.8)   +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    facet_wrap( ~ CellType) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'rho') +
    ylab('SenMayo Activity Score')
  
  pdf(height = 6, 
      width = 6,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/",
               "Spearman_SenMayo_HDS_byZoneMergedXiao.pdf"))
  print(SenMayo_HDS_byZone)
  dev.off()
  
  png(height = 6, 
      width = 6,
      units = 'in',
      res = 300,
      file = 
        paste0(outputPath,
               "Results/CellFateAnalysis/PNGs/",
               "Spearman_SenMayo_HDS_byZoneMergedXiao.png"))
  print(SenMayo_HDS_byZone)
  dev.off()
  
  SenMayo_HDS_conditionsMixed <-  
    ggplot(data = mergedDataFrame,
           aes(x = HDS,
              y = AUCell_SenMayo)) + 
    geom_point(size = 0.8,
               alpha = 0.8)   +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    stat_cor(method = 'spearman',
             cor.coef.name = 'rho')

  
  # 3. Thesholds for classification defined by quartiles
  quantilesHDSall <- 
    quantile(mergedObj@meta.data$HDS, 
             probs = c(0.25, 0.5, 0.75))
  
  quantilesSenMayoall <- 
    quantile(mergedObj@meta.data$AUCell_SenMayo, 
             probs = c(0.25, 0.5, 0.75,0.90))
  
  # floor and ceiling values for applying limits to the color scale
  # Squish values outside the limit to the nearest boundary
  floor_valueHDS <- quantile(mergedObj@meta.data$HDS, 0.01)  
  ceiling_valueHDS <- quantile(mergedObj@meta.data$HDS, 0.99)
  
  floor_valueSenMayo <- quantile(mergedObj@meta.data$AUCell_SenMayo, 0.01)  
  ceiling_valueSenMayo <- quantile(mergedObj@meta.data$AUCell_SenMayo, 0.99)
  
  umapSamples<-
    DimPlot(mergedObj, 
            group.by = 'condition',
            pt.size = 1,
            cells = sample(rownames(mergedObj@meta.data),
                           replace = FALSE)) + 
    scale_color_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height=unit(5,"mm"),
          legend.key.width =unit(5,"mm"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          text = element_text(size=18))
  
  png( height = 4.5, width = 4,
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/UmapXiao_Conditions.png") )
  print(umapSamples)
  dev.off()
  
  pdf( height = 4.5, width = 4,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/UmapXiao_Conditions.pdf") )
  print(umapSamples)
  dev.off()
  
  umapHDSmerged <- FeaturePlot(mergedObj,
                               features = "HDS",
                               cells = sample(rownames(mergedObj@meta.data),
                                              replace = FALSE),
                               pt.size = 1) + 
    theme(legend.position = "top", 
          legend.key.size = unit(15, "mm"),
          legend.key.height = unit(5,"mm"),
          text = element_text(size=20),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) +
    scale_colour_distiller(palette = "YlOrBr", 
                           direction = 1,
                           limits = c(floor_valueHDS, 
                                      ceiling_valueHDS),  
                           oob = scales::squish) + ggtitle("")
  
  umapSenMayomerged <- FeaturePlot(mergedObj,
                                   features = "AUCell_SenMayo",
                                   cells = sample(rownames(mergedObj@meta.data),
                                                  replace = FALSE),
                                   pt.size = 1) + 
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           limits = c(floor_valueSenMayo, 
                                      ceiling_valueSenMayo),  
                           oob = scales::squish) +
    theme(legend.position = "top", 
          legend.key.height=unit(5,"mm"),
          legend.key.width = unit(15, "mm"),
          text = element_text(size=20),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) + ggtitle("")
  
  
  densityPlotsHDS <- ggplot(data = mergedObj@meta.data,
                            aes(x = HDS, color = condition))+
    geom_density(lwd = 1)+ theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("HDS")
  
  densityPlotsSenMayo <- ggplot(data = mergedObj@meta.data,
                                aes(x = AUCell_SenMayo, color = condition))+
    geom_density(lwd = 1) + theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("SenMayo activity score")
  
  ggGrid1 <- cowplot::plot_grid(plotlist = list(umapHDSmerged, 
                                                umapSenMayomerged,
                                                densityPlotsHDS,
                                                densityPlotsSenMayo), 
                                rel_heights = c(1,0.35), nrow = 2) 
  
  png( height = 7, width = 8, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "UmapXiao_HDS_SenMayo_DensityPlots.png") )
  print(ggGrid1)
  
  dev.off()
  
  
  pdf( height = 7, width = 8, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "UmapXiao_HDS_SenMayo_DensityPlots.pdf") )
  print(ggGrid1)
  dev.off()
  
  
}
  
# Classify hepatocytes into cell fates based on HDS and SenMayo thresholds
  {
  
  # define column 'group' in metadata with cell fate  according to 
  # senMayo AUCell score and HDS 
  mergedObj@meta.data$CellID <- rownames(mergedObj@meta.data)
  mergedObj@meta.data$group <- NA
  
  # tag cells to cell fate by HDS and AUCell_SenMayo thresholds:
  # possible cell fates: healthy, damaged + senescent, and damaged + not senescent
  
  # threshold healthy: HDS <= 1. quartile 
  tempIndex <- mergedObj@meta.data$HDS <= quantilesHDSall[1]
  mergedObj@meta.data$group[tempIndex] <- 'Healthy'
  remove(tempIndex)
  
  # threshold: damaged + senescent
  # HDS > 3. quartiles & SenMayo > 3. quartile
  tempIndex <- 
    mergedObj@meta.data$AUCell_SenMayo > quantilesSenMayoall[3] &
    mergedObj@meta.data$HDS > quantilesHDSall[3]
  
  mergedObj@meta.data$group[tempIndex] <- 'DamagedSenescent'
  remove(tempIndex)
  
  # threshold: damaged + not senescent
  # HDS > 3. quartiles & SenMayo < 2.quartile  
  
  tempIndex <- 
    mergedObj@meta.data$AUCell_SenMayo <= quantilesSenMayoall[2] &
    mergedObj@meta.data$HDS > quantilesHDSall[3]
  mergedObj@meta.data$group[tempIndex] <- 'DamagedNotSenescent'
  remove(tempIndex)
  
  # tag NA's (cells that do not fulfill any of the criteria) as transition cells
  mergedObj@meta.data$group <- mergedObj@meta.data$group %>% replace_na('Transition')
  
  mergedObj@meta.data$group <- factor(mergedObj@meta.data$group, 
                                      ordered = TRUE,
                                      levels = c('Healthy', 
                                                 'Transition',
                                                 'DamagedNotSenescent',
                                                 'DamagedSenescent'),
                                      labels = c('Undamaged',
                                                 'Transition',
                                                 'Not-senescent-like, damaged',
                                                 'Senescent-like, damaged'))
  
  
  
  umapFate<-
    DimPlot(mergedObj, 
            group.by = 'group',
            pt.size = 1,
            cells = sample(rownames(mergedObj@meta.data))) + 
    scale_color_manual(values = c('Undamaged' = '#FDE725FF',
                                  'Transition' = '#7AD151FF',
                                  'Not-senescent-like, damaged' = '#414487FF',
                                  'Senescent-like, damaged' = '#440154FF'))+
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          legend.key.height=unit(5,"mm"),
          legend.text = element_text(size = 6))
  
  
  
  png( height = 4.5, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/UmapXiao_CellFateClassification.png") )
  print(umapFate)
  dev.off()
  
  mergedObj@meta.data$condition <- factor(mergedObj@meta.data$condition,
                                          ordered = TRUE,
                                          levels = c("3m NC","9m NC", 
                                                     "3m NASH", "9m NASH" ))
  
  umapSamples<-
    DimPlot(mergedObj, 
            group.by = 'condition',
            pt.size = 1,
            cells = sample(rownames(mergedObj@meta.data))) + 
    scale_color_colorblind() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          legend.key.height=unit(5,"mm"),
          text = element_text(size=12))
  
  png( height = 4.5, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/UmapXiao_Conditions.png") )
  print(umapSamples)
  dev.off()
  
  
  # Summarize the data: count the number of cells for each fate in each sample
  summarizedCellFates <- mergedObj@meta.data %>%  
    mutate(condition = factor(condition, 
                              levels = c('3m NC',
                                         '9m NC',
                                         '3m NASH',
                                         '9m NASH')),
           group = factor(group,
                          levels = c('Undamaged',
                                     'Transition',
                                     'Not-senescent-like, damaged',
                                     'Senescent-like, damaged'))) %>%
    group_by(condition, group) %>%
    summarize(Cell_Count = n(), .groups = "drop") %>%
    group_by(condition) %>% 
    mutate(proportion = Cell_Count / sum(Cell_Count))
  
  write.csv(summarizedCellFates[,1:4], file = paste0(outputPath,
                                                     "Results/CellFateAnalysis/CellFateCountsXiao.csv"))
  
  # Prepare a table-like data frame for annotations
  
  label_data <- summarizedCellFates %>% 
    group_by(condition) %>%
    mutate(n = sum(Cell_Count)) %>%
    dplyr::select(condition, group, Cell_Count,n) %>%
    pivot_wider(names_from = group, values_from = Cell_Count, values_fill = 0)
  
  # Convert the data frame to a long format for plotting as a table
  label_data_long <- label_data[,-2] %>%
    pivot_longer(
      cols = -condition,
      names_to = "group",
      values_to = "Cell_Count"
    ) %>% 
    mutate(group = factor(group,
                          levels = c('Undamaged',
                                     'Transition',
                                     'Not-senescent-like, damaged',
                                     'Senescent-like, damaged')))
  
  # Create a stacked bar plot
  stackedBarPlotCellFates <- 
    ggplot(summarizedCellFates, aes(x = condition, 
                                    y = proportion, 
                                    fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Not-senescent-like, damaged' = '#414487FF',
                                 'Senescent-like, damaged' = '#440154FF')) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      y = "Proportion of Hepatocytes (%)",
      x = "",
      title = "Normalized number of hepatocytes \n from each cell fate per sample "
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text (hjust = 1,
                                size = 18),
      plot.margin = margin(5,5,40,5),
      legend.position = "right",
      legend.text = element_text(size = 18)
    ) + geom_text(
      data = label_data,
      aes(x = condition, y = 1.05, label = n),  # Place labels above the bar
      # fontface = "bold",
      size = 7,
      inherit.aes = FALSE  # Prevent unwanted aesthetics from being inherited
    )
  
  png( height = 6, width = 10, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/Xiao_StackedBarplotCellFates_Percentages.png") )
  
  print(stackedBarPlotCellFates)
  
  dev.off()
  
  pdf( height = 6, width = 10, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/Xiao_StackedBarplotCellFates_Percentages.pdf") )
  
  print(stackedBarPlotCellFates)
  
  dev.off()
  
  
  
} 
  
  
# Find genes differentially expressed in different cell fates  
{
  # Find markers for each cell fate group comparison,
  
  # Normalize the data
  testData <- mergedObj
  testData <- PrepSCTFindMarkers(testData)
  Idents(testData) <- "group"

  
  # from the vignette : 
  # However, the p-values obtained from this analysis should be interpreted
  # with caution, because these tests treat each cell as an independent 
  # replicate and ignore inherent correlations between cells originating
  # from the same sample. Such analyses have been shown to find a large 
  # number of false positive associations, as has been demonstrated by 
  # Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, 
  # and others. Below, we show how pseudobulking can be used to account 
  # for such within-sample correlation.
  
  #
  
  # to find genes differentially expressed in each group of cells compared to 
  # all other cells 
  
  # differentially expressed genes 1 group vs rest of cells 
  
  allMarkers <- FindAllMarkers(testData) 
  allMarkers <- filter(allMarkers,
                       p_val_adj < 0.01 &  abs(avg_log2FC) >= 1)
  
  # differentially expressed genes in pairwise comparisons of all groups

  SenescenceDamageMarkers <- FindMarkers(testData,
                                         ident.1 =  'Senescent-like, damaged',
                                         ident.2 = 'Not-senescent-like, damaged')
  View(SenescenceDamageMarkers)
  
  SenescenceDamageMarkers$cluster <- rep('Senescent_vs_Not_senscent', 
                                         length(SenescenceDamageMarkers$p_val))
  
  SenescenceDamageMarkers$gene <- rownames(SenescenceDamageMarkers)
  
  # Are there positive markers intersecting with HDS or SenMayo? 
  
  write.csv(filter(SenescenceDamageMarkers,
                   p_val_adj < 0.05 &
                     avg_log2FC > 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1))),
            file = paste0(outputPath,
                          'Upregulated_Genes_Senescent_vs_NotSenescenct_DamagedHeps.csv'))
  
  positiveMarkers <- filter(SenescenceDamageMarkers,
                            p_val_adj < 0.05 &
                              avg_log2FC > 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))$gene
  
  intersect(positiveMarkers, senMayoGenes$Gene.murine.)
  
  intersect(positiveMarkers, HDAG$gene_symbol[1:42])
  
  
  # Are there negative markers intersecting with HDS or SenMayo? 
  
  write.csv( filter(SenescenceDamageMarkers,
                    p_val_adj < 0.05 &
                      avg_log2FC < -0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1))),
            file = paste0(outputPath,
                          'DownRegulated_Genes_Senescent_vs_NotSenescenct_DamagedHeps.csv'))
  
  
  negativeMarkers <- filter(SenescenceDamageMarkers,
                            p_val_adj < 0.05 &
                              avg_log2FC < -0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))$gene
  
  intersect(negativeMarkers, senMayoGenes$Gene.murine.)
  
  intersect(negativeMarkers, HDAG$gene_symbol[1:42])
  
  
  
  # NotSenescen_Healthy <- FindMarkers(testData,
  #                                     ident.1 = 'Not-senescent-like, damaged',
  #                                     ident.2 = 'Undamaged')
  # 
  # 
  # # 
  # # 
  # # HealthyvsDS <- FindMarkers(testData,
  # #                            ident.1 = 'Healthy',
  # #                            ident.2 = 'Damaged and senescent') 
  # #   filter(p_val_adj < 0.01 & abs(logfc.threshold) >= 0.5 )
  # # 
  # # 
  # # HealthyvsDNS <- FindMarkers(testData,
  # #                             ident.1 = 'Healthy',
  # #                             ident.2 = 'Damaged, not-senescent') 
  # #   filter(p_val_adj < 0.01 & abs(logfc.threshold) >= 0.5 )
  # # 
  # 
  # 
  # # PSEUDOBULK: 
  # 
  # # pseudobulk the counts based on zone and group 
  # pseudoTestData <- Seurat::AggregateExpression(testData, assays = "RNA",
  #                                               return.seurat = TRUE,
  #                                               group.by = c('CellType', 'group'))
  # 
  # tail(Cells(pseudoTestData))
  # 
  # Idents(pseudoTestData) <- pseudoTestData$group
  # 
  # bulkSenescenceDamageMarkers <- FindMarkers(object = pseudoTestData,
  #                                            ident.1 = "Damaged and senescent",
  #                                            ident.2 = "Damaged, not-senescent",
  #                                            test.use = "DESeq2") %>% arrange(p_val) 
  # 
  # View(filter(bulkSenescenceDamageMarkers, p_val < 0.01 & abs(avg_log2FC) > 0.5))


}

### II. ORA 

# II.1. GO term ORA of all markers per group (TIM's pipeline)
{
  source("hepatocyte-damage-score/Scripts/func_analysis.R") 
  
  # list of vectors 
  allMarkers_RobertGO.50 <- lapply(unique(allMarkers$cluster),
                                function(nname){
                                  # convert to entrez
                                  universe <- unique( entr2gName$entrezgene_id)
                                  universe <- as.character( universe[!is.na(universe)] )
                      
                                  # select genes passing significance
                                  # and foldchange threshold that belong each
                                  # group (== nname)
                                  gset <- filter(allMarkers,
                                                 p_val_adj < 0.01 & 
                                                   abs(avg_log2FC) >= 1 &
                                                   cluster == nname )$gene
                            
                                  # prepare gene set, convert ensembleIDs to entrezIDs
                                  geneset <- unique(entr2gName$entrezgene_id[match( gset,  entr2gName$external_gene_name)])
                                  geneset <- geneset[!is.na(geneset)]
                                  geneset <- as.character(geneset)
                                  
                                  print(length(gset))
                                  print(length(intersect(geneset,colnames(gomatrix))))
                                  
                                  # apply Robert's function that given a sparse
                                  # matrix of GO terms (columns = genes, 
                                  # rows = GO terms)
                                  # a geneset of interest and a background 
                                  # set of genes (universe)
                                  # carry out clustering with members diverging 
                                  # by at most cut_max genes, 
                                  # and do enrichment testing.
                                  # Note, multiplicity adjustment 
                                  # is performed for the representative terms only.
                                  RobertGO.50 <- tryCatch( sf.clusterGoByGeneset( gomatrix, geneset, universe,
                                                                                  min.genes= 3 , cut_max = 50  ),
                                                           error = function(e) NA )
                                  print(head(RobertGO.50$results))
                                  # if( !is.na(RobertGO.50)) RobertGO.50$results$Term <- goterms[ match(
                                  #   RobertGO.50$results$GO.ID, names(goterms))]
                                  return(RobertGO.50)
                                      })
  
  
   names(allMarkers_RobertGO.50) <- unique(allMarkers$cluster)
   
   # table of average Foldchanges per GO term per group 
   
  
  ### top functions per model
  # get a union of top N func in each dataset, select from primary terms only
  
  # Fisher test significance threshold
  thrsh <- 0.05
  goRes <- allMarkers_RobertGO.50
  
  topN_GOrb  <- Reduce(union, lapply(seq(goRes), function(ii, N = 10)
  {
    print(ii)
    datt <- goRes[[ii]]
    # select from primary terms only
    datt <- subset(datt$results,  Primary.Fisher.adj < thrsh)
    datt <- datt[ order(datt$Fisher) ,]
    # return N top primar terms
    return(datt$GO.ID[ 1:min(N, nrow(datt)) ])
  }))
  
  # Bar plot enriched top 10 only primary GO terms with adj. p-value < 0.05 
  GOrobert_barplot_clustered(iid = 1:4,
                             GOtoPlot = topN_GOrb, datGO = goRes)

  


  
  LFC_goTerms <- sapply(unique(allMarkers$cluster), function(jj){
    
    tempMarkers <- filter(allMarkers,
                          p_val_adj < 0.01 & 
                            abs(avg_log2FC) >= 1 &
                            cluster == jj)
    
    xList <- sapply(topN_GOrb, 
                    function(ii){
                      
                      genes_inGO <- go2gName[[ii]]
                      temp <- na.omit(tempMarkers[genes_inGO,])
                      x <- mean(temp$avg_log2FC)
                      # names(x) <- ii
                      
                      return(x)
                    })
    return(xList)
  } )
  
  LFC_goTerms <- as.data.frame(LFC_goTerms)
  colnames(LFC_goTerms) <- unique(allMarkers$cluster)
  
  
  
  # same but with fold changes 
  GOrobert_barplot_clustered(iid = 1:4,
                             GOtoPlot = topN_GOrb, 
                             datGO = goRes,
                             datLFC = LFC_goTerms
                            )
  
  
  
  
cnetplot()
  
  
}

# II.2. GO term enrinchment of differentially expressed genes
# when comparing damaged and senescent cells against damaged not senescent cells

{
  
  # background genes convert to entrez
  universe <- unique(entr2gName$entrezgene_id)
  universe <- as.character(universe[!is.na(universe)])
  
  # select genes passing significance
  # and foldchange threshold that belong each
  gset <- rownames(filter(SenescenceDamageMarkers,
                 p_val_adj <= 0.01 & 
                   abs(avg_log2FC) >= 0.5) )
  
  ggvenn(list("damaged-senescent-like vs. not-senscent-like" = gset,
              "42 HDS markers" = HDAG$gene_symbol[c(1:42)]), 
    c("damaged-senescent-like vs. not-senscent-like", "42 HDS markers"),
    fill_color = c("#56B4E9", "#009E73"),
    stroke_size = 0.5,
    set_name_size = 4
  )
  
  # prepare gene set, convert ensembleIDs to entrezIDs
  geneset <- unique(entr2gName$entrezgene_id[match(gset, entr2gName$external_gene_name)])
  geneset <- geneset[!is.na(geneset)]
  geneset <- as.character(geneset)
  
  print(length(gset))
  print(length(intersect(geneset,colnames(gomatrix))))

  SenDama_RobertGO.50 <- tryCatch( sf.clusterGoByGeneset( gomatrix, geneset, universe,
                                                  min.genes= 3 , cut_max = 5  ),
                           error = function(e) NA )
  print(head(SenDama_RobertGO.50$results))

  ### top functions per model
  # get a union of top N func in each dataset, select from primary terms only
  
  testAppend <- c(allMarkers_RobertGO.50, 
                  list(SenDama_RobertGO.50 ))
  names(testAppend)[[5]] <- 'DS vs. D-Not-S'
  
  # Fisher test significance threshold
  thrsh <- 0.05
  
  topN_GOrbAppend  <- Reduce(union, lapply(seq(testAppend), function(ii, N = 15)
  {
    print(ii)
    datt <- testAppend[[ii]]
    # select from primary terms only
    datt <- subset(datt$results,  Primary.Fisher.adj < thrsh)
    datt <- datt[ order(datt$Fisher) ,]
    # return N top primar terms
    return(datt$GO.ID[ 1:min(N, nrow(datt)) ])
  }))


  # 
  # Bar plot enriched top 15 only primary GO terms with adj. p-value < 0.05 
  GOrobert_barplot_clustered(GOtoPlot = topN_GOrbAppend,
                             iid = 1:5,
                             datGO = testAppend)
  
  
  
  tail(rbind.data.frame(allMarkers, SenescenceDamageMarkers))
  clusterNames <-unique(rbind.data.frame(allMarkers, SenescenceDamageMarkers)$cluster)
  LFC_goTermsAppend <- sapply(clusterNames, function(jj){
    
    tempMarkers <- filter(rbind.data.frame(allMarkers, 
                                           SenescenceDamageMarkers),
                          p_val_adj < 0.01 & 
                            abs(avg_log2FC) >= 0.5 &
                            cluster == jj)
    print(head(tempMarkers))
    
    xList <- sapply(topN_GOrbAppend, 
                    function(ii){
                      
                      genes_inGO <- go2gName[[ii]]
                      temp <- na.omit(tempMarkers[genes_inGO,])
                      x <- mean(temp$avg_log2FC)
                      return(x)
                    })
    return(xList)
  } )
  
  LFC_goTermsAppend <- as.data.frame(LFC_goTermsAppend)
  colnames(LFC_goTermsAppend) <- unique(clusterNames)
  
  
  
  # same but with fold changes 
  
  
  GOrobert_barplot_clustered(iid = 5,
                             GOtoPlot = topN_GOrbAppend, 
                             datGO = testAppend,
                             datLFC = LFC_goTermsAppend)
  
  
  
  png( height = 15, width = 15, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DSvsDnS_GOrobert_barplot_clustered_N15_all.png") )
  
  GOrobert_barplot_clustered(iid = 5,
                             GOtoPlot = topN_GOrbAppend, 
                             datGO = testAppend,
                             datLFC = LFC_goTermsAppend)
  dev.off()
  
  
  pdf( height = 15, width = 15, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_DSvsDnS_GOrobert_barplot_clustered_N15_all.pdf") )
  GOrobert_barplot_clustered(iid = 5,
                             GOtoPlot = topN_GOrbAppend, 
                             datGO = testAppend,
                             datLFC = LFC_goTermsAppend)
  dev.off()
  
  
  
  
  
  
  
}


# II.3. GO ORA (Paula's Version)
{
  
  library(clusterProfiler)
  
  # a.1. ORA GO (ALL) markers senescent-like-damaged hepatocytes
  

  egoPosBP <- enrichGO(gene = positiveMarkers,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  
  #with universe 
 
  egoPos_simplified <- simplify(
    egoPosBP,
    cutoff = 0.5,             # Similarity threshold (0.7 is a common default)
    by = "p.adjust",          # Retain the most significant term
    select_fun = min,         # Select term with the smallest p-value within each cluster
    measure = "Wang"          # Semantic similarity measure ("Wang" is common for GO terms)
  )
  
  heatplot(enrichplot::pairwise_termsim(egoPos_simplified), 
           showCategory=30,
           label_format = 60)
  
  png( height = 8, width = 10, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DS_Markers_ORA_GO_ALL_top30with05similarity.png") )
  
  
  mutate(egoPosBP, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 30, label_format = 60, font.size = 15) +
    facet_grid(~ONTOLOGY)
  
  dev.off()
  
 
   png( height = 8, width = 10, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DS_Markers_ORA_GO_ALL_simplified_top20with05similarity.png") )
  
  mutate(egoPos_simplified, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 30, label_format = 60, font.size = 15) +
    facet_grid(~ONTOLOGY)
  
   dev.off()
  # 
  # ## just cellular compartment ## WITH 
  # 
  # egoPosCC <- enrichGO(gene = filter(SenescenceDamageMarkers,
  #                                    p_val_adj < 0.05 &
  #                                      avg_log2FC > 0.5 &
  #                                      pct.1 > 0.3)$gene,
  #                      OrgDb         = org.Mm.eg.db,
  #                      keyType       = 'SYMBOL',
  #                      ont           = "CC",
  #                      pAdjustMethod = "BH",
  #                      pvalueCutoff  = 0.01,
  #                       qvalueCutoff  = 0.05) %>% simplify(
  #                        .,
  #                        cutoff = 0.5,
  #                        by = "p.adjust",
  #                        select_fun = min,
  #                        measure = "Wang")
  # 
  # png( height = 6, width = 6, 
  #      units = 'in',
  #      res = 300,
  #      file = paste0(outputPath,
  #                    "Results/CellFateAnalysis/PNGs/",
  #                    "Xiao_DS_Markers_ORA_GO_CC_simplified_top20with06similarity.png") )
  # 
  # mutate(egoPosCC, qscore = -log(p.adjust, base=10)) %>% 
  #   barplot(x = "qscore", showCategory = 10, label_format = 60, font.size = 12) 
  # 
  # dev.off()
  # 
  # heatplot(enrichplot::pairwise_termsim(egoPosCC), showCategory=20)
  
  # a.2. KEGG ORA
  
  kkPos <- enrichKEGG(gene = bitr(positiveMarkers, 
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Mm.eg.db)$ENTREZID,
                      keyType = 'ncbi-geneid',
                      organism     = 'mmu',
                      pvalueCutoff = 0.05)
  
  kkPos@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                   replacement = "", kkPos@result$Description, 
                                   fixed = TRUE)
  
  kkPos@result <- filter(kkPos@result, p.adjust < 0.05)
  
  
  png( height = 4, width = 6, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DS_Markers_ORA_KEGG_.png") )
  
  
  mutate(kkPos, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x = "qscore", showCategory = 20, label_format = 60, font.size = 12) +
    ggtitle('Sig. enriched KEGG Pathways in \n Damaged-Senescent-Like Hepatocyte Markers')
  
  dev.off()
  
  png( height = 4, width = 6, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DS_Markers_ORA_KEGG_HEATMAP.png") )

  
  heatplot(enrichplot::pairwise_termsim(kkPos), showCategory=13) 

  dev.off()
  
  toPlotPositive <- filter(SenescenceDamageMarkers,
                          p_val_adj < 0.05 &
                            avg_log2FC > 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  toPlotPositive[order(toPlotPositive$avg_log2FC, decreasing = FALSE),]
  
  # Plot distribution of the four genes 
   plotDSMarkers <- VlnPlot(testData,
                 features = toPlotPositive$gene[1:20], 
                 split.by = "group", 
                 stack = TRUE,
                 layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Not-senescent-like, damaged' = '#414487FF',
                                 'Senescent-like, damaged' = '#440154FF')) +
     theme(axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none")
   
   
   png( height = 4, width = 14, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DS_Markers_ViolinExpressionLevel_Top15.png") )
  
   plotDSMarkers
   dev.off()
   
   
   png( height = 4, width = 8, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DS_All_Markers_ViolinExpressionLevel.png") )
   
   VlnPlot(testData,
           features = c(filter(SenescenceDamageMarkers,
                               p_val_adj < 0.05 &
                                 avg_log2FC > 0.5 &
                                 pct.1 > 0.3)$gene), 
           split.by = "group", 
           stack = TRUE,
           layer = "data") + 
     scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                  'Transition' = '#7AD151FF',
                                  'Not-senescent-like, damaged' = '#414487FF',
                                  'Senescent-like, damaged' = '#440154FF')) +
     theme(axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none")
   
   dev.off()
 
   # b.1. ORA GO (ALL) markers not-senescent-like-damaged hepatocytes 
  
  
  egoNegALL <- enrichGO(gene = negativeMarkers,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
  
  egoNegALL_simplified <- simplify(egoNegALL, cutoff = 0.7, by = "p.adjust", 
                                  select_fun = min, measure = "Wang")
  
  
  png( height = 8, width = 12, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DNotS_Markers_ORA_GO_ALL_simplified_04similarityTop20.png") )
  
   mutate(egoNegALL_simplified, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory = 20, label_format = 65, font.size = 15) + 
     facet_grid(~ONTOLOGY)
   
   dev.off()
 
   
   
   png( height = 6, width = 12, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DNotS_Markers_ORA_GO_ALL_top30.png") )
   
   mutate(egoNegALL, qscore = -log(p.adjust, base=10)) %>% 
     barplot(x="qscore", showCategory = 30, label_format = 65, font.size = 15) + 
     facet_grid(~ONTOLOGY)
   
   dev.off()
 
  
  
  # b.2. ORA KEGG 
   
   kkNeg <- enrichKEGG(gene = bitr(negativeMarkers, 
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = org.Mm.eg.db)$ENTREZID,
                       keyType = 'ncbi-geneid',
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
   
   kkNeg@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                    replacement = "", kkNeg@result$Description, 
                                    fixed = TRUE)
   
   kkNeg@result <- filter(kkNeg@result, p.adjust < 0.05)
   
   
   png( height = 5, width = 8, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DnotS_Markers_ORA_KEGG_.png") )
   
   
   mutate(kkNeg, qscore = -log(p.adjust, base=10)) %>% 
     barplot(x = "qscore", showCategory = 20, label_format = 60, font.size = 12) +
     ggtitle('Sig. enriched KEGG Pathways in \n Damaged-Not-Senescent-Like \n Hepatocyte Markers')
   
   dev.off()
   
   png( height = 5, width = 6, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DS_Markers_ORA_KEGG_HEATMAP.png") )
   
   
   heatplot(enrichplot::pairwise_termsim(kkNeg), showCategory=13) 
   
   dev.off()
   toplotmarkers <- filter(SenescenceDamageMarkers,
          p_val_adj < 0.05 &
            avg_log2FC < -0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  toplotmarkers <- toplotmarkers[order(toplotmarkers$p_val_adj), ]
   
   plotDNSMarkers <- VlnPlot(testData,
                            features = toplotmarkers$gene[1:20], 
                            split.by = "group", 
                            stack = TRUE,
                            layer = "data") + 
     scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                  'Transition' = '#7AD151FF',
                                  'Not-senescent-like, damaged' = '#414487FF',
                                  'Senescent-like, damaged' = '#440154FF')) +
     theme(axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none") + ggtitle('Top 20 marker genes (Log2FC > 0.5 ) in damaged, \n not-senescent-like hepatocytes')
   
   
   png( height = 4, width = 15, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DNS_Top20_Markers_ViolinExpressionLevel.png") )
   
   plotDNSMarkers
   dev.off()
   
   plotDNSMarkers2 <- VlnPlot(testData,
                             features = toplotmarkers$gene[21:40], 
                             split.by = "group", 
                             stack = TRUE,
                             layer = "data") + 
     scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                  'Transition' = '#7AD151FF',
                                  'Not-senescent-like, damaged' = '#414487FF',
                                  'Senescent-like, damaged' = '#440154FF')) +
     theme(axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none") + ggtitle('Top 21-40 marker genes (Log2FC > 0.5 ) in damaged, \n not-senescent-like hepatocytes')
   
   
   png( height = 4, width = 15, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_DNS_Top16_30_Markers_ViolinExpressionLevel.png") )
   
   plotDNSMarkers2
   dev.off()
  
  
}

## III. GSEA 

# Do GSEA with ordered differentially expressed genes 
{
  
  
  # postive markers # 
  d <- filter(SenescenceDamageMarkers,
         avg_log2FC > 0 ) %>% 
    arrange(.,by_group = desc(avg_log2FC )) %>% 
    select(gene, avg_log2FC)
  
 mappingTemp <- bitr(d$gene, 
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)
 d[mappingTemp$SYMBOL,'gene'] <- mappingTemp$ENTREZID
  
  geneList <- d[,2]
  
  ## feature 2: named vector
  names(geneList) <- as.character(d[,1])
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  ## feature 3: decreasing orde
  
  gsePos <- gseGO(geneList     = geneList,
                  OrgDb        = org.Mm.eg.db,
                  ont          = "ALL",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE,
                  by = "fgsea",
                  scoreType = "pos")
  
  
     barplot(gsePos,x = "qvlue", showCategory = 10, 
             label_format = 60, font.size = 12) +
     ggtitle('')
  
  
  
  # use `showCategory` to select 
  # the displayed terms. It can be a number of a vector of terms.
  gsePosDotplot <- dotplot(gsePos, showCategory=30,
                           label_format = 45) + 
    ggtitle("GSEA up-regulated genes in senescent damaged cells")
  
  gsePos@result
  png( height = 16, width = 6, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DSvsDnS_gsePositiveMarkers.png") )
  
  gsePosDotplot
  
  dev.off()
  
  
  
  # negative markers # 
 
   d <- filter(SenescenceDamageMarkers,
              avg_log2FC < 0 ) %>% 
    arrange(.,by_group = avg_log2FC ) %>% 
    select(gene, avg_log2FC)
  
  mappingTemp <- bitr(d$gene, 
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)
  d[mappingTemp$SYMBOL,'gene'] <- mappingTemp$ENTREZID
  
  geneList <- abs(d[,2])
  
  ## feature 2: named vector
  names(geneList) <- as.character(d[,1])
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  ## feature 3: decreasing order 
  
  gseNeg <- gseGO(geneList     = geneList,
                  OrgDb        = org.Mm.eg.db,
                  ont          = "ALL",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  
  
  gseNegDotplot <- dotplot(gseNeg, showCategory = 30) + 
    ggtitle("GSEA up-regulated genes in not-senescent damaged cells")
  
  png( height = 16, width = 6, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DSvsDnS_gseNegativeMarkers.png") )
  
  gseNegDotplot 
  
  dev.off()
  
}

## IV. Plot expression of some interesting genes
# choose top genes from the single cell and bulk differential expression 
# analysis 
{  
  DEGs_up_SenDamaged_v_notSen <- 
    VlnPlot(testData,
          features = c("Igf1","Serpine1","Nrg1","C3","Vegfa","Igfbp4","Lcp1"), 
          split.by = "group", 
          stack = TRUE,
          layer = "data") + 
    scale_fill_manual(values = c('Healthy' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Damaged, not-senescent' = '#414487FF',
                                 'Damaged and senescent' = '#440154FF')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          title) + ggtitle("LFC > 0.5 in DS vs. DNS: genes overlap with SenMayo gene set")
  
  DEGs_down_Sen_vs_notSenBileAcidGenes <- 
    VlnPlot(testData,
            features = c(	"Cyp7a1","Cyp27a1","Sult2a8","Slco1b2","Ces1d", "Ces1f","Ces1c","Akr1c6","Pygl"), 
            split.by = "group", 
            stack = TRUE,
            layer = "data") + 
    scale_fill_manual(values =c('Undamaged' = '#FDE725FF',
                                'Transition' = '#7AD151FF',
                                'Not-senescent-like, damaged' = '#414487FF',
                                'Senescent-like, damaged' = '#440154FF')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          title) 
  
  
  png( height = 3.5, width = 8, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_DEGs_down_Sen_vs_notSen_BileAcid.png") )
  
  DEGs_down_Sen_vs_notSenBileAcidGenes
  
  dev.off()
  

  # Neg. DEGs in DNS |log2FC| > 2
  
  
  DNSplotNegative <- VlnPlot(testData,
                             features =  c('Grm8', 
                                           'Serpina1e','Tiam2',
                                           'Sult5a1', 'Selenbp2'), 
                             split.by = "group", 
                             stack = TRUE,
                             layer = "data",
                             log = TRUE
  ) +
    scale_fill_manual(values = c('Healthy' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Damaged, not-senescent' = '#414487FF',
                                 'Damaged and senescent' = '#440154FF'))
  
  png( height = 3, width = 9, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_DNS_cells_Abslog2FCover2_.png") )
  print(DNSplot)
  dev.off()
  
  DNSplotConditions <- VlnPlot(testData,
                               features =  c('Grm8', 
                                             'Serpina1e','Tiam2',
                                             'Sult5a1', 'Selenbp2'), 
                               group.by = "group",
                               split.by = "condition",
                               stack = TRUE,
                               layer = "data") 
  
  png( height = 4, width = 9, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_DNS_cells_Abslog2FCover2_Conditions.png") )
  print(DNSplotConditions)
  dev.off()
  
  # Pos and Neg. DEGs in Healthy |log2FC| > 2 (top 3 neg)
  
  Healthyplot <- VlnPlot(testData,
                         features =  c('Serpina1e','Selenbp2',
                                       'Apoa4', 'Rcan2',
                                       'Nrg1','Ephb2',
                                       'Slc22a3'), 
                         group.by = "group",
                         stack = TRUE,
                         layer = "data")
  
  HealthyplotConditions <- VlnPlot(testData,
                                   features =  c('Serpina1e','Selenbp2',
                                                 'Apoa4', 'Rcan2', 
                                                 'Nrg1','Ephb2',
                                                 'Slc22a3'), 
                                   group.by = "group",
                                   split.by = 'condition', 
                                   stack = TRUE,
                                   layer = "data")
  
  png( height = 4, width = 9, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_Healthy_cells_Abslog2FCover2.png") )
  print(Healthyplot)
  dev.off()
  
  png( height = 4, width = 9, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_Healthy_cells_Abslog2FCover2_Conditions.png") )
  print(HealthyplotConditions)
  dev.off()
  
  
  
  #Top 5 Positively DEGs in DS cells |log2FC| > 2
  
  upP <- VlnPlot(testData,
                 features =head(
                   filter(allMarkers, 
                          avg_log2FC >= 2 &
                            pct.1 >= 0.05 &
                            cluster == 'DamagedSenescent')$gene, 5), 
                 split.by = "group", 
                 stack = TRUE,
                 layer = "data") + 
    scale_fill_manual(values = c('Healthy' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Damaged, not-senescent' = '#414487FF',
                                 'Damaged and senescent' = '#440154FF'))  
  
  
  
  downP <- VlnPlot(testData,
                   features = head(
                     filter(allMarkers, 
                            avg_log2FC <= -2 &
                              pct.2 >= 0.05 &
                              cluster == 'DamagedSenescent')$gene, 5), 
                   split.by = "group", 
                   stack = TRUE,
                   layer = "data") + 
    scale_fill_manual(values = c('Healthy' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Damaged, not-senescent' = '#414487FF',
                                 'Damaged and senescent' = '#440154FF'))  
  ggtitle('Top sig. down-regulated genes in DS hepatocytes (log2FC =< -2)')
  
  
  
  png( height = 5, width = 10, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_DS_cells_log2FCover2_Upregulated.png") )
  print(upP)
  dev.off()
  
  
  png( height = 5, width = 10, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_DS_cells_log2FCover2_Downregulated.png") )
  print(downP)
  dev.off()
  
  
  png( height = 8, width = 9, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/DEGs_DS_cells_Abslog2FCover2.png") )
  print((upP/downP))
  dev.off()
  
  
  ### sames plots but separate cell by conditions 
  
  VlnPlot(testData,
          features =head(
            filter(allMarkers, 
                   avg_log2FC <= -2 &
                     pct.2 >= 0.05 &
                     cluster == 'DamagedSenescent')$gene), 
          group.by = "condition", 
          split.by = 'group', stack = TRUE,
          layer = "data") 
  
  VlnPlot(testData,
          features =head(
            filter(allMarkers, 
                   avg_log2FC >= 2 &
                     pct.1 >= 0.05 &
                     cluster == 'DamagedSenescent')$gene), 
          group.by = "condition", split.by = 'group', stack = TRUE,
          layer = "data") 
  
  
  
  # Comparison of raw counts, scaled/variance transformed 
  # and SCTransformed counts when plotting marker genes
  
  
  rawplot <- VlnPlot(testData,
                     features =head(
                       filter(allMarkers, 
                              avg_log2FC <= -2 &
                                pct.2 >= 0.05 &
                                cluster == 'DamagedSenescent')$gene), 
                     group.by = "group", stack = TRUE,
                     layer = "counts") + ggtitle('Raw Counts')
  
  
  scalePlot <- VlnPlot(testData,
                       features =head(
                         filter(allMarkers, 
                                avg_log2FC <= -2 &
                                  pct.2 >= 0.05 &
                                  cluster == 'DamagedSenescent')$gene), 
                       group.by = "group", stack = TRUE,
                       layer = "scale.data") + ggtitle('scaled counts')
  
  
  
  SCTplot <- VlnPlot(testData,
                     features =head(
                       filter(allMarkers, 
                              avg_log2FC <= -2 &
                                pct.2 >= 0.05 &
                                cluster == 'DamagedSenescent')$gene), 
                     group.by = "group", stack = TRUE,
                     layer = "data") + ggtitle('SCTransformed counts')
  
  patchPlots <- (rawplot/scalePlot/SCTplot) 
  
  
  png( height = 9, width = 12, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/CompareCountTransform_Top5DownregGenes_DamagedSenescentHep.png") )
  print(patchPlots)
  dev.off()
} 


# ROS stress markers per cell fate: do we see more in senescent hepatocytes? 
# how is the marker of these genes across HDS, particularly around the 
# HDS tipping point 

{
  
  # Gene set 1
 #  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
  
 Hallmark_ROS <- (read.csv(
   file = paste0('hepatocyte-damage-score/Data/Input/GeneSetsCellFates/',
                 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.v2024.1.Mm.tsv'),
   sep = '\t',
   header = TRUE
      ))
  
  Hallmark_ROS <- Hallmark_ROS[17,2]
  
  Hallmark_ROS <- unlist(strsplit(Hallmark_ROS,
                  split = ","))
  
    ROSmarkersViolin <- 
    VlnPlot(testData,
            features = Hallmark_ROS[1:5], 
            split.by = "group", 
            stack = TRUE,
            layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                 'Transition' = '#7AD151FF',
                                 'Not-senescent-like, damaged' = '#414487FF',
                                 'Senescent-like, damaged' = '#440154FF')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          title) + ggtitle("ROS")
    
    
    cells_rankingsAll <- AUCell_buildRankings(mergedObjExpr, 
                                              plotStats = TRUE)
    
    HallmarkROS_AUCell <- AUCell_calcAUC(Hallmark_ROS,
                                   cells_rankingsAll)
    
    if(identical(colnames(HallmarkROS_AUCell@assays@data$AUC), 
                 rownames(mergedObj@meta.data)) == TRUE){
      print('all good')
      mergedObj@meta.data$AUCell_HallmarkROS <- c(HallmarkROS_AUCell@assays@data$AUC)
    }
    
    VlnPlot(mergedObj,
            features = 'AUCell_HallmarkROS', 
            group.by = "group") + 
      scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                   'Transition' = '#7AD151FF',
                                   'Not-senescent-like, damaged' = '#414487FF',
                                   'Senescent-like, damaged' = '#440154FF')) +
      theme(axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            title) + ggtitle("Hallmark ROS Pathway AUCell Scores") +geom_violin(trim = TRUE) 
    
    
    ggplot(data = mergedObj@meta.data, 
           aes(x=group, y = AUCell_HallmarkROS) ) +
      geom_violin(trim = TRUE) +
      theme(axis.text.x = element_text(angle = 90)) +
      xlab("Cell Fates") 
    
    plot(mergedObj@meta.data$HDS, mergedObj@meta.data$AUCell_HallmarkROS)
    
    VlnPlot(mergedObj,
            features = 'AUCell_SenMayo', 
            group.by = "group") + 
      scale_fill_manual(values = c('Undamaged' = '#FDE725FF',
                                   'Transition' = '#7AD151FF',
                                   'Not-senescent-like, damaged' = '#414487FF',
                                   'Senescent-like, damaged' = '#440154FF')) +
      theme(axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            title) + ggtitle("Hallmark ROS Pathway AUCell Scores")
    
  
  
}