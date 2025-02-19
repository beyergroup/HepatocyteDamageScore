# try different thresholds
# 0.9 quantile

library(clusterProfiler)
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

{
  outputPath <- "hepatocyte-damage-score/Data/Output/"
  
  mergedObj<- readRDS(paste0(outputPath,
                        "MergedXiaoHDS_SenMayoHepatocytes.rds"))
  
  # define column 'group' in metadata with cell fate  according to 
  # senMayo AUCell score and HDS 
  mergedObj@meta.data$CellID <- rownames(mergedObj@meta.data)
  mergedObj@meta.data$group <- NA
  
  quantilesHDSall <- quantile(mergedObj@meta.data$HDS, probs = c(0.25, 0.5, 0.75))
  quantiles9NASH <- quantile(
    mergedObj@meta.data$AUCell_SenMayo[mergedObj@meta.data$condition == "9m NASH"],
    probs = c(0.25, 0.5, 0.75))
  
  # tag cells to cell fate by HDS and AUCell_SenMayo thresholds:
  # possible cell fates: healthy, damaged + senescent, and damaged + not senescent
  
  # threshold healthy: HDS <= 1. quartile 
  tempIndex <- mergedObj@meta.data$HDS <= quantilesHDSall[1]
  mergedObj@meta.data$group[tempIndex] <- 'Healthy'
  remove(tempIndex)
  
  # threshold: damaged + senescent
  # HDS > 3. quartiles & SenMayo > 3. quartile
  tempIndex <- 
    mergedObj@meta.data$AUCell_SenMayo > quantiles9NASH[3] &
    mergedObj@meta.data$HDS > quantilesHDSall[3]
  
  mergedObj@meta.data$group[tempIndex] <- 'DamagedSenescent'
  remove(tempIndex)
  
  # threshold: damaged + not senescent
  # HDS > 3. quartiles & SenMayo < 2.quartile  
  
  tempIndex <- 
    mergedObj@meta.data$AUCell_SenMayo <= quantiles9NASH[2] &
    mergedObj@meta.data$HDS >  quantilesHDSall[3]
  mergedObj@meta.data$group[tempIndex] <- 'DamagedNotSenescent'
  remove(tempIndex)
  
  # intermediate cell fate
  
  tempIndex <- 
    mergedObj@meta.data$AUCell_SenMayo > quantiles9NASH[2] &
    mergedObj@meta.data$AUCell_SenMayo <= quantiles9NASH[3] &
    mergedObj@meta.data$HDS >  quantilesHDSall[3]
  mergedObj@meta.data$group[tempIndex] <- 'unresolved'
  remove(tempIndex)
  
  # tag NA's (cells that do not fulfill any of the criteria) as transition cells
  mergedObj@meta.data$group <- mergedObj@meta.data$group %>% replace_na('Transition')
  
  mergedObj@meta.data$group <- factor(mergedObj@meta.data$group, 
                                      ordered = TRUE,
                                      levels = c('Healthy', 
                                                 'Transition',
                                                 'DamagedNotSenescent',
                                                 'unresolved',
                                                 'DamagedSenescent'),
                                      labels = c('Undamaged',
                                                 'Transition',
                                                 'Damaged-Not-Senescent',
                                                 'Damaged-Unresolved',
                                                 'Damaged-Senescent'))
  
  
  
  umapFate<-
    DimPlot(mergedObj, 
            group.by = 'group',
            pt.size = 1,
            cells = sample(rownames(mergedObj@meta.data),
                           replace = FALSE)) + 
    scale_color_manual(values = c('Undamaged' = '#CCBB44',
                                  'Transition' = '#228833',
                                  'Damaged-Not-Senescent' = '#4477AA',
                                  'Damaged-Unresolved' = '#EE6677',
                                  'Damaged-Senescent' = '#AA3377'))+
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          legend.key.height=unit(5,"mm"),
          legend.key.width = unit(5, "mm"),
          legend.text = element_text(size = 6))



  png( height = 4.5, width = 4,
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/UmapXiao_CellFateClassificationUpdatedThresholds.png") )
  print(umapFate)
  dev.off()
  
  pdf( height = 4.5, width = 4,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/UmapXiao_CellFateClassificationUpdatedThresholds.pdf") )
  print(umapFate)
  dev.off()
  
  mergedObj@meta.data$condition <- factor(mergedObj@meta.data$condition,
                                          ordered = TRUE,
                                          levels = c("3m NC","9m NC", 
                                                     "3m NASH", "9m NASH" ))

  
  # Summarize the data: count the number of cells for each fate in each sample
  summarizedCellFates <- mergedObj@meta.data %>%  
    dplyr::mutate(condition = factor(condition, 
                              levels = c('3m NC',
                                         '9m NC',
                                         '3m NASH',
                                         '9m NASH')),
           group = factor(group,
                          levels = c('Undamaged',
                                     'Transition',
                                     'Damaged-Not-Senescent',
                                     'Damaged-Unresolved',
                                     'Damaged-Senescent'))) %>%
    dplyr::group_by(condition, group) %>%
    dplyr::summarise(Cell_Count = n(), .groups = "drop") %>%
    dplyr::group_by(condition) %>% 
    dplyr::mutate(proportion = Cell_Count / sum(Cell_Count))
  
  #write.csv(summarizedCellFates[,1:4], file = paste0(outputPath,
   #                                                  "Results/CellFateAnalysis/CellFateCountsXiao.csv"))
  
  # Prepare a table-like data frame for annotations
  
  label_data <- summarizedCellFates %>% 
    dplyr::group_by(condition) %>%
    dplyr::mutate(n = sum(Cell_Count)) %>%
    dplyr::select(condition, group, Cell_Count, n) %>%
    pivot_wider(names_from = group, values_from = Cell_Count, values_fill = 0)
  
  # Convert the data frame to a long format for plotting as a table
  label_data_long <- label_data[,-2] %>%
    pivot_longer(
      cols = -condition,
      names_to = "group",
      values_to = "Cell_Count"
    ) %>% 
    dplyr::mutate(group = factor(group,
                          levels = c('Undamaged',
                                     'Transition',
                                     'Damaged-Not-Senescent',
                                     'Damaged-Unresolved',
                                     'Damaged-Senescent')))
  
  # Create a stacked bar plot
  stackedBarPlotCellFates <- 
    ggplot(summarizedCellFates, aes(x = condition, 
                                    y = proportion, 
                                    fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
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
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "mm"),
      legend.title = element_blank()
    ) + geom_text(
      data = label_data,
      aes(x = condition, y = 1.05, label = n),  # Place labels above the bar
      # fontface = "bold",
      size = 7,
      inherit.aes = FALSE  # Prevent unwanted aesthetics from being inherited
    )
  
  png( height = 6, width = 8, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/Xiao_StackedBarplotCellFates_Percentages_NewThreshold.png") )
  
  print(stackedBarPlotCellFates)
  
  dev.off()
  
  pdf( height = 6, width = 8, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/Xiao_StackedBarplotCellFates_Percentages_NewThreshold.pdf") )
  
  print(stackedBarPlotCellFates)
  
  dev.off()
  
  
  
} 

# 
{
  
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
  
  # allMarkers <- FindAllMarkers(testData) 
  # allMarkers <- filter(allMarkers,
  #                      p_val_adj < 0.01 &  abs(avg_log2FC) >= 1)
  # 
  # differentially expressed genes in pairwise comparisons of all groups
  
  SenescenceDamageMarkers <- FindMarkers(testData,
                                         ident.1 =  'Damaged-Senescent',
                                         ident.2 = 'Damaged-Not-Senescent')
  View(SenescenceDamageMarkers)

  
  SenescenceDamageMarkers$gene <- rownames(SenescenceDamageMarkers) 
  
  positiveMarkers <- filter(SenescenceDamageMarkers,
                            p_val_adj < 0.05 &
                              avg_log2FC > 0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  
  positiveMarkers <- positiveMarkers[order(positiveMarkers$avg_log2FC, decreasing = TRUE),]
  
  intersect(positiveMarkers$gene, senMayoGenes$Gene.murine.)
  
  intersect(positiveMarkers$gene, HDAG$gene_symbol[1:42])
  intersect(positiveMarkers$gene, senMayoGenes$Gene.murine.)
  
  venn.plot <- venn.diagram(
    x = list(Set1 = HDAG[1:42,"gene_symbol"], Set2 = senMayoGenes$Gene.murine.,
             Set3 = positiveMarkers$gene),
    category.names = c("HDS 42 genes","SenMayo genes", "Positive marker genes"),
    filename = NULL,  # Set to NULL to plot in R instead of saving to a file
    output = TRUE,
    main = "Intersection of gene sets"
  )
  dev.off()
  grid.draw(venn.plot)
  
  
  
  # NEGATIVE MARKERS
  
  negativeMarkers <- filter(SenescenceDamageMarkers,
                            p_val_adj < 0.05 &
                              avg_log2FC < -0.5 & ((pct.1 >= 0.1) | (pct.2 >= 0.1)))
  
  intersect(negativeMarkers$gene, senMayoGenes$Gene.murine.)
  
  intersect(negativeMarkers$gene, HDAG$gene_symbol[1:42])
  
  negativeMarkers <- negativeMarkers[order(negativeMarkers$avg_log2FC, 
                                           decreasing = FALSE),]
  
  head(formatC(positiveMarkers$p_val, format = "e", digits = 2))
  
  
}

# Save differentially expressed tables for supplementary table in paper 

{
  tempData <- data.frame('Gene Name' = positiveMarkers$gene,
                         'Adj. p value' = formatC(positiveMarkers$p_val_adj, 
                                                  format = "e",
                                                  digits = 2),
                         'Avg. log2FC' = round(positiveMarkers$avg_log2FC, 
                                               digits = 2),
                         'Percent group 1' = round(positiveMarkers$pct.1, 
                                                               digits = 2),
                         'percent group 2' = round(positiveMarkers$pct.2, 
                                                                digits = 2))

  write.csv( tempData,
             file = paste0(outputPath,
                           'PositiveMarkers_DSvsDNS_NewThreshold.csv'))
  
  remove(tempData)
  
  
  tempData <- data.frame('Gene Name' = negativeMarkers$gene,
                         'Adj. p value' = formatC(negativeMarkers$p_val_adj, 
                                                  format = "e",
                                                  digits = 2),
                         'Avg. log2FC' = round(negativeMarkers$avg_log2FC, 
                                               digits = 2),
                         'Percent group 1' = round(negativeMarkers$pct.1, 
                                                   digits = 2),
                         'percent group 2' = round(negativeMarkers$pct.2, 
                                                   digits = 2))
  
  write.csv( tempData,
             file = paste0(outputPath,
                           'NegativeMarkers_DSvsDNS_NewThreshold.csv'))
  remove(tempData)
  
}


# ORA of marker genes

{
  egoPosBP <- enrichGO(gene = positiveMarkers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  egoPosMF <- enrichGO(gene = positiveMarkers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  egoPosBP_simplified <- simplify(
    egoPosBP,
    cutoff = 0.7,             # Similarity threshold (0.7 is a common default)
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
}



# Violin Plots Up-regulated Genes

{
  plotposMarkersIntersect <- VlnPlot(testData,
                           features = intersect(positiveMarkers$gene, 
                                                senMayoGenes$Gene.murine.), 
                           split.by = "group", 
                           stack = TRUE,
                           layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 12),
          legend.position = "none")
  
  # Genes annotated positive regulation of ERK1/ERK2 cascade:  Nrg1, C3, Jun, Vegfa and Igf1 
  
  plotposMarkersErK<- VlnPlot(testData,
                                     features = c("Nrg1", "C3",
                                                  "Jun", "Vegfa",
                                                  "Igf1"), 
                                     split.by = "group", 
                                     stack = TRUE,
                                     layer = "data",
                              flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10,
                                     angle = 90))
  
  png( height = 4, width = 4, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_PositiveMarkers_ERK1_2_positiveRegulation.png") )
  
  print(plotposMarkersErK)
  
  dev.off()
  
  
  pdf( height = 4, width = 4, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_PositiveMarkers_ERK1_2_positiveRegulation.pdf") )
  
  print(plotposMarkersErK)
  
   dev.off()
   
   
  # Lipid metabolims 
   plotposLipid<- VlnPlot(testData,
                               features = c("Spp1", "Apoa4",
                                             "Scd1","Rtn4"), 
                               split.by = "group", 
                               stack = TRUE,
                               layer = "data",
                          flip = TRUE) + 
     scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                  'Transition' = '#228833',
                                  'Damaged-Not-Senescent' = '#4477AA',
                                  'Damaged-Unresolved' = '#EE6677',
                                  'Damaged-Senescent' = '#AA3377')) +
     theme(axis.title.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = "none",
           legend.direction = "horizontal",
           legend.text = element_text(size = 5,
                                      angle = 90))
   
   png( height = 4, width = 4, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_PositiveMarkers_LipidGOterms.png") )
   
   print(plotposLipid)
   
   dev.off()
   
   
   pdf( height = 4, width = 4, 
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/",
                      "Xiao_PositiveMarkers_LipidGOterms.pdf") )
   
   print(plotposLipid)
   
   dev.off()
   
   # GO wound healing 
   plotposWoundHealing<- VlnPlot(testData,
                          features = c("Serpine1", "Fgl1",
                                       "S100a10","Map3k1"), 
                          split.by = "group", 
                          stack = TRUE,
                          layer = "data",
                          flip = TRUE) + 
     scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                  'Transition' = '#228833',
                                  'Damaged-Not-Senescent' = '#4477AA',
                                  'Damaged-Unresolved' = '#EE6677',
                                  'Damaged-Senescent' = '#AA3377')) +
     theme(axis.title.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = "none",
           legend.direction = "horizontal",
           legend.text = element_text(size = 5,
                                      angle = 90))
   
   png( height = 4, width = 4, 
        units = 'in',
        res = 300,
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/PNGs/",
                      "Xiao_PositiveMarkers_WoundHealingGOterms.png") )
   
   print(plotposWoundHealing)
   
   dev.off()
   
   
   pdf( height = 4, width = 4, 
        file = paste0(outputPath,
                      "Results/CellFateAnalysis/",
                      "Xiao_PositiveMarkers_WoundHealingGOterms.pdf") )
   
   print(plotposWoundHealing)
   
   dev.off()
   
   
   
  
  
  
  
  
  
}

# ORA Down-regualated genes 
{
  egoNegBP <- enrichGO(gene = negativeMarkers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  egoNegMF <- enrichGO(gene = negativeMarkers$gene,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  egoNegBP_simplified <- simplify(
    egoPosBP,
    cutoff = 0.7,             # Similarity threshold (0.7 is a common default)
    by = "p.adjust",          # Retain the most significant term
    select_fun = min,         # Select term with the smallest p-value within each cluster
    measure = "Wang"          # Semantic similarity measure ("Wang" is common for GO terms)
  )
  
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
  
  
}

# violin plots downregulated genes
{
  # GO BP: lipid catabolic process 
  lipidCat <-c("Lipg","Akr1c6","Hpgd", "Akr1c20","Abhd6",
               "Gpcpd1","Plcg1", "Irs1","Ldlr", "Prdx6")
  
  # phospholidpid catabolic process
  phospholipidCat <- c("Abhd6","Gpcpd1","Plcg1","Ldlr","Prdx6")
  
  # GO BP terms: bile acid metabolic process
  # transport and secretion 
  bileAcidGo<- c("Cyp27a1","Cyp1a2","Sult2a8","Ces1d","Ces1c","Ces1f","Akr1d1")
  
  # bile acid and bile salt transport 
  bileAcidTransport <- c("Slco1a1","Slc10a1","Aqp9","Slco1b2")
  
  
  #  lipid catabolims genes (also include phospholipid cat.)
  
  plotNegLipidCat <- VlnPlot(testData,
                             features = lipidCat, 
                             split.by = "group", 
                             stack = TRUE,
                             layer = "data",
                             flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10,
                                     angle = 90))
  
 
   png( height = 10, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_NegativeMarkers_LipidCat.png") )
  
  print(plotNegLipidCat)
  
  dev.off()
  
  
  pdf( height = 10, width = 4, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_NegativeMarkers_LipCat.pdf") )
  
  print(plotNegLipidCat)
  
  dev.off()
  
  # Bile acid: includes the two rate limiting enzymes, and genes annotated with
  
  plotNegBileAcidMet <- VlnPlot(testData,
                                features = c("Cyp7a1",
                                             bileAcidGo,
                                             bileAcidTransport), 
                                split.by = "group", 
                                stack = TRUE,
                                layer = "data",
                                flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10,
                                     angle = 90))
  
  
  png( height = 15, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_NegativeMarkers_BileAcids.png") )
  
  print(plotNegBileAcidMet)
  
  dev.off()
  
  
  pdf( height = 15, width = 4, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_NegativeMarkers_BileAcids.pdf") )
  
  print(plotNegBileAcidMet)
  
  dev.off()
  
  # only transport gene 
  plotNegBileAcidTransport <- VlnPlot(testData,
                                features = bileAcidTransport, 
                                split.by = "group", 
                                stack = TRUE,
                                layer = "data",
                                flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 12),
          legend.position = "none")
  
  
  
  plotNegPhospholipid <- VlnPlot(testData,
                                features = phospholipidCat, 
                                split.by = "group", 
                                stack = TRUE,
                                layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 12),
          legend.position = "none")
  
  
  # amino acid metabolic process 
  
  aminoacidMet <- c("Glul","Prodh","Ahcy","Sardh",
  "Hgd","Aldh1a1","Ido2","Kyat3","Gclc","Cth","Kmo",
  "Prodh2","Bckdha","Gnmt")
  
  plotNegAminoAcidMet <- VlnPlot(testData,
                                 features = aminoacidMet, 
                                 split.by = "group", 
                                 stack = TRUE,
                                 layer = "data",
                                 flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10,
                                     angle = 90))
  
  png( height = 15, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_NegativeMarkers_AminoAcid.png") )
  
  print(plotNegAminoAcidMet)
  
  dev.off()
  
  
  pdf( height = 15, width = 4, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_NegativeMarkers_AminoAcid.pdf") )
  
  print(plotNegAminoAcidMet)
  
  dev.off()
  
  # GO BP term: xenobiotic metabolic process
  
  
  
  xenobioticGenes <- c("Cyp2c29","Cyp2c54","Cyp1a2","Cyp2c50",
    "Sult2a8","Cyp2d9","Ugt2b1","Cyp2a5","Cyp2c38",
    "Cyp2e1","Ahr","Cyp3a11","Pon3","Ugt2b5","Gsta3","Rorc",
    "Acsl1","Fmo1","Fmo5")
  
  xenobioticSelectedGenes <- c("Cyp2c54",
                               "Ahr",
                               "Cyp2c38",
                               "Rorc",
                               "Fmo1",
                               "Cyp2d9",
                               "Cyp2c29",
                               "Ugt2b1",
                               "Cyp2e1","Cyp3a11",
                               "Ugt2b5","Gsta3","Rorc",
                               "Acsl1")
  
  plotNegXenoSelect <- VlnPlot(testData,
                                 features = xenobioticSelectedGenes, 
                                 split.by = "group", 
                                 stack = TRUE,
                                 layer = "data",
                                 flip = TRUE) + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 10,
                                     angle = 90))
  
  png( height = 10, width = 4, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_NegativeMarkers_XenobioticMet.png") )
  
  print(plotNegXenoSelect)
  
  dev.off()
  
  
  pdf( height = 10, width = 4, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_NegativeMarkers_XenobioticMet.pdf") )
  
  print(plotNegXenoSelect)
  
  dev.off()
  
  
}


# selected upregulated genes
{
  
  # metadata column order 
  testData@meta.data$group <- factor(mergedObj@meta.data$group,
                                     ordered = TRUE,
                                     levels = c("Damaged-Senescent",
                                                'Damaged-Unresolved',
                                                "Damaged-Not-Senescent", 
                                                "Transition",
                                                "Undamaged" ))
  
  testData@meta.data$group <- factor(mergedObj@meta.data$group,
                                     ordered = TRUE,
                                     levels = c("Undamaged",
                                                "Transition",
                                                "Damaged-Not-Senescent",
                                                'Damaged-Unresolved',
                                                "Damaged-Senescent" ))
  

  Idents(testData) <- factor(Idents(testData), levels= c("Undamaged",
                                                         "Transition",
                                                         "Damaged-Not-Senescent",
                                                         'Damaged-Unresolved',
                                                         "Damaged-Senescent" ))
  

  #  constitutive androstane receptor (CAR) (NR1I3)
  
  selectedDownLong <-c("Cyp7a1","Cyp27a1","Sult2a8", "Nr1i3","Prdx6","Gpcpd1",
               "Irs1","Cyp1a2", "Cyp2c29", "Ahcy", "Sardh")
  
  selectedDownSmall <-c("Cyp7a1","Cyp27a1","Sult2a8", "Nr1i3","Prdx6","Gpcpd1",
                        "Irs1","Cyp1a2", "Cyp2c29", "Ahcy", "Sardh")
  
  
  plotSelectedDownLong <- VlnPlot(testData,
                                 features = selectedDownLong, 
                                 split.by = "group", 
                                 stack = TRUE,
                                 layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  png( height = 4, width = 10, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_NegativeMarkers_SelectedLong.png") )
  
  print(plotSelectedDownLong)
  
  dev.off()
  
  
  pdf( height = 4, width = 10, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_NegativeMarkers_SelectedLong.pdf") )
  
  print(plotSelectedDownLong)
  
  dev.off()
  
  # upregulated genes 
  
  upregulatedGenesSelected <- c("Spp1", "Rtn4", "Scd1", "Nrg1", "C3", "Jun","Vegfa", "Igf1")
  
  
  plotSelectedUp <- VlnPlot(testData,
                                  features = upregulatedGenesSelected, 
                                  split.by = "group", 
                                  stack = TRUE,
                                  layer = "data") + 
    scale_fill_manual(values = c('Undamaged' = '#CCBB44',
                                 'Transition' = '#228833',
                                 'Damaged-Not-Senescent' = '#4477AA',
                                 'Damaged-Unresolved' = '#EE6677',
                                 'Damaged-Senescent' = '#AA3377')) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  png( height = 4, width = 10, 
       units = 'in',
       res = 600,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "Xiao_PositiveMarkers_Selected.png") )
  
  print(plotSelectedUp)
  
  dev.off()
  
  
  pdf( height = 4, width = 10, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "Xiao_PositiveMarkers_Selected.pdf") )
  
  print(plotSelectedUp)
  
  dev.off()
  
}
