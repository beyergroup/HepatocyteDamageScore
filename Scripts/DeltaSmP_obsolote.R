# obsolete analysis cell fate delta smp

#  Read in pseudo-proliferation index genes 
# list originally human --> converted to mouse orthologs manually using MGI
{
  proliferationIndex <- read.csv(file = 'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/Pseudo-Proliferation-Index-Genes_ConvertedToMouseByPaula.csv',
                                 sep = ';',
                                 header = TRUE,
                                 stringsAsFactors = FALSE)
  
  proliferationIndex <- proliferationIndex[-c(1,2),]
  colnames(proliferationIndex) <- c('human_gene_names', 'UNIPROT_ID', 'mouse_gene' )
  
  # Do they intersect with HDS?
  
  intersect(HDAG$gene_symbol[1:42], proliferationIndex$mouse_gene)
  # no intersection
}


# Proliferation Index AUCell score calculation
{
  
  cells_rankingsAll <- AUCell_buildRankings(mergedObjExpr, 
                                            plotStats = TRUE)
  
  Prolcells_All <- AUCell_calcAUC(proliferationIndex$mouse_gene,
                                  cells_rankingsAll)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	1 (4% of 25)
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/AUCellExploreThresholds_ProlifIndex_MergedXiaoSamples.pdf") )
  
  AUCell_exploreThresholds(Prolcells_All, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  
  png( height = 6, width = 12, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/AUCellExploreThresholds_ProlifIndex_MergedXiaoSamples.png") )
  
  AUCell_exploreThresholds(Prolcells_All, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  if(identical(colnames(Prolcells_All@assays@data$AUC), 
               rownames(mergedObj@meta.data)) == TRUE){
    print('all good')
    mergedObj@meta.data$AUCell_ProlifIndex <- c(Prolcells_All@assays@data$AUC)
  }
  
  
}


# DaHep pre-cancerous state signature
{
  
  # daHep Markers (pre-cancerous cells)
  daHep <- read.csv(file = 
                      'hepatocyte-damage-score/Data/Input/GeneSetsCellFates/Carlessi_daHepMakers.csv')
  
  
  daHepUp <- daHep$X[daHep$avg_log2FC > 0 & abs(daHep$avg_log2FC) > 0.5 & daHep$p_val_adj < 0.05 ]
  
  daHepDown <- daHep$X[daHep$avg_log2FC < 0 & abs(daHep$avg_log2FC) > 0.5 & daHep$p_val_adj < 0.05]
  
  intersect(HDAG$gene_symbol[1:42], daHepDown)
  
  cells_daHepUp <- AUCell_calcAUC(daHepUp,
                                  cells_rankingsAll)
  
  # Genes in the gene sets NOT available in the dataset: 
  #   geneSet: 	11 (4% of 310)
  
  # Plot AUC histogram
  
  pdf( height = 6, width = 12, 
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/",
                     "AUCellExploreThresholds_DaHepUp_MergedXiaoSamples.pdf") 
  )
  
  AUCell_exploreThresholds(cells_daHepUp, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  
  png( height = 6, width = 12, 
       units = 'in',
       res = 300,
       file = paste0(outputPath,
                     "Results/CellFateAnalysis/PNGs/",
                     "AUCellExploreThresholds_DaHepUp_MergedXiaoSamples.png") 
  )
  
  AUCell_exploreThresholds(cells_daHepUp, 
                           plotHist = TRUE, 
                           assign = TRUE)
  dev.off()
  
  
  if(identical(colnames(cells_daHepUp@assays@data$AUC), 
               rownames(mergedObj@meta.data)) == TRUE){
    print('all good')
    mergedObj@meta.data$AUCell_daHepUp <- c(cells_daHepUp@assays@data$AUC)
  }
  
  FeaturePlot(mergedObj,
              features = 'AUCell_daHepUp')  
  
  
  
}



# Calculate delta AUCell SenMayo - AUCell Prolifertion Index (deltaSmP)
# assigne a cell fate to hepatocytes (either senscencent or proliferative)
# based con a deltaSmP_centered threshold 

{

mergedDataFrame <- mergedObj@meta.data


mergedDataFrame$CellType <- factor(mergedDataFrame$CellType,
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
  mutate( deltaSmP = AUCell_SenMayo - AUCell_ProlifIndex)
mergedDataFrame$deltaSmP_centered <- 
  mergedDataFrame$deltaSmP - mean(mergedDataFrame$deltaSmP)

# distribution of deltaSmP (centered)

png(height = 4,
    width = 6,
    units = 'in',
    res = 300,
    file =
      paste0(outputPath,
             "Results/CellFateAnalysis/PNGs/",
             "DeltaSmPcentered_Histogram.png"))

quantiles <- quantile(mergedDataFrame$deltaSmP_centered, probs = c(0.1,0.25, 0.5, 0.75,0.9))

ggplot(mergedDataFrame, aes(x = deltaSmP_centered)) + 
  geom_histogram(binwidth = 0.002) +
  geom_vline(aes(xintercept = quantiles[1]), color = 'red') +
  geom_vline(aes(xintercept = quantiles[5]), color = 'red') + 
  geom_text(aes(quantiles[1]), 
            y = Inf, 
            label = paste0('10%-Q: ',round(quantiles[1],2)), 
            vjust = 2, 
            hjust = 1,
            size = 4) +
  geom_text(aes(quantiles[5]), 
            y =Inf, 
            label = paste0('90%-Q: ',round(quantiles[5], 2)),
            vjust = 2, 
            hjust = 0,
            size = 4)

dev.off()

png(height = 4,
    width = 6,
    units = 'in',
    res = 300,
    file =
      paste0(outputPath,
             "Results/CellFateAnalysis/PNGs/",
             "DistributionAUCellSenMayo_Histogram.png"))

quantilesSenMayo <- quantile(mergedDataFrame$AUCell_SenMayo, probs = c(0.1,0.25, 0.5, 0.75,0.9))

ggplot(mergedDataFrame, aes(x = AUCell_SenMayo)) + 
  geom_histogram(binwidth = 0.002) +
  geom_vline(aes(xintercept = quantilesSenMayo[1]), color = 'red') +
  geom_vline(aes(xintercept = quantilesSenMayo[5]), color = 'red') + 
  geom_text(aes(quantilesSenMayo[1]), 
            y = Inf, 
            label = paste0('10%-Q: ',round(quantilesSenMayo[1],2)), 
            vjust = 2, 
            hjust = 1,
            size = 4) +
  geom_text(aes(quantilesSenMayo[5]), 
            y =Inf, 
            label = paste0('90%-Q: ',round(quantilesSenMayo[5], 2)),
            vjust = 2, 
            hjust = 0,
            size = 4)

dev.off()

png(height = 4,
    width = 6,
    units = 'in',
    res = 300,
    file =
      paste0(outputPath,
             "Results/CellFateAnalysis/PNGs/",
             "DistributionHDS_Histogram.png"))

quantilesHDS <- quantile(mergedDataFrame$HDS, probs = c(0.1,0.25, 0.5, 0.75,0.9))

ggplot(mergedDataFrame, aes(x = HDS)) + 
  geom_histogram(binwidth = 0.002) +
  geom_vline(aes(xintercept = quantilesHDS[1]), color = 'red') +
  geom_vline(aes(xintercept = quantilesHDS[5]), color = 'red') + 
  geom_text(aes(quantilesHDS[1]), 
            y = Inf, 
            label = paste0('10%-Q: ',round(quantilesHDS[1],2)), 
            vjust = 2, 
            hjust = 1,
            size = 4) +
  geom_text(aes(quantilesHDS[5]), 
            y =Inf, 
            label = paste0('90%-Q: ',round(quantilesHDS[5], 2)),
            vjust = 2, 
            hjust = 0,
            size = 4)

dev.off()


## categorize as senscence or proliferative cellfate
mergedDataFrame <- mergedDataFrame %>% mutate(CellFate = case_when(
  deltaSmP_centered >= quantiles[5] ~ 'S',
  deltaSmP_centered < quantiles[5] & deltaSmP_centered > quantiles[1] ~ 'undefined',
  deltaSmP_centered <= quantiles[1] ~ 'P'
))

# Catergorize by combination of HDS and SenMayo signal
mergedDataFrame <- mergedDataFrame %>% mutate(SenescenceDamageStatus = case_when(
  AUCell_SenMayo >= quantilesSenMayo[5] & HDS >= quantilesHDS[5] ~ 'HighHDS_HighSI',
  AUCell_SenMayo <= quantilesSenMayo[3] & HDS >= quantilesHDS[5] ~ 'HighHDS_LowSI',
  HDS <= quantilesHDS[3] ~ 'lowHDSroot'
))

mergedDataFrame$SenescenceDamageStatus[is.na(mergedDataFrame$SenescenceDamageStatus)] <- 'undefined'

if(identical(rownames(mergedDataFrame), 
             rownames(mergedObj@meta.data)) == TRUE){
  print('all good')
  
  mergedObj@meta.data$deltaSmP_centered <- mergedDataFrame$deltaSmP_centered
  mergedObj@meta.data$CellFate <- mergedDataFrame$CellFate
  mergedObj@meta.data$SenescenceDamageStatus <- mergedDataFrame$SenescenceDamageStatus
  
}





p1 <- (FeaturePlot(mergedObj,'HDS', min.cutoff = -0.15, max.cutoff = -0.025 )|FeaturePlot(mergedObj, 'deltaSmP_centered', min.cutoff = -0.05, max.cutoff = 0.05))/(DimPlot(mergedObj, group.by = c('CellType','condition')))
p2 <- (FeaturePlot(mergedObj,'HDS', min.cutoff = -0.15, max.cutoff = -0.025)|FeaturePlot(mergedObj, 'AUCell_SenMayo', min.cutoff = 0.02, max.cutoff = 0.07))/(DimPlot(mergedObj, group.by = c('CellType','condition')))
p3 <- (FeaturePlot(mergedObj,'HDS', min.cutoff = -0.15, max.cutoff = -0.025)|FeaturePlot(mergedObj, 'AUCell_ProlifIndex', max.cutoff = 0.04))/(DimPlot(mergedObj, group.by = c('CellType','condition')))
p4 <- DimPlot(mergedObj, group.by = c('CellFate','condition'), order = c('S','P', 'undefined'), alpha = 0.6)


png(height = 6,
    width = 8,
    units = 'in',
    res = 300,
    file = paste0(outputPath,
                  "Results/CellFateAnalysis/PNGs/",
                  "Umap_Merged_deltaSmpCentered.png"))
p1
dev.off()

png(height = 6,
    width = 8,
    units = 'in',
    res = 300,
    file = paste0(outputPath,
                  "Results/CellFateAnalysis/PNGs/",
                  "Umap_Merged_SenMayo.png"))
p2
dev.off()

png(height = 6,
    width = 8,
    units = 'in',
    res = 300,
    file = paste0(outputPath,
                  "Results/CellFateAnalysis/PNGs/",
                  "Umap_Merged_ProlifIndex.png"))
p3
dev.off()

png(height = 4,
    width = 8,
    units = 'in',
    res = 300,
    file = paste0(outputPath,
                  "Results/CellFateAnalysis/PNGs/",
                  "Umap_Merged_CellFate.png"))
p4
dev.off()
}

# 1. Distribution of deltaSmP across HDS is not random
# with increasing HDS, absolute delta SmP increases 

deltaSmP_HDS_bySample <-
  ggplot(data = mergedDataFrame,
         aes(x = HDS,
             y = deltaSmP_centered)) + 
  geom_point(size = 0.5)   +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  facet_wrap( ~ condition) +
  stat_cor(method = 'spearman',
           cor.coef.name = 'rho')

AbsDeltaSmP_HDS_bySample <-
  ggplot(data = mergedDataFrame,
         aes(x = HDS,
             y = abs(deltaSmP_centered))) + 
  geom_point(size = 0.5)   +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  facet_wrap( ~ condition) +
  stat_cor(method = 'spearman',
           cor.coef.name = 'rho')

pdf(height = 6, 
    width = 6,
    file = 
      paste0(outputPath,
             "Results/CellFateAnalysis/",
             "Spearman_deltaSmP_HDS_bySampleMergedXiao.pdf"))
print(deltaSmP_HDS_bySample)
dev.off()

png(height = 6, 
    width = 6,
    units = 'in',
    res = 300,
    file = 
      paste0(outputPath,
             "Results/CellFateAnalysis/PNGs/",
             "Spearman_deltaSmP_HDS_bySampleMergedXiao.png"))
print(deltaSmP_HDS_bySample)
dev.off()


pdf(height = 6, 
    width = 6,
    file = 
      paste0(outputPath,
             "Results/CellFateAnalysis/",
             "Spearman_AbsDeltaSmP_HDS_bySampleMergedXiao.pdf"))
print(AbsDeltaSmP_HDS_bySample)
dev.off()


png(height = 6, 
    width = 6,
    units = 'in',
    res = 300,
    file = 
      paste0(outputPath,
             "Results/CellFateAnalysis/PNGs/",
             "Spearman_AbsDeltaSmP_HDS_bySampleMergedXiao.png"))
print(AbsDeltaSmP_HDS_bySample)
dev.off()


# 2. Cell annotated as potentially ongoing senescent or proliferative 
# cell fate, the thresholds ware based on visual examination of 
# the distribution of centered deltaSmP

plotCategoryAssign <-
  ggplot(data = mergedDataFrame,
         aes(x = HDS,
             y = deltaSmP_centered,
             color = CellFate)) + 
  geom_point(size = 0.85)   +
  scale_color_manual(values = 
                       c('S' = 'orange', 
                         'undefined' = 'grey', 
                         'P' = 'skyblue')) +
  theme_bw() +
  facet_wrap( ~ condition) +
  theme(text = element_text(size = 20)) +
  geom_text(data = counts, 
            aes(x = label_x, y = label_y,
                label = paste(CellFate, 
                              count, 
                              sprintf("(%.1f%%)", 
                                      percentage), 
                              sep = " "),
                color = CellFate,
                fontface = "bold"),
            hjust = 0) + 
  theme(legend.position = "none")

pdf( height = 6, width = 6, 
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/centeredDeltaSmP_HDS_AnnotationCellFateMerged.pdf") )
print(plotCategoryAssign)
dev.off()

png( height = 6, width = 6, 
     units = 'in',
     res = 300,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/centeredDeltaSmP_HDS_AnnotationCellFateMerged.png") )
print(plotCategoryAssign)
dev.off()


# 3. Plot genes that sig. correlated with HDS in each group of cells 
# either proliferative or senscent group 
# genes are considered correlated if they have adj. p-value < 0.01 and 
# abs(ro) >= 0.2


genesPassCorr <- 
  ggplot(data = dfForPlot, aes(x = ro_Ses, y = ro_Prolif)) + 
  geom_point(size = 1.5, aes(color = GeneCategory)) + 
  theme_bw() +
  theme(text = element_text(size=20)) + 
  scale_color_manual( values = c('intersect' = 'pink', 
                                 'only proliferative' = 'skyblue',
                                 'only senescent' = 'orange',
                                 'not correlated' = 'grey')) +
  theme(
    legend.position = c(0.8, 0.2), 
    legend.background = element_rect(fill = "white", color = "black"),  
    legend.title = element_text(size = 12),  # Customize the legend title size
    legend.text = element_text(size = 10)# Customize the legend text size
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 5)  # Make the dots in the legend bigger
  )) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_fixed() + # makes plot canvas square
  labs(color = 'Gene significant in:',
       x = 'Senescent Cells Rho',
       y = 'Proliferative Cells Rho')



pdf( height = 6, width = 6, 
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/CorHDS_GenesCellFates_Xiao_rhoOver_0.2_adjPvalule0.01.pdf") )
print(genesPassCorr)
dev.off()

png( height = 6, width = 6, 
     units = 'in',
     res = 300,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/CorHDS_GenesCellFates_Xiao_rhoOver_0.2_adjPvalule0.01.png") )
print(genesPassCorr)
dev.off()







}
  
