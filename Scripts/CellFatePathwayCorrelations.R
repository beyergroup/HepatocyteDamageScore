
# PATHWAY ACTIVITY PER BIN (first attempt)
{
  
  # Load necessary libraries
  library(msigdbr)
  
  ########### AUCELL
  
  
  # Load necessary libraries
  
  # Assume data is a Seurat object
  # Convert reactome dataframe to a list of gene sets
  gene_sets <- split(msigdb$gene_symbol, msigdb$gs_name)
  
  # Calculate AUCell scores
  identical(rownames(mergedDataFrame), rownames(mergedObj@meta.data))
  
  
  
  
  # AUCell scores of all pathways in gene_sets 
  # output list of 
  cellCat <- c('S', 'P', 'undefined')
  binsNames <- 
    unique(mergedObj@meta.data$HDS_equalBinsCellFate[mergedObj@meta.data$CellFate == 'S'])
  
  AUCellScoresPerCellFate <- lapply(seq(cellCat), function(ii){
    tempSubset <- subset(mergedObj, 
                         subset = CellFate == cellCat[ii] )
    print('done with subsetting')
    cell_rankings <- 
      AUCell_buildRankings(tempSubset@assays$SCT@data, 
                           plotStats = TRUE)
    remove(tempSubset)
    
    aucell_scores <- AUCell_calcAUC(gene_sets, 
                                    cell_rankings)
    
    return(getAUC(aucell_scores))
    print('done')
    
  }
  )
  
  mergedDataFrame$cellBarcode <- rownames(mergedDataFrame)
  # calculate the correlation of pathway activity with HDS once of
  # each cell fate group
  # function returns a list of three correlation test results 
  ScoreCorrHDSperCellFate <- lapply(seq(cellCat), function(ii){
    temp <- AUCellScoresPerCellFate[[ii]] 
    
    HDScells <- mergedDataFrame %>%
      filter(cellBarcode %in% colnames(temp))
    
    if(identical(colnames(temp), HDScells$cellBarcode)){
      print('yay')
    }
    HDScells <- HDScells$HDS 
    head(HDScells)
    
    tempRes <-
      apply(temp, 1, function(x){
        cor.test(x = x, y = HDScells, method = 'spearman')
      })
    
    class(tempRes)
    
    return(tempRes)
    
  }) 
  
  # reformat results
  PathwayCorrelations <- lapply(seq(cellCat), function(ii){
    pvalues <- c(ScoreCorrHDSperCellFate[[ii]][[1]]$p.value)
    rhovalues <- c(ScoreCorrHDSperCellFate[[ii]][[1]]$estimate)
    
    for (i in 2:1600){
      ptoadd <- ScoreCorrHDSperCellFate[[ii]][[i]]$p.value
      pvalues <- c(pvalues, ptoadd)
      
      rhotoadd <- ScoreCorrHDSperCellFate[[ii]][[i]]$estimate
      rhovalues <- c(rhovalues, rhotoadd)
      
    } 
    
    # add ajusted p-value 
    
    PathwayCorrelations <- data.frame(
      'pathway' = names(ScoreCorrHDSperCellFate[[ii]]),
      'pvalue' = pvalues,
      'padjust' = p.adjust(pvalues, method = "BH"),
      'rho' = rhovalues,
      'CellFate' = rep(cellCat[ii], length(pvalues))) 
    
    return(PathwayCorrelations)
  })
  
  
  # pathway correlations with HDS in all groups 
  PathwayCorrelations <- rbind.data.frame(PathwayCorrelations[[1]],
                                          PathwayCorrelations[[2]],
                                          PathwayCorrelations[[3]])
  
  PathwayCorrelations <- PathwayCorrelations %>% 
    mutate( sigcor = ifelse(padjust < 0.01 & abs(rho) >= 0.5, 
                            'correlated', 
                            'not correlated')) 
  
  selectedPathways <- PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated']
  # PathwaysOnlyCorS <-PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated' &
  #                                                PathwayCorrelations$CellFate == 'S' ]
  # 
  # PathwaysOnlyCorP <-PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated' &
  #                                                  PathwayCorrelations$CellFate == 'P' ]
  # 
  avgAUCellScorePerBin <- lapply(seq(cellCat), function(ii){
    
    tempData <- mergedDataFrame[mergedDataFrame$CellFate == cellCat[ii], c('cellBarcode', 'HDS_equalBinsCellFate')]
    bins <- unique(tempData$HDS_equalBinsCellFate)
    i <- ii
    avgPerBin <- lapply(seq(bins), function(jj){
      barcodes <- tempData$cellBarcode[tempData$HDS_equalBinsCellFate == bins[jj]]
      tempAUC <- AUCellScoresPerCellFate[[i]][, barcodes]
      tempData <- data.frame('avgScore' = rowSums(tempAUC),
                             'pathway' = rownames(tempAUC))
      tempData <- tempData %>%
        filter(pathway %in% selectedPathways)
      return(tempData)
    })
    
    
    
    return(avgPerBin)
    
  }
  )
  
  View(avgAUCellScorePerBin[[1]][[1]])
  
  # reformat results
  avgAUCellScorePerBin <- lapply(seq(cellCat), function(ii){  
    
    temp <- reduce(avgAUCellScorePerBin[[ii]],
                   left_join, by = 'pathway')
    
    rownames(temp) <- temp$pathway
    
    temp <- temp[,-2]
    
    nbins <- ncol(temp)
    binNames <- unique(mergedDataFrame[mergedDataFrame$CellFate == cellCat[ii],
                                       'HDS_equalBinsCellFate'])
    j <- 1
    binColumn <- rep(binNames[j], length(rownames(temp)))
    print(binColumn)
    j <- j + 1
    for(j in 2:nbins){
      binColumn <- 
        c(binColumn,rep(binNames[j], 
                        length(rownames(temp))
        )
        )
    }
    print(binColumn)                  
    temp <- data.frame('pathway' = rep(rownames(temp), 
                                       nbins),
                       'activity_score' = unlist(temp),
                       'bin' = binColumn,
                       'CellFate' = rep(cellCat[ii], nrow(temp))
    )
    
    return(temp)
    
  })
  
  allData <- rbind(avgAUCellScorePerBin[[1]],
                   avgAUCellScorePerBin[[2]],
                   avgAUCellScorePerBin[[3]])
  
  ## HERE 
  
  # Create the ggplot
  
  plotTest <- ggplot(allData, aes(x = bin, y = pathway, color = activity_score)) +
    geom_point(size = 4) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(
      title = "Pathway Activity Over Time",
      x = "Timepoint",
      y = "Pathway Activity",
      color = "Activity Level"
    ) +
    theme_minimal() + theme(axis.text.x = 
                              element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~CellFate)
  
  pdf(height = 20, width = 18, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/TEST_1.pdf"))
  plotTest
  dev.off()
  
  
  # Add interferon pathway scores to Seurat object metadata
  for (pathway in ifn_pathways) {
    SenesceSeuratBin1@meta.data[[pathway]] <- aucell_scores@assays@data[[pathway]]
  }
  
  # UMAP plot with condition, cell type, and interferon pathways
  DimPlot(SenesceSeuratBin1, reduction = "umap", group.by = c("condition"))
  
  # Feature plots for the interferon pathways
  FeaturePlot(SenesceSeuratBin1, features = ifn_pathways, reduction = "umap", ncol = 2)
  
  
}

# Bin cells in equally sized bins by HDS REGARDLESS of group 
# this could lead to come bins having less member or more members of a group 
# a. we will do a pathway activity analysis on each bin of each group 

{
  
  #cell type categories
  # cellCat <- c('S', 'P', 'undefined')
  # 
  # binsN <- 10
  #   
  # quantilesBins <- quantile(mergedDataFrame$HDS,
  #                           probs = seq(0, 1, length.out = binsN +1))
  # 
  # labelsBins <- sapply(seq(1:(binsN)), function(ii){
  #     
  #     label <- paste0('(',round(quantilesBins[ii], 3),
  #                     ',',round(quantilesBins[ii+1],3),']')  
  #     return(label)
  #     })
  #   
  #   
  #   mergedDataFrame$HDS_equalBins <- cut(mergedDataFrame$HDS,
  #                                     breaks = quantilesBins,
  #                                     include.lowest = TRUE,
  #                                     labels = labelsBins)
  #   
  # 
  # colorsToUse <- c('orange','skyblue','grey')
  # 
  # lapply(seq(cellCat), function(ii){
  #   
  #   ggplot(data = mergedDataFrame[mergedDataFrame$CellFate == cellCat[ii],], 
  #          aes(x=HDS_equalBins, y = deltaSmP_centered) ) +
  #     geom_violin(trim = FALSE) +
  #     theme(axis.text.x = element_text(angle = 90)) +
  #     xlab("HDS (bins) ") +
  #     stat_summary(fun.y=median, geom="point", size=2, color=colorsToUse[ii])
  #   
  # })
  # 
  # table(mergedDataFrame$HDS_equalBins,
  #       mergedDataFrame$CellFate)
  # 
  # mergedDataFrame$HDS_equalBins <- factor(mergedDataFrame$HDS_equalBins)
  # 
  # ggplot(data = mergedDataFrame, 
  #        aes(x=HDS_equalBins, y = deltaSmP_centered) ) +
  #   geom_violin(trim = FALSE) +
  #   theme(axis.text.x = element_text(angle = 90)) +
  #   xlab("HDS (bins) ") +
  #   stat_summary(fun.data=mean_sdl, mult=1, 
  #                geom="pointrange", color="red")
  # 
  # # STORE bin information in Seurat object metadata
  # 
  # 
  # mergedDataFrame$cellBarcode <- rownames(mergedDataFrame)
  # mergedObj@meta.data$cellBarcode <- rownames(mergedObj@meta.data$cellBarcode)
  # 
  # mergedDataFrame <- mergedDataFrame %>%
  #   arrange(match(mergedDataFrame$cellBarcode, mergedObj@meta.data$cellBarcode))
  # 
  # if(identical(mergedDataFrame$cellBarcode, mergedObj@meta.data$cellBarcode) == TRUE){
  #  print('all good')
  #    mergedObj@meta.data$deltaSmP_centered <- mergedDataFrame[,'deltaSmP_centered']
  #   mergedObj@meta.data$CellFate <- mergedDataFrame[,'CellFate']
  #   mergedObj@meta.data$HDS_equalBins <- mergedDataFrame[,'HDS_equalBins']
  #   print('done')
  # }
  # 
}


# PATHWAY ACTIVITY PER BIN (second attempt)
{
  
  
  
  ########### AUCELL
  
  
  # Load necessary libraries
  
  # Assume data is a Seurat object
  # Convert reactome dataframe to a list of gene sets
  gene_sets <- split(msigdb$gene_symbol, msigdb$gs_name)
  
  
  
  # AUCell scores of all pathways in gene_sets 
  binsNames <- unique(mergedObj@meta.data$HDS_equalBins)
  
  ## !! BEWARE REALLY TIME INTENSIVE! 15 min?
  
  cell_rankings <-
    AUCell_buildRankings(mergedObjExpr, 
                         plotStats = TRUE)
  
  aucell_scores <- AUCell_calcAUC(gene_sets,
                                  cell_rankings)
  
  AUCellScores <- getAUC(aucell_scores)
  #remove(aucell_scores)
  
  
  
  # calculate the correlation of pathway activity with HDS once of
  # each cell fate group
  # function returns a list of three correlation test results 
  ScoreCorrHDSperCellFate <- lapply(seq(cellCat), function(ii){
    
    cellIDs <- mergedDataFrame$cellBarcode[mergedDataFrame$CellFate == cellCat[ii]]
    HDScells <- mergedDataFrame$HDS[mergedDataFrame$CellFate == cellCat[ii]]
    
    print(head(HDScells))
    
    tempRes <-
      apply(AUCellScores[,cellIDs], 1, function(x){
        cor.test(x = x, y = HDScells, method = 'spearman')
      })
    
    class(tempRes)
    
    return(tempRes)
    
  }) 
  
  # reformat results
  PathwayCorrelations <- lapply(seq(cellCat), function(ii){
    pvalues <- c(ScoreCorrHDSperCellFate[[ii]][[1]]$p.value)
    rhovalues <- c(ScoreCorrHDSperCellFate[[ii]][[1]]$estimate)
    
    for (i in 2:1600){
      ptoadd <- ScoreCorrHDSperCellFate[[ii]][[i]]$p.value
      pvalues <- c(pvalues, ptoadd)
      
      rhotoadd <- ScoreCorrHDSperCellFate[[ii]][[i]]$estimate
      rhovalues <- c(rhovalues, rhotoadd)
      
    } 
    
    # add ajusted p-value 
    
    PathwayCorrelations <- data.frame(
      'pathway' = names(ScoreCorrHDSperCellFate[[ii]]),
      'pvalue' = pvalues,
      'padjust' = p.adjust(pvalues, method = "BH"),
      'rho' = rhovalues,
      'CellFate' = rep(cellCat[ii], length(pvalues))) 
    
    return(PathwayCorrelations)
  })
  
  
  # pathway correlations with HDS in all groups 
  PathwayCorrelations <- rbind.data.frame(PathwayCorrelations[[1]],
                                          PathwayCorrelations[[2]],
                                          PathwayCorrelations[[3]])
  
  PathwayCorrelations <- PathwayCorrelations %>% 
    mutate( sigcor = ifelse(padjust < 0.01 & abs(rho) >= 0.5, 
                            'correlated', 
                            'not correlated')) 
  
  # # alternative
  # PathwayCorrelations <- PathwayCorrelations %>% 
  #   mutate( sigcorGroup = ifelse(padjust < 0.01 & abs(rho) >= 0.4 & CellFate == 'S', 
  #                           's_correlated', 
  #                           's_not_correlated')) 
  # 
  
  selectedPathways <- PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated']
  
  #PathwaysOnlyCorS <-PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated' &
  #                                                   PathwayCorrelations$CellFate == 'S']
  
  # PathwaysOnlyCorP <-PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated' &
  #                                                 PathwayCorrelations$CellFate == 'P' ]
  
  #PathwaysOnlyCorUndefined <-PathwayCorrelations$pathway[PathwayCorrelations$sigcor == 'correlated' &
  ##                                                PathwayCorrelations$CellFate == 'undefined' ]
  
  avgAUCellScorePerBin <- lapply(seq(cellCat), function(ii){
    
    tempData <- mergedDataFrame[mergedDataFrame$CellFate == cellCat[ii], c('cellBarcode', 'HDS_equalBins')]
    bins <- unique(tempData$HDS_equalBins)
    i <- ii
    avgPerBin <- lapply(seq(bins), function(jj){
      barcodes <- tempData$cellBarcode[tempData$HDS_equalBins == bins[jj]]
      tempAUC <- AUCellScores[, barcodes]
      tempData <- data.frame('avgScore' = rowMeans(tempAUC),
                             'pathway' = rownames(tempAUC))
      tempData <- tempData %>%
        filter(pathway %in% selectedPathways)
      return(tempData)
    })
    
    
    
    return(avgPerBin)
    
  }
  )
  
  
  # reformat results
  avgAUCellScorePerBin <- lapply(seq(cellCat), function(ii){  
    
    temp <- reduce(avgAUCellScorePerBin[[ii]],
                   left_join, by = 'pathway')
    
    rownames(temp) <- temp$pathway
    
    temp <- temp[,-2]
    
    nbins <- ncol(temp)
    binNames <- unique(mergedDataFrame[mergedDataFrame$CellFate == cellCat[ii],
                                       'HDS_equalBins'])
    j <- 1
    binColumn <- rep(binNames[j], length(rownames(temp)))
    print(binColumn)
    j <- j + 1
    for(j in 2:nbins){
      binColumn <- 
        c(binColumn,rep(binNames[j], 
                        length(rownames(temp))
        )
        )
    }
    print(binColumn)                  
    temp <- data.frame('pathway' = rep(rownames(temp), 
                                       nbins),
                       'activity_score' = unlist(temp),
                       'bin' = binColumn,
                       'CellFate' = rep(cellCat[ii], nrow(temp))
    )
    
    return(temp)
    
  })
  
  allData <- rbind(avgAUCellScorePerBin[[1]],
                   avgAUCellScorePerBin[[2]],
                   avgAUCellScorePerBin[[3]])
  
  allData$ScaledCentered_activity_score <- scale(allData$activity_score, 
                                                 center = TRUE, scale = TRUE)
  
  ## HERE 
  
  # Create the ggplot
  # 
  colorScaleBreaks <- seq(floor(min(allData$ScaledCentered_activity_score)),
                          ceiling(max(allData$ScaledCentered_activity_score)),
                          by = 0.2)
  
  plotTest <- ggplot(allData, aes(x = bin, y = pathway, 
                                  color = ScaledCentered_activity_score)) +
    geom_point(size = 4, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D",
                           breaks = colorScaleBreaks,
                           #breaks = seq(-1.8, 5, by = 0.4) 
    ) +
    labs(
      title = "Pathway activity sig. correlated \neither in cell fate group with\n |rho| >= 0.5 & p < 0.01",
      x = "HDS bins",
      y = "Pathway",
      color = "Mean Activity Score per Bin\n(centered and scaled)"
    ) +
    facet_grid(~ CellFate) +
    theme_classic() + theme(axis.text.x = 
                              element_text(angle = 90, vjust = 0.5, hjust=1),
                            legend.key.height = unit(1.5, 'cm'),
                            legend.key.width = unit(0.5, 'cm'))
  
  pdf(height = 8, width = 12, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/ScaledCentered_All_CellBinnedByHDS_PathwayCorrelationWithHDS_Rho0.5.pdf"))
  plotTest
  dev.off()
  
  png(height = 8, 
      width = 12,
      units = 'in',
      res = 300, 
      file = paste0(outputPath, 
                    'Results/CellFateAnalysis/PNGs/',
                    'scaledCentered_All_CellBinnedByHDS_PathwayCorrelationWithHDS_Rho0.5.png'))
  plotTest
  dev.off()
  
  
  
  temp <- mergedDataFrame[mergedDataFrame$CellFate == 'S',]
  temp <- mergedDataFrame
  temp$ReactomeCellCycle <- AUCellScores['REACTOME_CELL_CYCLE', temp$cellBarcode]
  temp$ReactomeCilumAssembly <- AUCellScores['REACTOME_CILIUM_ASSEMBLY', temp$cellBarcode]
  temp$ReactomeCHREBP <- AUCellScores['REACTOME_CHREBP_ACTIVATES_METABOLIC_GENE_EXPRESSION', temp$cellBarcode]
  temp$ReactomeRRNAProcessing <- AUCellScores['REACTOME_RRNA_PROCESSING', temp$cellBarcode]
  
  
  plotPathway <- ggplot(temp, aes(x = HDS, 
                                  y = ReactomeCellCycle, 
                                  color = ReactomeCellCycle)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D"
                           #breaks = seq(-1.8, 5, by = 0.4)
    ) +
    labs(
      title = "Pathway activity correlated in each cell fate (rho > 0.5)",
      x = "HDS bins",
      y = "Pathway",
      color = "AUCellScore"
    ) +
    theme_classic() + theme(axis.text.x = 
                              element_text(angle = 90, vjust = 0.5, hjust=1),
                            legend.key.height = unit(1.5, 'cm'),
                            legend.key.width = unit(0.5, 'cm'))
  
  plotPathway
  
  
  pathways <- unique(reactome_unique$gs_name)
  senescence_indices <- grep("SENESCENCE", 
                             pathways, 
                             ignore.case = TRUE)
  
  senscence_pathways <- pathways[senescence_indices]
  
  # Add interferon pathway scores to Seurat object metadata
  for (pathway in senscence_pathways) {
    
    temp <- mergedDataFrame
    temp$tempslot <- AUCellScores[pathway,
                                  mergedDataFrame$cellBarcode]
    
    
    gg <- ggplot(temp, aes(x = HDS_equalBins,
                           y = tempslot,
                           color = tempslot)) +
      geom_point(size = 2, alpha = 1) +  # Adjust size as needed
      scale_colour_viridis_b(option = "D") +
      labs(title = paste0(pathway, ' across HDS bins'),
           x = "HDS bins",
           y = pathway,
           color = "AUCellScore") +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.key.height = unit(1.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'))
    
    print(gg)
    
  }
  
  pathways <- unique(reactome_unique$gs_name)
  senescence_indices <- grep("SENESCENCE", 
                             pathways, 
                             ignore.case = TRUE)
  
  senscence_pathways <- c(pathways[senescence_indices], 'REACTOME_DNA_REPAIR',
                          'REACTOME_CELL_CYCLE_MITOTIC')
  
  # Add interferon pathway scores to Seurat object metadata
  for (pathway in senscence_pathways) {
    
    temp <- mergedDataFrame
    temp$tempslot <- AUCellScores[pathway,
                                  mergedDataFrame$cellBarcode]
    
    
    gg <- ggplot(temp, aes(x = HDS_equalBins,
                           y = tempslot,
                           color = tempslot)) +
      geom_boxplot(size = 1) +  # Adjust size as needed
      scale_colour_viridis_b(option = "D") +
      labs(title = paste0(pathway, ' across HDS bins'),
           x = "HDS bins",
           y = pathway,
           color = "AUCellScore") +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.key.height = unit(1.5, 'cm'),
            legend.key.width = unit(0.5, 'cm')) +
      facet_grid(~senescent)
    
    print(gg)
    
  }
  
  temp <- mergedDataFrame
  temp$tempslot <- AUCellScores['REACTOME_SUPPRESSION_OF_APOPTOSIS',
                                mergedDataFrame$cellBarcode]
  
  
  gg <- ggplot(temp[temp$CellFate == 'P',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'Apoptosis Supression in Proliferative Cells',
         x = "HDS bins",
         y = 'REACTOME_SUPPRESSION_OF_APOPTOSIS',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(gg)
  
  pp <- ggplot(temp[temp$CellFate == 'S',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'Apoptosis Supression in Senescent Cells',
         x = "HDS bins",
         y = 'REACTOME_SUPPRESSION_OF_APOPTOSIS',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(pp)
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Sen_CellBinnedByHDS_ApoptosisSuppression.pdf"))
  pp
  dev.off()
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Prolif_CellBinnedByHDS_ApoptosisSuppression.pdf"))
  gg
  dev.off()
  
  
  #### SAME BUT DIFFERENT PATHWAY
  
  temp <- mergedDataFrame
  temp$tempslot <- AUCellScores['REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS',
                                mergedDataFrame$cellBarcode]
  
  
  gg <- ggplot(temp[temp$CellFate == 'P',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'RIPK1_MEDIATED_REGULATED_NECROSIS in Proliferative Cells',
         x = "HDS bins",
         y = 'REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(gg)
  
  pp <- ggplot(temp[temp$CellFate == 'S',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'RIPK1_MEDIATED_REGULATED_NECROSIS in Senescent Cells',
         x = "HDS bins",
         y = 'REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(pp)
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Prolif_CellBinnedByHDS_RIPK1regulatednecrosis.pdf"))
  gg
  dev.off()
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Ses_CellBinnedByHDS_RIPK1regulatednecrosis.pdf"))
  pp
  dev.off()
  
  #### SAME BUT DIFFERENT PATHWAY
  
  temp <- mergedDataFrame
  temp$tempslot <- AUCellScores['REACTOME_ONCOGENE_INDUCED_SENESCENCE',
                                mergedDataFrame$cellBarcode]
  
  
  gg <- ggplot(temp[temp$CellFate == 'P',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'REACTOME_ONCOGENE_INDUCED_SENESCENCE in Proliferative Cells',
         x = "HDS bins",
         y = 'REACTOME_ONCOGENE_INDUCED_SENESCENCE',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(gg)
  
  pp <- ggplot(temp[temp$CellFate == 'S',], aes(x = HDS_equalBins,
                                                y = tempslot,
                                                color = tempslot)) +
    geom_point(size = 2, alpha = 1) +  # Adjust size as needed
    scale_colour_viridis_b(option = "D") +
    labs(title = 'REACTOME_ONCOGENE_INDUCED_SENESCENCE in Senescent Cells',
         x = "HDS bins",
         y = 'REACTOME_ONCOGENE_INDUCED_SENESCENCE',
         color = "AUCellScore") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.key.height = unit(1.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'))
  
  print(pp)
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Prolif_CellBinnedByHDS_REACTOME_ONCOGENE_INDUCED_SENESCENCE.pdf"))
  gg
  dev.off()
  
  pdf(height = 6, width = 8, 
      file = paste0(outputPath, 
                    "Results/CellFateAnalysis/Ses_CellBinnedByHDS_REACTOME_ONCOGENE_INDUCED_SENESCENCE.pdf"))
  pp
  dev.off()
  
  
  
  # UMAP plot with condition, cell type, and interferon pathways
  DimPlot(SenesceSeuratBin1, reduction = "umap", group.by = c("condition"))
  
  # Feature plots for the interferon pathways
  FeaturePlot(SenesceSeuratBin1, 
              features = ifn_pathways, 
              reduction = "umap", 
              ncol = 2)
  
}

# try a different labeling procedure:
# cells are considered senescent 
{
  
  # 
  mergedDataFrame$senescent <- rep('no', length(mergedDataFrame$orig.ident))
  mergedDataFrame$senescent[mergedDataFrame$AUCell_SenMayo >= 0.065] <- 'yes'
  
  ## 
  
  table(mergedDataFrame$senescent)
  
  ggplot(data = mergedDataFrame, 
         aes(x=HDS_equalBins, y = AUCell_SenMayo) ) +
    geom_violin(trim = FALSE) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("HDS (bins) ") +
    stat_summary(fun.data=mean_sdl, mult = 1, 
                 geom="pointrange", color="red")
  
  ggplot(data = mergedDataFrame, 
         aes(x=HDS_equalBins, y = AUCell_SenMayo, color = senescent) ) + geom_point()
  
  
  
}

