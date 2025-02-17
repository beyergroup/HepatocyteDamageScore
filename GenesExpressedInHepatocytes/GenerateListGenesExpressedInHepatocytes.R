library(Seurat)
library(sctransform)
library(ggplot2)
library(tidyverse)


## Load merged, filtered, transformed data (output from MergeSubsetToFindGenes.R)

outputPath <- 'hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/'

merged.mouse.lca.hep.nucSeq <- 
  readRDS(paste0(outputPath,
                 'mouseMergedHepatocytesnucSeqQCfitleredSCT.rds'))


# Detection Rate in hepatocytes (LFC between single nucleus and single cell)
foldchange <- FoldChange(merged.mouse.lca.hep.nucSeq,
                         ident.1 = 'nucSeq', 
                         group.by = 'typeSample',
                         assay = "SCT")

detection_rate_df <- data.frame('Genes' = row.names(foldchange), 
                                'detection_rate' = foldchange[,2],
                                'log2_dec_rate' = log2(foldchange[,2]+ 0.0001))

ggplot(data= detection_rate_df, aes(x=detection_rate)) + geom_histogram(binwidth = 0.005)
ggplot(data = detection_rate_df, aes(x = log2_dec_rate)) + geom_density()

summaryStats <- summary(detection_rate_df$log2_dec_rate, na.rm = TRUE)
quantilesLog2DetRate <- quantile(detection_rate_df$log2_dec_rate)

# Plotting theshold 
plotDetectionRate <- ggplot(data= detection_rate_df, aes(x=log2_dec_rate)) +
  geom_histogram(binwidth = 0.05) + 
  geom_vline(xintercept = c(quantilesLog2DetRate[2], 
                            quantilesLog2DetRate[3],
                            quantilesLog2DetRate[4]),
             colour = 'red') + 
  xlab('log2(detection rate of gene + 0.0001)') + 
  ylab('number of genes') + 
  geom_text(aes(x=quantilesLog2DetRate[2], 
                label="1. Quartile", 
                y=1100), 
            colour="red", 
            angle=90, 
            vjust = -1) +
  geom_text(aes(x=quantilesLog2DetRate[3], 
                label="2. Quartile", 
                y=1100), 
            colour="red", 
            angle=90, 
            vjust = -1) +
  geom_text(aes(x=quantilesLog2DetRate[4], 
                label="3. Quartile", 
                y=1100), 
            colour="red", 
            angle=90,
            vjust = -1)

# Save Output

# save gene list 

detection_rate_df <- 
  detection_rate_df[detection_rate_df$log2_dec_rate >  -5, ]

detection_rate_df <- 
  detection_rate_df[order(detection_rate_df$log2_dec_rate, decreasing = TRUE),]

rownames(detection_rate_df) <- 1:length(detection_rate_df$Genes)



# Save gene list

write.csv(detection_rate_df, 
          file = paste0(outputPath, 'GenesExpressedInHepatocytes.csv'))

# Save plots

ggsave(
  'plotDetectionRate.png',
  plot = plotDetectionRate,
  path = outputPath
)

# Save session info
write.table(toLatex(sessionInfo()), file = 'sessionInfo.txt', sep = '/')






