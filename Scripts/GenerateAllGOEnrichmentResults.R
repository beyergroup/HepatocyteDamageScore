# Functional characterization of HDS markers 
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(clusterProfiler)


HDAG <- read.csv(file = "hepatocyte-damage-score/Data/Output/HDAG.csv")
HDAG <- HDAG[(1:42),]
HDAG_neg42 <- HDAG[HDAG$direction_foldchange < 0,]
HDAG_pos42 <- HDAG[HDAG$direction_foldchange >0,]

egoALL <- enrichGO(gene = HDAG_neg42$gene_symbol,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.01)

write.csv(egoALL@result[], file = "hepatocyte-damage-score/Data/Output/Results/Negative42HDSMarkers_ClusterProfilerEnrichGOResults.csv")

egoALLpositive <- enrichGO(gene = HDAG_pos42$gene_symbol,
            OrgDb         = org.Mm.eg.db,
            keyType       = 'SYMBOL',
            ont           = "ALL",
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.01,
            qvalueCutoff  = 0.01)


write.csv(egoALL[egoALL@result$p.adjust < 0.01,], 
          file = "hepatocyte-damage-score/Data/Output/Results/Negative42HDSMarkers_ClusterProfilerEnrichGOResults.csv")

write.csv(egoALLpositive[egoALLpositive@result$p.adjust < 0.01,], 
          file = "hepatocyte-damage-score/Data/Output/Results/Positive42HDSMarkers_ClusterProfilerEnrichGOResults.csv")


#### GO enrichment of genes sig. up- or down-regulated in senescent damaged 
# hepatocytes compared to not-senescent-like damaged hepatocytes

# Upregulated Genes 
SenPosMarkers <- read.csv(file =
                            paste0("hepatocyte-damage-score/Data/Output/",
                                   'PositiveMarkers_DSvsDNS_NewThreshold.csv'))
SenPosMarkers <- SenPosMarkers$Gene.Name

egoSenPos <- enrichGO(gene = SenPosMarkers,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.01)

write.csv(egoSenPos[egoSenPos@result$p.adjust < 0.01,], 
          file = paste0("hepatocyte-damage-score/Data/Output/Results/",
                        "UpregulatedSenGenes_ClusterProfilerEnrichGOResults.csv"))


SenNegMarkers <- read.csv(file =
                            paste0("hepatocyte-damage-score/Data/Output/",
                                   'NegativeMarkers_DSvsDNS_NewThreshold.csv'))

SenNegMarkers <- SenNegMarkers$Gene.Name

egoSenNeg <- enrichGO(gene = SenNegMarkers,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.01)

write.csv(egoSenNeg[egoSenNeg@result$p.adjust < 0.01,], 
          file = paste0("hepatocyte-damage-score/Data/Output/Results/",
                        "DownregulatedSenGenes_ClusterProfilerEnrichGOResults.csv"))
