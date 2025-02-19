library(Seurat)

# script to create expression matrix of only nuqSeq counts annotaed as 
# Hepatocytes from two cohorts of the Liver Cell Atlas murine data set


merged.mouse.lca.hep.nucSeq <- 
  readRDS('hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/countTable_mouseMerged_Hepatocytes_nucSeq.rds')

merged.mouse.lca.hep.nucSeq[["percent.mt"]] <- 
  PercentageFeatureSet(merged.mouse.lca.hep.nucSeq,
                       pattern = "^mt-")

VlnPlot(merged.mouse.lca.hep.nucSeq, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
        ncol = 3, pt.size = F) 


VlnPlot(merged.mouse.lca.hep.nucSeq,
        features = c("percent.mt" ), pt.size = F, group.by = 'typeSample' ) 

 FeatureScatter(merged.mouse.lca.hep.nucSeq,
                feature1 = "nCount_RNA", 
                feature2 = "percent.mt")

merged.mouse.lca.hep.nucSeq <- subset(merged.mouse.lca.hep.nucSeq, 
                                      subset = nFeature_RNA > 200 &
                                        nFeature_RNA < 6000 &
                                        percent.mt <= 5 & 
                                        nCount_RNA < 30000 &
                                      typeSample == 'nucSeq')

unique(merged.mouse.lca.hep.nucSeq@meta.data$celltype)
# [1] "Hepatocytes"
# this is what we want!

# extract counts to a matrix
exprMatrices <- GetAssayData(object = merged.mouse.lca.hep.nucSeq , 
                             slot = 'counts')

temp <- merged.mouse.lca.hep.nucSeq@meta.data

annotations <- data.frame('Sample' = temp$sample , 
                          'cell_id' = rownames(temp),
                          'condition' = gsub(pattern = "_.*", 
                                             replacement = "", 
                                             rownames(temp)) )

exprMatrices <- as.matrix(exprMatrices)

saveRDS(exprMatrices, 
        'hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/ExprMatrix_QCfiltered_mouseMerged_Hepatocytes_nucSeq.rds')
write.csv(annotations, 
          'hepatocyte-damage-score/Data/Input/scRNAseq/LiverCellAtlas/ExprMatrix_QCfiltered_mouseMerged_Hepatocytes_nucSeq_annotations.csv')
