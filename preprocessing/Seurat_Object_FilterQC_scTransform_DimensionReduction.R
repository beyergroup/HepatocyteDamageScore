library(Seurat)
library(sctransform)

lca.mouseMerged <- readRDS(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/countTable_mouseMergedConditions_randomsubset.rds')


lca.mouseMerged[["percent.mt"]] <- PercentageFeatureSet(lca.mouseMerged, pattern = "^mt-")

lca.mouseMerged <- subset(lca.mouseMerged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5 & nCount_RNA < 15000)

lca.mouseMerged  <- SCTransform(lca.mouseMerged , verbose = FALSE)

lca.mouseMerged <- RunPCA(lca.mouseMerged, verbose = FALSE)
lca.mouseMerged  <- RunUMAP(lca.mouseMerged , dims = 1:17, verbose = FALSE)

lca.mouseMerged <- FindNeighbors(lca.mouseMerged, dims = 1:17)
lca.mouseMerged<- FindClusters(lca.mouseMerged, resolution = 0.5)



saveRDS(lca.mouseMerged, file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/countTable_mouseMergedConditions_randomsubset_sctransformed.rds')
