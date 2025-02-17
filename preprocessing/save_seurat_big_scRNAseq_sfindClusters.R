
library(Seurat)
library(sctransform)

lca.mouseStSt.Subset <- readRDS(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/countTable_mouseStSt_randomsubset_sctransformed.rds')

lca.mouseNafld.Subset <- readRDS(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/countTable_mouseNafld_randomsubset_sctransformed.rds')


lca.mouseStSt.Subset <- FindNeighbors(lca.mouseStSt.Subset, dims = 1:17)
lca.mouseStSt.Subset <- FindClusters(lca.mouseStSt.Subset, resolution = 0.4)

lca.mouseNafld.Subset <- FindNeighbors(lca.mouseNafld.Subset, dims = 1:16)
lca.mouseNafld.Subset <- FindClusters(lca.mouseNafld.Subset, resolution = 0.5)

saveRDS(lca.mouseStSt.Subset, file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/countTable_mouseStSt_randomsubset_sctransformed.rds')
saveRDS(lca.mouseNafld.Subset , file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/countTable_mouseNafld_randomsubset_sctransformed.rds')

