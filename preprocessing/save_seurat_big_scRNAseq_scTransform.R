library(Seurat)
library(sctransform)

#lca.mouseStSt.Subset <- readRDS(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/countTable_mouseStSt_randomsubset.rds')

lca.mouseNafld.Subset <- readRDS(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/countTable_mouseNafld_randomsubset.rds')

#annotation.mouseStSt.Subset <- 
 # read.csv(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/annot_mouseStStAll_randomsubset.csv')

annotation.mouseNafld.Subset <- 
  read.csv(file = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/annot_mouseNafldAll_randomsubset.csv')

#Idents(lca.mouseStSt.Subset) <- annotation.mouseStSt.Subset$annot
Idents(lca.mouseNafld.Subset) <- annotation.mouseNafld.Subset$annot

#lca.mouseStSt.Subset[["percent.mt"]] <- PercentageFeatureSet(lca.mouseStSt.Subset, 
 #                                                           pattern = "^mt-")

lca.mouseNafld.Subset[["percent.mt"]] <- PercentageFeatureSet(lca.mouseNafld.Subset, 
                                                             pattern = "^mt-")

#lca.mouseStSt.Subset <- subset(lca.mouseStSt.Subset, 
                              # subset = nFeature_RNA > 200 &
                               #  nFeature_RNA < 7000 &
                                # percent.mt < 10 & 
                                 #nCount_RNA < 50000)

lca.mouseNafld.Subset <- subset(lca.mouseNafld.Subset, 
                               subset = nFeature_RNA > 200 &
                                 nFeature_RNA < 3000 &
                                 percent.mt < 7.5 & 
                                 nCount_RNA < 15000)
#vars.to.regress = "percent.mt"was left out to see

#lca.mouseStSt.Subset <- SCTransform(lca.mouseStSt.Subset, 
#                                    verbose = FALSE)

lca.mouseNafld.Subset <- SCTransform(lca.mouseNafld.Subset, 
                                    verbose = FALSE)

#lca.mouseStSt.Subset <- RunPCA(lca.mouseStSt.Subset, verbose = FALSE)
lca.mouseNafld.Subset <- RunPCA(lca.mouseNafld.Subset, verbose = FALSE)
# dims = 17 because there are 17 cell types 
#lca.mouseStSt.Subset  <- RunUMAP(lca.mouseStSt.Subset , dims = 1:17, verbose = FALSE)

# dims = 16 because there are 16 cell types in atlas
lca.mouseNafld.Subset  <- RunUMAP(lca.mouseNafld.Subset , dims = 1:16, verbose = FALSE)


#saveRDS(lca.mouseStSt.Subset, file = 'countTable_mouseStSt_randomsubset_sctransformed.rds')
saveRDS(lca.mouseNafld.Subset , file = 'countTable_mouseNafld_randomsubset_sctransformed.rds')


