library(Seurat)
## Script to merge two data sets from two conditions
## DATA SET 1
path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/'
path_annotation = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseStStAll.csv'


annotation.mouseStSt <- read.csv(
  file = path_annotation)


counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
        cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
        features = paste0(path_data,'features.tsv.gz', sep = '' ),
        feature.column = 1) 

sample_indices <-sample(c(1:length(annotation.mouseStSt$cell)), 25000, 
                        replace = F )

annotation.mouseStSt<-  annotation.mouseStSt[sample_indices, ]

lca.mouseStSt <- CreateSeuratObject(counts = counts, 
                                    project = "liver_cell_atlas_stst", 
                                    min.cells = 3)

lca.mouseStSt <- subset(lca.mouseStSt, cells = annotation.mouseStSt[,'cell'] )

##DATA SET 2

path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/'
path_annotation = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseNafldAll.csv'

annotation.mouseNafld <- read.csv(
  file = path_annotation)


counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
        cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
        features = paste0(path_data,'features.tsv.gz', sep = '' ),
        feature.column = 1) 

sample_indices <-sample(c(1:length(annotation.mouseNafld$cell)), 25000, 
                        replace = F )

annotation.mouseNafld<-  annotation.mouseNafld[sample_indices, ]

lca.mouseNafld <- CreateSeuratObject(counts = counts, 
                                    project = "liver_cell_atlas_nafld", 
                                    min.cells = 3)

lca.mouseNafld <- subset(lca.mouseNafld, cells = annotation.mouseNafld[,'cell'] )


lca.StStNafld.combined <- merge(lca.mouseStSt, y = lca.mouseNafld, add.cell.ids = c("StSt", "Nafld"), project = "LiverCellAtlas")


# SAVE 
write.csv(annotation.mouseStSt, file = 'annot_mouseStSt_ForMerging_randomsubset.csv')
write.csv(annotation.mouseNafld, file = 'annot_mouseNafld_ForMerging_randomsubset.csv')
saveRDS(lca.StStNafld.combined, file = 'countTable_mouseMergedConditions_randomsubset.rds')
