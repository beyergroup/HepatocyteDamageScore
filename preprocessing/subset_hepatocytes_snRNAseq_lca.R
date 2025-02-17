library(Seurat)

### Script for the cluster
## StSt
path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseStSt/countTable_mouseStSt/'
path_annotation = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseStStAll.csv'


annotation.mouseStSt <- read.csv(
  file = path_annotation)


counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
                  cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
                  features = paste0(path_data,'features.tsv.gz', sep = '' ),
                  feature.column = 1)

lca.mouse.hep.nucSEq.stst <- CreateSeuratObject(counts = counts, 
                                   project = "liver_cell_atlas_stst",
                                   min.cells = 3)

cells_to_keep <- annotation.mouseStSt[annotation.mouseStSt$annot == 'Hepatocytes', 'cell']


#subset: keep only celltype == 'Hepatocytes' and typeSample == 'nucSeq'
lca.mouse.hep.nucSEq.stst <- subset(lca.mouse.hep.nucSEq.stst, cells = cells_to_keep )
saveRDS(lca.mouse.hep.nucSEq.stst, file = 'countTable_mouseStSt_Hepatocytes_nucSeq.rds')

## Nafld 

path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/'
path_annotation = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseNafldAll.csv'

annotation.mouseNafld <- read.csv(
  file = path_annotation)

counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
                  cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
                  features = paste0(path_data,'features.tsv.gz', sep = '' ),
                  feature.column = 1)

lca.mouse.hep.nucSEq.Nafld <- CreateSeuratObject(counts = counts, 
                                                project = "liver_cell_atlas_Nafld",
                                                min.cells = 3)

cells_to_keep <- annotation.mouseNafld[annotation.mouseNafld$annot == 'Hepatocytes', 'cell']


#subset: keep only celltype == 'Hepatocytes' and typeSample == 'nucSeq'
lca.mouse.hep.nucSEq.Nafld <- subset(lca.mouse.hep.nucSEq.Nafld, cells = cells_to_keep)
saveRDS(lca.mouse.hep.nucSEq.Nafld, file = 'countTable_mouseNafld_Hepatocytes_nucSeq.rds')
