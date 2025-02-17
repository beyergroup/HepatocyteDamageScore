### Filter Universal Damage Signature Gene List numerically
# script to run on cluster 
library(Seurat)
path_data =  '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/rawData_mouseNafld/countTable_mouseNafld/'
path_annotation = '/data/public/pungerav/liver_disease_score/data_disease_score/scRNAseq/Liver_Cell_Atlas/annot_mouseNafldAll.csv'

#annotation.mouseStSt <- read.csv(
#  file = 'scRNAseq/Liver Cell Atlas Data /annot_mouseStStAll.csv')

annotation.mouseStSt <- read.csv(
  file = path_annotation)


counts <- ReadMtx(mtx =  paste0(path_data,'matrix.mtx.gz', sep = '' ), 
        cells = paste0(path_data,'barcodes.tsv.gz', sep = '' ),
        features = paste0(path_data,'features.tsv.gz', sep = '' ),
        feature.column = 1) 

# plot 

## filtering 

## a) create subset small enough to 'play with' 
# using equivalent of sample function of Seutat -> randomly select 50 000 cells 

sample_indices <-sample(c(1:length(annotation.mouseStSt$cell)), 50000, 
                        replace = F )

annotation.mouseStSt.Subset <-  annotation.mouseStSt[sample_indices, ]

lca.mouseStSt <- CreateSeuratObject(counts = counts, 
                                    project = "liver_cell_atlas_nafld", 
                                    min.cells = 3)

lca.mouseStSt.Subset <- subset(lca.mouseStSt, 
                               cells = annotation.mouseStSt.Subset[,'cell'] )

write.csv(annotation.mouseStSt.Subset, file = 'annot_mouseNafldAll_randomsubset.csv')
saveRDS(lca.mouseStSt.Subset, file = 'countTable_mouseNafld_randomsubset.rds')


# scRNA-seq QC metric.
# lca.mouseStSt[["percent.mt"]] <- PercentageFeatureSet(lca.mouseStSt, pattern = "^mt-")
# jpeg(file="vln_plot.jpeg")
# VlnPlot(lca.mouseStSt, features = c("nFeature_RNA", "nCount_RNA","percent.mt" ), 
#         ncol = 3,)
# dev.off()
# 
# jpeg(file="plot1.jpeg")
# FeatureScatter(lca.mouseStSt, feature1 = "nCount_RNA", 
#                         feature2 = "nFeature_RNA" )
# dev.off()
# 
# jpeg(file="plot2.jpeg")
# FeatureScatter(lca.mouseStSt, feature1 = "nCount_RNA", 
#                         feature2 = "percent.mt")
# dev.off()


## b) calculate average gene expression across all hepatocytes for each gene
# should result in vector of legth n = number of genes, index = gene names

## c) calculate average gene expression across all cells excluding hepatocytes 
# for each gene. Should result in vector of length n = number of genes, index =
# gene names 

## d) create a subset of all cells (keep all columns) but only features in 
#  universal damage signature 
