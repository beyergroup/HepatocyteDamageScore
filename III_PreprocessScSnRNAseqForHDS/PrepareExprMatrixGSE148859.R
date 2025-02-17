# Expression Matrices from scRNAseq GSE148859 to Apply HDS 
# output: expression matrix including only hepatocyte counts
# and annonation file to those cells 
library(Seurat)

sce.combined <- 
  load("~/Desktop/HepatocyteDamageScore/scRNAseq/GSE148859/GSE148859_SingleCellExprement.Rdata")

subsetHepatocytes <- subset(x = sce.combined, 
                             idents = c("Hepatocytes 1",
                                        "Hepatocytes 2",
                                        "Hepatocytes 3",
                                        "Hepatocytes 4"))

remove(sce.combined) 

temp <- subsetHepatocytes@meta.data

scAnnotations <- data.frame('sample' = temp$orig.ident, 
                             'cellId' = temp$Barcode,
                             'condition' = gsub(pattern = 'HOHO', 'Tak1_HEP+D138', 
                                                temp$batchid)
                            )

sc_annotations$condition <- gsub(pattern = 'HO', 'Tak1_HEP',
                                 sc_annotations$condition) 
row.names(scAnnotations) <- row.names(temp)

write.csv(scAnnotations, 
          '~/Desktop/HepatocyteDamageScore/scRNAseq/GSE148859/GSE148859HepatocytesScAnnotations.csv')

exprMatrix <- GetAssayData(object = subsetHepatocytes, 
                           slot = 'counts')
saveRDS(exprMatrix, 
        file = '~/Desktop/HepatocyteDamageScore/scRNAseq/GSE148859/GSE148859HepatocytesExprMatrix.rds')


