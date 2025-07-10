##########################################################
### Convert mouse gene list to human gene list and vs. ###
##########################################################
# Tim's original function modified by Paula 

### function to convert Human to mouse genes
genes <- read.csv('hepatocyte-damage-score/Data/Output/HDAG.csv')
genes <- genes$gene_symbol

fun_homoTO.FROMmouse <- function(genes, TO = TRUE ){
  # if TO parameter is TRUE then Homo -> Mouse, if FALSE then Mouse -> Homo
  
  options(connectionObserver = NULL)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(Orthology.eg.db)
  
  if(isTRUE(TO)){
    egs <- mapIds(org.Hs.eg.db, genes, "ENTREZID","SYMBOL")
    mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
    mapped$MUS <- mapIds(org.Mm.eg.db, 
                         as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
    mapped$MUS.ens <- mapIds(org.Mm.eg.db, 
                             as.character(mapped$Mus.musculus), "ENSEMBL", "ENTREZID")
  } else {
    egs <- mapIds(org.Mm.eg.db, genes, "ENTREZID","SYMBOL")
    mapped <- select( Orthology.eg.db, egs, "Homo.sapiens", "Mus.musculus")
    mapped$HOMO <- mapIds(org.Hs.eg.db, 
                          as.character(mapped$Homo.sapiens), "SYMBOL", "ENTREZID")
    mapped$HOMO[ is.na(mapped$Homo.sapiens)] <- NA
     # toupper(rownames(mapped))[ is.na( mapped$Homo.sapiens ) ]
    mapped$HOMO[rownames(mapped)=="Cd59a"]<- "CD59"
    mapped$HOMO.ens <- mapIds(org.Hs.eg.db, 
                              as.character(mapped$HOMO), "ENSEMBL", "SYMBOL")
    
  }
  return(mapped)
}

mapped_MusToHom <- fun_homoTO.FROMmouse(genes, FALSE)
# !!!!! WARNING !!!!!!!!
# beware that genes in column "HumanEnsembl" marked as NA, are not identified
# correctly in "HumanGeneSymbol" column
TopHDSGenesConverted <- 
  data.frame("MouseGeneSymbol" = head(genes,42),
             "HumanEnsembl" = head(mapped_MusToHom$HOMO.ens, 42),
             "HumanGeneID" = head(unname(unlist(mapped_MusToHom$HOMO)), 42))
sum(which(is.na(mapped_MusToHom$Homo.sapiens)) <= 42)
# ==> 8 from 42 top genes not found 
# This will be identified looked into manually and decided what to do about them

allGenes <- data.frame(MouseGeneSymbol = genes,
                       HumanEnsembl = mapped_MusToHom$HOMO.ens,
                       HumanGeneID = unname(unlist(mapped_MusToHom$HOMO)))

write.csv(allGenes, 
            file = 
              'hepatocyte-damage-score/Data/Output/HDAGmappedToHumanGenes.csv')

write.csv(TopHDSGenesConverted, 
          file = 
            'hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenes.csv')
