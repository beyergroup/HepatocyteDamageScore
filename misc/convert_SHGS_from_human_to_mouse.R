# Convert human marker list to mouse marker list 

  fun_human_to_mouse <- function(genes, TO = TRUE ){
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

#SHGS 
shgs_genes_human <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025.rds")
shgs_genes_human <- shgs_genes_human$SHGS


shgs_genes_mouse <- fun_human_to_mouse(shgs_genes_human , TRUE)

# will add the missing ones manually at the end of list after looking them up 
shgs_genes_final <- unname(unlist(shgs_genes_mouse$MUS))
shgs_genes_final <- c(shgs_genes_final, c("Gbp2","Cxcl5","Cxcl3","Myl12a"))

saveRDS(shgs_genes_final, file ="/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds" )
