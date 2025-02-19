#=================================================================#
# ========== GO annotation and plotting ========= =
#=================================================================#
options(connectionObserver = NULL)
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")
source("hepatocyte-damage-score/Data/Input/justGO_ROBERT.R")
library("org.Mm.eg.db")
library("GO.db")
library("DBI")
library( "biomaRt")
library(stringr)
library(ggplot2)
# BiocManager::install(version = "3.18")
# devtools::install_version('dbplyr', version = '2.3.4')
library(dbplyr)

# This is script is mostly written by Tim 
# and adapted to apply on HDS

outputPath <- 'hepatocyte-damage-score/Data/Output/Results/'

# load biomart annotations
mart_mouse <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl")

entr2gName <- getBM( attributes=c('external_gene_name', 
                                  "entrezgene_id" ),  
                     mart = mart_mouse)

esmbID2gName <- getBM( attributes=c('ensembl_gene_id', 
                                    "external_gene_name" ),  
                       mart = mart_mouse)


# load a sparse binary matrix indicating membership 
# of genes (columns) in GO terms (rows)

gomatrix <- readRDS("hepatocyte-damage-score/Data/Input/gomatrix.31.08.23.rda")
  
### load GO term to ID, ID to ontology mappings   
goterms <- Term(GOTERM)
goont <- Ontology(GOTERM)
  

# load a gene background 

# Tim uses a global background with all mouse genes
tssGR <- readRDS("hepatocyte-damage-score/Data/Input/mouseTSS.GR.rda")
bckgr <- tssGR@ranges@NAMES

#My background, list of genes detected in hepatocyte snRNAseq data set
# bckgr <- read.csv('hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv')
# bckgr <- bckgr$Genes
  
  
# define background gene set = universe, extract the entrez gene IDs 
# that correspond to gene symbols 

universe <- unique(
  entr2gName$entrezgene_id[match(bckgr, entr2gName$external_gene_name)]
  ) 
universe <- as.character(universe[!is.na(universe)])
  
  
# load genes sets for over representation analysis  
  
genesDownSenes <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorSenNegative_adj_pvalue_0.01_rho_0.2.csv')
  )

genesDownSenes <- genesDownSenes$gene
  
genesDownProlif <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorProNegative_adj_pvalue_0.01_rho_0.2.csv'))

genesDownProlif <- genesDownProlif$gene

genesUpSenes <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorSenPositive_adj_pvalue_0.01_rho_0.2.csv')
  )

genesUpSenes <- genesUpSenes$gene

genesUpProlif <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorProPositive_adj_pvalue_0.01_rho_0.2.csv')
)
genesUpProlif <- genesUpProlif$gene 
  
genesIntersectDown <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorIntersectNegative_adj_pvalue_0.01_rho_0.2.csv')
    )
genesIntersectDown <- genesIntersectDown$gene
  
genesIntersectUp <- read.csv(
  paste0(outputPath,
         'CellFateAnalysis/genesCorIntersectPositive_adj_pvalue_0.01_rho_0.2.csv')
)
genesIntersectUp <- genesIntersectUp$gene

geneSet_list <- list(genesDownSenes, genesDownProlif, genesIntersectDown,
                     genesUpSenes, genesUpProlif, genesIntersectUp )
  
# run function on a list of gene sets = geneSet_list
  GOrobert_res <- lapply( seq(geneSet_list), function(ii)
    {
    print(ii)
    
    # define background gene set
     universe <- unique( entr2gName$entrezgene_id[match(  bckgr , entr2gName$external_gene_name)] ) 
     universe <- as.character( universe[!is.na(universe)] )
    
    
    # prepare gene set, convert ensembleIDs to entrezIDs
    geneset <-  geneSet_list[[ii]]
    geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
    geneset <- geneset[!is.na(geneset)]
    geneset <- as.character(geneset)
    print(length(intersect(geneset,colnames(gomatrix))))
    
    # apply Robert's function that given a sparse matrix of GO terms 
    # (columns = genes, rows = GO terms),
    # a geneset of interest, and a background set of genes (universe)
    # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
    # Note, multiplicity adjustment is performed for the representative terms only.
    RobertGO <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                       min.genes = 3, cut_max = 10 )
    
    return(RobertGO)
  })  
 
  
### select top N (maximum) significantly (qvalue <0.1) enriched "primary" terms per tested geneset
  
  names(GOrobert_res) <- c('NegativeSenscence',
                           'NegativeProliferative',
                           'NegativeIntersect',
                           'PositiveSenescence',
                           'PositiveProliferative',
                           'PositiveIntersect')
  
  #remove the elements of the list that do not have results,
  # i.e., the gene lists that do not have any enriched terms 
   GOrobert_res <- within(GOrobert_res, 
                          rm('PositiveSenescence',
                             'PositiveProliferative'))
  
  thrsh <- 0.1  #threshold for q-values to select primary terms
  topN_GOrb  <- Reduce(union , lapply(seq(GOrobert_res), function(ii, N = 35)
    {
    print(ii)
    datt <- GOrobert_res[[ii]]
    # select from primary terms only
     datt <- subset(datt$results,  Is.primary==TRUE & Primary.Fisher.adj < thrsh )
    #datt <- subset(datt$results,  Fisher < thrsh )
    datt <- datt[ order(datt$Fisher) ,]
    # return N top primar terms
    return( datt$GO.ID[ 1:min(N, nrow(datt)) ])
  }) )
  
  
### prepare data and make barplot for the selected terms
  {
    # make a table from a list of results
    toPlot <- Reduce( rbind , lapply( seq(GOrobert_res), function(ii)
      {
      XX <- GOrobert_res[[ii]]$results[ 
        match(  topN_GOrb, GOrobert_res[[ii]]$results$GO.ID), ]
      XX$dataset <- names(GOrobert_res)[ii]
      XX$GO.ID <- topN_GOrb
      return(XX)
    }))
    
    # make sure you habe all term names correct
    toPlot$Term <- goterms[ match( toPlot$GO.ID, names(goterms))]
    # make sure you habe all ontologies correct
    toPlot$Ontology <- goont[ match( toPlot$GO.ID, names(goont))]
    
    toPlot$log2Enrichment <- ifelse( toPlot$Fisher>0.05| 
                                        is.na(toPlot$log2Enrichment), NA,
                                      toPlot$log2Enrichment  )
    # 
    # trim term names, if longer than N charachters
     toPlot$Term <- sapply( sub( "__.*","",toPlot$Term) , str_trunc, 50)
     toPlot$Primary <- sapply( sub( "__.*","",toPlot$Primary) , str_trunc, 50)
    # # order columns same as tested gene sets

    toPlot$dataset <- factor(toPlot$dataset, 
                             levels = names(GOrobert_res))
    # 
    # plot!
   plotToSave <- ggplot(toPlot,
                        aes(x = -log10(Fisher),
                            y = reorder(Term,
                                      -log10(Fisher),
                                      FUN = mean,
                                      na.rm = TRUE),
                            fill = log2Enrichment)) + 
      scale_fill_binned(type = "viridis", 
                        na.value = "grey", 
                        n.breaks =6 ) +  
      ylab("GO terms") +
      geom_bar(stat = "identity", 
               colour="black") + 
      facet_grid(rows =  vars(Ontology), 
                 cols = vars(dataset),
                 scales = "free_y",
                 space = "free_y"
                 )+
      theme_bw() + theme( text = element_text(size = 20 ),
                          axis.text.x = element_text(size = 16))
    

    
  }
  
  pdf(height = 12, 
      width = 20, 
      file = paste0(outputPath, 
                    "CellFateAnalysis/",
                    "GoAnalysisBACKGROUD_allGenes_AllNegativeGenesets_MinGenes3CutMax10Top35_PrimaryTerms0.1.pdf"))
  plotToSave
  dev.off()
  
  png(height = 12, width = 20,
      units = 'in',
      res = 400,
      file = paste0(outputPath,
                    "CellFateAnalysis/PNGs/",
                    "GoAnalysisBACKGROUD_allGenes_AllNegativeGenesets_MinGenes3CutMax10Top35_PrimaryTerms0.1.png"))
  plotToSave
  dev.off() 
  
  ### 
  
  topN_GOrb  <- Reduce(union , lapply(4, function(ii, N = 35)
  {
    print(ii)
    datt <- GOrobert_res[[ii]]
    # select from primary terms only
    datt <- subset(datt$results,  Is.primary==TRUE & Primary.Fisher.adj < thrsh )
    #datt <- subset(datt$results,  Fisher < thrsh )
    datt <- datt[ order(datt$Fisher) ,]
    # return N top primar terms
    return( datt$GO.ID[ 1:min(N, nrow(datt)) ])
  }) )
  
  toPlot <- Reduce( rbind , lapply( 4, function(ii)
  {
    XX <- GOrobert_res[[ii]]$results[ 
      match(  topN_GOrb, GOrobert_res[[ii]]$results$GO.ID), ]
    XX$dataset <- names(GOrobert_res)[ii]
    XX$GO.ID <- topN_GOrb
    return(XX)
  }))
  
  
  toPlot$Term <- goterms[ match( toPlot$GO.ID, names(goterms))]
  # make sure you habe all ontologies correct
  toPlot$Ontology <- goont[ match( toPlot$GO.ID, names(goont))]
  
  toPlot$log2Enrichment <- ifelse( toPlot$Fisher>0.05| 
                                     is.na(toPlot$log2Enrichment), NA,
                                   toPlot$log2Enrichment  )
  # 
  # trim term names, if longer than N charachters
  toPlot$Term <- sapply( sub( "__.*","",toPlot$Term) , str_trunc, 50)
  toPlot$Primary <- sapply( sub( "__.*","",toPlot$Primary) , str_trunc, 50)
  # # order columns same as tested gene sets
  
  toPlot$dataset <- factor(toPlot$dataset, 
                           levels = names(GOrobert_res))
  # 
  # plot!
  plotToSave <- ggplot(toPlot,
                       aes(x = -log10(Fisher),
                           y = reorder(Term,
                                       -log10(Fisher),
                                       FUN = mean,
                                       na.rm = TRUE),
                           fill = log2Enrichment)) + 
    scale_fill_binned(type = "viridis", 
                      na.value = "grey", 
                      n.breaks =6 ) +  
    ylab("GO terms") +
    geom_bar(stat = "identity", 
             colour="black") + 
    facet_grid(rows =  vars(Ontology), 
               scales = "free_y",
               space = "free_y"
    )+
    theme_bw() + theme( text = element_text(size = 20 ),
                        axis.text.x = element_text(size = 16))
  
  pdf(height = 12, 
      width = 20, 
      file = paste0(outputPath, 
                    "CellFateAnalysis/",
                    "GoAnalysisBACKGROUD_allGenes_PositiveIntersect_MinGenes3CutMax10Top35_PrimaryTerms0.1.pdf"))
  plotToSave
  dev.off()
  
  png(height = 12, width = 20,
      units = 'in',
      res = 400,
      file = paste0(outputPath,
                    "CellFateAnalysis/PNGs/",
                    "GoAnalysisBACKGROUD_allGenes_PositiveIntersect_MinGenes3CutMax10Top35_PrimaryTerms0.1.png"))
  plotToSave
  dev.off() 
  