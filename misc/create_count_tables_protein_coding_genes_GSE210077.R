library(Seurat)
library(stringr)
library(GSEABase)
library(dplyr)
library(patchwork)
library(AUCell)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(biomaRt)

# In this data set, the authors have converted Ensembl Gene IDs to HGNC symbols, but have kept ENSEMBL Gene IDs for the sequences 
# that don't have corresponsding symbol. In some samples, they have kept ncbi sequence accession numbers instead. 
# In samples from donors 4 to 6 they have also kept a the suffix on the Ensebl Gene IDs for the version. 
# This is all messing up the FindMarkers() step at the end of the cell fate analysis. 
# I will now filter all rows that do not contain protein coding genes and make sure the HGNC symbols feature names only contain this nomenclature
# They are keeping many long non coding RNAs and pseudogenes that we do not need for this analysis. 

path_to_data <-  "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/human/Rahman_snRNAseq_GSE210077/"
sample_annotation <- read.csv(paste0(path_to_data, "GSE210077_sample_annotation.csv"))
cell_properties <- read.csv(paste0(path_to_data,"cell_properties_healthy_diseased_nucseq.csv"))


# Use Ensembl BioMart (human)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get mapping: Ensembl gene ID and biotype
protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

# Just keep the gene IDs
protein_coding_genes <- protein_coding_genes$ensembl_gene_id

# Retrieve the mapping between Ensembl IDs and HGNC symbols
gene_map_protein_coding <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = protein_coding_genes,  
  mart = ensembl
)

# Remove genes that don't have gene symbol name (883 from 22379)
gene_map_protein_coding <- gene_map_protein_coding[gene_map_protein_coding$hgnc_symbol != "", ]

# Name rows with Ensembl Gene ID 
# rownames(gene_map) <- gene_map$ensembl_gene_id

# Control for duplicates: there are no duplicates
length(gene_map_protein_coding[duplicated(gene_map_protein_coding$ensembl_gene_id), 1 ])
length(gene_map_protein_coding[duplicated(gene_map_protein_coding$hgnc_symbol), 1 ])
# There are 2911 gene symbols duplicated and no ensembl IDs duplicated

# We are only interest in hepatocytes, thefore, we extract the cell barcodes from the cells annotated as hepatocytes 
# by the authors of the publication 

cell_properties <-cell_properties %>% 
  dplyr::filter(., cell_type_final_healthy %in% c("Hep_1","Hep_2","Hep_3") | cell_type_final_injured %in% c("Central_Hep","Portal_Hep","IJ_1_Hep","IJ_2_Hep") )

#IDs have a suffix that has to be removed -> to format "AAACCCAAGACTACCT-1"
clean_ids <- gsub(".*_", "", cell_properties$index)
clean_ids<- stringr::str_replace(string = clean_ids, pattern = c(".1-0"), replacement = "-1")
clean_ids<- stringr::str_replace(string = clean_ids, pattern = c(".1-1"), replacement = "-1")

# create a list of count tables that are filtered by: 
# - only columns that have cell id corresponding to a hepatocyte
# - only rows that have at least 3 > RowSums 

cell_properties$cell_id  <- clean_ids 
sample_list <- list()

sample_i <- 1
# Initialize the Seurat object with the raw (non-normalized data).
for (sample_i in 1:length(sample_annotation$Donor)) {

  # Extract the cell barcodes to keep for sample i

  temp_cell_properties <- cell_properties %>% 
    dplyr::filter(., sample_id == sample_annotation$Sample_ID[[sample_i]])

  print(sample_annotation$Sample_ID[[sample_i]])

  temp_counts <- ReadMtx(
    mtx = paste0(path_to_data,"samples_analysed/", sample_annotation$file_name[sample_i], "_matrix.mtx.gz"),
  cells = paste0(path_to_data, "samples_analysed/", sample_annotation$file_name[sample_i], "_barcodes.tsv.gz"),
  features = paste0(path_to_data, "samples_analysed/", sample_annotation$file_name[sample_i], "_features.tsv.gz"),
  feature.column = 1)

  # filter cells first
  temp_counts <- temp_counts[,temp_cell_properties$cell_id]
  # Remove rows of features with zero counts across cells
  temp_counts <- temp_counts[rowSums(temp_counts) > 3, ]

  sample_list[[sample_i]] <- temp_counts

  remove(temp_counts)
   
}

# Filter to keep only protein coding genes and convert to HGNC Gene Symbol 
# for ach sample manually, because they are all different 

count_table <- as.data.frame(sample_list[[1]])

any(duplicated(rownames(count_table)))  # Should be FALSE

# Add gene symbols to the count table
count_table$ensembl_id <- rownames(count_table)
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's false!

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

print(dim(sample_list[[1]]))
sample_list[[1]] <- as.matrix(count_table)
print(dim(sample_list[[1]]))

remove(count_table)

# Sample Donor 2

count_table <- as.data.frame(sample_list[[2]])

any(duplicated(rownames(count_table)))  # Should be FALSE -> is false

# Add gene symbols to the count table
count_table$ensembl_id <- rownames(count_table)
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's false!

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

print(dim(sample_list[[2]]))
sample_list[[2]] <- as.matrix(count_table)
print(dim(sample_list[[2]]))

remove(count_table)

# Sample Donor 3

count_table <- as.data.frame(sample_list[[3]])

any(duplicated(rownames(count_table)))  # Should be FALSE -> is false

# Add gene symbols to the count table
count_table$ensembl_id <- rownames(count_table)
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's false!

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

print(dim(sample_list[[3]]))
sample_list[[3]] <- as.matrix(count_table)
print(dim(sample_list[[3]]))

# Sample Donor 4 (!!) Beware for the donors 4 to 6 the ENSEMBL gene ID's have a suffix

count_table <- as.data.frame(sample_list[[4]])
any(duplicated(sub("\\..*", "", rownames(count_table))))  # Should be FALSE -> is false

# Add gene symbols to the count table
count_table$ensembl_id <- sub("\\..*", "", rownames(count_table))
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's TRUE, we need to aggregate counts by gene symbol
sum(duplicated(count_table$hgnc_symbol)) # one gene 

# Which gene -> PINX1
count_table[duplicated(count_table$hgnc_symbol), c("ensembl_id", "hgnc_symbol")]
#ensembl_id hgnc_symbol
# 5957 ENSG00000254093       PINX1
# row with aggregsted values fro both ensebl ids of duplicated gene symbol
new_row <- count_table[count_table$hgnc_symbol == "PINX1",] %>% dplyr::select(-ensembl_id) %>% dplyr::select(-hgnc_symbol) %>% colSums()
# remove the rows of that gene, we will add the aggregated values later
rownames(count_table[count_table$hgnc_symbol == "PINX1",])
# [1] "5956" "5957"
count_table <-count_table %>% slice(., -c(5956,5957))

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

count_table<-rbind(count_table, new_row) 
rownames(count_table) <- c(rownames(count_table)[1:(dim(count_table)[1]-1)],"PINX1" )


print(dim(sample_list[[4]]))
sample_list[[4]] <- as.matrix(count_table)
print(dim(sample_list[[4]]))

# Sample Donor 5 (!!) Beware for the donors 4 to 6 the ENSEMBL gene ID's have a suffix

count_table <- as.data.frame(sample_list[[5]])
any(duplicated(sub("\\..*", "", rownames(count_table))))  # Should be FALSE -> is false

# Add gene symbols to the count table
count_table$ensembl_id <- sub("\\..*", "", rownames(count_table))
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's TRUE, we need to aggregate counts by gene symbol
sum(duplicated(count_table$hgnc_symbol)) # two genes 

# Which genes? -> PINX1 and POLR2J3
count_table[duplicated(count_table$hgnc_symbol), c("ensembl_id", "hgnc_symbol")]
#ensembl_id hgnc_symbol
# 5746 ENSG00000285437     POLR2J3
# 6017 ENSG00000254093       PINX1

# row with aggregated values fro both ensebl ids of duplicated gene symbol
new_row_POLR2J3 <- count_table[count_table$hgnc_symbol == "POLR2J3",] %>% dplyr::select(-ensembl_id) %>% dplyr::select(-hgnc_symbol) %>% colSums()
new_row_PINX1 <- count_table[count_table$hgnc_symbol == "PINX1",] %>% dplyr::select(-ensembl_id) %>% dplyr::select(-hgnc_symbol) %>% colSums()

# remove the rows of that gene, we will add the aggregated values later
rownames(count_table[count_table$hgnc_symbol == "POLR2J3",])
# [1] "5744" "5746"
rownames(count_table[count_table$hgnc_symbol == "PINX1",])
# [1] "6016" "6017"
count_table <-count_table %>% slice(., -c(5744,5746,6016,6017))

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

count_table<-rbind(count_table, new_row_PINX1, new_row_POLR2J3) 
rownames(count_table) <- c(rownames(count_table)[1:(dim(count_table)[1]-2)],"PINX1", "POLR2J3")


print(dim(sample_list[[5]]))
sample_list[[5]] <- as.matrix(count_table)
print(dim(sample_list[[5]]))

# Sample Donor 6 (!!) Beware for the donors 4 to 6 the ENSEMBL gene ID's have a suffix

count_table <- as.data.frame(sample_list[[6]])
any(duplicated(sub("\\..*", "", rownames(count_table))))  # Should be FALSE -> is false

# Add gene symbols to the count table
count_table$ensembl_id <- sub("\\..*", "", rownames(count_table))
# Filter out empty gene symbols
count_table <- count_table %>% filter(.,ensembl_id %in% gene_map_protein_coding$ensembl_gene_id)

any(duplicated(gene_map_protein_coding$ensembl_gene_id))  # Should be FALSE -> it's FALSE
count_table <- left_join(count_table, gene_map_protein_coding,
                                 by = c("ensembl_id" = "ensembl_gene_id"))

sum(is.na(count_table$hgnc_symbol)) # 0, so nothing to remove 
any(duplicated(count_table$hgnc_symbol)) # -> it's TRUE, we need to aggregate counts by gene symbol
sum(duplicated(count_table$hgnc_symbol)) # two genes 

# Which genes? -> PINX1 and POLR2J3
count_table[duplicated(count_table$hgnc_symbol), c("ensembl_id", "hgnc_symbol")]
#ensembl_id hgnc_symbol
# 5736 ENSG00000285437     POLR2J3
# 6001 ENSG00000254093       PINX1

# row with aggregated values fro both ensebl ids of duplicated gene symbol
new_row_POLR2J3 <- count_table[count_table$hgnc_symbol == "POLR2J3",] %>% dplyr::select(-ensembl_id) %>% dplyr::select(-hgnc_symbol) %>% colSums()
new_row_PINX1 <- count_table[count_table$hgnc_symbol == "PINX1",] %>% dplyr::select(-ensembl_id) %>% dplyr::select(-hgnc_symbol) %>% colSums()

# remove the rows of that gene, we will add the aggregated values later
rownames(count_table[count_table$hgnc_symbol == "POLR2J3",])
# [1] "5734" "5736"
rownames(count_table[count_table$hgnc_symbol == "PINX1",])
# [1] "6000" "6001"
count_table <-count_table %>% slice(., -c(5734,5736,6000,6001))

rownames(count_table) <- count_table$hgnc_symbol
count_table <- count_table %>% 
  dplyr::select(-ensembl_id) %>%
  dplyr::select(-hgnc_symbol)

count_table<-rbind(count_table, new_row_PINX1, new_row_POLR2J3) 
rownames(count_table) <- c(rownames(count_table)[1:(dim(count_table)[1]-2)],"PINX1", "POLR2J3")


print(dim(sample_list[[6]]))
sample_list[[6]] <- as.matrix(count_table)
print(dim(sample_list[[6]]))

saveRDS(object = sample_list, file = paste0(path_to_data, "count_tables_sample_list_.RDS"))

