# Create Count Matrix 
#set working directory to the one that contains all reads 
setwd(dir = '../GSE114261/')
run_list <- list.files(path = ".", pattern = '.txt', all.files = T,
           full.names = T)
#last one is not a run
run_list <- run_list[-9]
filename <- run_list[1]

# initiate count table with first run
temp <- read.delim(filename)
count_table <- data.frame('gene_symbol' = temp[,1], 'run1' = temp[,2])
names(count_table)[names(count_table) == 'run1'] <- run_list[1]

# add the other runs
i = 2
for (i in 2:(length(run_list))) {
  run <- run_list[i]
  temp <- read.delim(file = run)
  
  if (identical(count_table[,1],temp[,1]) == TRUE){
    count_table$temp <- temp[,2]
    names(count_table)[names(count_table) == 'temp'] <- run
    print('yay!')
  } else { print('error: different genes')}
  
}
annotation <- read.csv(file = 'SraRunTable.txt')
colnames(count_table) <- c('gene_symbol',annotation[1:8,"Run"])


library(biomaRt)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "mmusculus")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

#Ensemble IDs instead of Gene Symbols 

genes <- count_table$gene_symbol

#searchFilters(mart = ensembl, pattern = "symbol")
#searchFilters(mart = ensembl, pattern = "ensemble*")
G_list <- getBM( filter = 'ensembl_gene_id', 
                 attributes= c("ensembl_gene_id","mgi_symbol"), 
                 values = genes,mart= ensembl)

# regard that, when ensembl ID has not matching gene symbol in G_list, 
# the following command will assign NA as the name of the gene 
count_table$gene_symbol <- G_list$mgi_symbol[match( count_table$gene_symbol,
                                                    G_list$ensembl_gene_id 
)]
#this shows us how many genes could not get assigned a gene symbol
length(count_table$gene_symbol[is.na(count_table$gene_symbol)] )
# 3796 genes lost

# remove them Gene.Symbols == NAs 
count_table <- count_table[is.na(count_table$gene_symbol) == FALSE ,]



#does this match: YES -> bug fixed 
length( unique( G_list$mgi_symbol))
length( unique( count_table$gene_symbol))
length( count_table$gene_symbol)
#length( count_table$gene_symbol) > length( unique( count_table$gene_symbol))
# aggregate repeated genes

# this does not make sense 
count_table.ag.p <- aggregate(cbind(SRR7140599,SRR7140600,SRR7140601,SRR7140602,SRR7140603,SRR7140604, 
                                    SRR7140605,SRR7140606) ~ gene_symbol, 
                              data = count_table, FUN = mean )

#check if final length makes sense
length( unique( count_table.ag.p$gene_symbol))

# remove first row if empty
count_table.ag.p <- count_table.ag.p[-1,]
count_table <- count_table.ag.p

## save 
write.table(count_table, 
          'GSE114261_count_table', sep = ' ', 
          row.names = F, col.names = T)

remove(count_table)
remove(temp)


