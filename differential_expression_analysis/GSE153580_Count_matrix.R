# Create Count Matrix 
#set working directory to the one that contains all reads 
setwd(dir = '../GSE153580/')
# (!!) in this case I went manually through text files and removed header and tail comments 
run_list <- c('11.txt','12.txt','15.txt','17.txt','19.txt','20.txt','22.txt','23.txt',
              '24.txt','26.txt','27.txt','28.txt','32.txt','33.txt',
              '34.txt','36.txt','38.txt','39.txt')

# initiate count table with first run
temp <- read.delim(run_list[1], sep = '\t', header = FALSE)
count_table <- data.frame('gene_symbol' = temp[,1], 'run1' = temp[,2])
names(count_table)[names(count_table) == 'run1'] <- run_list[1]

# add the other runs
i = 2
for (i in 2:(length(run_list))) {
  run <- run_list[i]
  temp <- read.delim(file = run, sep = '\t', header = FALSE)
  
  if (identical(count_table[,1],temp[,1]) == TRUE){
    count_table$temp <- temp[,2]
    names(count_table)[names(count_table) == 'temp'] <- run
    print('yay!')
  } else { print('error: different genes')}
  
}

annotation <- read.csv('SraRunTable.txt')
colnames(count_table) <- c("gene_symbol", annotation$Run)

## replace Emseble IDs 
library(biomaRt)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "mmusculus")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

# CASE 1: Ensemble IDs instead of Gene Symbols 
tail(count_table,20)
#get rid of rows with p-value = NA 
#count_table <- count_table[is.na(count_table$pvalue) == FALSE ,]
genes <- count_table$gene_symbol
# for transcript ensembl IDs do this: 
genes <- sub( "\\..*","",genes)
count_table$gene_symbol <- sub('\\..*','', count_table$gene_symbol ) 
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

# remove them Gene.Symbols == NAs 
count_table <- count_table[is.na(count_table$gene_symbol) == FALSE ,]

#does this match: YES -> bug fixed 
length( unique( G_list$mgi_symbol))
length( unique( count_table$gene_symbol))

# this does not make sense 
count_table.ag.p <- aggregate(cbind(SRR12121750,SRR12121751,SRR12121752,SRR12121753,SRR12121754,SRR12121755,
SRR12121756, SRR12121757, SRR12121758, SRR12121759, SRR12121760, SRR12121761, SRR12121762,SRR12121763,SRR12121764,SRR12121765,SRR12121766,SRR12121767) ~ gene_symbol, 
                              data = count_table, FUN = mean )

#check if final length makes sense
length( unique( count_table.ag.p$gene_symbol))

# remove first row if empty
count_table.ag.p <- count_table.ag.p[-1,]
count_table.ag.p[,-1] <-round(count_table.ag.p[,-1],0)

count_table <- count_table.ag.p

write.table(count_table, file = 'GSE153580_count_table', sep = ' ', row.names = F,
            col.names = T)

remove(count_table)
remove(temp)
