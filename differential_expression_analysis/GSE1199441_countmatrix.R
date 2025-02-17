# Create Count Matrix 
#set working directory to the one that contains all reads 
setwd(dir = '../GSE119441/')
temp_count_table <- read.delim('GSE119441_genecounts.txt')
temp_count_table$gene_symbol <- temp_count_table$name

annotation <-  read.csv('SraRunTable.txt')

count_table <- cbind('gene_symbol' = temp_count_table$gene_symbol,
                               temp_count_table[,2:17]
                          )
count_table$gene_symbol <- temp_count_table$gene_symbol 
colnames(count_table) <- c('gene_symbol',annotation$Run)


length( unique( count_table$gene_symbol))
length( count_table$gene_symbol)

#length( count_table$gene_symbol) > length( unique( count_table$gene_symbol))
# aggregate repeated genes

# this does not make sense 
count_table.ag.p <- aggregate(cbind(SRR7782890,SRR7782891,SRR7782892,SRR7782893,
                                    SRR7782894,SRR7782895,SRR7782896,SRR7782897,
                                    SRR7782898,SRR7782899,SRR7782900,SRR7782901,
                                    SRR7782902,SRR7782903,SRR7782904,SRR7782905) ~ gene_symbol, 
                              data = count_table, FUN = mean )

#check if final length makes sense
length( unique( count_table.ag.p$gene_symbol))

# remove first row if empty
count_table <- count_table.ag.p

## save 
write.table(count_table, 
            'GSE119441_count_table', sep = ' ', 
            row.names = F, col.names = T)
