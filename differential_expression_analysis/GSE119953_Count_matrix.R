# Create Count Matrix 
#set working directory to the one that contains all reads
setwd(dir = 'GSE119953/')
run_list <- read.delim('List_of_files.txt',header = F)
annotation <- read.delim('SraRunTable.txt', sep = ',')

# initiate count table with first run
temp <- read.delim(run_list[1,], header = F)
#delete last 5 rows (they contain things like alignment not unique, etc)
temp <- head(temp,-5)
# initiate count table
count_table <- data.frame('gene_symbol' = temp[,1], 'run1' = temp[,2])
names(count_table)[names(count_table) == 'run1'] <- annotation$GEO_Accession..exp.[1]
remove(temp)

# add the other runs
i = 2
for (i in 2:(length(run_list[,1]))) {
  print(run_list[i,])
  temp <- read.delim(run_list[i,], header = F)
  temp <- head(temp,-5)
  
  if (identical(count_table[,1],temp[,1]) == TRUE){
    count_table$temp <- temp[,2]
    names(count_table)[names(count_table) == 'temp'] <- 
      annotation$GEO_Accession..exp.[i]
    print('yay!')
  } else { print('error: different genes')}
  
}



write.table(count_table, file = 'GSE119953_count_table', sep = ' ', row.names = F,
            col.names = T)
remove(count_table)
remove(temp)
