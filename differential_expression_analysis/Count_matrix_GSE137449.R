# Create Count Matrix 
#set working directory to the one that contains all reads 


setwd(dir = '../GSE137449/')

temp_males <- read.delim('GSE137449_Rawcounts_Males_Cyp2bnull_NASH_study.txt', 
                         sep = '\t')
# remove knock-outs
temp_males <- temp_males[,c(1:9)]

temp_females <- read.delim('GSE137449_Rawcounts_Females_Cyp2bnull_NASH_study.txt', 
                           sep = '\t')
# remove knock-outs
temp_females <- temp_females[,c(1:9)]

# add females, check if genes are in same order

if (identical(temp_males[,1],temp_females[,1]) == TRUE){
    count_table <- data.frame('gene_symbol' = temp_males[,1], 
                              temp_females[,c(2:9)], 
                              temp_males[,c(2:9)])
    print('yay!')
} else { print('error: different genes')}

length( unique( count_table$gene_symbol))
length( count_table$gene_symbol)

count_table.ag.p <- aggregate(cbind(WT_NDF_1,WT_NDF_2,WT_NDF_3,WT_NDF_4,WT_CDF_1,
                                    WT_CDF_2,WT_CDF_3,WT_CDF_4,WT_NDM_1,WT_NDM_2,
                                    WT_NDM_3,WT_NDM_4,WT_CDM_1,WT_CDM_2,WT_CDM_3,
                                    WT_CDM_4) ~ gene_symbol, 
                              data = count_table, FUN = mean )


write.table(count_table.ag.p, file = 'GSE137449_count_table', sep = ' ', row.names = F,
            col.names = T)
remove(count_table)
remove(count_table.ag.p)
remove(temp_females)
remove(temp_males)
