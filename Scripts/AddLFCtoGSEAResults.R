# how mant down/up regulated HDS genes in term 
library(stringr)
library(ggplot2)
res <- read.csv('~/Desktop/HepatocyteDamageScore/Results/HUDScharacterization/SelectedG_Profiler_Results.csv',
        sep = ';')
HUDS <-read.csv('~/Desktop/HepatocyteDamageScore/Universal Damage Signature/HUDS.csv')

genesInTermList <- lapply(seq(1,length(res$Genes)), function(ii){
  genesInTerm <- str_split_1(res$Genes[ii], pattern = ',')
  return(genesInTerm)
})

LFCinHUDSList <- lapply(seq(1,length(genesInTermList)), function(ii){
  LFCsInTerm <- HUDS$direction_foldchange[match(genesInTermList[[ii]], toupper(HUDS$gene_symbol))]
  return(LFCsInTerm)
  })

# add number of - or + to term
res$negLFC <- NA
res$posLFC <- NA
i <- 1

while (i <= 22) {
  res$negLFC[i] <- sum(LFCinHUDSList[[i]] == -1 )
  res$posLFC[i] <- sum(LFCinHUDSList[[i]] == 1 )
  i <- i + 1
}

write.csv(res, '~/Desktop/HepatocyteDamageScore/Results/HUDScharacterization/SelectedG_Profiler_Results_withLFCs.csv' )
