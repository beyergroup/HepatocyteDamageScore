# Generate List of Hepatocyte Damage Associated Genes #

# load our functions to generate list of genes associated with damage
source("~/Repositories/cell-damage-score/SharedFunctions.R")

# Input
# - all tables of DEGs per condition/treatment/disease
# - list of genes expressed in hepatocytes:'GenesExpressedInHepatocytes.csv'
# - function written by us damage_signature.func()
#


# load DEG tables
path <- 'hepatocyte-damage-score/Data/Input/DEGTables/'
studies <- c('GSE99010_CCl4_deg.csv',
             'GSE99010_WesternDiet_deg.csv',
             'GSE97234_AH_deg.csv',
             'GSE97234_ASH_deg.csv',
             'GSE83240_DEN_deg.csv',
             'GSE153580_STZ_deg.csv',
             'GSE148849_fastfooddiet_deg.csv',
             'GSE138419_AMLN_deg.csv',
             'GSE137449_CDAHFD_diet_deg.csv',
             'GSE135050_hfcfdiet_deg.csv',
             'GSE132040_young(6_9_12)vs_old(18_21_24)_deg.csv',
             'GSE123894_Fructose_only_DBA2Jstrain_deg.csv',
             'GSE119953_HFD_deg.csv',
             'GSE119953_DDC_deg.csv',
             'GSE119441_HFD_deg.csv',
             'GSE119441_PFOA_deg.csv',
             'GSE114261_STAM_deg.csv',
             'GSE111828_Acetaminophen_deg.csv'
)

DElist <- lapply( seq( studies ), function(ii){
  read.csv(paste0(path,studies[[ii]]))
})

names(DElist) <- basename(studies)

# load list of genes to be considered as possibly expressed in hepatocytes
geneFilter <- 
  read.csv(file = 
             'hepatocyte-damage-score/Data/Output/GenesExpressedInHepatocytes/GenesExpressedInHepatocytes.csv',
           row.names = 1)
geneFilter <- geneFilter$Genes

# create HDAG

HDAG <- damage_signature.func(DElist = DElist, 
                      geneFilter = geneFilter)

rownames(HDAG)<- 1:length(HDAG$gene_symbol)

write.csv(HDAG,'hepatocyte-damage-score/Data/Output/HDAG.csv')

