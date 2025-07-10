# Figure 3 & Supplementary Figure 6 - HDS panels

# Figure 3 Panel A: Xiao et al - human snRNAseq
{
# HDS of Human hepatocytes from GSE189600 (Xiao et. al, 2023)
# - input is raw data, human version of HDAG list, our functions
# - subset, then normalized counts for each sample and the merged all 
# normalized counts
# - then added annotation from authors to metadata and extracted counts from 
# cells annotated as hepatocytes (four types) 
# - then I calculated the weighted HDS of those hepatocytes and made vioin
# plots of the HDS distribution per condition and hepatocyte type


library(Seurat)
library(rtracklayer)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)


pathData = 
  'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/HumanData/'

outputPath = 'hepatocyte-damage-score/Data/Output/Results/'

samplesNames <- c('GSM6808755_Y21822A1', 'GSM6808756_Y31122A5',
                  'GSM6808757_Y31122A8', 'GSM6808758_Y21822A2',
                  'GSM6808759_Y31122A6', 'GSM6808760_Y31122A7')

# read raw counts 
countsA1 <- ReadMtx(mtx =  paste0(pathData, samplesNames[1],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[1],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[1],'_features.tsv'),
                    feature.column = 2) 

A1Seurat <- CreateSeuratObject(counts = countsA1, 
                               project = samplesNames[1],
                               min.cells = 3)

A1Seurat[["percent.mt"]] <- PercentageFeatureSet(A1Seurat, 
                                                 pattern = "^MT-")


A1Seurat <- subset(A1Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A1Seurat <- SCTransform(A1Seurat)


remove(countsA1)
gc()

###

countsA2 <- ReadMtx(mtx =  paste0(pathData, samplesNames[2],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[2],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[2],'_features.tsv'),
                    feature.column = 2) 

A2Seurat <- CreateSeuratObject(counts = countsA2, 
                               project = samplesNames[2],
                               min.cells = 3)

A2Seurat[["percent.mt"]] <- PercentageFeatureSet(A2Seurat, 
                                                 pattern = "^MT-")
A2Seurat <- subset(A2Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A2Seurat <- SCTransform(A2Seurat)

remove(countsA2)
gc()


## 

countsA3 <- ReadMtx(mtx =  paste0(pathData, samplesNames[3],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[3],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[3],'_features.tsv'),
                    feature.column = 2) 

A3Seurat <- CreateSeuratObject(counts = countsA3, 
                               project = samplesNames[3],
                               min.cells = 3)

A3Seurat[["percent.mt"]] <- PercentageFeatureSet(A3Seurat, 
                                                 pattern = "^MT-")

A3Seurat <- subset(A3Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A3Seurat <- SCTransform(A3Seurat)

remove(countsA3)
gc()

##

countsA4 <- ReadMtx(mtx =  paste0(pathData, samplesNames[4],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[4],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[4],'_features.tsv'),
                    feature.column = 2) 

A4Seurat <- CreateSeuratObject(counts = countsA4, 
                               project = samplesNames[4],
                               min.cells = 3)

A4Seurat[["percent.mt"]] <- PercentageFeatureSet(A4Seurat, 
                                                 pattern = "^MT-")

A4Seurat <- subset(A4Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A4Seurat <- SCTransform(A4Seurat)

remove(countsA4)
gc()
##
countsA5 <- ReadMtx(mtx =  paste0(pathData, samplesNames[5],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[5],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[5],'_features.tsv'),
                    feature.column = 2) 

A5Seurat <- CreateSeuratObject(counts = countsA5, 
                               project = samplesNames[5],
                               min.cells = 3)

A5Seurat[["percent.mt"]] <- PercentageFeatureSet(A5Seurat, 
                                                 pattern = "^MT-")

A5Seurat <- subset(A5Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A5Seurat <- SCTransform(A5Seurat)


remove(countsA5)
gc()

countsA6 <- ReadMtx(mtx =  paste0(pathData, samplesNames[6],'_matrix.mtx'), 
                    cells = paste0(pathData, samplesNames[6],'_barcodes.tsv'),
                    features = paste0(pathData, samplesNames[6],'_features.tsv'),
                    feature.column = 2) 

A6Seurat <- CreateSeuratObject(counts = countsA6, 
                               project = samplesNames[6],
                               min.cells = 3)

A6Seurat[["percent.mt"]] <- PercentageFeatureSet(A6Seurat, 
                                                 pattern = "^MT-")

A6Seurat <- subset(A6Seurat, subset = 
                     nFeature_RNA > 500 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)

A6Seurat <- SCTransform(A6Seurat)


remove(countsA6)
gc()

allSamples <- list(A1Seurat, A2Seurat, A3Seurat,
                   A4Seurat, A5Seurat, A6Seurat)

#If you want to merge the normalized data matrices as well as the raw count matrices
# pass merge.data = TRUE

mergedSamples <- merge(allSamples[[1]], y = c(allSamples[[2]],
                                              allSamples[[3]],
                                              allSamples[[4]],
                                              allSamples[[5]],
                                              allSamples[[6]]), 
                       add.cell.ids = c("Normal1", "Normal2", "Normal3",
                                        "NASH1", "NASH2", "NASH3"), 
                       project = "HumanXiaoEtAl",
                       merge.data = TRUE)

remove(allSamples)
remove(A1Seurat, A2Seurat, A3Seurat, A4Seurat, A5Seurat, A6Seurat)

metaDataAnnotations <- read.table(
  paste0(pathData,
         'scitranslmed.adc9653_data_files_s1_to_s6/humanSnRNAseqSeuratAnnotations.csv'), 
  sep = ';', dec = ',', header = TRUE)

mergedSamples@meta.data$cell_id <- rownames(mergedSamples@meta.data)


temp <- metaDataAnnotations
rownames(temp) <- temp$X
temp <- temp[intersect(temp$X, rownames(mergedSamples@meta.data)),]

mergedSamples@meta.data$cellBarcode <- NA
mergedSamples@meta.data$cellBarcode[match(temp$X, rownames(mergedSamples@meta.data))] <- 
  temp$X
mergedSamples@meta.data$CellType[match(temp$X, rownames(mergedSamples@meta.data))] <- 
  temp$CellCluster

# mark whick cells to filter out of Seurat object (cells that are not annotated)
mergedSamples@meta.data$deleteNA <- 
  rep("not na", length(mergedSamples@meta.data$CellType))
mergedSamples@meta.data$deleteNA[is.na(mergedSamples@meta.data$CellType)] <-
  "na"

mergedSamples@meta.data <- subset(mergedSamples@meta.data,
                                  subset = deleteNA == 'not na')

saveRDS(object = mergedSamples, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanSamplesMerged.rds')

mergedHepatocytes <- subset(x = mergedSamples, 
                            subset = CellType ==  "hPC-Hep" |
                              CellType == "hNASH-Hep"  |
                              CellType == "hPP-Hep" |
                              CellType == "hInt-Hep" )

saveRDS(object = mergedHepatocytes, 
        file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')
# Very suspicious: after subsetting the subset has as many cells as genes 

# calculate HDS with HDAG list mapped to human genes
mergedHepatocytes <- readRDS( file = 'hepatocyte-damage-score/Data/Output/GSE189600SeuratAllHumanHepMerged.rds')

humanHDAG <- 
  read.csv(
    file = 'hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesOnlyHumanGenesManual.csv',
    sep = ';')

# load functions to calculate HDS
source('SharedFunctions.R')

# create expression Matrix
ExprMatrix <- GetAssayData(mergedHepatocytes)

# dummy values, because rank should as in list
humanHDAG$mean_rank <- 1:42

l <- dim(ExprMatrix)[2]

HDS <- DS_calc.func(ExprMatrix[,1:floor((l/2))], DSignature = humanHDAG,
                    geneIDname = 'HumanGeneID'
                    
)
gc()

HDS <- c(HDS, DS_calc.func(ExprMatrix[,(floor(l/2)+1):l], DSignature = humanHDAG,
                           geneIDname = 'HumanGeneID'
) 
)

length(HDS) == l 

# Prepare HDS results for plotting

DataFrameHDS <- data.frame('HDS' = unname(HDS),
                           'CellID' = names(unlist(HDS)))

identical(DataFrameHDS$CellID, rownames(mergedHepatocytes@meta.data))

mergedHepatocytes@meta.data$sample <- 
  gsub('_.*', replacement = '', mergedHepatocytes@meta.data$cellBarcode)

DataFrameHDS$sample <- mergedHepatocytes@meta.data$sample

DataFrameHDS$condition <-gsub('[1-3]_.*', 
                              replacement = '', 
                              mergedHepatocytes@meta.data$cellBarcode)


DataFrameHDS$HepatocyteType <- mergedHepatocytes@meta.data$CellType


DataFrameHDS$condition <- factor(DataFrameHDS$condition,
                                 ordered = TRUE,
                                 levels = c('Normal', 'NASH'))

DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                      ordered = TRUE,
                                      levels = c('hPP-Hep', 'hInt-Hep',
                                                 'hPC-Hep','hNASH-Hep'))

# One-sided Wilcoxon Rank Sum Test 
x <- DataFrameHDS$HDS[DataFrameHDS$condition == "NASH"]
y <- DataFrameHDS$HDS[DataFrameHDS$condition == "Normal"]
wilcox.test(x = x, y = y, data = DataFrameHDS, alternative = "greater")
#data:  x and y
#W = 208473462, p-value < 2.2e-16

# Ploting 

# Functions for ploting statistics
median_IQR <- function(x){
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}


PlotHDSSamples <- 
  ggplot( DataFrameHDS,
          aes( x = sample,
               y = HDS,
               color = condition
          )) + 
  geom_violin(trim = TRUE, lwd = 1.5) + 
  scale_color_colorblind() +
  theme( text = element_text(size=25),
         plot.title = element_text(size=25),
         legend.position = 'right',
         axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5
         )) + guides(colour = guide_legend(override.aes = aes(label = ""),
                               title="condition"))+  labs(y = "HDS", x = '')  + 
  theme_classic() +
  stat_summary(geom = "linerange",
               fun.data = median_IQR, 
               size = 1.5 ,
               show.legend = FALSE, 
               position = position_dodge(0.95)) +
  stat_summary(
    fun.y = "median",
    size = 1,
    geom = "crossbar",
    show.legend = FALSE,
    width = 0.2,
    position = position_dodge(0.95))

}

# Figure 3Panel B & C: Watson et al - human snRNAseq
{
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
  
  # samples 4,5, and 6 have ENSEMBL IDs with a version suffix! 
  
  # count tables were pr
  
  path_to_data <-  "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/human/Rahman_snRNAseq_GSE210077/"
  sample_annotation <- read.csv(paste0(path_to_data, "GSE210077_sample_annotation.csv"))
  cell_properties <- read.csv(paste0(path_to_data,"cell_properties_healthy_diseased_nucseq.csv"))
  
  # Hepatocyte Damage Associated Genes and function for calculating HDS 
  humanHDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesManuallyModified",
    sep = '\t')
  
  # dummy values, because rank should as in list
  humanHDAG$mean_rank <- 1:42
  
  source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')
  
  # filter hepaotcytes
  cell_properties <-cell_properties %>% 
    dplyr::filter(., cell_type_final_healthy %in% c("Hep_1","Hep_2","Hep_3") | cell_type_final_injured %in% c("Central_Hep","Portal_Hep","IJ_1_Hep","IJ_2_Hep") )
  # Convert to format "AAACCCAAGACTACCT-1"
  clean_ids <- gsub(".*_", "", cell_properties$index)
  clean_ids<- stringr::str_replace(string = clean_ids, pattern = c(".1-0"), replacement = "-1")
  clean_ids<- stringr::str_replace(string = clean_ids, pattern = c(".1-1"), replacement = "-1")
  cell_properties$cell_id  <- clean_ids 
  
  # Create Seurat objects with count tables previously 
  # cleaned up in script "create_count_tables_protein_coding_genes_GSE210077.R"
  # stores in misc/
  count_table_list <- readRDS(paste0(path_to_data,"count_tables_sample_list_.RDS"))
  sample_list <- list()
  sample_i <- 1
  
  # Initialize the Seurat object with the raw (non-normalized data).
  for (sample_i in 1:length(sample_annotation$Donor)) {
    
    temp_cell_properties <- cell_properties %>% dplyr::filter(., sample_id == sample_annotation$Sample_ID[[sample_i]])
    
    print(sample_annotation$Sample_ID[[sample_i]])
    
    sample_list[[sample_i]] <- CreateSeuratObject(
      counts = count_table_list[[sample_i]],
      project = sample_annotation$Sample_ID[[sample_i]] ,
      min.cells = 3, min.features = 200 
    )
    
    print(identical(rownames(sample_list[[sample_i]]@meta.data), temp_cell_properties$cell_id))
    
    sample_list[[sample_i]] <- PercentageFeatureSet(sample_list[[sample_i]],  pattern = "^MT-", col.name = "percent.mt")
    
    metadata_to_add <- temp_cell_properties[rownames(sample_list[[sample_i]]@meta.data) %in% temp_cell_properties$cell_id, 
                                            c("Condition", "Sex", "Age", "cell_type_final_injured", "cell_type_final_healthy", "cell_id")]
    identical(metadata_to_add$cell_id, rownames(sample_list[[sample_i]]@meta.data) )
    sample_list[[sample_i]]@meta.data <- cbind(sample_list[[sample_i]]@meta.data, metadata_to_add )
    
    sample_list[[sample_i]]<- SCTransform(sample_list[[sample_i]], vars.to.regress = "percent.mt", verbose = TRUE)
    sample_list[[sample_i]] <- RunPCA(sample_list[[sample_i]], verbose = TRUE)
    sample_list[[sample_i]] <- RunUMAP(sample_list[[sample_i]], dims = 1:30, verbose = TRUE)
    sample_list[[sample_i]] <- FindNeighbors(sample_list[[sample_i]], dims = 1:30, verbose = FALSE)
    sample_list[[sample_i]] <- FindClusters(sample_list[[sample_i]], verbose = FALSE)
    
    
  }

  # saveRDS(object = sample_list, file = paste0(path_to_data, "sample_seurat_list.RDS" ))
  
  sample_list <- readRDS(file = paste0(path_to_data, "sample_seurat_list.RDS"))
  
  # Calculate HDS per sample
  DimPlot(sample_list[[6]], group.by = "cell_type_final_injured")
  
  HDS_list <- list()
  sample_i <- 1
  
  for(sample_i in 1:length(sample_annotation$Donor)){
    
    HDS_list[[sample_i]] <- DS_calc.func(exprMatrices = 
                                           GetAssayData(sample_list[[sample_i]],
                                                        assay = 'SCT',
                                                        layer = 'counts'),
                                         DSignature = humanHDAG,
                                         geneIDname = 'HumanGeneID')
    
    if(identical(names(HDS_list[[sample_i]]), 
                 rownames(sample_list[[sample_i]]@meta.data) ) == TRUE){
      sample_list[[sample_i]]@meta.data$HDS <- c(unlist(HDS_list[[sample_i]]))
    }
    
    sample_list[[sample_i]]
    
  }
  
  saveRDS(sample_list, file = paste0(path_to_data, "GSE21007_Seurat_objects_with_HDS.rds"))
  sample_list <- readRDS(file = paste0(path_to_data, "GSE21007_Seurat_objects_with_HDS.rds"))
  
  # Prepare data to plot
  
  data_to_plot <- rbind(
    sample_list[[1]]@meta.data,
    sample_list[[2]]@meta.data,
    sample_list[[3]]@meta.data,
    sample_list[[4]]@meta.data,
    sample_list[[5]]@meta.data,
    sample_list[[6]]@meta.data
  )
  
  
  data_to_plot <- data_to_plot  %>% 
    mutate(.,fibrosis_level = case_when(
      (orig.ident == "AM042"  | orig.ident == "AM061"| orig.ident == "AM048") ~ "Healthy",
      orig.ident == "AM031" ~ "Fibrosis_F4",
      orig.ident == "AM062" ~ "Fibrosis_F2",
      orig.ident == "AM072" ~ "Fibrosis_F3"))
  
  data_to_plot$Condition <- factor(data_to_plot$Condition, levels = c("Normal", "Disease"))
  data_to_plot$fibrosis_level <- factor(data_to_plot$fibrosis_level, levels = c("Healthy", "Fibrosis_F2", "Fibrosis_F3", "Fibrosis_F4"))
  
  
  ############
  # Plotting #
  ############
  
  # Functions for ploting statistics
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  ## Final Figures
  
  # Copy paste results to plot with violin plots
  x <- data_to_plot$HDS[data_to_plot$Condition == "Disease"]
  y <- data_to_plot$HDS[data_to_plot$Condition == "Normal"]
  wilcox.test(x = x, y = y, data = data_to_plot, alternative = "greater")
  
  
  ggplot(data_to_plot, aes(y = HDS, x = Condition, colour = Condition)) +
    geom_violin(trim = TRUE, lwd = 1.5) + 
    scale_color_colorblind() +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR, 
                 size = 1.5 ,
                 show.legend = FALSE, 
                 position = position_dodge(0.95)) +
    stat_summary(
      fun.y = "median",
      size = 1,
      geom = "crossbar",
      show.legend = FALSE,
      width = 0.2,
      position = position_dodge(0.95)) +
    theme_classic() +
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank()) + 
    annotate("text",
             label= "One-Sided Wilcoxon-Test,\n p-value < 2.2e-16", y = 0.07, x = 1)
  
  ggsave(filename = paste0(path_to_data,"final_figures/", 
                           "Violin_plots_HDS-Hepatocyes-six_samples.png"),units = "cm",
         height = 15, width = 15)
  ggsave(filename = paste0(path_to_data,"final_figures/", 
                           "Violin_plots_HDS-Hepatocyes-six_samples.pdf"),units = "cm",
         height = 15, width = 15)
  
  
  ggplot(data_to_plot, aes(x = HDS, color= fibrosis_level)) +
    geom_density(lwd = 1.5, trim = TRUE) + facet_grid(Age ~ .) + scale_color_colorblind() + theme_bw() +
    theme( text = element_text(size = 24)
    ) + theme_bw()
  
  ggsave(filename = paste0(path_to_data, "final_figures/",
                           "Density_plots_HDS-Hepatocyes-six_samples_facet_grid_age.png"),units = "cm",
         height = 15, width = 20)
  
  ggsave(filename = paste0(path_to_data, "final_figures/",
                           "Density_plots_HDS-Hepatocyes-six_samples_facet_grid_age.pdf"),units = "cm",
         height = 15, width = 20)
  
}

# Figure 3 Panel D & E 
# Supplementary Figure 6 B: human spatial RNA-seq 10xVisum - Liver Cell Atlas

# We analyze patients with more than 10% steatosis.
# Patient information available in original publication ( patients H35 and H37) 
# Run on work station
{
  ### Applying HDS on Spatial Transcriptomics ###
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(gridExtra)
  library(AUCell)
  library(GSEABase)
  library(ggpubr)
  library(hdf5r)
  library(arrow)
  library(biomaRt)
  
  
  # Use Ensembl BioMart (human)
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  
  protein_coding_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    filters = "biotype",
    values = "protein_coding",
    mart = mart
  )
  genes_to_keep <- protein_coding_genes$hgnc_symbol[protein_coding_genes$hgnc_symbol != ""]
  
  setwd("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/")
  
  path_to_data <-  "Data/human/Liver_Cell_Atlas_spatial_RNAseq_GSE192742/"
  sample_directories <- list.files(path_to_data)
  # remove first row, because it contains cells annotations
  sample_directories <- sample_directories[2:(length(sample_directories)-2)]
  sample_annotation <- read.csv(paste0(path_to_data,"cell_annotations.csv" ))
  
  # Manually annotated regions of interest per sample 
  roi_file_names <- list.files(paste0(path_to_data,"ROI_annotations_Martin_Kann/"))
  sample_rois_files <- list()
  # There are two files per sample, use rbind to create one per sample
  for (file in 1:(length(roi_file_names))) {
    sample_rois_files[[file]] <- read.delim(paste0(path_to_data,"ROI_annotations_Martin_Kann/",roi_file_names[file]), sep = ",")
    colnames(sample_rois_files[[file]]) <- c("Barcode", "Cell_Identity")
  }
  
  list_manual_cell_identities <- list()
  
  # Create annotation for H35 sample 1 
  sample_rois_files[[1]]$Cell_Identity[sample_rois_files[[1]]$Barcode %in% sample_rois_files[[2]]$Barcode ] <- "Steatotic_Areas"
  missing_indices <- setdiff(sample_rois_files[[2]]$Barcode,sample_rois_files[[1]]$Barcode)
  sample_rois_files[[1]] <- rbind(sample_rois_files[[1]], sample_rois_files[[2]][sample_rois_files[[2]]$Barcode %in% missing_indices ,])
  setdiff(sample_rois_files[[2]]$Barcode,sample_rois_files[[1]]$Barcode)
  table(sample_rois_files[[1]]$Cell_Identity)
  
  list_manual_cell_identities[[1]] <- sample_rois_files[[1]]
  
  # Create annotation for H35 sample 2 
  sample_rois_files[[4]]$Cell_Identity <- "Steatotic_Areas"
  sample_rois_files[[3]]$Cell_Identity[sample_rois_files[[3]]$Barcode %in% sample_rois_files[[4]]$Barcode ] <- "Steatotic_Areas"
  missing_indices <- setdiff(sample_rois_files[[4]]$Barcode,sample_rois_files[[3]]$Barcode)
  sample_rois_files[[3]] <- rbind(sample_rois_files[[3]], sample_rois_files[[4]][sample_rois_files[[4]]$Barcode %in% missing_indices ,])
  table(sample_rois_files[[3]]$Cell_Identity)
  
  list_manual_cell_identities[[2]] <- sample_rois_files[[3]]
  
  # H37 
  sample_rois_files[[6]]$Cell_Identity <- "Steatotic_Areas"
  sample_rois_files[[5]]$Cell_Identity[sample_rois_files[[5]]$Barcode %in% sample_rois_files[[6]]$Barcode ] <- "Steatotic_Areas"
  missing_indices <- setdiff(sample_rois_files[[6]]$Barcode,sample_rois_files[[5]]$Barcode)
  sample_rois_files[[5]] <- rbind(sample_rois_files[[5]], sample_rois_files[[6]][sample_rois_files[[6]]$Barcode %in% missing_indices ,])
  table(sample_rois_files[[5]]$Cell_Identity)
  
  list_manual_cell_identities[[3]] <- sample_rois_files[[5]]
  
  names(list_manual_cell_identities) <- c( "Human_H35_sample1","Human_H35_sample2","Human_H37")
  
  
  #### Preparation of data before HDS calculation 
  
  h35_sample_1 <- Load10X_Spatial(
    data.dir = paste0(path_to_data,"Human_H35_sample1","/"),
    assay = "Spatial",
    slice = "Human_H35_sample1")
  # add sample ID to meta data 
  h35_sample_1@meta.data$sample_name  <- "Human_H35_sample1"
  # Keep only protein coding genes 
  h35_sample_1 <- subset(h35_sample_1, features = genes_to_keep, cells = list_manual_cell_identities[["Human_H35_sample1"]]$Barcode )
  h35_sample_1@meta.data$region_status <- "no steatosis"
  h35_sample_1@meta.data$region_status[rownames(h35_sample_1@meta.data) %in% list_manual_cell_identities[["Human_H35_sample1"]]$Barcode[list_manual_cell_identities[["Human_H35_sample1"]]$Cell_Identity == "Steatotic_Areas"]] <- "steatosis"
  
  SpatialPlot(h35_sample_1 ,
              group.by = "region_status", 
              pt.size.factor = 3,
              crop = TRUE,
              alpha = 0.7)
  
  
  h35_sample_2 <- Load10X_Spatial(
    data.dir = paste0(path_to_data,"Human_H35_sample2","/"),
    assay = "Spatial",
    slice = "Human_H35_sample2")
  # add sample ID to meta data 
  h35_sample_2@meta.data$sample_name  <- "Human_H35_sample2"
  # Keep only protein coding genes 
  h35_sample_2 <- subset(h35_sample_2, features = genes_to_keep, cells = list_manual_cell_identities[["Human_H35_sample2"]]$Barcode )
  h35_sample_2@meta.data$region_status <- "no steatosis"
  h35_sample_2@meta.data$region_status[rownames(h35_sample_2@meta.data) %in% list_manual_cell_identities[["Human_H35_sample2"]]$Barcode[list_manual_cell_identities[["Human_H35_sample2"]]$Cell_Identity == "Steatotic_Areas"]] <- "steatosis"
  
  SpatialPlot(h35_sample_2 ,
              group.by = "region_status", 
              pt.size.factor = 3,
              crop = TRUE,
              alpha = 0.7)
  
  
  h37_sample<- Load10X_Spatial(
    data.dir = paste0(path_to_data, "Human_H37","/"),
    assay = "Spatial",
    slice =  "Human_H37")
  # add sample ID to meta data 
  h37_sample@meta.data$sample_name  <- "Human_H37"
  # Keep only protein coding genes 
  h37_sample <- subset(h37_sample, features = genes_to_keep, cells = list_manual_cell_identities[[ "Human_H37"]]$Barcode )
  h37_sample@meta.data$region_status <- "no steatosis"
  h37_sample@meta.data$region_status[rownames(h37_sample@meta.data) %in% list_manual_cell_identities[[ "Human_H37"]]$Barcode[list_manual_cell_identities[[ "Human_H37"]]$Cell_Identity == "Steatotic_Areas"]] <- "steatosis"
  
  SpatialPlot(h37_sample ,
              group.by = "region_status", 
              pt.size.factor = 3,
              crop = TRUE,
              alpha = 0.7)
  
  
  ### Annotation succesful
  sample_list <- list(h35_sample_1, h35_sample_2, h37_sample)
  names(sample_list) <- names(list_manual_cell_identities)
  #sample_annotation$spot <- gsub("_.*", "", sample_annotation$spot)
  
  # Hepatocyte Damage Associated Genes and function for calculating HDS 
  humanHDAG <- read.csv(
    file = "Data/Output/HDAGTop42mappedToHumanGenesManuallyModified",
    sep = '\t')
  # dummy values, because rank should as in list
  humanHDAG$mean_rank <- 1:42
  
  source('../SharedFunctions.R')
  
  # load all data into a list, each element is the Seurat object corresponsing to a sample
  sample_nCount_violin_plots <- list()
  
  sample_i <- 1 
  for(sample_i in 1:length(sample_list)){
    sample_nCount_violin_plots[[sample_i]] <- VlnPlot(
      sample_list[[sample_i]], 
      features = "nCount_Spatial", pt.size = 0) + 
      theme(axis.text = element_text(size = 10)) + theme(legend.position = "none")
    
  }
  
  all_nCount_violin_plots <- 
    sample_nCount_violin_plots[[1]] | sample_nCount_violin_plots[[2]] | sample_nCount_violin_plots[[3]] 
  
  all_nCount_violin_plots
  
  # In GEO Accession the sample IDs have a short file, which is what is used in the annotation file to identify the sample
  # Note: each sample name in samples_directories corresponds to an entry here unique(sample_annotation$sample)
  
  # sample_list[[1]]: Human H35 (sample1) --> JBO14 | steatotic
  # sample_list[[2]]: Human H35 (sample2) --> JBO15 | steatotic
  # sample_list[[3]]: Human H37 --> JBO19           | steatotic
  
  sample_annotation <- sample_annotation %>%
    mutate(., sample_name = case_when(
      sample == "JBO14" ~ "Human_H35_sample1",
      sample == "JBO15" ~ "Human_H35_sample2",
      sample == "JBO18" ~ "Human_H36",
      sample == "JBO19" ~ "Human_H37",
      sample == "JBO22" ~ "Human_H38"
    ))
  
  
  sample_i <- 1
  for(sample_i in 1:length(sample_list)){
    
    sample_list[[sample_i]]<-subset(sample_list[[sample_i]], 
                                    nFeature_Spatial > 250 )
    
    # print(paste0("After: ", length(rownames(sample_list[[sample_name]]@meta.data))))
    sample_list[[sample_i]] <- NormalizeData(sample_list[[sample_i]])
    # Calculate HDS 
    HDS_temp <- DS_calc.func(
      exprMatrices = GetAssayData(sample_list[[sample_i]], layer = "data"),
      DSignature = humanHDAG,
      geneIDname = "HumanGeneID")
    
    # Add to metadata 
    sample_list[[sample_i]]@meta.data$HDS <- unname(HDS_temp)
    
    gc()
    
    
  }
  
  
  # Prepare data for plotting
  data_to_plot <- rbind(
    sample_list[[1]]@meta.data,
    sample_list[[2]]@meta.data,
    sample_list[[3]]@meta.data)
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., patient = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "H35",
      sample_name == "Human_H37" ~ "H37"))
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., all_info = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "female, 59 y/o, 35-40% stet., no fib.",
      sample_name == "Human_H37" ~ "male, 58 y/o 70% stet., clear fib."))
  
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., steatosis = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "35-40%",
      sample_name == "Human_H37" ~ "70%"))
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., age = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "59 y/o",
      sample_name == "Human_H37" ~ "58 y/o"))
  
  data_to_plot <- data_to_plot %>%
    mutate(., fibrosis = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "no fibrosis",
      sample_name == "Human_H37" ~ "clear pericellular"
    ))
  
  data_to_plot <- data_to_plot  %>% 
    mutate(., sex = case_when(
      (sample_name == "Human_H35_sample1" | sample_name == "Human_H35_sample2") ~ "female",
      (sample_name == "Human_H37") ~ "male"))
  
  data_to_plot$region_status <- factor(data_to_plot$region_status, levels = c("no steatosis","steatosis" ))
  data_to_plot$steatosis <- factor(data_to_plot$steatosis, levels = c("35-40%", "70%"))
  data_to_plot$fibrosis <- factor(data_to_plot$fibrosis, levels = c("no fibrosis","clear pericellular"))
  
  # Statistical tests
  
  y <- data_to_plot$HDS[(data_to_plot$patient == "H37" | data_to_plot$patient == "H35") & data_to_plot$region_status == "no steatosis"]
  x <- data_to_plot$HDS[(data_to_plot$patient == "H37"| data_to_plot$patient == "H35") & data_to_plot$region_status == "steatosis"]
  wilcox.test(x,y, alternative = "greater")
  # p-value = 6.086e-14
  
  # Plotting 
  
  # Functions for ploting statistics
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  ggplot(data_to_plot, aes(y = HDS, x = region_status, color = region_status)) +
    geom_violin(trim = TRUE, lwd = 1) + 
    scale_color_colorblind() +
    stat_summary(geom = "linerange", fun.data = median_IQR, size = 1, show.legend = FALSE, position = position_dodge(0.95)) +
    stat_summary(fun.y = "median", geom = "crossbar", show.legend = FALSE, width = 0.2,position = position_dodge(0.95)) +
    theme_base() + theme( text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  #ggsave(filename = paste0(path_to_data, "/output/Violin_plots_Spot_Steatosis_Status.png"),units = "cm",
  #      height = 15, width = 15)
  
  # Plot HDS on H&E staining 
  
  # Choose a floor and ceiling thesholds for color range: 10%-tile and 90%-tile 
  floor_ceiling <- quantile(data_to_plot$HDS, probs= c(0.1, 0.9))
  #floor_ceiling <- c(min(data_to_plot$HDS), max(data_to_plot$HDS)) 
  sample_i <- 1
  list_spatial_plots <- list()
  
  for (sample_i in 1:length(sample_list)) {
    
    list_spatial_plots[[sample_i]] <- 
      SpatialFeaturePlot(sample_list[[sample_i]], features = "HDS", pt.size.factor = 1,crop = FALSE, image.alpha = 1) +
      scale_fill_distiller(palette = "YlOrBr", direction = 1, limits = floor_ceiling,name = 'HDS', oob = squish) +
      theme(legend.key.size = unit(2, 'cm'), legend.position = "bottom" , text = element_text(size = 20)) + 
      ggtitle(names(sample_list)[sample_i])
    
    # ggsave(filename = paste0(path_to_data, "/output/",names(sample_list)[sample_i], "_Spatial_HDS.png"), 
    #       units = "cm", height = 30, width = 30)
    
    
  }
  # 
  # medianlines_groups <- data_to_plot %>%
  #   group_by(sample_name, region_status) %>%
  #   summarise(mean_HDS = median(HDS, na.rm = TRUE))
  # # Calculates mean of each distribution
  # ggplot(data_to_plot, aes(x = HDS, color = region_status)) +
  #   geom_density(lwd = 1.5, trim = TRUE) + facet_grid(sample_name ~ .) +
  #   theme( text = element_text(size = 8))  + geom_vline(aes(xintercept= mean_HDS, color = region_status), medianlines_groups, lwd = 1) +
  #   scale_color_colorblind() + theme_base() 
  # 
  # ggsave(filename = paste0(path_to_data, "/output/","DensityPlot_PerSample_with_Median.png"), units = "cm", height = 18, width = 25)
  # ggsave(filename = paste0(path_to_data, "/output/","DensityPlot_PerSample_with_Median.pdf"), units = "cm", height = 18, width = 25)
  # 
  # spatial plot with steatotic zones colored
  list_spatial_plots_rois <- list()
  for (sample_i in 1:length(sample_list)) {
    
    list_spatial_plots_rois[[sample_i]] <- 
      SpatialPlot(sample_list[[sample_i]],
                  group.by = "region_status", 
                  pt.size.factor = 1,
                  crop = FALSE,
                  alpha = 0.7) +
      theme(
        legend.key.size = unit(2, 'cm'), 
        legend.position = "bottom" ,
        text = element_text(size = 20)) + 
      ggtitle(names(sample_list)[sample_i])
    
    #ggsave(filename = paste0(path_to_data, "/output/",names(sample_list)[sample_i], "_Spot_Manual_Identity_Spatial_HDS.png"), 
    #      units = "cm", height = 30, width = 30)
    
    
  }
  
  # Figure for paper: 
  
  # we will only consider patients with more than 10 %, as in in the original publication these are considered diseased
  # and to be able to compare enough steatotic regions within the same patient sample
  
  data_to_plot_figure <- data_to_plot %>% dplyr::filter(., patient == c("H35", "H37"))
  
  ggplot(data_to_plot_figure, aes(y = HDS, x = patient, color = region_status)) +
    geom_violin(trim = TRUE, lwd = 1) + 
    scale_color_colorblind() +
    stat_summary(geom = "linerange", fun.data = median_IQR, size = 1, show.legend = FALSE, position = position_dodge(0.95)) +
    stat_summary(fun.y = "median", geom = "crossbar", show.legend = FALSE, width = 0.2,position = position_dodge(0.95)) +
    theme_base() + theme( text = element_text(size = 20))
  
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Violin_plots_HDS_perPatient_SteatoticSpot_H35_H37.png"),units = "cm",
  #        height = 15, width = 20)
  # 
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Violin_plots_HDS_perPatient_SteatoticSpot_H35_H37.pdf"),units = "cm",
  #        height = 15, width = 20)
  
  # Spatial plots with HDS overlay
  
  floor_ceiling <- quantile(data_to_plot_figure$HDS, probs= c(0.1, 0.9))
  
  SpatialFeaturePlot(sample_list[["Human_H35_sample1"]], 
                     features = "HDS", 
                     pt.size.factor = 1,
                     crop = FALSE, 
                     image.alpha = 1) +
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1, 
                         limits = floor_ceiling, 
                         name = 'HDS', 
                         oob = squish) +
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" , text = element_text(size = 20)) + 
    ggtitle("Patient H35 (slice 1)")
  
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H35_sample1", "_Spatial_HDS.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H35_sample1", "_Spatial_HDS.pdf"), 
  #        units = "cm", height = 30, width = 30)
  
  SpatialFeaturePlot(sample_list[["Human_H35_sample2"]], 
                     features = "HDS", 
                     pt.size.factor = 1,
                     crop = FALSE, 
                     image.alpha = 1) +
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1, 
                         limits = floor_ceiling, 
                         name = 'HDS', 
                         oob = squish) +
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" , text = element_text(size = 20)) + 
    ggtitle("Patient H35 (slice 2)")
  
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H35_sample2", "_Spatial_HDS.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H35_sample2", "_Spatial_HDS.pdf"), 
  #        units = "cm", height = 30, width = 30)
  
  
  SpatialFeaturePlot(sample_list[["Human_H37"]], 
                     features = "HDS", 
                     pt.size.factor = 1,
                     crop = FALSE, 
                     image.alpha = 1) +
    scale_fill_distiller(palette = "YlOrBr",
                         direction = 1, 
                         limits = floor_ceiling, 
                         name = 'HDS', 
                         oob = squish) +
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" , text = element_text(size = 20)) + 
    ggtitle("Patient H37")
  
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H37", "_Spatial_HDS.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, 
  #                          "/output/final_figures/","Human_H37", "_Spatial_HDS.pdf"), 
  #        units = "cm", height = 30, width = 30)
  
  
  
  
  # Classification Plots for supplementary data
  colors_regions <- pal_brewer("qual")(2)
  names(colors_regions) <- unique(sample_list[["Human_H35_sample1"]]@meta.data$region_status)
  
  SpatialPlot(sample_list[["Human_H35_sample1"]],
              group.by = "region_status", 
              pt.size.factor = 1.5,
              crop = FALSE,
              alpha = 1, 
              cols = colors_regions) + 
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" ,
          text = element_text(size = 20)) + 
    ggtitle("Patient H35 (slice 1)")
  
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H35_sample1",
  #                          "Region_classification.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H35_sample1",
  #                          "Region_classification.pdf"), 
  #        units = "cm", height = 30, width = 30)
  
  SpatialPlot(sample_list[["Human_H35_sample2"]],
              group.by = "region_status", 
              pt.size.factor = 1.5,
              crop = FALSE,
              alpha = 1, 
              cols = colors_regions) + 
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" ,
          text = element_text(size = 20)) + 
    ggtitle("Patient H35 (slice 2)")
  
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H35_sample2",
  #                          "Region_classification.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H35_sample2",
  #                          "Region_classification.pdf"), 
  #        units = "cm", height = 30, width = 30)
  
  SpatialPlot(sample_list[["Human_H37"]],
              group.by = "region_status", 
              pt.size.factor = 1.5,
              crop = FALSE,
              alpha = 1, 
              cols = colors_regions) + 
    theme(legend.key.size = unit(2, 'cm'), 
          legend.position = "bottom" ,
          text = element_text(size = 20)) + 
    ggtitle("Patient H37")
  
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H37",
  #                          "Region_classification.png"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # ggsave(filename = paste0(path_to_data, "/output/final_figures/",
  #                          "Human_H37",
  #                          "Region_classification.pdf"), 
  #        units = "cm", height = 30, width = 30)
  # 
  # 
  
  
  
  
}

# Figure 3 Panel F: HDS is not dependent on hepatocyte zonation, Xiao data

{
  outputPath <- "hepatocyte-damage-score/Data/Output/"
  dfHDS <- readRDS(file = paste0(outputPath, 
                                 'Xiao2023HDSsnRNAseqAllHep.rds'))
  
  dfHDS$sample <- ordered(dfHDS$sample, 
                          levels = c("NC 3 months", 
                                     "NASH 3 months", 
                                     "NC 9 months", 
                                     "NASH 9 months"))
  
  my_comparisons <- list( c("NC 3 months", "NASH 3 months"), 
                          c("NC 9 months", "NASH 9 months"), 
                          c("NC 3 months", "NC 9 months"),
                          c("NASH 3 months","NASH 9 months" ))
  
  
  HDSbyHepType <- 
    ggplot(dfHDS[dfHDS$HepatocyteType != 'mNASH-Hep1' & dfHDS$HepatocyteType != 'mNASH-Hep2', ],
           aes( x = sample,
                y = HDS,
                color = condition
           )) + 
    geom_violin(trim = TRUE, 
                lwd = 1.5) + 
    ylab("HDS") +
    scale_color_colorblind() +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR,
                 size = 1.5 ,
                 show.legend = FALSE,
                 position = position_dodge(0.95)) +
    stat_summary( fun.y = "median",  
                  size = 1, 
                  geom = "crossbar",
                  show.legend = FALSE,
                  width = 0.2,
                  position = position_dodge(0.95)) +
    theme_bw() + 
    theme( text = element_text(size = 24) ,
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, 
                                      hjust = 1),
           legend.title=element_blank()) +
    facet_wrap(~HepatocyteType) 
   # stat_compare_means(comparisons = my_comparisons) +
    #stat_summary(
    #  fun.data = stat_box_data,
    #  geom = "text",
    #  position = position_dodge(1),
    #  size = 3
    #)
  
  #pdf( height = 8, width = 12, 
 #      file = paste0(outputPath,
  #                   "Results/FiguresManuscript/HDSXiao2023HDSbyHepZone.pdf")) 
  print(HDSbyHepType)
  #dev.off()
  
}

# Supplementary Figure 6 Panel A: HDS is not dependent on hepatocyte zonation,
# Carlessi data 
{
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(patchwork)
  library(scales)
  
  # optional: load RDS with SCTransformed counts with HDS and SHGS calculated 
  # seurat_objects <- readRDS(paste0(path_to_data, "seurat_object_SCT_HDS_and_Senescence_scores.rds"))
  
  path_to_data <- "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/mouse/carlessi_snRNAseq_GSE200366/"
  setwd(path_to_data)
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  #SHGS 
  shgs_genes_mouse <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds")
  
  # load functions to calculate HDS
  source('/cellfile/datapublic/pungerav/cell-damage-score/SharedFunctions.R')
  
  seurat_object <- readRDS("Hep.rds")
  View(seurat_object@meta.data)
  
  DimPlot(seurat_object)
  
  seurat_object <- SCTransform(seurat_object)
  
  counts <- GetAssayData(seurat_object,
                         assay = 'SCT',
                         layer = 'counts')
  
  ## Calculate HDS
  HDS_calc <- DS_calc.func(exprMatrices = counts, DSignature = HDAG)
  
  if(identical(rownames(seurat_object@meta.data), names(HDS_calc)) == TRUE){
    seurat_object@meta.data$HDS <- unname(HDS_calc)
    print("Matching IDs!")
  }
  
  seurat_object@meta.data <- seurat_object@meta.data %>% 
    mutate(zone = case_when(
      cell_type == "Zone_1_Hep" ~ "Periportal",
      cell_type == "Zone_2_Hep" ~ "Intermediate",
      cell_type == "Zone_3_Hep" ~ "Pericentral",
      cell_type == "daHep" ~ "damaged Hepatocytes" 
    ))
  
  seurat_object@meta.data$zone <- factor(seurat_object@meta.data$zone, levels = c("Periportal", "Intermediate","Pericentral", "damaged Hepatocytes"), ordered = TRUE)
  
  
  test_periportal <- seurat_object@meta.data[seurat_object@meta.data$zone == "Periportal",]
  pairwise.wilcox.test(test_periportal$HDS,
                       g = test_periportal$Group, alternative = "greater", p.adjust.method = "none")
  
  test_intermediate <- seurat_object@meta.data[seurat_object@meta.data$zone == "Intermediate",]
  pairwise.wilcox.test(test_intermediate$HDS,
                       g = test_intermediate$Group, alternative = "greater", p.adjust.method = "none")
  
  test_pericentral<- seurat_object@meta.data[seurat_object@meta.data$zone == "Pericentral",]
  pairwise.wilcox.test(test_pericentral$HDS,
                       g = test_pericentral$Group, alternative = "greater", p.adjust.method = "none")
  
  
  
  
  
  median_IQR <- function(x){
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  
  
  ggplot( seurat_object@meta.data[seurat_object@meta.data$zone != "damaged Hepatocytes",],
          aes( x = Group,
               y = HDS,
               color = Group)) + 
    geom_violin(trim =TRUE, lwd = 1.5) + 
    scale_color_colorblind() +
    theme( text = element_text(size=25),
           plot.title = element_text(size=24),
           legend.position = 'bottom',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "HDS", x = 'condition') +
    facet_wrap(~zone) +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR, 
                 size = 1.5 ,
                 show.legend = FALSE, 
                 position = position_dodge(0.95)) +
    stat_summary(
      fun.y = "median",
      size = 1,
      geom = "crossbar",
      show.legend = FALSE,
      width = 0.2,
      position = position_dodge(0.95)) + theme_classic()
  
  ggsave(filename = paste0(path_to_data,"Violin_plots_HDS__per_zoneall_conditions.png"),
         units = "cm",height = 12, width = 24) 
  
  ggsave(filename = paste0(path_to_data,"Violin_plots_HDS__per_zoneall_conditions.pdf"),
         units = "cm",height = 12, width = 24) 
  
  ggplot( seurat_object@meta.data[seurat_object@meta.data$zone != "damaged Hepatocytes" &
                                    seurat_object@meta.data$Group != "TAA",],
          aes( x = Group,
               y = HDS,
               color = Group)) + 
    geom_violin(trim =TRUE, lwd = 1.5) + 
    scale_color_colorblind() +
    theme( text = element_text(size=25),
           plot.title = element_text(size=24),
           legend.position = 'bottom',
           axis.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = aes(label = ""),
                                 title="Health Status")) +
    labs(y = "HDS", x = 'condition') +
    facet_wrap(~zone) +
    stat_summary(geom = "linerange",
                 fun.data = median_IQR, 
                 size = 1.5 ,
                 show.legend = FALSE, 
                 position = position_dodge(0.95)) +
    stat_summary(
      fun.y = "median",
      size = 1,
      geom = "crossbar",
      show.legend = FALSE,
      width = 0.2,
      position = position_dodge(0.95)) + theme_classic()
  
  ggsave(filename = paste0(path_to_data,"Violin_plots_HDS__per_zone_no_TAA.png"),
         units = "cm",height = 12, width = 24) 
  
  ggsave(filename = paste0(path_to_data,"Violin_plots_HDS__per_zone_no_TAA.pdf"),
         units = "cm",height = 12, width = 24) 
  
}

# Figure 3 Panel G,TSCAN - 3 month NASH mouse  & 
# Supplementary Figure 6 Panel C: - rest of samples (3m NC, 9m NC, 9m NASH)
{
  # # HDS of mouse hepatocytes from GSE189600 (Xiao et. al, 2023)
  # - input is raw data, human version of HDAG list, our functions
  # - added annotation from authors to metadata and extracted counts from 
  # cells annotated as hepatocytes (four types) 
  # - subset,SCTransformed counts 
  # cell filtering as described in supplementary data from Xiaoo et al. 
  # nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 2 for mice 
  
  library(Seurat)
  library(rtracklayer)
  library(sctransform)
  library(ggplot2)
  library(tidyverse)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(glmGamPoi)
  library(TSCAN)
  library(scater)
  library(scales)
  
  # hepatocyte damage associated genes
  HDAG <- 
    read.csv(
      file = 'hepatocyte-damage-score/Data/Output/HDAG.csv')
  
  # load functions to calculate HDS
  source('SharedFunctions.R')
  
  pathData = 
    'hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/'
  
  outputPath = 'hepatocyte-damage-score/Data/Output/'
  
  FileNames <- list.files(path = pathData)
  
  metaDataAnnotations <- read.table('hepatocyte-damage-score/Data/Input/scRNAseq/GSE189600/MouseData/scitranslmed.adc9653_data_file_s1.csv', 
                                    sep = ';', dec = ',', header = TRUE)
  
  metaDataAnnotations$cellBarcode <- 
    gsub('.*_', replacement = '', metaDataAnnotations$X)
  metaDataAnnotations$condition <- 
    gsub('_.*', replacement = '', metaDataAnnotations$X)
  
  
  # I. Read raw data & process & calculate all values: HDS, Pseudotime, etc
  {
    # read raw counts 
    counts3moNC <- ReadMtx(mtx =  paste0(pathData, FileNames[3]), 
                           cells = paste0(pathData, FileNames[1]),
                           features = paste0(pathData, FileNames[2]),
                           feature.column = 2) 
    
    Seurat3moNC <- CreateSeuratObject(counts = counts3moNC, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = FileNames[1]),
                                      min.cells = 3)
    
    class(Seurat3moNC)
    
    remove(counts3moNC)
    
    Seurat3moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNC, 
                                                        pattern = "^mt-")
    
    
    Seurat3moNC  <- subset(Seurat3moNC , 
                           subset = 
                             nFeature_RNA > 500 &
                             nFeature_RNA < 8000 &
                             percent.mt < 2)
    
    temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC3mo',]
    
    rownames(temp) <- temp$cellBarcode
    
    temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNC@meta.data)),]
    
    Seurat3moNC@meta.data$cellBarcode <- NA
    
    Seurat3moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- 
      temp$cellBarcode
    
    Seurat3moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- 
      temp$CellCluster
    
    Seurat3moNC@meta.data$deleteNA <- 
      rep("not na", length(Seurat3moNC@meta.data$CellType))
    Seurat3moNC@meta.data$deleteNA[is.na(Seurat3moNC@meta.data$CellType)] <-
      "na"
    
    Seurat3moNC <- subset(Seurat3moNC,
                          subset = deleteNA == 'not na')
    
    
    Seurat3moNC <- subset(Seurat3moNC ,
                          subset = 
                            CellType ==  "PC-Hep" |
                            CellType == "mNASH-Hep1"  |
                            CellType == "mNASH-Hep2" |
                            CellType == "PP-Hep" |
                            CellType == "Int-Hep" )
    
    Seurat3moNC <- SCTransform(Seurat3moNC)
    
    
    ## Calculate HDS
    HDSSeurat3moNC <- DS_calc.func(exprMatrices = 
                                     GetAssayData(Seurat3moNC,
                                                  assay = 'SCT',
                                                  layer = 'counts'),
                                   DSignature = HDAG)
    
    HDSSeurat3moNCRNA <- DS_calc.func(exprMatrices = 
                                        GetAssayData(Seurat3moNC,
                                                     assay = 'RNA',
                                                     layer = 'counts'),
                                      DSignature = HDAG)
    
    Seurat3moNC@meta.data$condition <- '3m NC'
    Seurat3moNC@meta.data$HDS <- unname(HDSSeurat3moNC)
    Seurat3moNC@meta.data$HDSrawCounts <- unname(HDSSeurat3moNCRNA)
    
    
    # dimension reduction & PT trajectory inference
    Seurat3moNC <- RunPCA(Seurat3moNC)
    Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:20)
    Seurat3moNC <- FindNeighbors(Seurat3moNC, dims = 1:20)
    Seurat3moNC <- FindClusters(Seurat3moNC, 
                                res = 0.3)
    DimPlot(Seurat3moNC, label = TRUE)
    
    scExp3mo <- as.SingleCellExperiment(Seurat3moNC,
                                        assay = 'SCT')
    colLabels(scExp3mo) <- Seurat3moNC@meta.data$SCT_snn_res.0.3
    
    pseudo.mnn3mo  <- TSCAN::quickPseudotime( scExp3mo,
                                              use.dimred = "PCA", 
                                              dist.method = 'mnn')
    plot(pseudo.mnn3mo$mst)
    mnn.pseudo3mo <- 
      averagePseudotime(pseudo.mnn3mo$ordering)
    
    remove(temp)
    gc()
    
    
    # read raw counts 9moNC
    counts9moNC <- ReadMtx(mtx =  paste0(pathData, FileNames[6]), 
                           cells = paste0(pathData, FileNames[4]),
                           features = paste0(pathData, FileNames[5]),
                           feature.column = 2) 
    
    Seurat9moNC <- CreateSeuratObject(counts = counts9moNC, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = FileNames[4]),
                                      min.cells = 3)
    
    Seurat9moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNC, 
                                                        pattern = "^mt-")
    
    Seurat9moNC  <- subset(Seurat9moNC , 
                           subset = 
                             nFeature_RNA > 500 &
                             nFeature_RNA < 8000 &
                             percent.mt < 2)
    
    remove(counts9moNC)
    
    temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC9mo',]
    rownames(temp) <- temp$cellBarcode
    temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNC@meta.data)),]
    
    Seurat9moNC@meta.data$cellBarcode <- NA
    Seurat9moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- 
      temp$cellBarcode
    Seurat9moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- 
      temp$CellCluster
    Seurat9moNC@meta.data$deleteNA <- 
      rep("not na", length(Seurat9moNC@meta.data$CellType))
    Seurat9moNC@meta.data$deleteNA[is.na(Seurat9moNC@meta.data$CellType)] <-
      "na"
    
    Seurat9moNC <- subset(Seurat9moNC,
                          subset = deleteNA == 'not na')
    
    
    Seurat9moNC <- subset(Seurat9moNC ,
                          subset = 
                            CellType ==  "PC-Hep" |
                            CellType == "mNASH-Hep1"  |
                            CellType == "mNASH-Hep2" |
                            CellType == "PP-Hep" |
                            CellType == "Int-Hep" )
    
    Seurat9moNC <- SCTransform(Seurat9moNC)
    
    # calculate HDS:
    
    HDSSeurat9moNC <- DS_calc.func(exprMatrices = 
                                     GetAssayData(Seurat9moNC,
                                                  assay = 'SCT',
                                                  layer = 'counts'),
                                   DSignature = HDAG)
    
    HDSSeurat9moNCRNA <- DS_calc.func(exprMatrices = 
                                        GetAssayData(Seurat9moNC,
                                                     assay = 'RNA',
                                                     layer = 'counts'),
                                      DSignature = HDAG)
    
    Seurat9moNC@meta.data$condition <- '9m NC'
    Seurat9moNC@meta.data$HDS <- unname(HDSSeurat9moNC)
    Seurat9moNC@meta.data$HDSrawCounts <- unname(HDSSeurat9moNCRNA)
    
    
    
    # dimension reduction & PT trajectory inference
    Seurat9moNC <- RunPCA(Seurat9moNC)
    Seurat9moNC <- RunUMAP(Seurat9moNC, dims = 1:20)
    Seurat9moNC <- FindNeighbors(Seurat9moNC, dims = 1:20)
    Seurat9moNC <- FindClusters(Seurat9moNC,
                                resolution = 0.2)
    DimPlot(Seurat9moNC, label = TRUE)
    
    scExp9mo <- as.SingleCellExperiment(Seurat9moNC,
                                        assay = 'SCT')
    colLabels(scExp9mo) <- Seurat9moNC@meta.data$SCT_snn_res.0.2
    
    pseudo.mnn9mo  <- TSCAN::quickPseudotime( scExp9mo,
                                              use.dimred = "PCA", 
                                              dist.method = 'mnn')
    plot(pseudo.mnn9mo$mst)
    mnn.pseudo9mo <- 
      averagePseudotime(pseudo.mnn9mo$ordering)
    
    remove(temp)
    
    # 3 month NASH
    
    counts3moNASH <- ReadMtx(mtx =  paste0(pathData, FileNames[9]), 
                             cells = paste0(pathData, FileNames[7]),
                             features = paste0(pathData, FileNames[8]),
                             feature.column = 2) 
    
    Seurat3moNASH <- CreateSeuratObject(counts = counts3moNASH, 
                                        project = gsub('_.*', 
                                                       replacement = '', 
                                                       x = FileNames[9]),
                                        min.cells = 3)
    
    Seurat3moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNASH, 
                                                          pattern = "^mt-")
    
    Seurat3moNASH  <- subset(Seurat3moNASH , 
                             subset = 
                               nFeature_RNA > 500 &
                               nFeature_RNA < 8000 &
                               percent.mt < 2)
    
    remove(counts3moNASH)
    
    temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH3mo',]
    rownames(temp) <- temp$cellBarcode
    temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNASH@meta.data)),]
    
    Seurat3moNASH@meta.data$cellBarcode <- NA
    Seurat3moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- 
      temp$cellBarcode
    Seurat3moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- 
      temp$CellCluster
    Seurat3moNASH@meta.data$deleteNA <- 
      rep("not na", length(Seurat3moNASH@meta.data$CellType))
    Seurat3moNASH@meta.data$deleteNA[is.na(Seurat3moNASH@meta.data$CellType)] <-
      "na"
    
    Seurat3moNASH <- subset(Seurat3moNASH,
                            subset = deleteNA == 'not na')
    
    
    Seurat3moNASH <- subset(Seurat3moNASH ,
                            subset = 
                              CellType ==  "PC-Hep" |
                              CellType == "mNASH-Hep1"  |
                              CellType == "mNASH-Hep2" |
                              CellType == "PP-Hep" |
                              CellType == "Int-Hep" )
    Seurat3moNASH <- SCTransform(Seurat3moNASH)
    
    ## Calculate HDS
    HDSSeurat3moNASH <- DS_calc.func(exprMatrices = 
                                       GetAssayData(Seurat3moNASH,
                                                    assay = 'SCT',
                                                    layer = 'counts'),
                                     DSignature = HDAG)
    
    HDSSeurat3moNASHRNA <- DS_calc.func(exprMatrices = 
                                          GetAssayData(Seurat3moNASH,
                                                       assay = 'RNA',
                                                       layer = 'counts'),
                                        DSignature = HDAG)
    
    Seurat3moNASH@meta.data$condition <- '3m NASH'
    Seurat3moNASH@meta.data$HDS <- unname(HDSSeurat3moNASH)
    Seurat3moNASH@meta.data$HDSrawCounts <- unname(HDSSeurat3moNASHRNA)
    
    remove(temp)
    gc()
    
    #dimension reduction & PT trajectory inference
    Seurat3moNASH <- RunPCA(Seurat3moNASH)
    Seurat3moNASH <- RunUMAP(Seurat3moNASH, dims = 1:14)
    Seurat3moNASH <- FindNeighbors(Seurat3moNASH, dims = 1:14)
    # adjusted resolution to make it more comparable to annotation and other metadata
    Seurat3moNASH <- FindClusters(Seurat3moNASH, resolution = 0.1)
    DimPlot(Seurat3moNASH, label = TRUE)
    
    scExp3moNASH <- as.SingleCellExperiment(Seurat3moNASH,
                                            assay = 'SCT')
    colLabels(scExp3moNASH) <- Seurat3moNASH@meta.data$SCT_snn_res.0.1
    
    pseudo.mnn3moNASH  <- TSCAN::quickPseudotime( scExp3moNASH,
                                                  use.dimred = "PCA", 
                                                  dist.method = 'mnn')
    plot(pseudo.mnn3moNASH$mst)
    mnn.pseudo3moNASH <- 
      averagePseudotime(pseudo.mnn3moNASH$ordering)
    
    
    # 9 month NASH
    
    counts9moNASH <- ReadMtx(mtx =  paste0(pathData, FileNames[12]), 
                             cells = paste0(pathData, FileNames[10]),
                             features = paste0(pathData, FileNames[11]),
                             feature.column = 2) 
    
    Seurat9moNASH <- CreateSeuratObject(counts = counts9moNASH, 
                                        project = gsub('_.*', 
                                                       replacement = '', 
                                                       x = FileNames[12]),
                                        min.cells = 3)
    
    Seurat9moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNASH, 
                                                          pattern = "^mt-")
    
    Seurat9moNASH  <- subset(Seurat9moNASH , 
                             subset = 
                               nFeature_RNA > 500 &
                               nFeature_RNA < 8000 &
                               percent.mt < 2)
    
    temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH9mo',]
    rownames(temp) <- temp$cellBarcode
    temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNASH@meta.data)),]
    
    Seurat9moNASH@meta.data$cellBarcode <- NA
    Seurat9moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
      temp$cellBarcode
    Seurat9moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
      temp$CellCluster
    Seurat9moNASH@meta.data$deleteNA <- 
      rep("not na", length(Seurat9moNASH@meta.data$CellType))
    Seurat9moNASH@meta.data$deleteNA[is.na(Seurat9moNASH@meta.data$CellType)] <-
      "na"
    
    Seurat9moNASH <- subset(Seurat9moNASH,
                            subset = deleteNA == 'not na')
    
    
    Seurat9moNASH <- subset(Seurat9moNASH ,
                            subset = 
                              CellType ==  "PC-Hep" |
                              CellType == "mNASH-Hep1"  |
                              CellType == "mNASH-Hep2" |
                              CellType == "PP-Hep" |
                              CellType == "Int-Hep" )
    
    remove(counts9moNASH, temp)
    
    Seurat9moNASH <- SCTransform(Seurat9moNASH)
    
    ## calculate HDS
    HDSSeurat9moNASH <- DS_calc.func(exprMatrices = 
                                       GetAssayData(Seurat9moNASH,
                                                    assay = 'SCT',
                                                    layer = 'counts'),
                                     DSignature = HDAG)
    
    HDSSeurat9moNASHRNA <- DS_calc.func(exprMatrices = 
                                          GetAssayData(Seurat9moNASH,
                                                       assay = 'RNA',
                                                       layer = 'counts'),
                                        DSignature = HDAG)
    
    
    Seurat9moNASH@meta.data$condition <- '9m NASH'
    Seurat9moNASH@meta.data$HDS <- unname(HDSSeurat9moNASH)
    Seurat9moNASH@meta.data$HDSrawCounts <- unname(HDSSeurat9moNASHRNA)
    
    # dimension reduction & PT trajectory inference
    Seurat9moNASH <- RunPCA(Seurat9moNASH)
    Seurat9moNASH <- RunUMAP(Seurat9moNASH, dims = 1:14)
    Seurat9moNASH <- FindNeighbors(Seurat9moNASH, dims = 1:14)
    Seurat9moNASH <- FindClusters(Seurat9moNASH, resolution = 0.3)
    DimPlot(Seurat9moNASH, label = TRUE)
    
    scExp9moNASH <- as.SingleCellExperiment(Seurat9moNASH,
                                            assay = 'SCT')
    colLabels(scExp9moNASH) <- Seurat9moNASH@meta.data$SCT_snn_res.0.3
    
    pseudo.mnn9moNASH  <- TSCAN::quickPseudotime( scExp9moNASH,
                                                  use.dimred = "PCA", 
                                                  dist.method = 'mnn')
    plot(pseudo.mnn9moNASH$mst)
    mnn.pseudo9moNASH <- 
      averagePseudotime(pseudo.mnn9moNASH$ordering)
    
    PTresultsList <- list(PTvalues = list('3moNC' = mnn.pseudo3mo, 
                                          '9moNC'= mnn.pseudo9mo, 
                                          '3moNASH'= mnn.pseudo3moNASH, 
                                          '9moNASH' = mnn.pseudo9moNASH),
                          PTordering = list('3moNC' =pseudo.mnn3mo,
                                            '9moNC'= pseudo.mnn9mo,
                                            '3moNASH'= pseudo.mnn3moNASH,
                                            '9moNASH'= pseudo.mnn9moNASH),
                          scExpObjects = list('3moNC' = scExp3mo,
                                              '9moNC'= scExp9mo,
                                              '3moNASH'= scExp3moNASH,
                                              '9moNASH'= scExp9moNASH))
    
   # saveRDS(PTresultsList, 
    #        file = paste0(outputPath, 'Results/Xiao2023_PseudoTimeTSCAN.rds'))
    
    
    }
  
  # II. Store all values in dataframe for easy plotting with ggplot
  {
    ## bring all HDS from all samples to one dataframe
    # sanity check: meta data matches HDS 
    identical(names(unlist(HDSSeurat3moNC)),
              rownames(Seurat3moNC@meta.data))
    identical(names(unlist(HDSSeurat9moNC)),
              rownames(Seurat9moNC@meta.data))
    identical(names(unlist(HDSSeurat3moNASH)),
              rownames(Seurat3moNASH@meta.data))
    identical(names(unlist(HDSSeurat9moNASH)),
              rownames(Seurat9moNASH@meta.data))
    # all TRUE: good to go
    DataFrameHDS <- data.frame('HDS' = unname(HDSSeurat3moNC),
                               'CellID' = names(unlist(HDSSeurat3moNC)))
    
    DataFrameHDS$condition <- rep('NC', length(DataFrameHDS$HDS))
    DataFrameHDS$age <- rep('3 months', length(DataFrameHDS$HDS))
    DataFrameHDS$sample <- rep('NC 3 months', length(DataFrameHDS$HDS))
    identical(DataFrameHDS$CellID, rownames(Seurat3moNC@meta.data))
    DataFrameHDS$HepatocyteType <- Seurat3moNC@meta.data$CellType
    
    
    DataFrameHDS <- rbind.data.frame(DataFrameHDS,
                                     data.frame('HDS' = unname(HDSSeurat9moNC),
                                                'CellID' = names(unlist(HDSSeurat9moNC)),
                                                'condition' = rep('NC', 
                                                                  length(HDSSeurat9moNC)),
                                                'age' = rep('9 months', 
                                                            length(HDSSeurat9moNC)),
                                                'sample' = rep('NC 9 months', 
                                                               length(HDSSeurat9moNC)),
                                                'HepatocyteType' = Seurat9moNC@meta.data$CellType),
                                     
                                     data.frame('HDS' = unname(HDSSeurat3moNASH),
                                                'CellID' = names(unlist(HDSSeurat3moNASH)),
                                                'condition' = rep('NASH', 
                                                                  length(HDSSeurat3moNASH)),
                                                'age' = rep('3 months', 
                                                            length(HDSSeurat3moNASH)),
                                                'sample' = rep('NASH 3 months', 
                                                               length(HDSSeurat3moNASH)),
                                                'HepatocyteType' = Seurat3moNASH@meta.data$CellType),
                                     
                                     data.frame('HDS' = unname(HDSSeurat9moNASH),
                                                'CellID' = names(unlist(HDSSeurat9moNASH)),
                                                'condition' = rep('NASH', 
                                                                  length(HDSSeurat9moNASH)),
                                                'age' = rep('9 months', 
                                                            length(HDSSeurat9moNASH)),
                                                'sample' = rep('NASH 9 months', 
                                                               length(HDSSeurat9moNASH)),
                                                'HepatocyteType' = Seurat9moNASH@meta.data$CellType)
                                     
                                     
    )
    
    Before <- DataFrameHDS
    
    DataFrameHDS$condition <- factor(DataFrameHDS$condition,
                                     ordered = TRUE,
                                     levels = c("NC","NASH" ))
    
    DataFrameHDS$age <- factor(DataFrameHDS$age,
                               ordered = TRUE, 
                               levels = c("3 months",
                                          "9 months"))
    
    DataFrameHDS$sample <- factor(DataFrameHDS$sample,
                                  ordered = TRUE, 
                                  levels = c("NC 3 months",
                                             "NC 9 months",
                                             "NASH 3 months",
                                             "NASH 9 months"))
    
    
    DataFrameHDS$HepatocyteType <- factor(DataFrameHDS$HepatocyteType,
                                          ordered = TRUE,
                                          levels = c("PP-Hep",
                                                     "Int-Hep",
                                                     "PC-Hep",
                                                     "mNASH-Hep1",
                                                     "mNASH-Hep2"))
    
    
    remove(Before)
    #saveRDS(DataFrameHDS, file = paste0(outputPath, 'Xiao2023HDSsnRNAseqAllHep.rds'))
    
  }
  
  ## III. Plot UMAPs with CellType annotation, Pseudotime and HDS
  {
    # 1. Plot Pseudotime Inferece Results
    
    UmapPT3moNC <- plotUMAP(scExp3mo,
                            colour_by=I(mnn.pseudo3mo),
                            point_size = 1,
                            point_alpha = 0.7) +
      scale_color_viridis(option = "B", 
                          direction = -1
                          # limits = c(min(I(mnn.pseudo3mo)),
                          #           max(I(mnn.pseudo3mo)))
      ) +
      geom_line(data = pseudo.mnn3mo$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, group=edge)) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "PT") 
    
    UmapPT9moNC <- plotUMAP(scExp9mo,
                            colour_by=I(mnn.pseudo9mo),
                            point_size = 1,
                            point_alpha = 0.7) +
      scale_color_viridis(option = "B", 
                          direction = -1
                          # limits = c(min(I(mnn.pseudo9mo)),
                          #             max(I(mnn.pseudo9mo)))
      ) +
      geom_line(data = pseudo.mnn9mo$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, group=edge)) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "PT") 
    
    UmapPT3moNASH <- plotUMAP(scExp3moNASH,
                              colour_by=I(mnn.pseudo3moNASH),
                              point_size = 1,
                              point_alpha = 0.7) +
      scale_color_viridis(option = "B", 
                          direction = -1
                          # limits = c(min(I(mnn.pseudo3moNASH)),
                          #           max(I(mnn.pseudo3moNASH)))
      ) +
      geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, group=edge)) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "PT")
    
    UmapPT9moNASH <- plotUMAP(scExp9moNASH,
                              colour_by=I(mnn.pseudo9moNASH),
                              point_size = 1,
                              point_alpha = 0.7) +
      scale_color_viridis(option = "B", 
                          direction = -1
                          #  limits = c(min(I(mnn.pseudo9moNASH)),
                          #             max(I(mnn.pseudo9moNASH)))
      ) +
      geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, group=edge)) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "PT") 
    
    ## Plot HDS: 
    # color palette will use min and max values HDS from each sample as 
    # reference points --> each UMAP needs it's own HDS legend. 
    # this is necessary cause the color resolution is needed to compare to pseudo
    # time trajectory 
    
    # Find global color boarders: 1%-5% top and bottom 
    summary(DataFrameHDS)
    length(DataFrameHDS$sample)
    
    limitBottom1percent <- DataFrameHDS[order(DataFrameHDS$HDS),'HDS'][round((15825)*0.05)]
    limitTop1percent <- DataFrameHDS[order(DataFrameHDS$HDS, decreasing = TRUE),'HDS'][round((15825)*0.01)]
    
    UmapHDS3mo <-plotUMAP(scExp3mo,
                          colour_by ='HDS',
                          point_size = 0.5,
                          point_alpha = 1) +
      geom_line(data = pseudo.mnn3mo$connected$UMAP, 
                mapping=aes(x=umap_1, 
                            y=umap_2, 
                            group=edge)) +
      # scale_color_viridis(option = "turbo", 
      #                     direction = 1) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1,
                            limits = c(limitBottom1percent,
                                       limitTop1percent),
                            oob = squish
      ) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "HDS") 
    
    UmapHDS9mo <-plotUMAP(scExp9mo,
                          colour_by ='HDS',
                          point_size = 0.5,
                          point_alpha = 1) +
      geom_line(data = pseudo.mnn9mo$connected$UMAP, 
                mapping=aes(x=umap_1, 
                            y=umap_2, 
                            group=edge)) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1,
                            limits = c(limitBottom1percent,
                                       limitTop1percent),
                            oob = squish) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "HDS") 
    
    UmapHDS3moNASH <-plotUMAP(scExp3moNASH,
                              colour_by ='HDS',
                              point_size = 0.5,
                              point_alpha = 1
    ) +
      geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, 
                            y=umap_2, 
                            group=edge)) +
      # scale_color_viridis(option = "turbo", direction = 1,
      #                     limits = c(limitBottom1percent,
      #                                limitTop1percent)) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1,
                            limits = c(limitBottom1percent,
                                       limitTop1percent),
                            oob = squish) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "HDS") 
    
    UmapHDS9moNASH <-plotUMAP(scExp9moNASH,
                              colour_by ='HDS',
                              point_size = 0.5,
                              point_alpha = 1) +
      geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, 
                            y=umap_2, 
                            group=edge)) +
      # scale_color_viridis(option = "turbo", direction = 1,
      #                     limits = c(limitBottom1percent,
      #                                limitTop1percent),
      #                     oob = squish) +
      scale_color_distiller(palette = "YlOrBr",
                            direction = 1,
                            limits = c(limitBottom1percent,
                                       limitTop1percent),
                            oob = squish) +
      theme(legend.position = 'top',
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16)) +
      labs(color = "HDS") 
    
    ##### Plot Hepatocyte Zone Annotation
    
    UMAPct3moNC <- plotUMAP(scExp3mo, 
                            colour_by='CellType',
                            point_size = 0.5,
                            point_alpha = 0.7) +
      geom_line(data = pseudo.mnn3mo$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, 
                            group=edge)) +
      scale_color_colorblind() +
      theme(legend.position = 'none',
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            axis.ticks = element_blank(),
            axis.text = element_blank()
      ) 
    
    UMAPct9moNC <- plotUMAP(scExp9mo, 
                            colour_by='CellType',
                            point_size = 0.5,
                            point_alpha = 0.7) +
      geom_line(data = pseudo.mnn9mo$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, 
                            group=edge)) +
      scale_color_colorblind() +
      theme(legend.position = 'none',
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            axis.ticks = element_blank(),
            axis.text = element_blank()
      ) 
    
    
    
    UMAPct3moNASH <- plotUMAP(scExp3moNASH, 
                              colour_by='CellType',
                              point_size = 0.5,
                              point_alpha = 0.7) +
      geom_line(data = pseudo.mnn3moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, 
                            group=edge)) +
      scale_color_colorblind() +
      theme(legend.position = 'none',
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            axis.ticks = element_blank(),
            axis.text = element_blank()
      ) 
    
    
    UMAPct9moNASH <- plotUMAP(scExp9moNASH, 
                              colour_by='CellType',
                              point_size = 0.5,
                              point_alpha = 0.7) +
      geom_line(data = pseudo.mnn9moNASH$connected$UMAP, 
                mapping=aes(x=umap_1, y=umap_2, 
                            group=edge)) +
      scale_color_colorblind() +
      theme(legend.position = 'none',
            legend.key.size = unit(1.5, 'cm'),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            axis.ticks = element_blank(),
            axis.text = element_blank()
      ) 
    
    patchWorkCTs <- (UMAPct3moNC | UMAPct9moNC | UMAPct3moNASH | UMAPct9moNASH )
    patchWorkHDS <- (UmapHDS3mo | UmapHDS9mo | UmapHDS3moNASH | UmapHDS9moNASH)
    patchWorkPT <- (UmapPT3moNC | UmapPT9moNC | UmapPT3moNASH | UmapPT9moNASH)
    
    allPatch <- (patchWorkCTs / patchWorkHDS / patchWorkPT ) + theme(legend.position = 'none')
  
    
    # For the manuscript:
    
    Umaps3moNC <- (UMAPct3moNC | UmapHDS3mo | UmapPT3moNC)
    Umaps9moNC <- (UMAPct9moNC | UmapHDS9mo | UmapPT9moNC)
    Umaps3moNASH <- (UMAPct3moNASH | UmapHDS3moNASH | UmapPT3moNASH)
    Umaps9moNASH <- (UMAPct9moNASH | UmapHDS9moNASH | UmapPT9moNASH)

    
  }
  
}

# Supplementary Figure 6 Panels: D
{
  # Pseudotime trajectory inference on Xiao et al. NASH mouse hepatocytes 
  # 
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
  library(tidyverse)
  library(rtracklayer)
  library(Matrix)
  library(scater)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(scales)
  
  pathData <- "Data/mouse/Xiao_snRNAseq_GSE189600/"
  outputPath <- "Data/mouse/output/trajectory_inference/"
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "Data/Output/HDAG.csv",
    sep = ',')
  # load HDS function
  source('SharedFunctions.R')
  
  # we need only the files from the 9 month NASH mouse
  file_names <- list.files(path = pathData)
  
  
  metaDataAnnotations <- read.table("Data/mouse/annotation_files/Xiao_snRNAseq_GSE189600_scitranslmed.adc9653_data_file_s1.csv",
                                    sep = ';', 
                                    dec = ',', 
                                    header = TRUE)
  
  metaDataAnnotations$cellBarcode <- 
    gsub('.*_', replacement = '', metaDataAnnotations$X)
  metaDataAnnotations$condition <- 
    gsub('_.*', replacement = '', metaDataAnnotations$X)
  
  ###
  # 3-month old mouse fed Control Diet 
  ###
  
  counts3moNC <- ReadMtx(mtx =  paste0(pathData, file_names[3]), 
                         cells = paste0(pathData, file_names[1]),
                         features = paste0(pathData, file_names[2]),
                         feature.column = 2) 
  
  Seurat3moNC <- CreateSeuratObject(counts = counts3moNC, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = file_names[1]),
                                    min.cells = 3)
  
  remove(counts3moNC)
  
  Seurat3moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNC, pattern = "^mt-")
  
  
  Seurat3moNC <- subset(Seurat3moNC, subset = 
                          nFeature_RNA > 500 &
                          nFeature_RNA < 8000 &
                          percent.mt < 2 )
  
  temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC3mo',]
  rownames(temp) <- temp$cellBarcode
  temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNC@meta.data)),]
  
  Seurat3moNC@meta.data$cellBarcode <- NA
  Seurat3moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- temp$cellBarcode
  Seurat3moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNC@meta.data))] <- temp$CellCluster
  Seurat3moNC@meta.data$deleteNA <- rep("not na", length(Seurat3moNC@meta.data$CellType))
  Seurat3moNC@meta.data$deleteNA[is.na(Seurat3moNC@meta.data$CellType)] <-"na"
  Seurat3moNC <- subset(Seurat3moNC,subset = deleteNA == 'not na')
  
  Seurat3moNC <- subset(Seurat3moNC ,
                        subset = CellType ==  "PC-Hep" | CellType == "mNASH-Hep1"  | CellType == "mNASH-Hep2" |
                          CellType == "PP-Hep" | CellType == "Int-Hep" )
  
  Seurat3moNC <- SCTransform(Seurat3moNC)
  Seurat3moNC_HDS <- DS_calc.func(exprMatrices = GetAssayData(Seurat3moNC,assay = "SCT", layer = "counts"), DSignature = HDAG)
  
  identical(names(unlist(Seurat3moNC_HDS)), rownames(Seurat3moNC@meta.data))
  # TRUE, so:
  Seurat3moNC@meta.data$HDS <- unname(Seurat3moNC_HDS)
  
  Seurat3moNC <- RunPCA(Seurat3moNC,feautres = VariableFeatures(Seurat3moNC))
  DimPlot(Seurat3moNC , reduction = "pca", group.by = 'CellType')
  ggsave('XiaoMouse3moNC_PC1andPC2_HepatocyteZones.png', path = outputPath, width = 16, height = 12, unit = "cm")
  DimPlot(Seurat3moNC , reduction = "pca",dims = c(3,4), group.by = 'CellType')                      
  Seurat3moNC <- RunUMAP(Seurat3moNC, dims = 1:14)
  Seurat3moNC  <- FindNeighbors(Seurat3moNC , dims = 1:14)
  Seurat3moNC <- FindClusters(Seurat3moNC, resolution = 0.5)
  
  # checking out the cluster labels
  head(Idents(Seurat3moNC), 5)
  DimPlot(Seurat3moNC, reduction = "umap", group.by = 'CellType')

  # Trajectory inference with Monocle3 
  ## Step 1: Normalize and pre-process the data
  monocle_3moNC <- as.cell_data_set(Seurat3moNC, assay = "RNA")
  monocle_3moNC <- preprocess_cds(monocle_3moNC, num_dim = 100)
  ## Step 2: Reduce the dimensions using UMAP
  monocle_3moNC <- reduce_dimension(monocle_3moNC )
  ## Step 3: Cluster the cells
  monocle_3moNC <- cluster_cells(monocle_3moNC )
  ## Step 5: Learn a graph
  monocle_3moNC <- learn_graph(monocle_3moNC )
  # Plot monocle cluster and graph colored by zone
  plot_zone_3moNC <- plot_cells(monocle_3moNC, label_groups_by_cluster = FALSE, color_cells_by = "CellType", 
                                label_cell_groups = TRUE, cell_size = 1 ) + scale_color_colorblind() + 
    ggtitle("3-month-old Normal Chow") + theme(legend.position = "bottom")
  
  plot_zone_3moNC

  ## Step 6: Order cells
  monocle_3moNC <- order_cells(monocle_3moNC)
  
  plot_pseudotime_3moNC <- plot_cells(monocle_3moNC, color_cells_by = "pseudotime", cell_size = 1 ) +theme(legend.position = "bottom")
  
  plot_HDS_monocle3moNC <- plot_cells(monocle_3moNC, color_cells_by = "HDS", cell_size = 1 )  + 
    scale_color_distiller(palette = "YlOrBr",direction = 1) +
    theme(legend.position = 'top',
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)) +
    labs(color = "HDS") + theme(legend.position = "bottom")
  
  panel3moNC <- (plot_zone_3moNC | plot_pseudotime_3moNC | plot_HDS_monocle3moNC)
  panel3moNC 
  

  
  ### 
  # 9-month NC
  ###
  
  
  counts9moNC <- ReadMtx(mtx =  paste0(pathData, file_names[6]), 
                         cells = paste0(pathData, file_names[4]),
                         features = paste0(pathData, file_names[5]),
                         feature.column = 2) 
  
  Seurat9moNC <- CreateSeuratObject(counts = counts9moNC, 
                                    project = gsub('_.*', 
                                                   replacement = '', 
                                                   x = file_names[4]),
                                    min.cells = 3)
  
  remove(counts9moNC)
  
  Seurat9moNC[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNC, pattern = "^mt-")
  
  
  Seurat9moNC <- subset(Seurat9moNC, subset = 
                          nFeature_RNA > 500 &
                          nFeature_RNA < 8000 &
                          percent.mt < 2 )
  
  temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NC9mo',]
  rownames(temp) <- temp$cellBarcode
  temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNC@meta.data)),]
  
  Seurat9moNC@meta.data$cellBarcode <- NA
  Seurat9moNC@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- temp$cellBarcode
  Seurat9moNC@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNC@meta.data))] <- temp$CellCluster
  Seurat9moNC@meta.data$deleteNA <- rep("not na", length(Seurat9moNC@meta.data$CellType))
  Seurat9moNC@meta.data$deleteNA[is.na(Seurat9moNC@meta.data$CellType)] <-"na"
  Seurat9moNC <- subset(Seurat9moNC,subset = deleteNA == 'not na')
  
  Seurat9moNC <- subset(Seurat9moNC ,
                        subset = CellType ==  "PC-Hep" | CellType == "mNASH-Hep1"  | CellType == "mNASH-Hep2" |
                          CellType == "PP-Hep" | CellType == "Int-Hep" )
  
  Seurat9moNC <- SCTransform(Seurat9moNC)
  Seurat9moNC_HDS <- DS_calc.func(exprMatrices = GetAssayData(Seurat9moNC,assay = "SCT", layer = "counts"), DSignature = HDAG)
  
  identical(names(unlist(Seurat9moNC_HDS)), rownames(Seurat9moNC@meta.data))
  # TRUE, so:
  Seurat9moNC@meta.data$HDS <- unname(Seurat9moNC_HDS )
  
  Seurat9moNC <- RunPCA(Seurat9moNC,feautres = VariableFeatures(Seurat9moNC))
  DimPlot(Seurat9moNC , reduction = "pca", group.by = 'CellType')
  ggsave('XiaoMouse9moNC_PC1andPC2_HepatocyteZones.png', path = outputPath, width = 16, height = 12, unit = "cm")
  DimPlot(Seurat9moNC , reduction = "pca",dims = c(3,4), group.by = 'CellType')                      
  Seurat9moNC <- RunUMAP(Seurat9moNC, dims = 1:14)
  Seurat9moNC  <- FindNeighbors(Seurat9moNC , dims = 1:14)
  Seurat9moNC <- FindClusters(Seurat9moNC, resolution = 0.5)
  
  # checking out the cluster labels
  head(Idents(Seurat9moNC), 5)
  DimPlot(Seurat9moNC, reduction = "umap", group.by = 'CellType')

  
  # Trajectory inference with Monocle3 
  ## Step 1: Normalize and pre-process the data
  monocle_9moNC <- as.cell_data_set(Seurat9moNC, assay = "RNA")
  monocle_9moNC <- preprocess_cds(monocle_9moNC, num_dim = 100)
  ## Step 2: Reduce the dimensions using UMAP
  monocle_9moNC <- reduce_dimension(monocle_9moNC)
  ## Step 3: Cluster the cells
  monocle_9moNC <- cluster_cells(monocle_9moNC)
  ## Step 5: Learn a graph
  monocle_9moNC <- learn_graph(monocle_9moNC)
  
  # Plot monocle cluster and graph colored by zone
  plot_zone_9moNC <- plot_cells(monocle_9moNC, label_groups_by_cluster = FALSE, color_cells_by = "CellType", 
                                label_cell_groups = TRUE, cell_size = 1 ) + scale_color_colorblind() + 
    ggtitle("9-month-old Normal Chow") + theme(legend.position = "bottom")
  
  plot_zone_9moNC 
  
  ## Step 6: Order cells
  monocle_9moNC  <- order_cells(monocle_9moNC)
  
  plot_pseudotime_9moNC <- plot_cells(monocle_9moNC, color_cells_by = "pseudotime", cell_size = 1 ) +theme(legend.position = "bottom")
  
  plot_HDS_monocle9moNC <- plot_cells(monocle_9moNC, color_cells_by = "HDS", cell_size = 1 )  + 
    scale_color_distiller(palette = "YlOrBr",direction = 1) +
    theme(legend.position = 'top',
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)) +
    labs(color = "HDS") + theme(legend.position = "bottom")
  
  panel9moNC <- (plot_zone_9moNC | plot_pseudotime_9moNC | plot_HDS_monocle9moNC)
  panel9moNC 
  

  
  
  
  ###
  # 3-month NASH
  ###
  
  counts3moNASH <- ReadMtx(mtx =  paste0(pathData, file_names[9]), 
                           cells = paste0(pathData, file_names[7]),
                           features = paste0(pathData, file_names[8]),
                           feature.column = 2) 
  
  Seurat3moNASH <- CreateSeuratObject(counts = counts3moNASH, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = file_names[7]),
                                      min.cells = 3)
  
  remove(counts3moNASH)
  
  Seurat3moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat3moNASH, pattern = "^mt-")
  
  
  Seurat3moNASH <- subset(Seurat3moNASH, subset = 
                            nFeature_RNA > 500 &
                            nFeature_RNA < 8000 &
                            percent.mt < 2 )
  
  temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH3mo',]
  rownames(temp) <- temp$cellBarcode
  temp <- temp[intersect(temp$cellBarcode, rownames(Seurat3moNASH@meta.data)),]
  
  Seurat3moNASH@meta.data$cellBarcode <- NA
  Seurat3moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- temp$cellBarcode
  Seurat3moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat3moNASH@meta.data))] <- temp$CellCluster
  Seurat3moNASH@meta.data$deleteNA <- rep("not na", length(Seurat3moNASH@meta.data$CellType))
  Seurat3moNASH@meta.data$deleteNA[is.na(Seurat3moNASH@meta.data$CellType)] <-"na"
  Seurat3moNASH <- subset(Seurat3moNASH,subset = deleteNA == 'not na')
  
  Seurat3moNASH <- subset(Seurat3moNASH ,
                          subset = CellType ==  "PC-Hep" | CellType == "mNASH-Hep1"  | CellType == "mNASH-Hep2" |
                            CellType == "PP-Hep" | CellType == "Int-Hep" )
  
  Seurat3moNASH <- SCTransform(Seurat3moNASH)
  Seurat3moNASH_HDS <- DS_calc.func(exprMatrices = 
                                      GetAssayData(Seurat3moNASH,assay = "SCT", layer = "counts"), DSignature = HDAG)
  
  identical(names(unlist(Seurat3moNASH_HDS)), rownames(Seurat3moNASH@meta.data))
  # TRUE, so:
  Seurat3moNASH@meta.data$HDS <- unname(Seurat3moNASH_HDS)
  
  Seurat3moNASH <- RunPCA(Seurat3moNASH,feautres = VariableFeatures(Seurat3moNASH))
  DimPlot(Seurat3moNASH, reduction = "pca", group.by = 'CellType')
  DimPlot(Seurat3moNASH , reduction = "pca",dims = c(3,4), group.by = 'CellType')                      
  Seurat3moNASH <- RunUMAP(Seurat3moNASH, dims = 1:14)
  Seurat3moNASH <- FindNeighbors(Seurat3moNASH, dims = 1:14)
  Seurat3moNASH <- FindClusters(Seurat3moNASH, resolution = 0.3)
  
  # checking out the cluster labels
  head(Idents(Seurat3moNASH), 5)
  DimPlot(Seurat3moNASH, reduction = "umap", group.by = 'CellType')

  
  # Trajectory inference with Monocle3 
  ## Step 1: Normalize and pre-process the data
  monocle_3moNASH <- as.cell_data_set(Seurat3moNASH, assay = "RNA")
  monocle_3moNASH  <- preprocess_cds(monocle_3moNASH, num_dim = 100)
  ## Step 2: Reduce the dimensions using UMAP
  monocle_3moNASH <- reduce_dimension(monocle_3moNASH)
  ## Step 3: Cluster the cells
  monocle_3moNASH <- cluster_cells(monocle_3moNASH)
  ## Step 5: Learn a graph
  monocle_3moNASH <- learn_graph(monocle_3moNASH)
  
  # Plot monocle cluster and graph colored by zone
  plot_zone_3moNASH <- plot_cells(monocle_3moNASH, label_groups_by_cluster = FALSE, color_cells_by = "CellType", 
                                  label_cell_groups = TRUE, cell_size = 1 ) + scale_color_colorblind() + 
    ggtitle("3-month-old NASH") + theme(legend.position = "bottom")
  
  plot_zone_3moNASH

  ## Step 6: Order cells
  monocle_3moNASH  <- order_cells(monocle_3moNASH)
  
  plot_pseudotime_3moNASH <- plot_cells(monocle_3moNASH, color_cells_by = "pseudotime", cell_size = 1 ) +theme(legend.position = "bottom")
  
  plot_HDS_monocle3moNASH <- plot_cells(monocle_3moNASH, color_cells_by = "HDS", cell_size = 1 )  + 
    scale_color_distiller(palette = "YlOrBr",direction = 1) +
    theme(legend.position = 'top',
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)) +
    labs(color = "HDS") + theme(legend.position = "bottom")
  
  panel3moNASH <- (plot_zone_3moNASH| plot_pseudotime_3moNASH | plot_HDS_monocle3moNASH)
  panel3moNASH
  
  ####
  #  9-month NASH
  ####
  counts9moNASH <- ReadMtx(mtx =  paste0(pathData, file_names[12]), 
                           cells = paste0(pathData, file_names[10]),
                           features = paste0(pathData, file_names[11]),
                           feature.column = 2) 
  
  Seurat9moNASH <- CreateSeuratObject(counts = counts9moNASH, 
                                      project = gsub('_.*', 
                                                     replacement = '', 
                                                     x = file_names[10]),
                                      min.cells = 3)
  
  remove(counts9moNASH)
  
  Seurat9moNASH[["percent.mt"]] <- PercentageFeatureSet(Seurat9moNASH , 
                                                        pattern = "^mt-")
  
  
  Seurat9moNASH  <- subset(Seurat9moNASH , 
                           subset = 
                             nFeature_RNA > 500 &
                             nFeature_RNA < 8000 &
                             percent.mt < 2 )
  
  temp <- metaDataAnnotations[metaDataAnnotations$condition == 'NASH9mo',]
  rownames(temp) <- temp$cellBarcode
  temp <- temp[intersect(temp$cellBarcode, rownames(Seurat9moNASH@meta.data)),]
  
  Seurat9moNASH@meta.data$cellBarcode <- NA
  Seurat9moNASH@meta.data$cellBarcode[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
    temp$cellBarcode
  Seurat9moNASH@meta.data$CellType[match(temp$cellBarcode, rownames(Seurat9moNASH@meta.data))] <- 
    temp$CellCluster
  Seurat9moNASH@meta.data$deleteNA <- 
    rep("not na", length(Seurat9moNASH@meta.data$CellType))
  Seurat9moNASH@meta.data$deleteNA[is.na(Seurat9moNASH@meta.data$CellType)] <-
    "na"
  
  Seurat9moNASH <- subset(Seurat9moNASH,
                          subset = deleteNA == 'not na')
  
  
  Seurat9moNASH <- subset(Seurat9moNASH ,
                          subset = 
                            CellType ==  "PC-Hep" |
                            CellType == "mNASH-Hep1"  |
                            CellType == "mNASH-Hep2" |
                            CellType == "PP-Hep" |
                            CellType == "Int-Hep" )
  
  Seurat9moNASH<- SCTransform(Seurat9moNASH)
  
  Seurat9moNASH_HDS <- DS_calc.func(exprMatrices = GetAssayData(Seurat9moNASH,
                                                                assay = "SCT", layer = "counts"), DSignature = HDAG)
  
  identical(names(Seurat9moNASH_HDS), rownames(Seurat9moNASH@meta.data))
  # TRUE, so:
  
  Seurat9moNASH@meta.data$HDS <- unname(Seurat9moNASH_HDS)
  
  
  # Normalization & Dimension reduction & Clustering for Psuedotime 
  Seurat9moNASH <- RunPCA(Seurat9moNASH,
                          feautres =
                            VariableFeatures(Seurat9moNASH))
  
  DimPlot(Seurat9moNASH , reduction = "pca", group.by = 'CellType')
  DimPlot(Seurat9moNASH , reduction = "pca",dims = c(3,4), group.by = 'CellType')
  Seurat9moNASH <- RunUMAP(Seurat9moNASH, dims = 1:14)
  Seurat9moNASH  <- FindNeighbors(Seurat9moNASH , dims = 1:14)
  Seurat9moNASH <- FindClusters(Seurat9moNASH, 
                                resolution = 0.3)
  
  # checking out the cluster labels
  head(Idents(Seurat9moNASH), 5)
  
  DimPlot(Seurat9moNASH, reduction = "umap", group.by = 'CellType')
  
  ggsave('XiaoMouse9moNASH_Umap_colored_hepatocyte_zones.png', path = outputPath, width = 16, height = 12, unit = "cm")
  
  DimPlot(Seurat9moNASH, reduction = "umap", group.by = 'seurat_clusters')
  
  # Monocle3
  ## Step 1: Normalize and pre-process the data
  monocle_9moNASH <- as.cell_data_set(Seurat9moNASH, assay = "RNA")
  monocle_9moNASH  <- preprocess_cds(monocle_9moNASH, num_dim = 100)
  ## Step 2: Reduce the dimensions using UMAP
  monocle_9moNASH <- reduce_dimension(monocle_9moNASH)
  ## Step 3: Cluster the cells
  monocle_9moNASH <- cluster_cells(monocle_9moNASH)
  ## Step 5: Learn a graph
  monocle_9moNASH  <- learn_graph(monocle_9moNASH )
  
  # Plot monocle cluster and graph colored by zone
  plot_zone_9moNASH <- plot_cells(monocle_9moNASH, label_groups_by_cluster = FALSE, color_cells_by = "CellType", 
                                  label_cell_groups = TRUE, cell_size = 1 ) + scale_color_colorblind() + 
    ggtitle("9-month-old NASH") + theme(legend.position = "bottom")
  
  plot_zone_9moNASH
  ggsave('XiaoMouse9moNASH_Umap_Monocle3_colored_hepatocyte_zones.png', path = outputPath, width = 16, height = 12, unit = "cm")
  
  ## Step 6: Order cells
  monocle_9moNASH  <- order_cells(monocle_9moNASH)
  
  plot_pseudotime_9moNASH <- plot_cells(monocle_9moNASH, color_cells_by = "pseudotime", cell_size = 1 ) + 
    theme(legend.position = "bottom")
  
  plot_HDS_monocle9moNASH <- plot_cells(monocle_9moNASH, color_cells_by = "HDS", cell_size = 1 )  + 
    scale_color_distiller(palette = "YlOrBr",direction = 1) +
    theme(legend.position = 'top',
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)) +
    labs(color = "HDS") + theme(legend.position = "bottom")
  
  panel9moNASH <- (plot_zone_9moNASH| plot_pseudotime_9moNASH | plot_HDS_monocle9moNASH)
  panel9moNASH
  

  
}


