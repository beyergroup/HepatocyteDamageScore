# Cell state transition in damaged hepatocytes 

# Panel A, B & C:
#  Analysis of Carlesi snRNA-seq data set 
# - need output from script: ..... that creates the following object:
# seurat_object_SCT_HDS_and_Senescence_scores.rds
# - Cell state assignment
# - IQR and relationship between HDS 
# - Heteroskedasticity test
# FindMarkers() senescent-damaged vs. not-senescent-damaged-cells

# PANEL A:
{
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(patchwork)
  library(scales)
  library(lmtest)
  
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
  
  
  
  
  umap_conditions <- DimPlot(seurat_object, group.by = "Group",
                             pt.size = 1 , cells = sample(rownames(seurat_object@meta.data))) +
    scale_colour_colorblind()+
    theme(legend.position = "bottom",
          legend.key.height = unit(0.5,"cm"),
          legend.key.width = unit(1, "cm"),
          text = element_text(size = 16),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  
  # separate zones
  # Zone_1 -> Periportal hepatocytes
  # Zone_2 -> intermediate
  # Zone_3 ->  Pericentral 
  seurat_object@meta.data <- seurat_object@meta.data %>% 
    mutate(zone = case_when(
      cell_type == "Zone_1_Hep" ~ "Periportal",
      cell_type == "Zone_2_Hep" ~ "Intermediate",
      cell_type == "Zone_3_Hep" ~ "Pericentral",
      cell_type == "daHep" ~ "damaged Hepatocytes" 
    ))
  
  seurat_object@meta.data$zone <- factor(
    seurat_object@meta.data$zone,
    levels = c("Periportal", "Intermediate","Pericentral", "damaged Hepatocytes"), 
    ordered = TRUE)
  
  
  # SHGS - Senescencent Hepatocytes Gene Signature 
  
  cells_rankings <- AUCell_buildRankings(counts)
  shgs_score <- AUCell::AUCell_calcAUC(shgs_genes_mouse, cells_rankings)
  
  
  if(identical(colnames(shgs_score@assays@data$AUC), rownames(seurat_object@meta.data)) == TRUE) {
    print("Matching IDs!")
    seurat_object@meta.data$aucell_shgs <- c(shgs_score@assays@data$AUC)
  }
  
  # Compute thresholds for cell-fate assignment 
  quantiles_hds <- quantile(seurat_object@meta.data$HDS, probs = c(0.25, 0.75))
  quantiles_shgs <- quantile(seurat_object@meta.data$aucell_shgs[seurat_object@meta.data$Group == "TAA"], probs = c(0.25, 0.75))
  
  #> quantiles_hds
  #25%         75% 
  #  -0.15723161 -0.04037954 
  #  > quantiles_shgs
  #         25%        75% 
  #  0.01360055 0.03118557 
  
  # Add or overwrite the 'group' column based on the conditions
  seurat_object@meta.data <-  seurat_object@meta.data %>%
    mutate(cellfate = case_when(
      HDS <= quantiles_hds[1]  ~ "undamaged",
      HDS > quantiles_hds[1] & HDS <= quantiles_hds[2] ~ "transition",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[2]  ~ "damaged_senescent",
      HDS > quantiles_hds[2] & aucell_shgs <= quantiles_shgs[1] ~ "damaged_not_senescent",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[1] & aucell_shgs <= quantiles_shgs[2] ~ "damaged_unresolved",
      TRUE ~ NA_character_  # optional: assign NA if no condition is met
    ))
  
  # Control number of cells assigned per cell fate
  table( seurat_object@meta.data$cellfate)
  # Control if any cells have not been assigned a cell fate 
  sum(is.na( seurat_object@meta.data$cellfate))
  
  # Cell fates as ordered factors 
  seurat_object@meta.data$cellfate <- factor( seurat_object@meta.data$cellfate, 
                                              ordered = TRUE,
                                              levels = c('undamaged', 
                                                         'transition',
                                                         'damaged_not_senescent',
                                                         'damaged_unresolved',
                                                         'damaged_senescent'),
                                              labels = c('undamaged',
                                                         'transition',
                                                         'damaged-not-senescent',
                                                         'damaged-unresolved',
                                                         'damaged-senescent'))
  
  
  # Summarize the data: count the number of cells for each fate in each sample
  summarized_cellfates <-  seurat_object@meta.data %>%
    group_by(Group, cellfate) %>%
    summarize(Cell_Count = n(), .groups = "drop") %>%
    group_by(Group) %>% 
    mutate(proportion = Cell_Count / sum(Cell_Count))
  
  
  # Prepare a table-like data frame for annotations
  label_data <- summarized_cellfates %>% 
    group_by(Group) %>%
    mutate(n = sum(Cell_Count)) %>%
    dplyr::select(Group, cellfate, Cell_Count,n) %>%
    pivot_wider(names_from = cellfate, values_from = Cell_Count, values_fill = 0)
  
  # Convert the data frame to a long format for plotting as a table
  label_data_long <- label_data[,-2] %>%
    pivot_longer(
      cols = -Group,
      names_to = "cellfate",
      values_to = "Cell_Count"
    ) %>% 
    mutate(cellfate = factor(cellfate,
                             levels = c('undamaged',
                                        'transition',
                                        'damaged-not-senescent',
                                        'damaged-unresolved',
                                        'damaged-senescent')))
  
  
  
  ### PLOTS 
  
  # 1. UMAPS: conditions, HDS, SHGS, and cellfate assignment
  
  # floor and ceiling values for applying limits to the color scale
  # Squish values outside the limit to the nearest boundary
  floor_valueHDS <- quantile( seurat_object@meta.data$HDS, 0.01)  
  ceiling_valueHDS <- quantile( seurat_object@meta.data$HDS, 0.99)
  
  floor_valueSHGS <- quantile( seurat_object@meta.data$aucell_shgs, 0.01)  
  ceiling_valueSHGS <- quantile( seurat_object@meta.data$aucell_shgs, 0.99)
  
  set.seed(42)
  umapHDSmerged <- FeaturePlot( seurat_object, features = "HDS",
                                cells = sample(rownames( seurat_object@meta.data),
                                               replace = FALSE), pt.size = 1) + 
    theme(legend.position = "top", legend.key.size = unit(15, "mm"),
          legend.key.height = unit(5,"mm"),
          text = element_text(size=20),
          axis.title=element_blank(),axis.text=element_blank(),
          axis.ticks=element_blank()) +
    scale_colour_distiller(palette = "YlOrBr", 
                           direction = 1,
                           limits = c(floor_valueHDS, 
                                      ceiling_valueHDS),  
                           oob = scales::squish) + ggtitle("")
  
  umapSHGSmerged <- FeaturePlot( seurat_object,
                                 features = "aucell_shgs",
                                 cells = sample(rownames( seurat_object@meta.data),
                                                replace = FALSE),
                                 pt.size = 1) + 
    scale_colour_distiller(palette = "YlOrBr",
                           direction = 1,
                           limits = c(floor_valueSHGS, 
                                      ceiling_valueSHGS),  
                           oob = scales::squish) +
    theme(legend.position = "top", 
          legend.key.height=unit(5,"mm"),
          legend.key.width = unit(15, "mm"),
          text = element_text(size=20),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) + ggtitle("")
  
  
  densityPlotsHDS <- ggplot(data =  seurat_object@meta.data,
                            aes(x = HDS, color = Group))+
    geom_density(lwd = 1, trim = TRUE)+ theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("HDS")
  
  densityPlotsSHGS <- ggplot(data =  seurat_object@meta.data,
                             aes(x = aucell_shgs, color = Group))+
    geom_density(lwd = 1, trim = TRUE) + theme_classic() +
    scale_colour_colorblind() +
    theme(legend.position = "bottom",
          legend.key.height = unit(5,"mm"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=18)) +
    xlab("SHGS")
  
  grid_umaps <- cowplot::plot_grid(plotlist = list(umapHDSmerged, 
                                                   umapSHGSmerged,
                                                   densityPlotsHDS,
                                                   densityPlotsSHGS), 
                                   rel_heights = c(1,0.35), nrow = 2) 
  
  set.seed(42)
  umap_conditions <- DimPlot(seurat_object, group.by = 'Group', 
                             pt.size = 1, 
                             cells = sample(rownames( seurat_object@meta.data))) + 
    scale_color_colorblind() + 
    theme(legend.position = "bottom", axis.text.y = element_blank(),
          axis.text.x = element_blank(),axis.text = element_blank(), 
          axis.title = element_blank(), axis.ticks = element_blank(), 
          plot.title = element_blank(), legend.key.height=unit(5,"mm"),
          text = element_text(size=12)) 
  
  
}

# PANEL B: 
{
  #######
  # Plot SHGS activity vs HDS in Healthy, CDE and TAA
  ########
  
  data_to_plot <- seurat_object@meta.data
  dataSW <- data_to_plot
  dataSW <- dataSW[order(dataSW$HDS),] 
  
  # Initialize variables
  window_size <- 2000
  medians <- numeric()  # Store medians
  q1s <- numeric()
  q3s <- numeric() # Store IQRs
  positions <- numeric() # Store the middle position of the window
  
  # Perform sliding window calculations
  for (i in 1:(nrow(dataSW) - window_size + 1)) {
    # Get the current window
    window <- dataSW[i:(i + window_size - 1), ]
    
    # Calculate statistics for SenMayo
    medians <- c(medians, median(window$aucell_shgs))
    q1s <- c(q1s, quantile(window$aucell_shgs, 0.25))
    q3s <- c(q3s, quantile(window$aucell_shgs, 0.75))
    
    # Store the middle position of the window for plotting
    positions <- c(positions, mean(window$HDS))
  }
  
  # Combine results into a data frame
  resultsSWsmoothMedian <- data.frame(
    Position = positions,
    Median = medians,
    Q1 = q1s,
    Q3 = q3s
  )
  
  resultsSWsmoothMedian$IQR <- (abs(resultsSWsmoothMedian$Q3)- 
                                  abs(resultsSWsmoothMedian$Q1))
  
  
  IQR_MeanHDS <-  ggplot(resultsSWsmoothMedian, aes(x = Position, y = IQR)) + 
    geom_line(size = 1.5) + xlab("Mean HDS (2000 cells sliding window)") + 
    theme_classic() + theme(axis.text = element_text(size = 18),
                            text = element_text(size = 18)) + xlim(c(min(dataSW$HDS),max(dataSW$HDS)))
  
  shgs_hds_glm <- ggplot(data = seurat_object@meta.data, aes(x = HDS, y = aucell_shgs)) +
    geom_point(size =1,alpha = 0.8, color = "grey") + 
    stat_smooth(fill = "orange", color = "black", lwd = 1.5) + theme_classic() + 
    theme(axis.text = element_text(size = 18),text = element_text(size = 18)) + 
    xlim(c(min(seurat_object@meta.data$HDS),max(seurat_object@meta.data$HDS)))
  
  plot_all_glm  <-cowplot::plot_grid(plotlist = list(shgs_hds_glm,IQR_MeanHDS), 
                                     rel_heights = c(0.75,0.25), nrow = 2) + 
    ggtitle("Hepatocytes from all conditions")
  
  
  plot_all_glm
  
  # Heteroskedasticity test
  
  # 1.a. Using Residual Analysis
  # You can fit a linear model of one variable against the other and then test 
  # the residuals for heteroskedasticity:
  
  # Fit a linear model
  fit <- lm(aucell_shgs ~ HDS, data =  seurat_object@meta.data)
  # Perform Breusch-Pagan test for heteroskedasticity
  bptest(fit)
  
  # studentized Breusch-Pagan test
  #data:  fit
  #BP = 774.74, df = 1, p-value < 2.2e-16
}

#PANEL C: Cell state assginment
{
  
  umap_cellfate <- DimPlot( seurat_object, 
                            group.by = 'cellfate', pt.size = 1, cells = 
                              sample(rownames( seurat_object@meta.data), replace = FALSE)) + 
    scale_color_manual(values = c('undamaged' =  '#CCBB44',
                                  'transition' = '#228833',
                                  'damaged-not-senescent' = '#4477AA',
                                  'damaged-unresolved' = '#EE6677',
                                  'damaged-senescent' = '#AA3377'))+
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          legend.key.height=unit(5,"mm"),
          legend.text = element_text(size = 6))
  
  umap_cellfate
  
  
  
  # Stacked barplots of cellfate assginment
  
  # Create a stacked bar plot
  stackedBarPlotCellFates <- ggplot(summarized_cellfates, aes(x = Group, y = proportion, fill = cellfate)) +
    geom_bar(stat = "identity", position = "stack") + 
    scale_fill_manual(values = c('undamaged' =  '#CCBB44', 'transition' = '#228833', 'damaged-not-senescent' = '#4477AA',
                                 'damaged-unresolved' = '#EE6677','damaged-senescent' = '#AA3377')) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      y = "Proportion of Hepatocytes (%)",
      x = "",
      #  title = "Normalized number of hepatocytes \n from each cell fate per sample "
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text (hjust = 1,
                                size = 18),
      plot.margin = margin(5,5,40,5),
      legend.position = "right",
      legend.text = element_text(size = 18)
    ) + geom_text(
      data = label_data,
      aes(x = Group, y = 1.05, label = n),  # Place labels above the bar
      # fontface = "bold",
      size = 7,
      inherit.aes = FALSE  # Prevent unwanted aesthetics from being inherited
    )
  
  seurat_object <- PrepSCTFindMarkers( seurat_object)
  Idents( seurat_object) <- "cellfate"
  
  
  senescence_damage_markers <- FindMarkers( seurat_object, ident.1 =  'damaged-senescent',
                                            ident.2 = 'damaged-not-senescent')
  
  
  
}


# Panels D and E 
# Selected differentially expressed genes across cell states
{
  
  selected_downregulated_ds_markers <- c("Cyp7a1","Cyp27a1","Sult2a8","Abcb11","Nudt7", "Scp2", "Lpin1", "Igf1", "Sardh", "Cyp2e1","Hao1")
  selected_upregulated_ds_markers <- c("Rtn4","Iqgap1","Bcl2l1" )
  
  
  plot_all_selected_downregulated_markers <- VlnPlot(seurat_object, 
                                                     features = selected_downregulated_ds_markers, 
                                                     split.by = "cellfate", stack = TRUE, layer = "data") +
    scale_fill_manual(values = c('undamaged' =  '#CCBB44',
                                 'transition' = '#228833', 'damaged-not-senescent' = '#4477AA',
                                 'damaged-unresolved' = '#EE6677', 'damaged-senescent' = '#AA3377')) + 
    theme(axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 5))
  
  
  selected_positive_markers <- VlnPlot(seurat_object, 
                                       features = selected_upregulated_ds_markers, 
                                       split.by = "cellfate", stack = TRUE, layer = "data") +
    scale_fill_manual(values = c('undamaged' =  '#CCBB44',
                                 'transition' = '#228833', 'damaged-not-senescent' = '#4477AA',
                                 'damaged-unresolved' = '#EE6677', 'damaged-senescent' = '#AA3377')) + 
    theme(axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 5))
  
  selected_positive_markers
  
}

# Panel F: cell state assignment Wunderlich data
{
  # Applying HDS to Thomas scRNAseq data
  library(Seurat)
  library(AUCell)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(gridExtra)
  
  path_to_data <- "Data/mouse/wunderlich_NIK_NASH_CDA/"
  data <- readRDS(paste0(path_to_data,"nik_samples_annotated.rds"))
  data <- subset(data, subset = class_annotation == "Hepatocyte")
  data <- subset(data, subset = (Diet == "Nash" | Diet == "Cdaa") & Genetic == "NIK-FL")
  metadata_hepatocytes <- data@meta.data
  data <- SCTransform(data)
  
  # hepatocyte damage associated genes
  HDAG <- read.csv(
    file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAG.csv",
    sep = ',')
  
  #SHGS 
  shgs_genes_mouse <- readRDS("/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/misc/SHGS_Kuo_Du_2025_converted_to_mouse.rds")
  
  # load functions to calculate HDS
  source('SharedFunctions.R')
  
  
  HDS <- DS_calc.func( exprMatrices = GetAssayData(data,assay = 'SCT', layer = 'counts'), 
                       DSignature = HDAG )
  
  if(identical(rownames(data@meta.data), names(HDS)) == TRUE){
    data@meta.data$HDS <- unname(HDS)
    print("Matching IDs!")}
  
  
  # SHGS - Senescencent Hepatocytes Gene Signature 
  cells_rankings <- AUCell_buildRankings(GetAssayData(data,assay = 'SCT', layer = 'counts'))
  shgs_score <- AUCell::AUCell_calcAUC(shgs_genes_mouse, cells_rankings)
  
  png( height = 6, width = 12, units = 'in',res = 300,
       file = paste0(path_to_data,"shgs_aucell_histogram.png"))
  AUCell_exploreThresholds(shgs_score, plotHist = TRUE, assign = TRUE)
  dev.off()
  dev.off()
  
  if(identical(colnames(shgs_score@assays@data$AUC), rownames(data@meta.data)) == TRUE) {
    print("Matching IDs!")
    data@meta.data$aucell_shgs <- c(shgs_score@assays@data$AUC)}
  
  data@meta.data$Diet <- factor(data@meta.data$Diet,
                                ordered = TRUE,
                                levels = c("Nash", "Cdaa"))
  
  # Ploting 
  
  # use the same thresholds used in the analysis of the Carlessi data (Panel C)
  quantiles_hds <- c(-0.15723161, -0.04037954)
  quantiles_shgs <- c(0.01360055, 0.03118557)
  
  
  # Add or overwrite the 'group' column based on the conditions
  data@meta.data <- data@meta.data %>%
    mutate(cellfate = case_when(
      HDS <= quantiles_hds[1]  ~ "undamaged",
      HDS > quantiles_hds[1] & HDS <= quantiles_hds[2] ~ "transition",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[2]  ~ "damaged_senescent",
      HDS > quantiles_hds[2] & aucell_shgs <= quantiles_shgs[1] ~ "damaged_not_senescent",
      HDS > quantiles_hds[2] & aucell_shgs > quantiles_shgs[1] & aucell_shgs <= quantiles_shgs[2] ~ "damaged_unresolved",
      TRUE ~ NA_character_  # optional: assign NA if no condition is met
    ))
  
  # Control number of cells assigned per cell fate
  table(data@meta.data$cellfate)
  # Control if any cells have not been assigned a cell fate 
  sum(is.na(data@meta.data$cellfate))
  
  # Cell fates as ordered factors 
  data@meta.data$cellfate <- factor(data@meta.data$cellfate, 
                                    ordered = TRUE,
                                    levels = c('undamaged', 
                                               'transition',
                                               'damaged_not_senescent',
                                               'damaged_unresolved',
                                               'damaged_senescent'),
                                    labels = c('undamaged',
                                               'transition',
                                               'damaged-not-senescent',
                                               'damaged-unresolved',
                                               'damaged-senescent'))
  
  
  # Summarize the data: count the number of cells for each fate in each sample
  summarized_cellfates <- data@meta.data %>%
    group_by(Diet, cellfate) %>%
    summarize(Cell_Count = n(), .groups = "drop") %>%
    group_by(Diet) %>% 
    mutate(proportion = Cell_Count / sum(Cell_Count))
  
  
  # Prepare a table-like data frame for annotations
  label_data <- summarized_cellfates %>% 
    group_by(Diet) %>%
    mutate(n = sum(Cell_Count)) %>%
    dplyr::select(Diet, cellfate, Cell_Count,n) %>%
    pivot_wider(names_from = cellfate, values_from = Cell_Count, values_fill = 0)
  
  # Convert the data frame to a long format for plotting as a table
  label_data_long <- label_data[,-2] %>%
    pivot_longer(
      cols = -Diet,
      names_to = "cellfate",
      values_to = "Cell_Count"
    ) %>% 
    mutate(cellfate = factor(cellfate,
                             levels = c('undamaged',
                                        'transition',
                                        'damaged-not-senescent',
                                        'damaged-unresolved',
                                        'damaged-senescent')))
  
  
  # Create a stacked bar plot
  stackedBarPlotCellFates <- ggplot(summarized_cellfates, aes(x = Diet, y = proportion, fill = cellfate)) +
    geom_bar(stat = "identity", position = "stack") + 
    scale_fill_manual(values = c('undamaged' =  '#CCBB44', 'transition' = '#228833', 'damaged-not-senescent' = '#4477AA',
                                 'damaged-unresolved' = '#EE6677','damaged-senescent' = '#AA3377')) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      y = "Proportion of Hepatocytes (%)",
      x = "",
      #  title = "Normalized number of hepatocytes \n from each cell fate per sample "
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text (hjust = 1,
                                size = 18),
      plot.margin = margin(5,5,40,5),
      legend.position = "right",
      legend.text = element_text(size = 18)
    ) + geom_text(
      data = label_data,
      aes(x = Diet, y = 1.05, label = n),  # Place labels above the bar
      # fontface = "bold",
      size = 7,
      inherit.aes = FALSE  # Prevent unwanted aesthetics from being inherited
    )
  
  stackedBarPlotCellFates
  
}

