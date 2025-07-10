# Supplementary Figure 8 Panel A:
{
  # aligmnet of models 
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  
  
  # Load processed 3 snRNA-seq datasets
  # From code chunk for FIGURE 6 PANEL A:
  merged_carlessi <- readRDS("seurat_object_SCT_HDS_and_Senescence_scores.rds")
  
  # 36 weeks old mice --> bad quality data, strong batch effects, don't include in any analysis 
  # From code chunk for SUPPLEMENTARY FIGURE 1 C:
  liver_cell_atlas <- readRDS( "/liver_cell_atlas/Seurat_Object_OnlyHepatocytes_SelectedSamples_HDS.rds")
  liver_cell_atlas <- subset(liver_cell_atlas, subset = age == "24 weeks")
  
  #From code chunk SUPPLEMENTARY FIGURE 12 C:
  xiao_merged <- readRDS("Xiao_merged_seurat_SCT_HDS_and_SHGS.rds")
  
  
  View(merged_carlessi@meta.data)
  View(liver_cell_atlas@meta.data)
  View(xiao_merged@meta.data)
  
  all_data <- data.frame("dataset" = rep("carlessi", length(merged_carlessi@meta.data$orig.ident)),
                         "sample_id" = merged_carlessi@meta.data$SampleID,
                         "orig_condition" = merged_carlessi@meta.data$Group,
                         "HDS" =  merged_carlessi@meta.data$HDS )
  
  all_data <- rbind.data.frame(all_data,
                               data.frame("dataset" = rep("liver_cell_atlas", length(liver_cell_atlas@meta.data$orig.ident)),
                                          "sample_id" = liver_cell_atlas@meta.data$sample,
                                          "orig_condition" = liver_cell_atlas@meta.data$diet,
                                          "HDS" = liver_cell_atlas@meta.data$HDS),
                               data.frame("dataset" = rep("xiao", length(xiao_merged@meta.data$orig.ident)),
                                          "sample_id" = xiao_merged@meta.data$orig.ident,
                                          "orig_condition" = xiao_merged@meta.data$condition,
                                          "HDS" = xiao_merged@meta.data$HDS))
  
  unique(all_data$orig_condition)
  
  # Add or overwrite the 'group' column based on the conditions
  all_data <- all_data %>%
    mutate(condition_summary = case_when( 
      ( orig_condition == "Healthy" | orig_condition == "SD" | orig_condition ==  "3m NC" | orig_condition ==  "9m NC" ) ~ "control",
      (orig_condition == "CDE" | orig_condition == "TAA" | orig_condition == "WD"| orig_condition == "3m NASH" | orig_condition == "9m NASH") ~ "experiment"
    ))
  
  all_data$condition_summary <- factor(all_data$condition_summary, levels = c("control", "experiment"), ordered = TRUE)
  
  ggplot(data = all_data,
         aes(x = HDS, group = sample_id, color = condition_summary),alpha = 0.3)+
    geom_line(stat = "density", adjust = 1.2, linewidth = 1.2, trim = TRUE, alpha = 0.5)+ theme_minimal()+
    scale_colour_colorblind()+
    theme(legend.position = "none",
          axis.ticks.y = element_blank()) 
  
  
  # HDS test if control samples are significantly better aligned than experimental samples for Suppl. Figure 8 b
  {
    HDSdat <- all_data
    
    # modify metadata
    HDSdat$sample <-  paste(HDSdat$sample_id , HDSdat$condition_summary, sep = "__")
    HDSdat$group <- paste(HDSdat$dataset , HDSdat$condition_summary, sep = "_")
    
    ## test on group level
    combinations <- combn( unique( HDSdat$group) , 2, simplify = FALSE)
    
    
    HDS.kstest <- Reduce( rbind ,lapply(seq(combinations), function(ii)
    {
      print(ii)
      pair <- combinations[[ii]]
      x <- HDSdat$HDS[HDSdat$group==pair[1] ]
      y <- HDSdat$HDS[HDSdat$group==pair[2] ]
      test <- tryCatch(  ks.test(x, y),  error = function(e) NA)
      if( length(test)==1) {
        sstat <- NA
        ppval <- NA
      } else {
        sstat <- test$statistic
        ppval <- test$p.value
      }
      data.frame(
        var1 = pair[1],
        var2 = pair[2],
        gtypeVar1 = sub(".*_","",pair[1]),
        gtypeVar2 = sub(".*_","",pair[2]),
        
        statistic = sstat,
        p_value = ppval,
        stringsAsFactors = FALSE
      ) } ) )
    
    
    wilcox.test( HDS.kstest$statistic[ HDS.kstest$gtypeVar1=="control" &
                                         HDS.kstest$gtypeVar2=="control"   ] ,
                 HDS.kstest$statistic[HDS.kstest$gtypeVar1!="control" &
                                        HDS.kstest$gtypeVar2!="control"],
                 alternative =  "less")
    
    
    ### test on sample level
    
    combinations.S <- combn( unique( HDSdat$sample), 2, simplify = FALSE)
    HDS.kstest.S <- Reduce( rbind ,lapply(seq(combinations.S), function(ii)
    {
      pairS <- combinations.S[[ii]]
      x <- HDSdat$HDS[HDSdat$sample==pairS[1] ]
      y <- HDSdat$HDS[HDSdat$sample==pairS[2] ]
      
      test <- tryCatch(  ks.test(x, y),  error = function(e) NA)
      if( length(test)==1 ) {
        sstat <- NA
        ppval <- NA
      } else {
        sstat <- test$statistic
        ppval <- test$p.value
      }
      data.frame(
        var1 = pairS[1],
        var2 = pairS[2],
        gtypeVar1 = sub(".*__","",pairS[1]),
        gtypeVar2 = sub(".*__","",pairS[2]),
        
        statistic = sstat,
        p_value = ppval,
        stringsAsFactors = FALSE
      ) }
    ) )
    
    
    wilcox.test( HDS.kstest.S$statistic[ HDS.kstest.S$gtypeVar1=="control" &
                                           HDS.kstest.S$gtypeVar2=="control"   ] ,
                 HDS.kstest.S$statistic[HDS.kstest.S$gtypeVar1!="control" &
                                          HDS.kstest.S$gtypeVar2!="control"],
                 alternative =  "less")
    
    
    
  }
  
}