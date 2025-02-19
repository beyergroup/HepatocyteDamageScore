## 1. Testing for Heteroskedasticity
# Heteroskedasticity refers to non-constant variance in the data, which is 
# what you're looking for. To test for it:

mergedObj <- readRDS(
  file = paste0(outputPath,
                "MergedXiaoHDS_SenMayoHepatocytes.rds"))

mergedDataFrame <- mergedObj@meta.data


mergedDataFrame$CellType <- factor(mergedDataFrame$CellCluster,
                                   ordered = TRUE,
                                   levels = c("PP-Hep",
                                              "Int-Hep",
                                              "PC-Hep",
                                              "mNASH-Hep1",
                                              "mNASH-Hep2"))

mergedDataFrame$condition <- factor(mergedDataFrame$condition,
                                    ordered = TRUE,
                                    levels = c("3m NC","9m NC", 
                                               "3m NASH", "9m NASH" ))



 
# 1.a. Using Residual Analysis
# You can fit a linear model of one variable against the other and then test 
# the residuals for heteroskedasticity:

# Fit a linear model
fit <- lm(AUCell_SenMayo ~ HDS, data = mergedDataFrame)

# Perform Breusch-Pagan test for heteroskedasticity
library(lmtest)
bptest(fit)

# studentized Breusch-Pagan test
# 
# data:  fit
# BP = 553.95, df = 1, p-value < 2.2e-16

# we can reject the null hypothesis that says homoscedasticity is present,
# meaning heteroskedastisticity is present 


# 1.b. Using Levene’s Test
# Levene's test checks if the variance of one variable differs across 
# groups of another variable:

# assumptions of the test: 
# 1. independent observations --> not really true cause 
# single cell data is not
# 2. test variable has a metric scale (it is true, I think )


library(car)
leveneTest( AUCell_SenMayo ~ cut(HDS, breaks = 10), data = mergedDataFrame)
leveneTest(AUCell_SenMayo ~ cut(HDS, breaks = 10), data = mergedDataFrame)
car::leveneTest(AUCell_SenMayo ~ cut(HDS, breaks = 10), data = mergedDataFrame)


# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
# group     9  63.411 < 2.2e-16 ***
#       15815

# H0 can be rejected, groups have significantly different variances



# 2. Smooth-Median with Interquartile Ranges
# This approach is more descriptive and less assumption-driven than formal 
# statistical tests. It’s useful for visualizing trends in variance across 
# the range of the second variable.

# smoothed_data <- mergedDataFrame[,c('AUCell_SenMayo','HDS')] %>%
#   mutate(bin = cut(HDS, breaks = 15)) %>% # Create 20 bins of phenotype2
#   group_by(bin) %>%
#   summarize(
#     bin_center = median(HDS),
#     median = median(AUCell_SenMayo),
#     Q1 = quantile(AUCell_SenMayo, 0.25),
#     Q3 = quantile(AUCell_SenMayo, 0.75)
#   )  
# 
# # Plot the smooth-median with IQR
# smoothPlot <- ggplot(smoothed_data, aes(x = bin_center)) +
#   geom_line(aes(y = median), color = "black", size = 1) + # Median line
#   geom_ribbon(aes(ymin = Q1, ymax = Q3), fill = "grey", alpha = 0.7) + # IQR ribbon
#   labs(
#     x = "HDS (Median of Binned Centers)",
#     y = "SenMayo Activity",
#     title = "Smooth Median with IQR"
#   ) + theme_bw()
# 
# png( height = 4, width = 6, 
#      units = 'in',
#      res = 300,
#      file = paste0(outputPath,
#                    "Results/CellFateAnalysis/PNGs/",
#                    "Xiao_SmoothMedianwithIQR_SenMayo_HDS_Xiao_mergedheps.png") )
# 
# smoothPlot
# 
# dev.off()

# 3. Smooth Median with IQR, sliding window approach 

# Sort data by Phenotype1
dataSW <- mergedDataFrame[,c('AUCell_SenMayo','HDS')]
dataSW <- dataSW[order(dataSW$HDS),] 

# Initialize variables
window_size <- 1000
medians <- numeric()  # Store medians
q1s <- numeric()
q3s <- numeric() # Store IQRs
positions <- numeric() # Store the middle position of the window

# Perform sliding window calculations
for (i in 1:(nrow(dataSW) - window_size + 1)) {
  # Get the current window
  window <- dataSW[i:(i + window_size - 1), ]
  
  # Calculate statistics for SenMayo
  medians <- c(medians, median(window$AUCell_SenMayo))
  q1s <- c(q1s, quantile(window$AUCell_SenMayo, 0.25))
  q3s <- c(q3s, quantile(window$AUCell_SenMayo, 0.75))
  
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


resultsSWsmoothMedian$IQR <- (abs(resultsSWsmoothMedian$Q1)+ 
                                abs(resultsSWsmoothMedian$Q3))

SWsmoothPlot <- ggplot(resultsSWsmoothMedian,
                       aes(x = Position)) +
  geom_line(aes(x = Position,
                y = Median,
                color = "Median"),
            color = "black", size = 1) +
  geom_ribbon(aes(ymin = Q1,
                  ymax = Q3),
              fill = "grey", alpha = 0.7) + # IQR ribbon
  labs(
    x = "HDS (1000 cells sliding window)",
    y = "Median SenMayo Activity per window",
    title = "Sliding Window Smooth Median with IQR"
  ) + theme_classic() + theme(text = element_text(size=15))

SWsmoothPlotWithDots <- ggplot(data = dataSW) + 
  geom_point(aes(x = HDS, y = AUCell_SenMayo),
             color = "grey", 
             alpha = 0.5, 
             size = 1.5
  ) + geom_line(data = resultsSWsmoothMedian,
                             aes(x = Position, 
                                 y = Median, 
                                 color = "Median"), 
                             color = "black", size = 1.5) + 
  geom_ribbon(data = resultsSWsmoothMedian,
              aes(x = Position,
                  ymin = Q1,
                  ymax = Q3), 
              fill = "orange", alpha = 0.5) + # IQR ribbon
  labs(
    x = "Mean HDS (1000 Cells Sliding Window)",
    y = "Median SenMayo Activity per Window",
  ) + theme_classic() + 
  theme(axis.text = element_text(size = 18),
        text = element_text(size=18)) +
  xlim(c(min(dataSW$HDS),max(dataSW$HDS)))


IQR_MeanHDS <-  ggplot(resultsSWsmoothMedian, aes(x = Position, y = IQR)) + 
  geom_line(size = 1.5) + xlab("Mean HDS (1000 cells sliding window)") + 
  theme_classic() + theme(axis.text = element_text(size = 18),
        text = element_text(size = 18)) + xlim(c(min(dataSW$HDS),max(dataSW$HDS)))

plotGrid1 <-cowplot::plot_grid(plotlist = list(SWsmoothPlotWithDots,
                                               IQR_MeanHDS), 
                                             rel_heights = c(0.75,0.25), nrow = 2)


png( height = 8, width = 6, 
     units = 'in',
     res = 600,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/",
                   "Grid_IGR_SlidingWindowsSmoothMedianIQRWithCell_SenMayo_HDS_Xiao_mergedheps.png") )

plotGrid1

dev.off()

pdf( height = 8, width = 6, 
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/",
                   "Grid_IGR_SlidingWindowsSmoothMedianIQRWithCell_SenMayo_HDS_Xiao_mergedheps.pdf") )

plotGrid1

dev.off()





## save plots separately 

png( height = 5, width = 7, 
     units = 'in',
     res = 600,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/",
                   "Xiao_SlidingWindowsSmoothMedianIQR_SenMayo_HDS_Xiao_mergedheps.png") )

SWsmoothPlot 

dev.off()

png( height = 5, width = 7, 
     units = 'in',
     res = 600,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/",
                   "Xiao_SlidingWindowsSmoothMedianIQR_SenMayo_HDS_Xiao_mergedheps_withCELLS.png") )

SWsmoothPlotWithDots

dev.off()


png( height = 5, width = 4, 
     units = 'in',
     res = 300,
     file = paste0(outputPath,
                   "Results/CellFateAnalysis/PNGs/",
                   "Xiao_IQR_slidingWindow_MeanHDS.png") )

ggplot(resultsSWsmoothMedian, aes(x = Position, y = IQR)) + 
  geom_line() + xlab("Mean HDS (1000 cells sliding window)") + theme_bw()



dev.off()



