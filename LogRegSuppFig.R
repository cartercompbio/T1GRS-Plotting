
# Replacement panels for Fig3
setwd('~/Code/T1D/')

# load in new delong, make auc plot, make density plot

# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)

# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)
linear <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_linear.txt", sep="\t", header=TRUE)
linear <- linear[,c('CV_ID','Probability')]
colnames(linear)<-c('CV_ID','LogReg')
ALL<-merge(ALL,linear,by='CV_ID')

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

# find FP and FN from ALL
# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- as.numeric(ALL$disease)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$DISEASE, predictor = ALL$Probability)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$LogReg)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")

# Print the results
print(delong_test)

# Create data frames for ggplot
# pROC stores sensitivities and specificities in reverse order, so we can just reverse them
prob_df <- data.frame(
  FPR = 1 - prob_roc$specificities,
  TPR = prob_roc$sensitivities,
  Model = "T1GRS"
)

grs2_df <- data.frame(
  FPR = 1 - grs2_roc$specificities,
  TPR = grs2_roc$sensitivities,
  Model = "LogReg"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(size = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="Discovery: Full Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("LogReg AUC = ", round(grs2_auc, 3)), 
           color = "goldenrod", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "LogReg" = "goldenrod"))

ggsave('figs_may/fig3/auc_T1GRS_vs_LogReg_full.pdf',width = 6,height =5)










# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_HLA_YESPCS_catboost.txt", sep="\t", header=TRUE)
linear <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_linear_MHCONLY.txt", sep="\t", header=TRUE)
linear <- linear[,c('CV_ID','Probability')]
colnames(linear)<-c('CV_ID','LogReg')
ALL<-merge(ALL,linear,by='CV_ID')

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

# find FP and FN from ALL
# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- as.numeric(ALL$disease)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$DISEASE, predictor = ALL$Probability)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$LogReg)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")

# Print the results
print(delong_test)

# Create data frames for ggplot
# pROC stores sensitivities and specificities in reverse order, so we can just reverse them
prob_df <- data.frame(
  FPR = 1 - prob_roc$specificities,
  TPR = prob_roc$sensitivities,
  Model = "T1GRS"
)

grs2_df <- data.frame(
  FPR = 1 - grs2_roc$specificities,
  TPR = grs2_roc$sensitivities,
  Model = "LogReg"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(size = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="Discovery: MHC Only",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("LogReg AUC = ", round(grs2_auc, 3)), 
           color = "goldenrod", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "LogReg" = "goldenrod"))

ggsave('figs_may/fig3/auc_T1GRS_vs_LogReg_MHCONLY.pdf',width = 6,height =5)










# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_NonHLA_YESPCS_catboost.txt", sep="\t", header=TRUE)
linear <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_linear_NonMHC.txt", sep="\t", header=TRUE)
linear <- linear[,c('CV_ID','Probability')]
colnames(linear)<-c('CV_ID','LogReg')
ALL<-merge(ALL,linear,by='CV_ID')

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

# find FP and FN from ALL
# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- as.numeric(ALL$disease)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$DISEASE, predictor = ALL$Probability)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$LogReg)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")

# Print the results
print(delong_test)

# Create data frames for ggplot
# pROC stores sensitivities and specificities in reverse order, so we can just reverse them
prob_df <- data.frame(
  FPR = 1 - prob_roc$specificities,
  TPR = prob_roc$sensitivities,
  Model = "T1GRS"
)

grs2_df <- data.frame(
  FPR = 1 - grs2_roc$specificities,
  TPR = grs2_roc$sensitivities,
  Model = "LogReg"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(size = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="Discovery: Non-MHC",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("LogReg AUC = ", round(grs2_auc, 3)), 
           color = "goldenrod", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "LogReg" = "goldenrod"))

ggsave('figs_may/fig3/auc_T1GRS_vs_LogReg_Non-MHC.pdf',width = 6,height =5)
