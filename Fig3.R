
# Replacement panels for Fig3
setwd('~/Code/T1D/')

# load in new delong, make auc plot, make density plot

# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)

# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)

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
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$Total_sum)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
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
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/auc_full.pdf',width = 6,height =5)


# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_NonHLA_YESPCS_catboost.txt", sep="\t", header=TRUE)

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

ALL$disease <- as.numeric(ALL$disease)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$DISEASE, predictor = ALL$Probability)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$sum_All5_nonHLA)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="Discovery: Non-MHC Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/auc_nonHLA.pdf',width = 6,height =5)


# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_HLA_YESPCS_catboost.txt", sep="\t", header=TRUE)

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

ALL$disease <- as.numeric(ALL$disease)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$disease, predictor = ALL$Probability)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$Total_sum_HLA)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="Discovery: MHC Only Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", '0.920'), # manually rounded, was giving errors
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = "Delong Test P < 0.0001", 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/auc_HLA.pdf',width = 6,height =5)




###############################
# NPOD DELONG ##
###############################

# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)

# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_T1GRS_Scores.tsv", sep="\t", header=TRUE)
GRS2 <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_GRS2_Scores.tsv", sep="\t", header=TRUE)
ALL<-merge(ALL,GRS2,on='ID')

# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- ifelse(ALL$Condition=='T1D',2,1)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$disease, predictor = ALL$T1GRS_ALL)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$Total_sum)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="NPOD: Full Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = paste0("Delong Test P < ",round(delong_test$p.value,4)), 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/NPOD_auc_full.pdf',width = 6,height =5)

# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_T1GRS_Scores.tsv", sep="\t", header=TRUE)
GRS2 <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_GRS2_Scores.tsv", sep="\t", header=TRUE)
ALL<-merge(ALL,GRS2,on='ID')

# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- ifelse(ALL$Condition=='T1D',2,1)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$disease, predictor = ALL$T1GRS_MHC)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$Total_sum_HLA)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="NPOD: MHC Only Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = paste0("Delong Test P < ",round(delong_test$p.value,4)), 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/NPOD_auc_MHC.pdf',width = 6,height =5)






# Read in your data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_T1GRS_Scores.tsv", sep="\t", header=TRUE)
GRS2 <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_GRS2_Scores.tsv", sep="\t", header=TRUE)
ALL<-merge(ALL,GRS2,on='ID')

# Ensure the outcome is a factor with correct levels (if needed)
# Assuming 'disease' is coded as 0/1 or a factor with levels "control"/"case"
# If it's numeric (0,1), that’s typically fine for pROC.
ALL$disease <- ifelse(ALL$Condition=='T1D',2,1)

# Compute ROC curves
prob_roc <- pROC::roc(response = ALL$disease, predictor = ALL$T1GRS_nonMHC)
grs2_roc <- pROC::roc(response = ALL$disease, predictor = ALL$sum_All5_nonMHC)

# Extract AUC values
prob_auc <- pROC::auc(prob_roc)
grs2_auc <- pROC::auc(grs2_roc)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc)
grs2_ci <- pROC::ci.auc(grs2_roc)
print(prob_ci)
print(grs2_ci)

# Perform DeLong's test to compare the two ROC curves
delong_test <- pROC::roc.test(prob_roc, grs2_roc, method = "delong")
print(delong_test$p.value)
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
  Model = "GRS2"
)

# Combine into one data frame
roc_df <- bind_rows(prob_df, grs2_df)

# Plot with ggplot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_path(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),# Add only bottom x-axis and left y-axis lines
                                        axis.line.x = element_line(color = "black", size = 0.5),
                                        axis.line.y = element_line(color = "black", size = 0.5)) +
  labs(title="NPOD: Non-MHC Model",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  # Annotate AUC values
  annotate("text", x = 0.65, y = 0.15, 
           label = paste0("T1GRS AUC = ", round(prob_auc, 3)), 
           color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, 
           label = paste0("GRS2 AUC = ", round(grs2_auc, 3)), 
           color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, 
           label = paste0("Delong Test P < ",round(delong_test$p.value,4)), 
           color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

ggsave('figs_may/fig3/NPOD_auc_NonMHC.pdf',width = 6,height =5)



###############################
#### Density plots ###
###############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Read the data
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)
# Ensure 'disease' is a factor or numeric 0/1
ALL$disease<-ALL$disease-1
ALL$disease <- factor(ALL$disease, levels = c(0,1), labels = c("No Disease", "Disease"))
ALL <- ALL[!is.na(ALL$disease),]
#ALL <- ALL[ALL$disease=='Disease',]

min_val <- min(ALL$Total_sum, na.rm = TRUE)
max_val <- max(ALL$Total_sum, na.rm = TRUE)
# Create a scaled version of Total_sum from 0 to 1
ALL$Total_sum_scaled <- (ALL$Total_sum - min_val) / (max_val - min_val)

# Use Pro as the ambiguity metric as specified
ALL$Ambiguity <- ALL$Pro

# Transform data into long format for the first two plots
ALL_long <- ALL %>%
  pivot_longer(cols = c("Probability", "Total_sum_scaled"), 
               names_to = "Model", 
               values_to = "Score") %>%
  mutate(Model = ifelse(Model == "Probability", "T1GRS", "GRS2"))

ALL_long$Model<-factor(ALL_long$Model,levels=c('T1GRS','GRS2'))
# Create density plots for T1GRS and GRS2
model_plot <- ggplot(ALL_long, aes(x = Score, fill = disease)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Model, ncol = 1) +
  scale_fill_manual(values = c("Disease" = "goldenrod", "No Disease" = "forestgreen")) +
  theme_minimal(base_size = 14) +
  labs(x = "Score",
       y = "Density",
       fill = "Disease Status") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  )

model_plot

# Calculate density overlap as a measure of ambiguity
# First, create a sequence of values across the score range
score_seq <- seq(0, 1, length.out = 100)

# Function to calculate overlap between two densities
calculate_overlap <- function(values1, values2, score_sequence) {
  # Estimate densities
  d1 <- density(values1, from = 0, to = 1, n = 100)
  d2 <- density(values2, from = 0, to = 1, n = 100)
  
  # Extract density values at score sequence points
  y1 <- approx(d1$x, d1$y, score_sequence)$y
  y2 <- approx(d2$x, d2$y, score_sequence)$y
  
  # Calculate minimum of the two densities at each point (overlap)
  overlap <- mapply(min, y1, y2)
  
  # Return data frame with scores and overlap values
  data.frame(Score = score_sequence, Overlap = overlap)
}

# Calculate overlap for T1GRS and GRS2
t1grs_overlap <- calculate_overlap(
  ALL$Probability[ALL$disease == "Disease"],
  ALL$Probability[ALL$disease == "No Disease"],
  score_seq
)
t1grs_overlap$`Model Asessed` <- "T1GRS"

grs2_overlap <- calculate_overlap(
  ALL$Total_sum_scaled[ALL$disease == "Disease"],
  ALL$Total_sum_scaled[ALL$disease == "No Disease"],
  score_seq
)
grs2_overlap$`Model Asessed` <- "GRS2"

# Combine overlap data
overlap_data <- rbind(t1grs_overlap, grs2_overlap)

# Create overlap plot
overlap_plot <- ggplot(overlap_data, aes(x = Score, y = Overlap, color = `Model Asessed`, fill = `Model Asessed`)) +
  geom_area(alpha = 0.5, position = "identity") +  # Add filled area
  geom_line(size = 1,show.legend = F) +  # Keep the line on top
  theme_minimal(base_size = 14) + scale_y_continuous(breaks = c(0, 1, 2), limits = c(0, 2)) +
  labs(title = "Prediction Overlap", 
       x = "Score", 
       y = "Overlap",
       color = "Model Asessed",
       fill = "Model Asessed") +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick")) +
  scale_fill_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick")) +  # Same colors for fill
  #scale_x_continuous(expand = expansion(add = c(1.5, 0))) +
  theme(
    title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=12),axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)
  )

library(gridtext)
# Arrange the plots into a single figure with a common title
combined_plot <- grid.arrange(
  model_plot, 
  overlap_plot,
  nrow = 2,
  heights = c(2, 1)
  #top = textGrob("T1GRS, GRS2, and Ambiguity Distributions", 
  #               gp = gpar(fontsize = 16, font = 2))
)

# Save the combined plot
ggsave("/Users/tjsears/Code/T1D/figs_april/fig3/density_with_ambiguity.pdf", 
       combined_plot, width = 8, height = 6)









##################################
# Reclassification Analysis ######
##################################


# First read the data and create the basic variables
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)
ALL <- ALL[!is.na(ALL$disease),]
ALL$disease <- factor(ALL$disease, levels = c(1,2), labels = c("No Disease", "Disease"))
clinical_data <- read.table('clinical_data/Fixed_RenalGeno_HLAfam_GENIE_UK_DCCT_GoKIND_T1DGC_phenotypes.txt', sep='\t', header=T)
ALL <- merge(clinical_data, ALL, by.x='FID', by.y='CV_ID')

df <- ALL
df$HLA_simple <- df$HLA_Status

# Get optimal cutpoints
library(cutpointr)
x <- cutpointr(df, Total_sum, disease)$optimal_cutpoint
df$GRS2_Assignment <- ifelse(df$Total_sum >= x, 'GRS2_disease', 'GRS2_no_disease')

y <- cutpointr(df, Probability, disease)$optimal_cutpoint
df$T1GRS_Assignment <- ifelse(df$Probability >= y, 'T1GRS_disease', 'T1GRS_no_disease')

# write out table of FIDs and accuracy
temp<-df[,c('FID','T1GRS_Assignment','DISEASE')]
write.table(temp,'/Users/tjsears/Code/T1D/data2025/accuracy.txt',sep='\t',col.names =T,quote=F)

# Compute confusion matrix values based on assignments
df$T1GRS_FP <- ifelse(df$T1GRS_Assignment == "T1GRS_disease" & df$DISEASE == 1, 1, 0)
df$T1GRS_FN <- ifelse(df$T1GRS_Assignment == "T1GRS_no_disease" & df$DISEASE == 2, 1, 0)
df$GRS2_FP <- ifelse(df$GRS2_Assignment == "GRS2_disease" & df$DISEASE == 1, 1, 0)
df$GRS2_FN <- ifelse(df$GRS2_Assignment == "GRS2_no_disease" & df$DISEASE == 2, 1, 0)

# Calculate summary statistics by HLA status
library(dplyr)
fp_summary <- df %>%
  group_by(HLA_Status) %>%
  summarize(
    T1GRS_FP_count = sum(T1GRS_FP),
    GRS2_FP_count = sum(GRS2_FP),
    Total_nonDisease = sum(DISEASE == 1)
  ) %>%
  mutate(
    T1GRS_FP_rate = T1GRS_FP_count / Total_nonDisease,
    GRS2_FP_rate = GRS2_FP_count / Total_nonDisease
  )

# For False Negatives
fn_summary <- df %>%
  group_by(HLA_Status) %>%
  summarize(
    T1GRS_FN_count = sum(T1GRS_FN),
    GRS2_FN_count = sum(GRS2_FN),
    Total_Disease = sum(DISEASE == 2)
  ) %>%
  mutate(
    T1GRS_FN_rate = T1GRS_FN_count / Total_Disease,
    GRS2_FN_rate = GRS2_FN_count / Total_Disease
  )

# Create new dataframes with differences between models
fp_diff <- fp_summary %>%
  mutate(Difference =  GRS2_FP_rate - T1GRS_FP_rate) %>%
  select(HLA_Status, Difference)

fn_diff <- fn_summary %>%
  mutate(Difference = GRS2_FN_rate - T1GRS_FN_rate) %>%
  select(HLA_Status, Difference)

# Get unique HLA groups
hla_groups <- unique(df$HLA_Status)

# Corrected Chi-squared tests for FALSE POSITIVES comparison
fp_p_values <- vector("numeric", length(hla_groups))
names(fp_p_values) <- hla_groups

for (i in 1:length(hla_groups)) {
  hla <- hla_groups[i]
  subdata <- df[df$HLA_Status == hla & df$DISEASE == 1, ]  # Only non-disease cases
  
  # Skip groups with too few samples
  if(nrow(subdata) < 10) {
    fp_p_values[i] <- NA
    next
  }
  
  # McNemar's test is more appropriate for comparing two classifiers on the same dataset
  # Create 2x2 table for matched pairs: rows=T1GRS (0/1), cols=GRS2 (0/1)
  contingency <- matrix(c(
    sum(subdata$T1GRS_FP == 0 & subdata$GRS2_FP == 0),  # Both correct
    sum(subdata$T1GRS_FP == 0 & subdata$GRS2_FP == 1),  # Only GRS2 wrong
    sum(subdata$T1GRS_FP == 1 & subdata$GRS2_FP == 0),  # Only T1GRS wrong
    sum(subdata$T1GRS_FP == 1 & subdata$GRS2_FP == 1)   # Both wrong
  ), nrow = 2)
  
  # Run McNemar's test
  if(sum(contingency) > 0 && (contingency[1,2] + contingency[2,1]) > 0) {
    test_result <- mcnemar.test(contingency, correct = FALSE)
    fp_p_values[i] <- test_result$p.value
  } else {
    fp_p_values[i] <- NA
  }
}

# Corrected Chi-squared tests for FALSE NEGATIVES comparison
fn_p_values <- vector("numeric", length(hla_groups))
names(fn_p_values) <- hla_groups

for (i in 1:length(hla_groups)) {
  hla <- hla_groups[i]
  subdata <- df[df$HLA_Status == hla & df$DISEASE == 2, ]  # Only disease cases
  
  # Skip groups with too few samples
  if(nrow(subdata) < 10) {
    fn_p_values[i] <- NA
    next
  }
  
  # McNemar's test for matched pairs
  contingency <- matrix(c(
    sum(subdata$T1GRS_FN == 0 & subdata$GRS2_FN == 0),  # Both correct
    sum(subdata$T1GRS_FN == 0 & subdata$GRS2_FN == 1),  # Only GRS2 wrong
    sum(subdata$T1GRS_FN == 1 & subdata$GRS2_FN == 0),  # Only T1GRS wrong
    sum(subdata$T1GRS_FN == 1 & subdata$GRS2_FN == 1)   # Both wrong
  ), nrow = 2)
  
  # Run test
  if(sum(contingency) > 0 && (contingency[1,2] + contingency[2,1]) > 0) {
    test_result <- mcnemar.test(contingency, correct = FALSE)
    fn_p_values[i] <- test_result$p.value
  } else {
    fn_p_values[i] <- NA
  }
}

# Format p-values for display
format_pvalue <- function(p) {
  if(is.na(p)) return("")
  if(p < 0.001) return("p < 0.001")
  if(p < 0.01) return(paste0("p = ", sprintf("%.3f", p)))
  if(p < 0.05) return(paste0("p = ", sprintf("%.2f", p)))
  return(paste0("p = ", sprintf("%.2f", p)))
}

fp_formatted_pvalues <- sapply(fp_p_values, format_pvalue)
fn_formatted_pvalues <- sapply(fn_p_values, format_pvalue)

# Get the number of HLA groups for positioning labels
n_groups <- length(hla_groups)

# Create factors with ordered levels for horizontal ordering
fp_diff$HLA_Status_ordered <- factor(fp_diff$HLA_Status, 
                                     levels = fp_diff$HLA_Status[order(fp_diff$Difference)])

fn_diff$HLA_Status_ordered <- factor(fn_diff$HLA_Status, 
                                     levels = fn_diff$HLA_Status[order(fn_diff$Difference)])




# PLOT

# FALSE POSITIVE DIFFERENCE PLOT - Flipped horizontally with corrected annotations
library(ggplot2)
fp_diff_plot <- ggplot(fp_diff, aes(x = HLA_Status_ordered, y = Difference, fill = Difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "slateblue", "FALSE" = "firebrick"),
                    labels = c("TRUE" = "T1GRS higher", "FALSE" = "GRS2 higher"),
                    name = "Comparison") +
  labs(
    title = "Difference in False Positive Rate (GRS2 - T1GRS)",
    x = "HLA Status",
    y = "Difference in False Positive Rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(face = 'bold'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5)  # Remove legend since we'll add our own annotations
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at 0
  annotate("segment", y = -0.01, yend = -0.04, x = 0.05, xend = 0.05, 
           arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  annotate("segment", y = 0.01, yend = 0.04, x = 0.05, xend = 0.05, 
           arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  annotate("text", y = -0.025, x = 0.3, label = "T1GRS Worse", size = 3.5) +
  annotate("text", y = 0.025, x = 0.3, label = "T1GRS Better", size = 3.5) +
  coord_flip() +  ylim(-0.05,0.1) + #ylim(-0.1,10)# Flip coordinates to make bars horizontal 
  scale_x_discrete(limits = rev(levels(fn_diff$HLA_Status_ordered)),expand = expansion(add = c(1.5, 0)))  # Reverse order so largest diff at top

# Add FP p-value annotations - corrected for horizontal orientatio
for (hla in rev(levels(fn_diff$HLA_Status_ordered))) {
  
  # Find position in the ordered factor
  y_pos <- which(rev(levels(fn_diff$HLA_Status_ordered) == hla))
  
  # Calculate position for annotation
  diff_val <- fp_diff$Difference[fp_diff$HLA_Status == hla]
  if(is.na(diff_val)) next
  
  x_direction <- ifelse(diff_val > 0, 1, -1)
  x_pos <- diff_val
  
  # Add annotation
  fp_diff_plot <- fp_diff_plot + 
    annotate("text", x = y_pos, y = (x_direction*0.01)+x_pos, 
             label = fp_formatted_pvalues[hla], 
             size = 3.5)
}

# FALSE NEGATIVE DIFFERENCE PLOT - Flipped horizontally with corrected annotations
fn_diff_plot <- ggplot(fn_diff, aes(x = HLA_Status_ordered, y = Difference, fill = Difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "slateblue", "FALSE" = "firebrick"),
                    labels = c("TRUE" = "T1GRS higher", "FALSE" = "GRS2 higher"),
                    name = "Comparison") +
  labs(
    title = "Difference in False Negative Rate (GRS2 - T1GRS)",
    x = "HLA Status",
    y = "Difference in False Negative Rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(face = 'bold'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5) # Remove legend since we'll add our own annotations
    #plot.margin = margin(l = 50, r = 10, t = 10, b = 10, unit = "pt")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at 0
  # Add arrows and labels
  annotate("segment", y = -0.01, yend = -0.04, x = 0.05, xend = 0.05, 
           arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  annotate("segment", y = 0.01, yend = 0.04, x = 0.05, xend = 0.05, 
           arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  annotate("text", y = -0.025, x = 0.3, label = "T1GRS Worse", size = 3.5) +
  annotate("text", y = 0.025, x = 0.3, label = "T1GRS Better", size = 3.5) +
  coord_flip() +  ylim(-0.05,0.1) + #ylim(-0.1,10)# Flip coordinates to make bars horizontal 
  scale_x_discrete(limits = rev(levels(fn_diff$HLA_Status_ordered)),expand = expansion(add = c(1.5, 0)))  # Reverse order so largest diff at top

# Add FN p-value annotations - corrected for horizontal orientation
for (i in 1:length(hla_groups)) {
  hla <- hla_groups[i]
  if(is.na(fn_p_values[i]) || fn_formatted_pvalues[i] == "") next  # Skip empty p-values
  
  # Find position in the ordered factor
  y_pos <- which(levels(fn_diff$HLA_Status_ordered) == hla)
  
  # Calculate position for annotation
  diff_val <- fn_diff$Difference[fn_diff$HLA_Status == hla]
  if(is.na(diff_val)) next
  
  x_direction <- ifelse(diff_val > 0, 1.1, -1.1)
  x_pos <- diff_val #+ (abs(diff_val) * x_direction)
  
  # Add annotation
  fn_diff_plot <- fn_diff_plot + 
    annotate("text", x = 7-y_pos, y = x_pos+0.01, 
             label = fn_formatted_pvalues[hla], 
             size = 3.5)
}

# Arrange both plots vertically
library(gridExtra)
combined_plot <- grid.arrange(fn_diff_plot, fp_diff_plot, ncol = 1)

# Save the combined plot
ggsave('figs_may/fig3/FPFN_difference_horizontal.pdf', width = 7, height = 11, plot = combined_plot)


###################
### Sankey Plot ###
###################

# Sankey plot showing where incorrect and unidentified pts ended up!

# label pts as unidentified or incorrect
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)

# read in extra metadata
meta<-read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta)<-c('FID',	'IID',	'SEX',	'DISEASE',	'PC1',	'PC2',	'PC3',	'PC4',	'T1DGC',	'DCCT',	'GENIE_UK',	'GoKIND')

ALL <- merge(ALL,meta,by.x='CV_ID',by.y='FID',all.y = T)

df<-ALL
df$HLA_simple<-df$HLA_Status

x<-cutpointr(df,Total_sum,DISEASE,na.rm=T)$optimal_cutpoint
df$GRS2_Assignment<-ifelse(df$Total_sum>=x,2,1)

y<-cutpointr(df,Probability,DISEASE)$optimal_cutpoint
df$T1GRS_Assignment<-ifelse(df$Probability>=y,2,1)

df$Incorrect_Missing<-ifelse(df$GRS2_Assignment!=df$DISEASE,'GRS2 Incorrect','Correct')
#df$Incorrect_Missing[df$T1GRS_Assignment!=df$DISEASE]<-'T1GRS Incorrect'
#df$Incorrect_Missing[df$T1GRS_Assignment!=df$DISEASE&df$GRS2_Assignment!=df$DISEASE]<-'Both Incorrect'

df$Incorrect_Missing[is.na(df$GRS2_Assignment)]<-'GRS2 Unable\nto Calculate'

df$T1GRS<-'T1GRS'

df$T1GRS_outcome<-'TP'
df$T1GRS_outcome[df$T1GRS_Assignment==1&df$DISEASE==1]<-'TN'
df$T1GRS_outcome[df$T1GRS_Assignment==1&df$DISEASE==2]<-'FN'
df$T1GRS_outcome[df$T1GRS_Assignment==2&df$DISEASE==1]<-'FP'

df<-df[df$Incorrect_Missing!='Correct',]
df$total_recovered<-ifelse(df$T1GRS_outcome=='TP'|df$T1GRS_outcome=='TN','Patient\nRecovered','Incorrect\nClassification')
df$Incorrect_Missing<-factor(df$Incorrect_Missing,levels=c('Both Incorrect','T1GRS Incorrect','GRS2 Incorrect','GRS2 Unable\nto Calculate'))

# sankey plot for df, with columns Incorrect_Missing -> T1GRS -> T1GRS_outcome. 
# color incorrect_mising levels red and grey, T1GRS level slateblue, T1GRS_outcome 4 other colors that make sense

# Load required packages
library(ggplot2)
library(ggalluvial)
library(dplyr)

# Prepare data for the Sankey diagram
# Convert from wide to long format (assuming df is your dataset)

sankey_data <- df#[,c('Incorrect_Missing','T1GRS','T1GRS_outcome')]
# Define custom colors
custom_colors <- c(
  # Colors for Incorrect_Missing levels
  "Both Incorrect" = "#AA78AA",
  "T1GRS Incorrect" = "#B0A8E2",
  "GRS2 Incorrect" = "#FF6666",  # Red
  "GRS2 Unable\nto Calculate" = "#AAAAAA",    # Grey
  
  # Color for T1GRS level
  "T1GRS" = "#B0A8E2",  # Slateblue (replace with your actual T1GRS levels)

  "TP" = "#7FB3D5",  # Soft blue-teal
  "TN" = "#76D7C4",  # Soft teal-green
  "FP" = "#F7DC6F",  # Soft yellow
  "FN" = "#F5B041",   # Soft orange
  
  "Incorrect\nClassification" = "#FFD791",
  "Patient\nRecovered" = "#A7C4C2"
)

# Create the Sankey plot
ggplot(sankey_data, 
       aes(axis1 = Incorrect_Missing, axis2=T1GRS,axis3 = T1GRS_outcome,axis4=total_recovered)) +
  geom_alluvium(
    aes(fill = Incorrect_Missing),
    width = 0.45, knot.pos = 0.3,spline_degree = 1#color = "grey88",alpha=0.5
  ) +
  geom_stratum(
    width = 0.45, 
    aes(fill = after_stat(stratum), 
        alpha = ifelse(grepl("Incorrect|Missing", after_stat(stratum)), 1, 
                       ifelse(grepl("Correct|Present", after_stat(stratum)), 0.7, 1))),
    color = "black",  # Change from grey to black for more defined borders
    size = 0.5,gap = 0.05
  ) +
  geom_text(
    stat = "stratum", 
    aes(label = paste0(after_stat(stratum), "\n(n=", after_stat(count), ")"),
        alpha = ifelse(grepl("Incorrect|Missing", after_stat(stratum)), 1, 0.8)), 
    size = 3,
  ) +
  scale_x_discrete(limits = c("Initial Status",'Re-Run through T1GRS', "T1GRS Final Classification",'Final Status')) +
  scale_fill_manual(values = custom_colors) +
  scale_alpha_identity() +
  labs(title = "Patients Recovered by T1GRS",
       x = NULL, y = "Count") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5, 5, 5, 5)
  )

# Save the plot
ggsave('figs_april/fig3/T1GRS_sankey.pdf', width = 8, height = 5)


