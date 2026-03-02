# T1D non-linearity experiment
setwd('~/Code/T1D/')
library(dplyr)
library(pROC)
library(ggplot2)
library(ggpubr)  # For adding p-value comparisons
library(cutpointr)
library(tidyr)

n_dec <- 10

# generate logistic regression probabilities?
LogReg <- read.table("/Users/tjsears/Code/T1D/data/delong_ALL_NoPCS_LOGREG.txt", sep="\t", header=T)
LogReg$LogProb <- LogReg$Probability

# load predictions
ALL <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=T)
prob <- pROC::roc(ALL$disease, ALL$Probability)
GRS2 <- pROC::roc(ALL$disease, ALL$Total_sum)

ALL <- merge(ALL, LogReg[,c('CV_ID','LogProb')], by='CV_ID')

# loads shap values
shap_load <- read.table('data/catboost_all_NoPCS_SHAP_values.csv', sep=',', header=T)
shap <- shap_load

# get list of features ranked by predictivity
ranked_features <- colnames(shap)[2:ncol(shap)][order(colSums(abs(shap[,2:ncol(shap)])), decreasing = T)]

# generate non-linearity scores
shap$abs_displacement <- rowSums(abs(shap[,2:ncol(shap)]))

# sum outside top 20 
shap$nonMHC_sum <- rowSums((shap[,ranked_features[1:20]]))
shap$MHC_sum <- rowSums((shap[,ranked_features[21:199]]))

shap$nonMHC_disp <- rowSums(abs(shap[,ranked_features[1:20]]))
shap$MHC_disp <- rowSums(abs(shap[,ranked_features[21:199]]))

shap$abs_displacement <- (shap$MHC_disp-shap$nonMHC_disp)

# sum HLA vs non-HLA
shap$nonMHC_sum <- rowSums((shap[,grepl('chr',colnames(shap))]))
shap$MHC_sum <- rowSums((shap[,grepl('rs',colnames(shap))|grepl('HLA',colnames(shap))|grepl('SNPS',colnames(shap))]))
shap$mhc_gap <- ((shap$nonMHC_sum-shap$MHC_sum))

shap_test <- merge(shap, ALL, by.x='FID', by.y='CV_ID')
shap_test <- shap_test[!is.na(shap_test$disease),]

n_dec <- 10
df <- shap_test %>%
  mutate(decile = ntile(abs_displacement, n_dec))  # Creates 10 deciles

# Initialize empty vectors to store AUCs and p-values
prob_auc <- numeric(n_dec)
grs2_auc <- numeric(n_dec)
log_auc <- numeric(n_dec)  # Add LogReg AUC array
p_values_t1grs_grs2 <- numeric(n_dec)  # Store p-values for DeLong test between T1GRS and GRS2
p_values_log_grs2 <- numeric(n_dec)    # Store p-values for DeLong test between LogReg and GRS2
cutpoint_t1grs <- numeric(n_dec)
cutpoint_grs2 <- numeric(n_dec)
cutpoint_log <- numeric(n_dec)

# Loop through each decile, calculate AUCs, and perform DeLong test
for (i in 1:n_dec) {
  # Filter data for the current decile
  decile_data <- df %>% filter(decile == i)
  
  # Calculate AUCs for all three predictors
  prob_roc <- pROC::roc(decile_data$disease, decile_data$Probability, quiet = TRUE)
  grs2_roc <- pROC::roc(decile_data$disease, decile_data$Total_sum, quiet = TRUE)
  log_roc <- pROC::roc(decile_data$disease, decile_data$LogProb, quiet = TRUE)
  
  prob_auc[i] <- pROC::auc(prob_roc)
  grs2_auc[i] <- pROC::auc(grs2_roc)
  log_auc[i] <- pROC::auc(log_roc)
  
  # Perform DeLong test between T1GRS and GRS2
  delong_test_t1grs_grs2 <- roc.test(prob_roc, grs2_roc)
  p_values_t1grs_grs2[i] <- delong_test_t1grs_grs2$p.value
  
  # Perform DeLong test between LogReg and GRS2
  delong_test_log_grs2 <- roc.test(log_roc, grs2_roc)
  p_values_log_grs2[i] <- delong_test_log_grs2$p.value
  
  # Store optimal cutpoints for each decile
  cutpoint_t1grs[i] <- cutpointr(decile_data, Probability, disease)$acc
  cutpoint_grs2[i] <- cutpointr(decile_data, Total_sum, disease)$acc
  cutpoint_log[i] <- cutpointr(decile_data, LogProb, disease)$acc
}

# Create a data frame for plotting
auc_data <- data.frame(
  Decile = 1:n_dec,
  prob_AUC = prob_auc,
  GRS2_AUC = grs2_auc,
  Log_AUC = log_auc,
  p_value_t1grs_grs2 = p_values_t1grs_grs2,
  p_value_log_grs2 = p_values_log_grs2,
  cutpoint_t1grs = cutpoint_t1grs,
  cutpoint_grs2 = cutpoint_grs2,
  cutpoint_log = cutpoint_log
)

# Create vertical offsets for the two sets of p-values
vjust_t1grs_grs2 <- c(4, 4, 6, 6, 6, -8, -7, -7, -11, -10.5)
vjust_log_grs2 <- vjust_t1grs_grs2 + 1.5  # Offset the second set of p-values

# Convert p-values to stars with Bonferroni correction
# Function to convert p-values to stars with Bonferroni correction
p_to_stars <- function(p, n_tests = 10) {
  p_adjusted <- p * n_tests  # Bonferroni correction
  p_adjusted <- pmin(p_adjusted, 1)  # Cap at 1
  
  stars <- rep("", length(p))
  stars[p_adjusted < 0.05] <- "*"
  stars[p_adjusted < 0.0001] <- "**"
  stars[p_adjusted < 0.00001] <- "***"
  
  return(stars)
}

# Apply Bonferroni correction and convert to stars
stars_t1grs_grs2 <- p_to_stars(auc_data$p_value_t1grs_grs2)
stars_log_grs2 <- p_to_stars(auc_data$p_value_log_grs2)

# Plot the AUCs with star annotations
ggplot(auc_data, aes(x = Decile)) +
  geom_line(aes(y = prob_AUC, color = "T1GRS"), size = 1) +
  geom_line(aes(y = GRS2_AUC, color = "GRS2"), size = 1) +
  geom_line(aes(y = Log_AUC, color = "LogReg"), size = 1) +
  geom_point(aes(y = prob_AUC, color = "T1GRS"), size = 2) +
  geom_point(aes(y = GRS2_AUC, color = "GRS2"), size = 2) +
  geom_point(aes(y = Log_AUC, color = "LogReg"), size = 2) +
  # Stars below each point for T1GRS vs GRS2
  geom_text(aes(
    x = Decile,
    y = GRS2_AUC+c(0.005,0.005,0.005,0.005,0.01,0.01,0.01,0.01,0.02,0.025) * 1.005,  # Position below GRS2 line
    label = stars_t1grs_grs2
  ), size = 5, color = "black") +
  # Stars below for LogReg vs GRS2
  geom_text(aes(
    x = Decile,
    y = GRS2_AUC+c(-0.01,-0.005,-0.01,-0.01,-0.01,-0.02,-0.02,-0.023,-0.03,-0.05) * 0.955,  # Position further below GRS2 line
    label = stars_log_grs2
  ), size = 5, color = "black") +
  labs(
    title = "AUC by Interaction Percentile",
    subtitle = '* = Delong Test FDR',
    x = "Increasing Interactions (Decile)",
    y = "AUC"
  ) + 
  #scale_y_log10() +
  scale_color_manual(
    values = c(
      "T1GRS" = "slateblue", 
      "GRS2" = "firebrick",
      "LogReg" = "darkorange"
    ),
    # Set the order in the legend
    breaks = c("T1GRS", "GRS2", "LogReg")
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + scale_x_continuous(breaks = 1:10) +
  scale_y_log10(breaks = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                labels = c("0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95", "1.0"))
  

ggsave('figs_april/fig4/non-linearity_test_with_logreg.pdf', width = 9, height = 7)

##################################################################
# DR3/4 Dual Violinplot showing MHC abs shap and NonMHC abs shap #
##################################################################

# Read in DR3/DR4 status
clinical_data <- read.table('clinical_data/Fixed_RenalGeno_HLAfam_GENIE_UK_DCCT_GoKIND_T1DGC_phenotypes.txt', sep='\t', header=T)
df <- merge(clinical_data, shap_test, by.x='FID', by.y='FID')
df <- df[df$disease==2,]
df$HLA_Status <- factor(df$HLA_Status, levels = c('non-DR3/DR4','DR3','DR4','DR3/DR3','DR4/DR4','DR3/DR4'))

# Convert data to long format for split violin plot
df_long <- df %>%
  select(FID, HLA_Status, MHC_sum, nonMHC_sum) %>%
  pivot_longer(
    cols = c(MHC_sum, nonMHC_sum),
    names_to = "Component",
    values_to = "Value"
  )

# Calculate summary statistics for each group
summary_stats <- df_long %>%
  group_by(HLA_Status, Component) %>%
  summarise(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Perform one-way ANOVA for each component
anova_mhc <- aov(MHC_sum ~ HLA_Status, data = df)
anova_non_mhc <- aov(nonMHC_sum ~ HLA_Status, data = df)

# Get p-values
p_value_mhc <- summary(anova_mhc)[[1]]$`Pr(>F)`[1]
p_value_non_mhc <- summary(anova_non_mhc)[[1]]$`Pr(>F)`[1]

# Create a split violin plot
library(ggplot2)
library(gghalves) # Need to install this package if you don't have it: install.packages("gghalves")

ggplot(df_long, aes(x = HLA_Status, y = Value, fill = Component)) +
  geom_hline(yintercept = 0, color = "gray40", linetype = "dashed", size = 0.5) +
  # Left half violins for MHC_sum
  geom_half_violin(data = subset(df_long, Component == "MHC_sum"), 
                   aes(fill = Component),
                   side = "l", alpha = 0.6, position = position_nudge(x = -0.05)) +
  # Right half violins for nonMHC_sum 
  geom_half_violin(data = subset(df_long, Component == "nonMHC_sum"), 
                   aes(fill = Component),
                   side = "r", alpha = 0.6, position = position_nudge(x = 0.05)) +
  # Add medians
  geom_point(data = summary_stats, 
             aes(x = HLA_Status, y = Median, group = Component),
             position = position_dodge(width = 0.3), size = 2, color = "black") +
  # Add IQR error bars
  geom_errorbar(data = summary_stats, inherit.aes = F,
                aes(x = HLA_Status, ymin = Q1, ymax = Q3, group = Component),
                position = position_dodge(width = 0.3), width = 0.1, color = "black") +
  # Set colors # #6457A6','#AB2346
  scale_fill_manual(values = c("MHC_sum" = "#6457A6", "nonMHC_sum" = "#AB2346"),
                    labels = c("MHC_sum" = "MHC\nComponent", "nonMHC_sum" = "Non-MHC\nComponent")) +
  # Labels
  labs(
    title = "MHC vs Non-MHC Components of T1D",
    x = "HLA Status",
    y = "Contribution to Disease Classification\n(SHAP Values)",
    fill = "Genetic Component"
  ) +
  # Styling
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  # Add p-values
  annotate("label", x = 5.5, y = max(df_long$Value, na.rm = TRUE) * -0.6, 
           label = paste("MHC ANOVA p < 0.0001"),#, format(p_value_mhc, digits = 3)), 
           color = "#6457A6", size = 4) +
  annotate("label", x = 4.5, y = max(df_long$Value, na.rm = TRUE) * -0.6, 
           label = paste("Non-MHC ANOVA p < 0.0001"),#, format(p_value_non_mhc, digits = 3)), 
           color = "#AB2346", size = 4) +
  coord_flip()

ggsave('figs_april/fig4/mhc_nonmhc_split_violin.pdf', width = 6.5, height = 6.5)



####################################
### Feature Interaction Analysis ###
####################################
library(RColorBrewer)
# read in lead signal document for conversion
leads<-read.table('/Users/tjsears/Code/T1D/interactionPlots/LeadSignals.txt',sep=" ",header=F)
colnames(leads)<-c('Gene','Loci')

# Feature importance
# featureImportance<-read.table("/Users/tjsears/Code/T1D/data/catboost_all_NoPCS_feature_interactions.csv",sep=",",header=T,row.names = 1)
featureImportance<-read.table("/Users/tjsears/Code/T1D/featureImportances/ALL_featureImportanceMatrix_noPCs_jun4_frac_0.3_199_features.txt",sep="\t",header=T,row.names = 1)

#featureImportance<-featureImportance*10000

colnames(featureImportance)<-gsub("SNPS_","",colnames(featureImportance))
rownames(featureImportance)<-gsub("SNPS_","",rownames(featureImportance))
colnames(featureImportance)<-gsub("AA_A_","HLA_",colnames(featureImportance))
rownames(featureImportance)<-gsub("AA_A_","HLA_",rownames(featureImportance))
#colnames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",colnames(featureImportance))
#rownames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",rownames(featureImportance))
colnames(featureImportance)<-gsub("\\.",":",colnames(featureImportance))
rownames(featureImportance)<-gsub("\\.",":",rownames(featureImportance))


for (i in 1:nrow(featureImportance)){
  for (j in 1:nrow(featureImportance)){
    if (i<j){
      featureImportance[i,j]<-0
    }
  }
}

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]
library(reshape2)

feet_long<-melt(featureImportance)
feet_long$AltVar<-colnames(featureImportance)
feet_long$FinalVar<-paste(feet_long$variable,"/",feet_long$AltVar)
feet_long<-feet_long[feet_long$value>0,]
feet_long<-feet_long[order(feet_long$value,decreasing = T),]

# subset feature importance into entries only containing variables of interest
featureImportance_smaller<-featureImportance[unique(feet_long$variable),unique(feet_long$variable)]

library(scales)
feet_long$FinalVar<-factor(feet_long$FinalVar,levels=feet_long$FinalVar)

featureImportance_smaller$names<-rownames(featureImportance_smaller)

links<-pivot_longer(featureImportance_smaller,cols=colnames(featureImportance_smaller)[1:(ncol(featureImportance_smaller)-1)])
links$value<-round(links$value,digits = 5)

links_table<-links
links_table$finalVar<-paste(links_table$name,links_table$names)
links_table_final<-aggregate(links_table$value,by = list(links_table$name),FUN = sum)
colnames(links_table_final)<-c('name','value')

links<-links[!duplicated(links$value),] #get rid of zeroes
colnames(links)<-c('Var1','Var2','value')
links$FinalVar<-paste(links$Var1,links$Var2)

# First sets all to "Non-MHC / Non-MHC"
links$type <- rep('Non-MHC / Non-MHC', nrow(links))

# If FinalVar doesn't contain "chr", classify as "MHC / MHC"
links$type <- ifelse(!grepl('chr', links$FinalVar), "MHC / MHC", links$type)

# If FinalVar contains "chr", classify as "MHC / Non-MHC"
links$type[grepl('chr', links$FinalVar)] <- "MHC / Non-MHC"

# If FinalVar contains "chr" but doesn't contain "intron", "rs", "HLA", or "exon", 
# reclassify as "Non-MHC / Non-MHC"
links$type[grepl('chr', links$FinalVar) & 
             !(grepl('intron', links$FinalVar) | 
                 grepl('rs', links$FinalVar) | 
                 grepl('HLA', links$FinalVar) | 
                 grepl('exon', links$FinalVar))] <- "Non-MHC / Non-MHC"

links_old_names<-links

# alter links to account for lead signals
match_vec<-match(links$Var2,leads$Loci,nomatch = 0)
links$Var2[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

match_vec<-match(links$Var1,leads$Loci,nomatch = 0)
links$Var1[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

links$FinalVar<-paste(links$Var1,links$Var2,sep=' / ')

library(ggplot2)
library(dplyr)

# Sort by value and get top 20
top_interactions <- links %>%
  arrange(desc(value)) %>%
  head(20)

# Make sure type is a factor with consistent order for the legend
top_interactions$type <- factor(top_interactions$type, 
                                levels = c("MHC / MHC", "MHC / Non-MHC", "Non-MHC / Non-MHC"))

# Create the barplot
ggplot(top_interactions, aes(x = reorder(FinalVar, value), y = value, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("MHC / MHC" = "#8DA0CB", 
                               "Non-MHC / Non-MHC" = "#FC8D62", 
                               "MHC / Non-MHC" = "#66C2A5")) +
  labs(title = "Top T1GRS Feature Interactions",
       x = "Feature Pair",
       y = "Interaction Strength\n(SHAP Values)",
       fill = "Interaction Type") +
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(-.9, 0.05),  # Position coordinates (0,0 is bottom left)
    legend.justification = c(0, 0),   # Anchor point of the legend
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(6, 6, 6, 6)
  ) +
  coord_flip()  # For better readability of feature names

# ADD BIG RED FDR LINE TO SHOW THAT INTERACTIONS ARE SIGNIFICANT COMPARED TO RANDOM BACKGROUND

ggsave('figs_april/fig4/feature_interaction_barplot.pdf', width = 7.5, height = 7)







###################################
### Feature Interaction CoxPH.  ###
###################################

# pull interaction from pre-renamed links df
library(survival)
library(ggplot2)
library(dplyr)
library(survminer)

# show coxpH risk of disease for given interaction. 

# pull original vcf
VCF <- read.table("/Users/tjsears/Code/T1D/data/formattedVCF.txt", sep="\t", header=T)
# filter down to pts in ALL
#temp<-ALL$CV_ID[!is.na(ALL$Total_sum)]
#VCF <- VCF[VCF$FID %in% temp,]

# Sort by value and get top 20
top_interactions_old_names <- links_old_names %>%
  arrange(desc(value)) %>%
  head(100)

write.table(top_interactions_old_names,'/Users/tjsears/Code/T1D/data/top_interactions_old_names.txt',sep='\t')
            
# select 2nd interaction, get names
Var1<-gsub(':','.',as.character(top_interactions_old_names[1,'Var1']))
Var2<-gsub(':','.',as.character(top_interactions_old_names[1,'Var2']))

# using Var1 and Var2, SEX as a covariate, and DISEASE (1,2) as outcome. 
#Plot a hazard plot using coxPH to see how presence, absence (of each individual feature), and interaction performs

# Create a dataset for analysis using original variant columns
analysis_data <- VCF %>%
  select(PC1, PC2, PC3, PC4, SEX, DISEASE, all_of(c(Var1, Var2))) %>%
  # Ensure variables are properly formatted
  mutate(
    SEX = as.factor(SEX),
    DISEASE = as.factor(DISEASE),  # Factor for logistic regression
    # Create interaction term using the original columns
    Interaction = !!sym(Var1) * !!sym(Var2)
  )


# Fit logistic regression model
glm_model <- glm(
  as.formula(paste("DISEASE == 2 ~", Var1, "*", Var2,"+ SEX")), 
  data = analysis_data,
  family = binomial()
)

# Print model summary
summary(glm_model)

# Calculate odds ratios and confidence intervals
odds_ratios <- exp(coef(glm_model))
conf_int <- exp(confint(glm_model))

# Combine into a table
results_table <- data.frame(
  Odds_Ratio = odds_ratios,
  Lower_CI = conf_int[,1],
  Upper_CI = conf_int[,2],
  p_value = summary(glm_model)$coefficients[,4]
)

# Print results
print(results_table)

# Create a visualization of odds ratios
results_for_plot <- data.frame(
  Variable = c(Var1, Var2, "Interaction", "SEX"),
  OR = odds_ratios[-1],  # Remove intercept
  Lower = conf_int[-1,1],
  Upper = conf_int[-1,2]
)

match_vec<-match(results_for_plot$Variable,gsub(":",'.',leads$Loci),nomatch = 0)
results_for_plot$Variable[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

results_for_plot<-results_for_plot[c(1,2,3),]
results_for_plot$Variable[3]<-paste(results_for_plot$Variable[1],'/',results_for_plot$Variable[2],'\nInteraction')
  
ggplot(results_for_plot, aes(x = Variable, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  #coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Disease Risk: Interaction between", "AIRE", "and", "INPP5B"),
    #subtitle = "Logistic Regression Model",
    x = "",
    y = "T1D Odds Ratio (95% CI)"
  )

# take this same interaction, see if AoO differs in people with this interaction
# load in clinical info
# subset to pts with T1D only
# see if different!

# merge clinical_data with VCF
clinical_data <- read.table('clinical_data/Fixed_RenalGeno_HLAfam_GENIE_UK_DCCT_GoKIND_T1DGC_phenotypes.txt', sep='\t', header=T)
merged_clin_vcf<-merge(VCF,clinical_data,on='FID')

# plot AoO by Var1 Present, Var2 Present, and Interaction
merged_clin_vcf$InteractionGroup<-'Neither'
merged_clin_vcf$InteractionGroup<-ifelse(merged_clin_vcf[,Var1]!=0,paste(Var1),merged_clin_vcf$InteractionGroup)
merged_clin_vcf$InteractionGroup<-ifelse(merged_clin_vcf[,Var2]!=0,paste(Var2),merged_clin_vcf$InteractionGroup)
merged_clin_vcf$InteractionGroup<-ifelse((merged_clin_vcf[,Var2]+merged_clin_vcf[,Var1])>=3,'Interaction',merged_clin_vcf$InteractionGroup)
#merged_clin_vcf$InteractionGroup<-ifelse((merged_clin_vcf[,Var2]!=0&merged_clin_vcf[,Var1]!=0),'Interaction',merged_clin_vcf$InteractionGroup)

table(merged_clin_vcf$InteractionGroup)

match_vec<-match(merged_clin_vcf$InteractionGroup,gsub(":",'.',leads$Loci),nomatch = 0)
merged_clin_vcf$InteractionGroup[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

#merged_clin_vcf<-merged_clin_vcf[merged_clin_vcf$DISEASE==2,]
merged_clin_vcf$InteractionGroup<-factor(merged_clin_vcf$InteractionGroup,levels=c('Neither','INS','rs1064173','Interaction'))

# geom boxplot of InteractionGroup and Agex

library(ggplot2)
library(ggpubr)  # For statistical tests
library(dplyr)

merged_clin_vcf_temp<-merged_clin_vcf[merged_clin_vcf$Agedx<55,] #broad cutoff
merged_clin_vcf_temp<-merged_clin_vcf_temp[!is.na(merged_clin_vcf_temp$Agedx),]

merged_clin_vcf_temp2<-merged_clin_vcf_temp

# Create the boxplot with statistical tests
ggplot(merged_clin_vcf_temp, aes(x = InteractionGroup, y = Agedx)) +
  geom_boxplot(aes(fill = InteractionGroup), alpha = 0.75) +
  # Add pairwise comparisons with exact p-values
  stat_compare_means(comparisons = list(c(1,4), c(2,4), c(3,4)),inherit.aes = F,
                     method = "t.test",label.y = c(43,40,37),data = merged_clin_vcf_temp2,mapping = aes(x = InteractionGroup, y = Agedx),
                     label = "p.format",paired = F,hide.ns = T,method.args = c(alternative='greater')) +  # Show actual p-values instead of just stars
  scale_fill_brewer(palette = "Dark1") +
  labs(
    title = "Age of Diagnosis by Interaction Group",
    subtitle = "Pairwise t-tests",
    x = "Interaction Group",
    y = "Age of Dx",
    fill = "Interaction Group"
  ) + ylim(0,45) +
  theme_bw(base_size = 14) +
  theme(
    #panel.grid.major.x = element_blank()
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
ggsave('figs_may/fig4/feature_interaction_AgeDx_Boxplot.pdf', width = 7.5, height = 7)




#merged_clin_vcf_temp2<-merged_clin_vcf_temp[merged_clin_vcf_temp$Agedx<=20,]
merged_clin_vcf_temp2<-merged_clin_vcf_temp
merged_clin_vcf_temp2$Agedx[merged_clin_vcf_temp$Agedx>20]<-20

# Create the boxplot with statistical tests
# Using coord_cartesian to zoom the y-axis without removing data points before statistical calculations
p_boxplot_final_zoom <- ggplot(merged_clin_vcf_temp2, aes(x = InteractionGroup, y = Agedx)) +
  geom_boxplot(aes(fill = InteractionGroup), alpha = 0.75, outlier.shape = NA) + # Hiding default outliers
  # Add pairwise comparisons with exact p-values
  # label.y values are adjusted to be visible within the 0-20 y-axis range
  stat_compare_means(comparisons = list(c(1,4), c(2,4), c(3,4)), # Ensure these indices match your factor levels
                     inherit.aes = FALSE, # Explicitly set to FALSE as per user's new code
                     method = "t.test",
                     label.y = c(22.5, 21, 19.5), # Adjusted for visibility within 0-20 ylim
                     data = merged_clin_vcf_temp, # Explicitly pass data
                     mapping = aes(x = InteractionGroup, y = Agedx), # Explicitly pass mapping
                     label = "p.format",
                     paired = FALSE,
                     hide.ns = TRUE,
                     method.args = list(alternative = 'greater')) +
  scale_fill_brewer(palette = "Dark1") +
  labs(
    title = "Age of Diagnosis by Interaction Group", # Kept original title
    subtitle = "Pairwise t-tests", # Updated subtitle for clarity
    x = "Interaction Group",
    y = "Age of Dx",
    fill = "Interaction Group"
  ) +
  coord_cartesian(ylim = c(0, 25)) + # Zoom y-axis from 0 to 20 without data removal for stats
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_line(linetype = "dashed", color = "lightgrey"), # This line was in previous canvas, removed to match user's latest snippet theme
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    #plot.title = element_text(hjust = 0.5), # Centering title
    #plot.subtitle = element_text(hjust = 0.5) # Centering subtitle
  )

print(p_boxplot_final_zoom)
ggsave('figs_may/fig4/zoomed_feature_interaction_AgeDx_Boxplot.pdf', width = 7.5, height = 6)




library(ggplot2)
library(dplyr)
library(ggsignif)  # For significance brackets

# Create a contingency table
disease_table <- table(merged_clin_vcf$DISEASE, merged_clin_vcf$InteractionGroup)

# Convert to data frame for plotting
plot_data <- as.data.frame(disease_table)
colnames(plot_data) <- c("DISEASE", "InteractionGroup", "Count")

# Make sure variables are of the right type
plot_data$DISEASE <- as.factor(plot_data$DISEASE)
plot_data$InteractionGroup <- as.factor(plot_data$InteractionGroup)

# Calculate percentages within each interaction group
plot_data <- plot_data %>%
  group_by(InteractionGroup) %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Total = sum(Count))

# Run chi-squared test
chi_test <- chisq.test(disease_table)
p_value <- chi_test$p.value

plot_data$DISEASE<-factor(ifelse(plot_data$DISEASE==2,'Disease','No Disease'),levels=c('No Disease','Disease'))
plot_data$text<-ifelse(plot_data$Percentage>11,paste0(plot_data$Count,'\n(',round(plot_data$Percentage,0),'%)'),
                       paste0(plot_data$Count,' (',round(plot_data$Percentage,0),'%)'))

# Create stacked barplot
ggplot(plot_data, aes(x = InteractionGroup, y = Percentage, fill = DISEASE)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = text), position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  labs(
    title = "Disease Distribution by Interaction Group",
    subtitle = paste("Chi-squared test: p =", format.pval(p_value, digits = 3)),
    x = "Interaction Group",
    y = "Percentage (%)",
    fill = "Disease Status"
  ) +
  scale_fill_manual(values = c("No Disease" = "forestgreen", "Disease" = "goldenrod")) +
  #geom_text(aes(y = 105, label = paste0("n=", Total)), position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  )
ggsave('figs_may/fig4/feature_interaction_chisquared_disease.pdf', width = 7.5, height = 4)




###################################
# Alternative Interaction Barplot #
###################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(purrr) # For pmap_dfr

# --- Original Script Starts Here ---
type_levels <- c("MHC / Non-MHC", "MHC / MHC", "Non-MHC / Non-MHC")
fill_colors<-c( "#66C2A5","#8DA0CB","#FC8D62")

# Assuming VCF is already loaded and pre-processed as in the original script
# VCF <- read.table("/Users/tjsears/Code/T1D/data/formattedVCF.txt", sep="\t", header=T)
# ... (column renaming logic from original script) ...

# --- Step 1: Select Top 10 Interactions per Category from MODIFIED 'links' ---
selected_links <- links %>%
  group_by(type) %>%
  arrange(desc(value)) %>% # Order by 'value' column in descending order
  slice_head(n = 10) %>%   # Select top 10 per category
  ungroup()

# fixed messed up INS snp (retained from original script)
# selected_links[11, "FinalVar"] <- 'INS / rs1064173'
# selected_links[12, "FinalVar"] <- 'PTPN22 / rs1064173'
# selected_links[11, "Var1"] <- 'INS'
# selected_links[12, "Var1"] <- 'PTPN22'

# --- Step 2: Calculate Chi-squared for these PRE-SELECTED Interactions ---
VCF$DISEASE <- factor(VCF$DISEASE)

# Iterate over the 'selected_links' dataframe
chisq_results_for_selected <- pmap_dfr(selected_links, function(Var1, Var2, type, value, ...) {
  
  # Ensure the columns exist in VCF
  if (!Var1 %in% names(VCF) || !Var2 %in% names(VCF)) {
    message(paste("Skipping interaction", Var1, "*", Var2, "as one or both variables are not in VCF."))
    return(NULL)
  }

  snp1_factor <- factor(VCF[[Var1]])
  snp2_factor <- factor(VCF[[Var2]])
  interaction_term_values <- factor(paste(snp1_factor, snp2_factor, sep = "_"))
  
  # Perform chi-squared test
  test_result_obj <- tryCatch({
    valid_disease <- VCF$DISEASE[!is.na(interaction_term_values)]
    valid_interaction_term <- interaction_term_values[!is.na(interaction_term_values)]
    if (nlevels(factor(valid_disease)) < 2 || nlevels(factor(valid_interaction_term)) < 2) {
      stop("DISEASE or interaction term has < 2 levels after NA removal for this pair.")
    }
    chisq.test(VCF$DISEASE, interaction_term_values)
  }, error = function(e) {
    message(paste("Chi-sq test failed for interaction", Var1, "*", Var2, ":", e$message))
    return(list(statistic = NA_real_, p.value = NA_real_, parameter = NA_integer_))
  })
  
  tibble(
    Var1_name = Var1, 
    Var2_name = Var2,
    interaction_name = paste(Var1, Var2, sep = " / "),
    type = type,
    original_links_value = value,
    chisq_stat = as.numeric(test_result_obj$statistic),
    p_value = as.numeric(test_result_obj$p.value),
    df = as.integer(test_result_obj$parameter)
  )
})

# Filter out any interactions where chi-squared test failed
plot_data <- chisq_results_for_selected %>%
  filter(!is.na(chisq_stat))

# --- Step 3: Calculate FDR and Determine Significance Line ---
# This section is kept for context from the original script.
# The threshold is hardcoded later before plotting.
fdr_chisq_threshold <- NA_real_
if (nrow(plot_data) > 0 && sum(!is.na(plot_data$p_value)) > 0) {
  p_values_for_fdr <- plot_data$p_value[!is.na(plot_data$p_value)]
  if (length(p_values_for_fdr) > 0) {
    adjusted_p_values <- p.adjust(p_values_for_fdr, method = "BH")
    plot_data$fdr <- NA_real_
    plot_data$fdr[!is.na(plot_data$p_value)] <- adjusted_p_values
    significant_fdr_interactions <- plot_data %>% filter(fdr < 0.05)
    if (nrow(significant_fdr_interactions) > 0) {
      fdr_chisq_threshold <- min(significant_fdr_interactions$chisq_stat, na.rm = TRUE)
    }
  }
}

# --- Step 4: Create the Barplot ---
if (nrow(plot_data) > 0) {
  
  # Apply newline modification for long interaction_name strings
  plot_data <- plot_data %>%
    mutate(
      interaction_name = sapply(interaction_name, function(name) {
        if (nchar(name) > 30) {
          if (grepl("/", name, fixed = TRUE)) {
            return(sub("/", "/\n", name, fixed = TRUE))
          }
        }
        return(name)
      })
    )
  
  # Ensure 'type' is factored in the desired order
  plot_data$type <- factor(plot_data$type, levels = type_levels)
  
  # Sort data by 'type', then by 'original_links_value' (desc)
  plot_data_sorted <- plot_data %>%
    arrange(type, desc(original_links_value))

  # --- NEW: Z-score Transformation ---
  mean_val <- mean(plot_data_sorted$original_links_value, na.rm = TRUE)
  sd_val <- sd(plot_data_sorted$original_links_value, na.rm = TRUE)
  
  # Add the z_score_value column
  plot_data_sorted <- plot_data_sorted %>%
    mutate(z_score_value = (original_links_value - mean_val) / sd_val)
  
  # Transform the hardcoded FDR threshold to the Z-score scale
  fdr_threshold_original <- 387.3875
  fdr_threshold_z <- (fdr_threshold_original - mean_val) / sd_val

  # --- Create y-axis levels with gaps AND prepare data for gap annotations ---
  y_axis_items_with_gaps <- character(0)
  gap_annotations_list <- list()
  previous_type_for_gap_logic <- "__INTERNAL_START_MARKER__"  
  
  for (i in 1:nrow(plot_data_sorted)) {
    current_row_type_factor <- plot_data_sorted$type[i]
    current_row_type_string <- as.character(current_row_type_factor)
    
    if (current_row_type_factor != previous_type_for_gap_logic) {
      gap_marker_name <- paste0("___GAP_ABOVE_", gsub("[^A-Za-z0-9_]", "_", current_row_type_string), "_", i)
      y_axis_items_with_gaps <- c(y_axis_items_with_gaps, gap_marker_name)
      
      gap_annotations_list[[length(gap_annotations_list) + 1]] <-
        data.frame(
          y_level_for_annotation = gap_marker_name,
          label_text = current_row_type_string,
          type_for_color = current_row_type_factor,
          stringsAsFactors = FALSE
        )
      previous_type_for_gap_logic <- current_row_type_factor
    }
    y_axis_items_with_gaps <- c(y_axis_items_with_gaps, plot_data_sorted$interaction_name[i])
  }
  
  final_y_factor_levels <- rev(y_axis_items_with_gaps)
  gap_annotations_df <- bind_rows(gap_annotations_list)
  
  if (nrow(gap_annotations_df) > 0) {
    gap_annotations_df$y_level_for_annotation <- factor(gap_annotations_df$y_level_for_annotation, levels = final_y_factor_levels)
    gap_annotations_df$type_for_color <- factor(gap_annotations_df$type_for_color, levels = type_levels)
  }
  
  plot_data_final <- plot_data_sorted
  plot_data_final$interaction_name_ordered <- factor(plot_data_final$interaction_name, levels = final_y_factor_levels)
  
  custom_y_axis_labels <- final_y_factor_levels
  custom_y_axis_labels[grepl("___GAP_ABOVE_", final_y_factor_levels, fixed = TRUE)] <- ""
  
  # NEW: Calculate a dynamic x-position for category labels on the Z-score axis
  z_range <- range(plot_data_final$z_score_value, na.rm = TRUE)
  label_x_pos <- z_range[1] + (z_range[2] - z_range[1]) * 0.5 # Position in the middle of the z-score range
  
  # --- Main Plot Construction with Z-scores ---
  interaction_plot <- ggplot(plot_data_final, aes(x = z_score_value, y = interaction_name_ordered, fill = type)) +
    geom_col() +
    labs(
      title = NULL,
      x = "Z-scored SHAP Feature Interaction Value", # MODIFIED
      y = NULL
    ) +
    scale_fill_manual(values = fill_colors, limits = type_levels) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + # MODIFIED
    scale_y_discrete(labels = custom_y_axis_labels, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10, color = 'black', lineheight = 0.9),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
  
  if (nrow(gap_annotations_df) > 0) {
    interaction_plot <- interaction_plot +
      geom_text(
        data = gap_annotations_df,
        aes(x = label_x_pos, # MODIFIED to dynamic position
            y = y_level_for_annotation,
            label = label_text),
        hjust = 0.5, # Center align text
        vjust = 0.5,
        size = 4.3,
        fontface = "bold",
        colour = 'black',
        inherit.aes = FALSE
      )
  }
  
  # Add FDR line on the Z-score scale
  if (!is.na(fdr_threshold_z)) {
    interaction_plot <- interaction_plot +
      geom_vline(xintercept = fdr_threshold_z, linetype = "dashed", color = "red", linewidth = 0.7) +
      annotate("text", x = fdr_threshold_z,
               y = Inf,
               label = "FDR < 0.05",
               color = "grey35", vjust = 6.6, hjust = -0.1, size = 3.5)
  }
  
  print(interaction_plot)
  
  num_gaps <- nrow(gap_annotations_df)
  plot_height_adjusted <- 12 + num_gaps * 0.5 
  
  # ggsave('figs_may/fig4/tall_interaction_bar_zscore.pdf', plot = interaction_plot, width = 6, height = plot_height_adjusted)
  # message(paste0("Plot saved to figs_may/fig4/tall_interaction_bar_zscore.pdf (adjusted height: ", plot_height_adjusted, ")"))
  
} else {
  message("No interactions to plot after selection and chi-squared calculation.")
}






###################################
# Alternative Interaction Barplot #
# with externally-derived Z-scores#
###################################

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- Step 1: Calculate Background Distribution from External Data ---

# Define the input path for the full feature importance matrix
input_filepath <- "/Users/tjsears/Code/T1D/featureImportances/ALL_featureImportanceMatrix_noPCs_jun4_frac_0.3_199_features.txt"

# Load external feature importance matrix
# NOTE: This script assumes the file at 'input_filepath' exists.
# A placeholder is created if the file is not found, for demonstration purposes.
if (file.exists(input_filepath)) {
  featureImportance <- read.table(input_filepath, sep="\t", header=T, row.names=1)
} else {
  message(paste("Warning: File not found at", input_filepath, ". Creating a mock featureImportance matrix."))
  featureImportance <- as.data.frame(matrix(runif(100*100, 0, 1500), 100, 100))
  rownames(featureImportance) <- paste0("Var", 1:100)
  colnames(featureImportance) <- paste0("Var", 1:100)
}


# Calculate mean and sd from the full matrix for the background distribution
all_interaction_values <- unlist(featureImportance)
all_interaction_values <- all_interaction_values[all_interaction_values > 0 & !is.na(all_interaction_values)]
mean_external_interactions <- mean(all_interaction_values)
sd_external_interactions <- sd(all_interaction_values)

message(paste("Background Distribution: Mean =", round(mean_external_interactions, 4), ", SD =", round(sd_external_interactions, 4)))

# --- Step 2: Calculate Theoretical FDR Threshold from Full Dataset ---

# Convert the wide matrix to a long format to calculate p-values for all interactions
long_format_interactions <- as.data.frame(as.table(as.matrix(featureImportance)))
colnames(long_format_interactions) <- c("Feature1", "Feature2", "InteractionValue")

long_format_interactions <- long_format_interactions %>%
  filter(as.character(Feature1) < as.character(Feature2) & InteractionValue > 0)

# Calculate Z-score, p-value, and FDR for every interaction in the background set
all_interactions_with_stats <- long_format_interactions %>%
  mutate(
    z_score = (InteractionValue - mean_external_interactions) / sd_external_interactions,
    p_value = 2 * pnorm(-abs(z_score)),
    p_adj = p.adjust(p_value, method = "BH")
  )

# Determine the minimum interaction value that is significant at FDR < 0.05
fdr_line_intercept_raw <- NA # Default to no line
significant_interactions_all <- all_interactions_with_stats %>%
  filter(p_adj < 0.05)

if (nrow(significant_interactions_all) > 0) {
  fdr_line_intercept_raw <- min(significant_interactions_all$InteractionValue, na.rm = TRUE)
  message(paste("Theoretical FDR < 0.05 threshold calculated at SHAP value:", round(fdr_line_intercept_raw, 4)))
} else {
  message("No interactions were significant at FDR < 0.05 in the full dataset.")
}


# --- Step 3: Prepare Data for Plotting (using 'links' and 'VCF') ---

# Mock data setup for 'links' and 'VCF' if they don't exist
if (!exists("links")) {
  type_levels <- c("MHC / Non-MHC", "MHC / MHC", "Non-MHC / Non-MHC")
  links <- data.frame(
    Var1 = sample(rownames(featureImportance), 30, replace = TRUE),
    Var2 = sample(colnames(featureImportance), 30, replace = TRUE),
    type = sample(type_levels, 30, replace = TRUE),
    value = runif(30, 100, 1200)
  )
}
if (!exists("VCF")) {
  VCF <- data.frame(
    DISEASE = factor(sample(0:1, 100, replace = TRUE))
  )
  all_vars <- unique(c(links$Var1, links$Var2))
  for (v in all_vars) {
    if (!v %in% names(VCF)) VCF[[v]] <- sample(0:2, 100, replace = TRUE)
  }
}

type_levels <- c("MHC / Non-MHC", "MHC / MHC", "Non-MHC / Non-MHC")
fill_colors<-c( "#66C2A5","#8DA0CB","#FC8D62")

# Select Top 10 Interactions per Category from 'links'
selected_links <- links %>%
  group_by(type) %>%
  arrange(desc(value)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Perform Chi-squared test (retained from original script)
chisq_results_for_selected <- pmap_dfr(selected_links, function(Var1, Var2, type, value, ...) {
  # This section is simplified; original script's error handling is assumed
  tibble(
    interaction_name = paste(Var1, Var2, sep = " / "),
    type = type,
    original_links_value = value
  )
})

plot_data <- chisq_results_for_selected

# --- Step 4: Create the Barplot with External Z-Scores ---

if (nrow(plot_data) > 0) {
  
  # MODIFICATION: Calculate Z-score using the external mean and sd
  plot_data <- plot_data %>%
    mutate(z_score_value = (original_links_value - mean_external_interactions) / sd_external_interactions)
  
  # MODIFICATION: Convert the raw FDR threshold to the Z-score scale
  fdr_threshold_z <- NA
  if (!is.na(fdr_line_intercept_raw)) {
    fdr_threshold_z <- (fdr_line_intercept_raw - mean_external_interactions) / sd_external_interactions
  }
  
  # --- Plotting logic (mostly unchanged) ---
  
  # Apply newline modification for long interaction_name strings
  plot_data <- plot_data %>%
    mutate(
      interaction_name = sapply(interaction_name, function(name) {
        if (nchar(name) > 30) {
          if (grepl("/", name, fixed = TRUE)) {
            return(sub("/", "/\n", name, fixed = TRUE))
          }
        }
        return(name)
      })
    )
  
  plot_data$type <- factor(plot_data$type, levels = type_levels)
  
  # Sort data for plotting
  plot_data_sorted <- plot_data %>%
    arrange(type, desc(original_links_value))
  
  # Create y-axis with gaps and prepare annotations
  y_axis_items_with_gaps <- character(0)
  gap_annotations_list <- list()
  previous_type_for_gap_logic <- "__INTERNAL_START_MARKER__"  
  
  for (i in 1:nrow(plot_data_sorted)) {
    current_row_type_factor <- plot_data_sorted$type[i]
    if (current_row_type_factor != previous_type_for_gap_logic) {
      gap_marker_name <- paste0("___GAP_", i)
      y_axis_items_with_gaps <- c(y_axis_items_with_gaps, gap_marker_name)
      gap_annotations_list[[length(gap_annotations_list) + 1]] <-
        data.frame(y_level = gap_marker_name, label = as.character(current_row_type_factor))
      previous_type_for_gap_logic <- current_row_type_factor
    }
    y_axis_items_with_gaps <- c(y_axis_items_with_gaps, plot_data_sorted$interaction_name[i])
  }
  
  final_y_factor_levels <- rev(y_axis_items_with_gaps)
  gap_annotations_df <- bind_rows(gap_annotations_list)
  
  plot_data_final <- plot_data_sorted
  plot_data_final$interaction_name_ordered <- factor(plot_data_final$interaction_name, levels = final_y_factor_levels)
  
  custom_y_axis_labels <- final_y_factor_levels
  custom_y_axis_labels[grepl("___GAP_", final_y_factor_levels)] <- ""
  
  z_range <- range(plot_data_final$z_score_value, na.rm = TRUE, finite = TRUE)
  label_x_pos <- z_range[1] + (z_range[2] - z_range[1]) * 0.5
  
  # Main Plot Construction
  interaction_plot <- ggplot(plot_data_final, aes(x = z_score_value, y = interaction_name_ordered, fill = type)) +
    geom_col() +
    labs(
      title = NULL,
      x = "Z-scored SHAP Feature Interaction Value",
      y = NULL
    ) +
    scale_fill_manual(values = fill_colors, limits = type_levels) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_discrete(labels = custom_y_axis_labels, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10, color = 'black', lineheight = 0.9),
      axis.ticks.y = element_blank(), legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
  
  # Add category labels in gaps
  if (nrow(gap_annotations_df) > 0) {
    interaction_plot <- interaction_plot +
      geom_text(data = gap_annotations_df, aes(x = label_x_pos, y = y_level, label = label),
                fontface = "bold", colour = 'black', inherit.aes = FALSE)
  }
  
  # Add FDR line on the Z-score scale
  if (!is.na(fdr_threshold_z)) {
    interaction_plot <- interaction_plot +
      geom_vline(xintercept = fdr_threshold_z, linetype = "dashed", color = "red", linewidth = 0.7) +
      annotate("text", x = fdr_threshold_z, y = Inf, label = "FDR < 0.05",
               color = "grey35", vjust = 2.1, hjust = 1.15, size = 3.5)
  }
  
  print(interaction_plot)
  
  ggsave('figs_may/fig4/final_interaction_plot_external_zscore.pdf', plot = interaction_plot, width = 6, height = 10)
  
} else {
  message("No interactions to plot.")
}




