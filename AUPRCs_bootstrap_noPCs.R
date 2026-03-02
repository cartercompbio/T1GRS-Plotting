# Bootstrap CI AUPRS

# Replacement panels for Fig3
setwd('~/Code/T1D/') # Make sure this path is correct for your environment

# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)
library(PRROC) # Added for Precision-Recall curves


#############
# Discovery #
#############

# --- Full Model ---

# Read in your data
ALL_full <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_NoPCS_catboost.txt", sep="\t", header=TRUE)

# read in extra metadata
meta_full <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_full) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')

ALL_full <- merge(ALL_full, meta_full, by.x='CV_ID', by.y='FID', all.y = TRUE) # Retain all meta rows

# Ensure the outcome is numeric and handle potential NAs from merge or in predictors
ALL_full$DISEASE <- as.numeric(as.character(ALL_full$DISEASE)) # Ensure it's numeric 0/1
# Assuming DISEASE might be 1 for control, 2 for case from PLINK format, convert to 0 and 1
# Or if it's already 0/1, this is fine. If it's 1/2, adjust accordingly:
# Example: ALL_full$DISEASE[ALL_full$DISEASE == 1] <- 0 # Control
# ALL_full$DISEASE[ALL_full$DISEASE == 2] <- 1 # Case
# For this script, we'll assume your 'DISEASE' column after merge is 0 for controls and 1 for cases.
# If it's different (e.g., 1 for controls, 2 for cases), you MUST adjust this.
# For pROC, it's good practice to make it a factor for clarity if levels are not 0/1
# ALL_full$DISEASE_factor <- factor(ALL_full$DISEASE, levels = c(0, 1), labels = c("Control", "Case"))

# Remove rows with NA in predictors or outcome, as they can't be used
ALL_full <- ALL_full[complete.cases(ALL_full$DISEASE, ALL_full$Probability, ALL_full$Total_sum), ]


# ROC Curves
prob_roc_full <- pROC::roc(response = ALL_full$DISEASE, predictor = ALL_full$Probability, quiet = TRUE)
grs2_roc_full <- pROC::roc(response = ALL_full$DISEASE, predictor = ALL_full$Total_sum, quiet = TRUE)

prob_auc_full <- pROC::auc(prob_roc_full)
grs2_auc_full <- pROC::auc(grs2_roc_full)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc_full)
grs2_ci <- pROC::ci.auc(grs2_roc_full)
print(prob_ci)
print(grs2_ci)

delong_test_full <- pROC::roc.test(prob_roc_full, grs2_roc_full, method = "delong")
print("DeLong Test for Full Model ROC:")
print(delong_test_full)

prob_df_full <- data.frame(FPR = 1 - prob_roc_full$specificities, TPR = prob_roc_full$sensitivities, Model = "T1GRS")
grs2_df_full <- data.frame(FPR = 1 - grs2_roc_full$specificities, TPR = grs2_roc_full$sensitivities, Model = "GRS2")
roc_df_full <- bind_rows(prob_df_full, grs2_df_full)

# Plot ROC for Full Model
plot_roc_full <- ggplot(roc_df_full, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        axis.line.x = element_line(color = "black", linewidth = 0.5),
                                        axis.line.y = element_line(color = "black", linewidth = 0.5)) +
  labs(title="Discovery: Full Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_full, 3)), color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_full, 3)), color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_full$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_full$p.value, 4)))), color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
print(plot_roc_full)
#ggsave('figs_may/fig3/auc_full.pdf', plot = plot_roc_full, width = 6, height = 5)



t1grs_data <- ALL_full %>% filter(!is.na(Probability))
grs2_data <- ALL_full %>% filter(!is.na(Total_sum))

# --- BOOTSTRAP FOR AUPRC CONFIDENCE INTERVALS ---
n_bootstraps <- 100
t1grs_auprc_boot <- numeric(n_bootstraps)
grs2_auprc_boot <- numeric(n_bootstraps)

cat(paste("Starting bootstrapping for", n_bootstraps, "iterations...\n"))

for (i in 1:n_bootstraps) {
  # --- T1GRS Bootstrap ---
  t1grs_boot_sample <- t1grs_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(t1grs_boot_sample$DISEASE)) == 2) {
    pr_t1grs_boot <- pr.curve(scores.class0 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 2],
                              scores.class1 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 1],
                              curve = FALSE)
    t1grs_auprc_boot[i] <- pr_t1grs_boot$auc.integral
  } else {
    t1grs_auprc_boot[i] <- NA # Mark as invalid if only one class is present
  }
  
  # --- GRS2 Bootstrap ---
  grs2_boot_sample <- grs2_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(grs2_boot_sample$DISEASE)) == 2) {
    pr_grs2_boot <- pr.curve(scores.class0 = grs2_boot_sample$Total_sum[grs2_boot_sample$DISEASE == 2],
                             scores.class1 = grs2_boot_sample$Total_sum[grs2_boot_sample$DISEASE == 1],
                             curve = FALSE)
    grs2_auprc_boot[i] <- pr_grs2_boot$auc.integral
  } else {
    grs2_auprc_boot[i] <- NA
  }
}

cat("Bootstrap complete.\n")

# Calculate the 95% confidence intervals from the bootstrap results
t1grs_ci <- quantile(t1grs_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
grs2_ci <- quantile(grs2_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
# --- END BOOTSTRAP SECTION ---


# --- AUPRC CURVE PLOTTING (using original NA-filtered datasets) ---

# Calculate PR curves on the clean, NA-filtered datasets
pr_t1grs <- pr.curve(scores.class0 = t1grs_data$Probability[t1grs_data$DISEASE == 2], 
                     scores.class1 = t1grs_data$Probability[t1grs_data$DISEASE == 1], curve = TRUE)
auprc_t1grs <- pr_t1grs$auc.integral

pr_grs2 <- pr.curve(scores.class0 = grs2_data$Total_sum[grs2_data$DISEASE == 2], 
                    scores.class1 = grs2_data$Total_sum[grs2_data$DISEASE == 1], curve = TRUE)
auprc_grs2 <- pr_grs2$auc.integral

# Create formatted text strings with CIs for the plot labels
t1grs_label_pr <- paste0("T1GRS AUPRC = ", round(auprc_t1grs, 3), " (", round(t1grs_ci[1], 3), "-", round(t1grs_ci[2], 3), ")")
grs2_label_pr <- paste0("GRS2 AUPRC = ", round(auprc_grs2, 3), " (", round(grs2_ci[1], 3), "-", round(grs2_ci[2], 3), ")")

# Print results to console
print(t1grs_label_pr)
print(grs2_label_pr)

# Create data frames for ggplot
pr_t1grs_df <- data.frame(Recall = pr_t1grs$curve[, 1], Precision = pr_t1grs$curve[, 2], Model = "T1GRS")
pr_grs2_df <- data.frame(Recall = pr_grs2$curve[, 1], Precision = pr_grs2$curve[, 2], Model = "GRS2")
pr_df <- bind_rows(pr_t1grs_df, pr_grs2_df)

# Calculate baseline for PR curve using the GRS2 data (or choose one consistently)
pr_baseline <- sum(grs2_data$T1D == 2) / nrow(grs2_data)

# Plot PR curves with ggplot
pr_plot_title <- "All of Us: Full Model\nAUPRC Curve"
plot_pr <- ggplot(pr_df, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = pr_baseline, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 16) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)) +
  labs(title = pr_plot_title,
       x = "Recall",
       y = "Precision",
       color = "Model") +
  ylim(0, 1) +
  annotate("text", x = 0.65, y = 0.25, label = t1grs_label_pr, color = "slateblue", size = 4) +
  annotate("text", x = 0.65, y = 0.20, label = grs2_label_pr, color = "firebrick", size = 4) +
  annotate("text", x = 0.65, y = 0.15, label = paste0("Baseline = ", round(pr_baseline, 3)), color = "grey50", size = 4) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))

print(plot_pr)







# --- Non-MHC Model ---
ALL_nonhla <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_NonHLA_NoPCS_catboost.txt", sep="\t", header=TRUE)
meta_nonhla <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_nonhla) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')
ALL_nonhla <- merge(ALL_nonhla, meta_nonhla, by.x='CV_ID', by.y='FID', all.y = TRUE)
ALL_nonhla$DISEASE <- as.numeric(as.character(ALL_nonhla$DISEASE))
ALL_nonhla <- ALL_nonhla[complete.cases(ALL_nonhla$DISEASE, ALL_nonhla$Probability, ALL_nonhla$sum_All5_nonHLA), ]

prob_roc_nonhla <- pROC::roc(response = ALL_nonhla$DISEASE, predictor = ALL_nonhla$Probability, quiet = TRUE)
grs2_roc_nonhla <- pROC::roc(response = ALL_nonhla$DISEASE, predictor = ALL_nonhla$sum_All5_nonHLA, quiet = TRUE)

prob_auc_nonhla <- pROC::auc(prob_roc_nonhla)
grs2_auc_nonhla <- pROC::auc(grs2_roc_nonhla)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_roc_nonhla)
grs2_ci <- pROC::ci.auc(grs2_auc_nonhla)
print(prob_ci)
print(grs2_ci)

delong_test_nonhla <- pROC::roc.test(prob_roc_nonhla, grs2_roc_nonhla, method = "delong")
print("DeLong Test for Non-MHC Model ROC:")
print(delong_test_nonhla)

prob_df_nonhla <- data.frame(FPR = 1 - prob_roc_nonhla$specificities, TPR = prob_roc_nonhla$sensitivities, Model = "T1GRS")
grs2_df_nonhla <- data.frame(FPR = 1 - grs2_roc_nonhla$specificities, TPR = grs2_roc_nonhla$sensitivities, Model = "GRS2")
roc_df_nonhla <- bind_rows(prob_df_nonhla, grs2_df_nonhla)

plot_roc_nonhla <- ggplot(roc_df_nonhla, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        axis.line.x = element_line(color = "black", linewidth = 0.5),
                                        axis.line.y = element_line(color = "black", linewidth = 0.5)) +
  labs(title="Discovery: Non-MHC Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_nonhla, 3)), color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_nonhla, 3)), color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_nonhla$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_nonhla$p.value, 4)))), color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
print(plot_roc_nonhla)
#ggsave('figs_may/fig3/auc_nonHLA.pdf', plot = plot_roc_nonhla, width = 6, height = 5)



t1grs_data <- ALL_nonhla %>% filter(!is.na(Probability))
grs2_data <- ALL_nonhla %>% filter(!is.na(sum_All5_nonHLA))

# --- BOOTSTRAP FOR AUPRC CONFIDENCE INTERVALS ---
n_bootstraps <- 100
t1grs_auprc_boot <- numeric(n_bootstraps)
grs2_auprc_boot <- numeric(n_bootstraps)

cat(paste("Starting bootstrapping for", n_bootstraps, "iterations...\n"))

for (i in 1:n_bootstraps) {
  # --- T1GRS Bootstrap ---
  t1grs_boot_sample <- t1grs_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(t1grs_boot_sample$DISEASE)) == 2) {
    pr_t1grs_boot <- pr.curve(scores.class0 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 2],
                              scores.class1 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 1],
                              curve = FALSE)
    t1grs_auprc_boot[i] <- pr_t1grs_boot$auc.integral
  } else {
    t1grs_auprc_boot[i] <- NA # Mark as invalid if only one class is present
  }
  
  # --- GRS2 Bootstrap ---
  grs2_boot_sample <- grs2_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(grs2_boot_sample$DISEASE)) == 2) {
    pr_grs2_boot <- pr.curve(scores.class0 = grs2_boot_sample$sum_All5_nonHLA[grs2_boot_sample$DISEASE == 2],
                             scores.class1 = grs2_boot_sample$sum_All5_nonHLA[grs2_boot_sample$DISEASE == 1],
                             curve = FALSE)
    grs2_auprc_boot[i] <- pr_grs2_boot$auc.integral
  } else {
    grs2_auprc_boot[i] <- NA
  }
}

cat("Bootstrap complete.\n")

# Calculate the 95% confidence intervals from the bootstrap results
t1grs_ci <- quantile(t1grs_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
grs2_ci <- quantile(grs2_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
# --- END BOOTSTRAP SECTION ---


# --- AUPRC CURVE PLOTTING (using original NA-filtered datasets) ---

# Calculate PR curves on the clean, NA-filtered datasets
pr_t1grs <- pr.curve(scores.class0 = t1grs_data$Probability[t1grs_data$DISEASE == 2], 
                     scores.class1 = t1grs_data$Probability[t1grs_data$DISEASE == 1], curve = TRUE)
auprc_t1grs <- pr_t1grs$auc.integral

pr_grs2 <- pr.curve(scores.class0 = grs2_data$sum_All5_nonHLA[grs2_data$DISEASE == 2], 
                    scores.class1 = grs2_data$sum_All5_nonHLA[grs2_data$DISEASE == 1], curve = TRUE)
auprc_grs2 <- pr_grs2$auc.integral

# Create formatted text strings with CIs for the plot labels
t1grs_label_pr <- paste0("T1GRS AUPRC = ", round(auprc_t1grs, 3), " (", round(t1grs_ci[1], 3), "-", round(t1grs_ci[2], 3), ")")
grs2_label_pr <- paste0("GRS2 AUPRC = ", round(auprc_grs2, 3), " (", round(grs2_ci[1], 3), "-", round(grs2_ci[2], 3), ")")

# Print results to console
print(t1grs_label_pr)
print(grs2_label_pr)



# --- MHC Only Model ---
ALL_hla <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_HLA_NoPCS_catboost.txt", sep="\t", header=TRUE)
meta_hla <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_hla) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')
ALL_hla <- merge(ALL_hla, meta_hla, by.x='CV_ID', by.y='FID', all.y = TRUE)
ALL_hla$DISEASE <- as.numeric(as.character(ALL_hla$DISEASE)) # Make sure this column is what you expect for disease status
ALL_hla <- ALL_hla[complete.cases(ALL_hla$DISEASE, ALL_hla$Probability, ALL_hla$Total_sum_HLA), ]

prob_roc_hla <- pROC::roc(response = ALL_hla$DISEASE, predictor = ALL_hla$Probability, quiet = TRUE)
grs2_roc_hla <- pROC::roc(response = ALL_hla$DISEASE, predictor = ALL_hla$Total_sum_HLA, quiet = TRUE)

prob_auc_hla <- pROC::auc(prob_roc_hla)
grs2_auc_hla <- pROC::auc(grs2_roc_hla)

# --- NEW: Extract AUC CIs and format them ---
prob_ci <- pROC::ci.auc(prob_auc_hla)
grs2_ci <- pROC::ci.auc(grs2_auc_hla)
print(prob_ci)
print(grs2_ci)

delong_test_hla <- pROC::roc.test(prob_roc_hla, grs2_roc_hla, method = "delong")
print("DeLong Test for MHC Only Model ROC:")
print(delong_test_hla)

prob_df_hla <- data.frame(FPR = 1 - prob_roc_hla$specificities, TPR = prob_roc_hla$sensitivities, Model = "T1GRS")
grs2_df_hla <- data.frame(FPR = 1 - grs2_roc_hla$specificities, TPR = grs2_roc_hla$sensitivities, Model = "GRS2")
roc_df_hla <- bind_rows(prob_df_hla, grs2_df_hla)

plot_roc_hla <- ggplot(roc_df_hla, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") +
  theme_minimal(base_size = 16) + theme(legend.title = element_blank(),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        axis.line.x = element_line(color = "black", linewidth = 0.5),
                                        axis.line.y = element_line(color = "black", linewidth = 0.5)) +
  labs(title="Discovery: MHC Only Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_hla, 3)), color = "slateblue", size = 5) +
  annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_hla, 3)), color = "firebrick", size = 5) +
  annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_hla$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_hla$p.value, 4)))), color = "black", size = 5) +
  scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
print(plot_roc_hla)
#ggsave('figs_may/fig3/auc_HLA.pdf', plot = plot_roc_hla, width = 6, height = 5)



t1grs_data <- ALL_hla %>% filter(!is.na(Probability))
grs2_data <- ALL_hla %>% filter(!is.na(Total_sum))

# --- BOOTSTRAP FOR AUPRC CONFIDENCE INTERVALS ---
n_bootstraps <- 100
t1grs_auprc_boot <- numeric(n_bootstraps)
grs2_auprc_boot <- numeric(n_bootstraps)

cat(paste("Starting bootstrapping for", n_bootstraps, "iterations...\n"))

for (i in 1:n_bootstraps) {
  # --- T1GRS Bootstrap ---
  t1grs_boot_sample <- t1grs_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(t1grs_boot_sample$DISEASE)) == 2) {
    pr_t1grs_boot <- pr.curve(scores.class0 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 2],
                              scores.class1 = t1grs_boot_sample$Probability[t1grs_boot_sample$DISEASE == 1],
                              curve = FALSE)
    t1grs_auprc_boot[i] <- pr_t1grs_boot$auc.integral
  } else {
    t1grs_auprc_boot[i] <- NA # Mark as invalid if only one class is present
  }
  
  # --- GRS2 Bootstrap ---
  grs2_boot_sample <- grs2_data %>% sample_n(size = nrow(.), replace = TRUE)
  if (length(unique(grs2_boot_sample$DISEASE)) == 2) {
    pr_grs2_boot <- pr.curve(scores.class0 = grs2_boot_sample$Total_sum_HLA[grs2_boot_sample$DISEASE == 2],
                             scores.class1 = grs2_boot_sample$Total_sum_HLA[grs2_boot_sample$DISEASE == 1],
                             curve = FALSE)
    grs2_auprc_boot[i] <- pr_grs2_boot$auc.integral
  } else {
    grs2_auprc_boot[i] <- NA
  }
}

cat("Bootstrap complete.\n")

# Calculate the 95% confidence intervals from the bootstrap results
t1grs_ci <- quantile(t1grs_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
grs2_ci <- quantile(grs2_auprc_boot, probs = c(0.025, 0.975), na.rm = TRUE)
# --- END BOOTSTRAP SECTION ---


# --- AUPRC CURVE PLOTTING (using original NA-filtered datasets) ---

# Calculate PR curves on the clean, NA-filtered datasets
pr_t1grs <- pr.curve(scores.class0 = t1grs_data$Probability[t1grs_data$DISEASE == 2], 
                     scores.class1 = t1grs_data$Probability[t1grs_data$DISEASE == 1], curve = TRUE)
auprc_t1grs <- pr_t1grs$auc.integral

pr_grs2 <- pr.curve(scores.class0 = grs2_data$Total_sum_HLA[grs2_data$DISEASE == 2], 
                    scores.class1 = grs2_data$Total_sum_HLA[grs2_data$DISEASE == 1], curve = TRUE)
auprc_grs2 <- pr_grs2$auc.integral

# Create formatted text strings with CIs for the plot labels
t1grs_label_pr <- paste0("T1GRS AUPRC = ", round(auprc_t1grs, 3), " (", round(t1grs_ci[1], 3), "-", round(t1grs_ci[2], 3), ")")
grs2_label_pr <- paste0("GRS2 AUPRC = ", round(auprc_grs2, 3), " (", round(grs2_ci[1], 3), "-", round(grs2_ci[2], 3), ")")

# Print results to console
print(t1grs_label_pr)
print(grs2_label_pr)



