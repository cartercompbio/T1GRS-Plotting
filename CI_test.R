# Replacement panels for Fig3
setwd('~/Code/T1D/') # Make sure this path is correct for your environment

# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)
library(PRROC)

# --- Initialize list to store all results ---
all_results <- list()

# --- Helper function for AUPRC with Bootstrap CI ---
calculate_auprc_ci <- function(scores_class0, scores_class1, n_boot = 200, seed = 123, curve = FALSE) {
  set.seed(seed) # for reproducibility
  
  # Observed AUPRC
  # Ensure there are scores for both classes
  if (length(scores_class0) == 0 || length(scores_class1) == 0) {
    warning("One or both classes have zero scores. Returning NA for AUPRC and CI.")
    return(list(auprc = NA, ci_low = NA, ci_high = NA, curve_data = NULL))
  }
  
  # Ensure there are enough distinct values if curve = TRUE
  # For just AUC, PRROC might be more lenient, but curve generation needs diversity.
  min_distinct_for_curve <- if (curve) 2 else 1
  
  if (length(unique(scores_class0)) < min_distinct_for_curve || length(unique(scores_class1)) < min_distinct_for_curve) {
    # If curve is not requested, and we only have one distinct value but multiple samples, PRROC might still give an AUC.
    # However, bootstrapping might still fail if all resampled values are identical.
    if (!curve && (length(unique(scores_class0)) >= 1 && length(unique(scores_class1)) >= 1)) {
      # Allow calculation of observed AUC if possible
    } else {
      warning("One or both classes have fewer than required distinct score values. Returning NA for AUPRC and CI.")
      return(list(auprc = NA, ci_low = NA, ci_high = NA, curve_data = if(curve) NA else NULL))
    }
  }
  
  
  pr_observed_obj <- PRROC::pr.curve(scores.class0 = scores_class0,
                                     scores.class1 = scores_class1,
                                     curve = curve)
  auprc_observed <- pr_observed_obj$auc.integral
  curve_data_observed <- if(curve) pr_observed_obj$curve else NULL
  
  boot_auprcs <- numeric(n_boot)
  for (i in 1:n_boot) {
    # Handle cases where original classes might be very small
    boot_scores_class0 <- if(length(scores_class0) > 0) sample(scores_class0, length(scores_class0), replace = TRUE) else numeric(0)
    boot_scores_class1 <- if(length(scores_class1) > 0) sample(scores_class1, length(scores_class1), replace = TRUE) else numeric(0)
    
    if (length(boot_scores_class0) == 0 || length(boot_scores_class1) == 0 ||
        length(unique(boot_scores_class0)) < 1 || length(unique(boot_scores_class1)) < 1 ) { # Need at least one unique value
      boot_auprcs[i] <- NA
      next
    }
    
    # Suppress warnings from pr.curve during bootstrap if desired, though errors are more critical
    pr_boot <- try(PRROC::pr.curve(scores.class0 = boot_scores_class0,
                                   scores.class1 = boot_scores_class1,
                                   curve = FALSE), silent = TRUE) # Only need AUC for bootstrap
    if (inherits(pr_boot, "try-error")) {
      boot_auprcs[i] <- NA
    } else {
      boot_auprcs[i] <- pr_boot$auc.integral
    }
  }
  
  valid_boot_auprcs <- boot_auprcs[!is.na(boot_auprcs)]
  if (length(valid_boot_auprcs) < n_boot * 0.8) {
    warning(paste("More than 20% of bootstrap iterations failed to compute AUPRC. N valid strap:", length(valid_boot_auprcs)))
  }
  
  if (length(valid_boot_auprcs) > 1) { # Need at least 2 for quantile
    ci_low <- quantile(valid_boot_auprcs, probs = 0.025, na.rm = TRUE)
    ci_high <- quantile(valid_boot_auprcs, probs = 0.975, na.rm = TRUE)
  } else {
    ci_low <- NA
    ci_high <- NA
    warning("Not enough valid bootstrap replicates to calculate CI.")
  }
  
  return(list(auprc = auprc_observed, ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high), curve_data = curve_data_observed))
}


# --- Section: Discovery - Full Model ---
cat("\n\n--- Results for Discovery: Full Model ---\n")
all_results[["Discovery_Full_Model"]] <- list()

ALL_full <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_ALL_YESPCS_catboost.txt", sep="\t", header=TRUE)
meta_full <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_full) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')
ALL_full <- merge(ALL_full, meta_full, by.x='CV_ID', by.y='FID', all.y = TRUE)
ALL_full$DISEASE <- as.numeric(as.character(ALL_full$DISEASE)) # PED: Expect 1=Control, 2=Case
ALL_full <- ALL_full[complete.cases(ALL_full$DISEASE, ALL_full$Probability, ALL_full$Total_sum), ]

# ROC Curves
prob_roc_full <- pROC::roc(response = ALL_full$DISEASE, predictor = ALL_full$Probability, quiet = TRUE, levels=c(1,2), direction="<") # Control=1, Case=2
grs2_roc_full <- pROC::roc(response = ALL_full$DISEASE, predictor = ALL_full$Total_sum, quiet = TRUE, levels=c(1,2), direction="<")   # Control=1, Case=2

prob_auc_full_val <- pROC::auc(prob_roc_full)
grs2_auc_full_val <- pROC::auc(grs2_roc_full)
prob_ci_auc_full <- pROC::ci.auc(prob_roc_full)
grs2_ci_auc_full <- pROC::ci.auc(grs2_roc_full)

all_results[["Discovery_Full_Model"]][["T1GRS_AUC"]] <- list(auc=prob_auc_full_val, ci_low=prob_ci_auc_full[1], ci_high=prob_ci_auc_full[3])
all_results[["Discovery_Full_Model"]][["GRS2_AUC"]] <- list(auc=grs2_auc_full_val, ci_low=grs2_ci_auc_full[1], ci_high=grs2_ci_auc_full[3])

cat("T1GRS AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", prob_auc_full_val, prob_ci_auc_full[1], prob_ci_auc_full[3]))
cat("GRS2 AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", grs2_auc_full_val, grs2_ci_auc_full[1], grs2_ci_auc_full[3]))

delong_test_full <- pROC::roc.test(prob_roc_full, grs2_roc_full, method = "delong")
all_results[["Discovery_Full_Model"]][["DeLong_P_value_AUC"]] <- delong_test_full$p.value
cat("DeLong Test for ROC comparison (T1GRS vs GRS2):\n")
print(delong_test_full)

# PR Curves
auprc_t1grs_full_obj <- calculate_auprc_ci(scores_class0 = ALL_full$Probability[ALL_full$DISEASE == 1], # Controls
                                           scores_class1 = ALL_full$Probability[ALL_full$DISEASE == 2], # Cases
                                           curve = TRUE, n_boot=100) # User can adjust n_boot
auprc_grs2_full_obj <- calculate_auprc_ci(scores_class0 = ALL_full$Total_sum[ALL_full$DISEASE == 1],   # Controls
                                          scores_class1 = ALL_full$Total_sum[ALL_full$DISEASE == 2],   # Cases
                                          curve = TRUE, n_boot=100)

all_results[["Discovery_Full_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_full_obj
all_results[["Discovery_Full_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_full_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_full_obj$auprc, auprc_t1grs_full_obj$ci_low, auprc_t1grs_full_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_full_obj$auprc, auprc_grs2_full_obj$ci_low, auprc_grs2_full_obj$ci_high))

# Plotting (ggsave commented out)
prob_df_full <- data.frame(FPR = 1 - prob_roc_full$specificities, TPR = prob_roc_full$sensitivities, Model = "T1GRS")
grs2_df_full <- data.frame(FPR = 1 - grs2_roc_full$specificities, TPR = grs2_roc_full$sensitivities, Model = "GRS2")
roc_df_full <- bind_rows(prob_df_full, grs2_df_full)
plot_roc_full <- ggplot(roc_df_full, aes(x = FPR, y = TPR, color = Model)) + geom_line(linewidth = 0.8) + geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: Full Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_full_val, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_full_val, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_full$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_full$p.value, 4)))), color = "black", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
print(plot_roc_full)
#ggsave('figs_may/fig3/auc_full.pdf', plot = plot_roc_full, width = 6, height = 5)

if(!is.null(auprc_t1grs_full_obj$curve_data) && !is.null(auprc_grs2_full_obj$curve_data)){
  pr_prob_df_full <- data.frame(Recall = auprc_t1grs_full_obj$curve_data[,1], Precision = auprc_t1grs_full_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_full <- data.frame(Recall = auprc_grs2_full_obj$curve_data[,1], Precision = auprc_grs2_full_obj$curve_data[,2], Model = "GRS2")
  pr_df_full <- bind_rows(pr_prob_df_full, pr_grs2_df_full)
  baseline_pr_full <- sum(ALL_full$DISEASE == 2) / nrow(ALL_full) # Case is 2
  plot_pr_full <- ggplot(pr_df_full, aes(x = Recall, y = Precision, color = Model)) + geom_line(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_full, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: Full Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0, 1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_full_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_full_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_full, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_full)
  #ggsave('figs_may/fig3/auprc_full.pdf', plot = plot_pr_full, width = 6, height = 5)
}


# --- Section: Discovery - Non-MHC Model ---
cat("\n\n--- Results for Discovery: Non-MHC Model ---\n")
all_results[["Discovery_NonMHC_Model"]] <- list()

ALL_nonhla <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_NonHLA_YESPCS_catboost.txt", sep="\t", header=TRUE)
meta_nonhla <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_nonhla) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')
ALL_nonhla <- merge(ALL_nonhla, meta_nonhla, by.x='CV_ID', by.y='FID', all.y = TRUE)
ALL_nonhla$DISEASE <- as.numeric(as.character(ALL_nonhla$DISEASE)) # PED: Expect 1=Control, 2=Case
ALL_nonhla <- ALL_nonhla[complete.cases(ALL_nonhla$DISEASE, ALL_nonhla$Probability, ALL_nonhla$sum_All5_nonHLA), ]

prob_roc_nonhla <- pROC::roc(response = ALL_nonhla$DISEASE, predictor = ALL_nonhla$Probability, quiet = TRUE, levels=c(1,2), direction="<")
grs2_roc_nonhla <- pROC::roc(response = ALL_nonhla$DISEASE, predictor = ALL_nonhla$sum_All5_nonHLA, quiet = TRUE, levels=c(1,2), direction="<")

prob_auc_nonhla_val <- pROC::auc(prob_roc_nonhla)
grs2_auc_nonhla_val <- pROC::auc(grs2_roc_nonhla)
prob_ci_auc_nonhla <- pROC::ci.auc(prob_roc_nonhla)
grs2_ci_auc_nonhla <- pROC::ci.auc(grs2_roc_nonhla)

all_results[["Discovery_NonMHC_Model"]][["T1GRS_AUC"]] <- list(auc=prob_auc_nonhla_val, ci_low=prob_ci_auc_nonhla[1], ci_high=prob_ci_auc_nonhla[3])
all_results[["Discovery_NonMHC_Model"]][["GRS2_AUC"]] <- list(auc=grs2_auc_nonhla_val, ci_low=grs2_ci_auc_nonhla[1], ci_high=grs2_ci_auc_nonhla[3])

cat("T1GRS AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", prob_auc_nonhla_val, prob_ci_auc_nonhla[1], prob_ci_auc_nonhla[3]))
cat("GRS2 AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", grs2_auc_nonhla_val, grs2_ci_auc_nonhla[1], grs2_ci_auc_nonhla[3]))

delong_test_nonhla <- pROC::roc.test(prob_roc_nonhla, grs2_roc_nonhla, method = "delong")
all_results[["Discovery_NonMHC_Model"]][["DeLong_P_value_AUC"]] <- delong_test_nonhla$p.value
cat("DeLong Test for ROC comparison (T1GRS vs GRS2):\n")
print(delong_test_nonhla)

auprc_t1grs_nonhla_obj <- calculate_auprc_ci(scores_class0 = ALL_nonhla$Probability[ALL_nonhla$DISEASE == 1],
                                             scores_class1 = ALL_nonhla$Probability[ALL_nonhla$DISEASE == 2],
                                             curve = TRUE, n_boot=100)
auprc_grs2_nonhla_obj <- calculate_auprc_ci(scores_class0 = ALL_nonhla$sum_All5_nonHLA[ALL_nonhla$DISEASE == 1],
                                            scores_class1 = ALL_nonhla$sum_All5_nonHLA[ALL_nonhla$DISEASE == 2],
                                            curve = TRUE, n_boot=100)

all_results[["Discovery_NonMHC_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_nonhla_obj
all_results[["Discovery_NonMHC_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_nonhla_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_nonhla_obj$auprc, auprc_t1grs_nonhla_obj$ci_low, auprc_t1grs_nonhla_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_nonhla_obj$auprc, auprc_grs2_nonhla_obj$ci_low, auprc_grs2_nonhla_obj$ci_high))

# Plotting (ggsave commented out)
# ... (plotting code for Non-MHC similar to Full Model, using _nonhla variables) ...
if(!is.null(auprc_t1grs_nonhla_obj$curve_data) && !is.null(auprc_grs2_nonhla_obj$curve_data)){
  # ROC plot (similar to full model, adapt variable names)
  prob_df_nonhla <- data.frame(FPR = 1 - prob_roc_nonhla$specificities, TPR = prob_roc_nonhla$sensitivities, Model = "T1GRS")
  grs2_df_nonhla <- data.frame(FPR = 1 - grs2_roc_nonhla$specificities, TPR = grs2_roc_nonhla$sensitivities, Model = "GRS2")
  roc_df_nonhla <- bind_rows(prob_df_nonhla, grs2_df_nonhla)
  plot_roc_nonhla <- ggplot(roc_df_nonhla, aes(x = FPR, y = TPR, color = Model)) + geom_line(linewidth = 0.8) + geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: Non-MHC Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_nonhla_val, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_nonhla_val, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_nonhla$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_nonhla$p.value, 4)))), color = "black", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_roc_nonhla)
  # #ggsave('figs_may/fig3/auc_nonHLA.pdf', plot = plot_roc_nonhla, width = 6, height = 5)
  
  # PR plot (similar to full model, adapt variable names)
  pr_prob_df_nonhla <- data.frame(Recall = auprc_t1grs_nonhla_obj$curve_data[,1], Precision = auprc_t1grs_nonhla_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_nonhla <- data.frame(Recall = auprc_grs2_nonhla_obj$curve_data[,1], Precision = auprc_grs2_nonhla_obj$curve_data[,2], Model = "GRS2")
  pr_df_nonhla <- bind_rows(pr_prob_df_nonhla, pr_grs2_df_nonhla)
  baseline_pr_nonhla <- sum(ALL_nonhla$DISEASE == 2) / nrow(ALL_nonhla)
  plot_pr_nonhla <- ggplot(pr_df_nonhla, aes(x = Recall, y = Precision, color = Model)) + geom_line(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_nonhla, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: Non-MHC Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0,1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_nonhla_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_nonhla_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_nonhla, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_nonhla)
  # #ggsave('figs_may/fig3/auprc_nonHLA.pdf', plot = plot_pr_nonhla, width = 6, height = 5)
}


# --- Section: Discovery - MHC Only Model ---
cat("\n\n--- Results for Discovery: MHC Only Model ---\n")
all_results[["Discovery_MHCOnly_Model"]] <- list()

ALL_hla <- read.table("/Users/tjsears/Code/T1D/data2025/catboost_update/delong_HLA_YESPCS_catboost.txt", sep="\t", header=TRUE)
meta_hla <- read.table("/Users/tjsears/Code/T1D/data2025/meta_all5_cat_cov.ped", sep=" ", header=F)
colnames(meta_hla) <- c('FID', 'IID', 'SEX', 'DISEASE', 'PC1', 'PC2', 'PC3', 'PC4', 'T1DGC', 'DCCT', 'GENIE_UK', 'GoKIND')
ALL_hla <- merge(ALL_hla, meta_hla, by.x='CV_ID', by.y='FID', all.y = TRUE)
ALL_hla$DISEASE <- as.numeric(as.character(ALL_hla$DISEASE)) # PED: Expect 1=Control, 2=Case
ALL_hla <- ALL_hla[complete.cases(ALL_hla$DISEASE, ALL_hla$Probability, ALL_hla$Total_sum_HLA), ]

prob_roc_hla <- pROC::roc(response = ALL_hla$DISEASE, predictor = ALL_hla$Probability, quiet = TRUE, levels=c(1,2), direction="<")
grs2_roc_hla <- pROC::roc(response = ALL_hla$DISEASE, predictor = ALL_hla$Total_sum_HLA, quiet = TRUE, levels=c(1,2), direction="<")

prob_auc_hla_val <- pROC::auc(prob_roc_hla)
grs2_auc_hla_val <- pROC::auc(grs2_roc_hla)
prob_ci_auc_hla <- pROC::ci.auc(prob_roc_hla)
grs2_ci_auc_hla <- pROC::ci.auc(grs2_roc_hla)

all_results[["Discovery_MHCOnly_Model"]][["T1GRS_AUC"]] <- list(auc=prob_auc_hla_val, ci_low=prob_ci_auc_hla[1], ci_high=prob_ci_auc_hla[3])
all_results[["Discovery_MHCOnly_Model"]][["GRS2_AUC"]] <- list(auc=grs2_auc_hla_val, ci_low=grs2_ci_auc_hla[1], ci_high=grs2_ci_auc_hla[3])

cat("T1GRS AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", prob_auc_hla_val, prob_ci_auc_hla[1], prob_ci_auc_hla[3]))
cat("GRS2 AUC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", grs2_auc_hla_val, grs2_ci_auc_hla[1], grs2_ci_auc_hla[3]))

delong_test_hla <- pROC::roc.test(prob_roc_hla, grs2_roc_hla, method = "delong")
all_results[["Discovery_MHCOnly_Model"]][["DeLong_P_value_AUC"]] <- delong_test_hla$p.value
cat("DeLong Test for ROC comparison (T1GRS vs GRS2):\n")
print(delong_test_hla)

auprc_t1grs_hla_obj <- calculate_auprc_ci(scores_class0 = ALL_hla$Probability[ALL_hla$DISEASE == 1],
                                          scores_class1 = ALL_hla$Probability[ALL_hla$DISEASE == 2],
                                          curve = TRUE, n_boot=100)
auprc_grs2_hla_obj <- calculate_auprc_ci(scores_class0 = ALL_hla$Total_sum_HLA[ALL_hla$DISEASE == 1],
                                         scores_class1 = ALL_hla$Total_sum_HLA[ALL_hla$DISEASE == 2],
                                         curve = TRUE, n_boot=100)

all_results[["Discovery_MHCOnly_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_hla_obj
all_results[["Discovery_MHCOnly_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_hla_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_hla_obj$auprc, auprc_t1grs_hla_obj$ci_low, auprc_t1grs_hla_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_hla_obj$auprc, auprc_grs2_hla_obj$ci_low, auprc_grs2_hla_obj$ci_high))

# Plotting (ggsave commented out)
# ... (plotting code for MHC Only similar to Full Model, using _hla variables) ...
if(!is.null(auprc_t1grs_hla_obj$curve_data) && !is.null(auprc_grs2_hla_obj$curve_data)){
  # ROC plot (similar to full model, adapt variable names)
  prob_df_hla <- data.frame(FPR = 1 - prob_roc_hla$specificities, TPR = prob_roc_hla$sensitivities, Model = "T1GRS")
  grs2_df_hla <- data.frame(FPR = 1 - grs2_roc_hla$specificities, TPR = grs2_roc_hla$sensitivities, Model = "GRS2")
  roc_df_hla <- bind_rows(prob_df_hla, grs2_df_hla)
  plot_roc_hla <- ggplot(roc_df_hla, aes(x = FPR, y = TPR, color = Model)) + geom_line(linewidth = 0.8) + geom_abline(intercept=0, slope=1, linetype="dashed", color = "grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: MHC Only Model", x = "1 - Specificity", y = "Sensitivity", color = "Model") + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUC = ", round(prob_auc_hla_val, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUC = ", round(grs2_auc_hla_val, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("DeLong P ", ifelse(delong_test_hla$p.value < 0.0001, "< 0.0001", paste0("= ",round(delong_test_hla$p.value, 4)))), color = "black", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_roc_hla)
  # #ggsave('figs_may/fig3/auc_HLA.pdf', plot = plot_roc_hla, width = 6, height = 5)
  
  # PR plot (similar to full model, adapt variable names)
  pr_prob_df_hla <- data.frame(Recall = auprc_t1grs_hla_obj$curve_data[,1], Precision = auprc_t1grs_hla_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_hla <- data.frame(Recall = auprc_grs2_hla_obj$curve_data[,1], Precision = auprc_grs2_hla_obj$curve_data[,2], Model = "GRS2")
  pr_df_hla <- bind_rows(pr_prob_df_hla, pr_grs2_df_hla)
  baseline_pr_hla <- sum(ALL_hla$DISEASE == 2) / nrow(ALL_hla)
  plot_pr_hla <- ggplot(pr_df_hla, aes(x = Recall, y = Precision, color = Model)) + geom_line(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_hla, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="Discovery: MHC Only Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0,1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_hla_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_hla_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_hla, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_hla)
  # #ggsave('figs_may/fig3/auprc_HLA.pdf', plot = plot_pr_hla, width = 6, height = 5)
}


# --- Section: NPOD AUPRC ---
cat("\n\n############################\n# NPOD AUPRC Results #\n############################\n")
all_results[["NPOD_AUPRC"]] <- list()

# Data for NPOD
ALL_npod_orig <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_T1GRS_Scores.tsv", sep="\t", header=TRUE)
GRS2_npod_orig <- read.table("/Users/tjsears/Code/T1D/data2025/Sets1_2_GRS2_Scores.tsv", sep="\t", header=TRUE)
ALL_npod_merged <- merge(ALL_npod_orig, GRS2_npod_orig, by='ID') # by='ID' or by.x/by.y if ID columns are named differently after read

# --- NPOD: Full Model AUPRC ---
cat("\n--- NPOD AUPRC: Full Model ---\n")
all_results[["NPOD_AUPRC"]][["Full_Model"]] <- list()

ALL_npod_full <- ALL_npod_merged
ALL_npod_full$DISEASE <- ifelse(ALL_npod_full$Condition == 'T1D', 2, 1) # Case (T1D) = 2, Control = 1
ALL_npod_full$DISEASE <- as.numeric(as.character(ALL_npod_full$DISEASE))
ALL_npod_full <- ALL_npod_full[complete.cases(ALL_npod_full$DISEASE, ALL_npod_full$T1GRS_ALL, ALL_npod_full$Total_sum), ]

auprc_t1grs_npod_full_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_full$T1GRS_ALL[ALL_npod_full$DISEASE == 1], # Controls
                                                scores_class1 = ALL_npod_full$T1GRS_ALL[ALL_npod_full$DISEASE == 2], # Cases
                                                curve = TRUE, n_boot=100)
auprc_grs2_npod_full_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_full$Total_sum[ALL_npod_full$DISEASE == 1],   # Controls
                                               scores_class1 = ALL_npod_full$Total_sum[ALL_npod_full$DISEASE == 2],   # Cases
                                               curve = TRUE, n_boot=100)

all_results[["NPOD_AUPRC"]][["Full_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_npod_full_obj
all_results[["NPOD_AUPRC"]][["Full_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_npod_full_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_npod_full_obj$auprc, auprc_t1grs_npod_full_obj$ci_low, auprc_t1grs_npod_full_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_npod_full_obj$auprc, auprc_grs2_npod_full_obj$ci_low, auprc_grs2_npod_full_obj$ci_high))

# Plotting (ggsave commented out)
if(!is.null(auprc_t1grs_npod_full_obj$curve_data) && !is.null(auprc_grs2_npod_full_obj$curve_data)){
  pr_prob_df_npod_full <- data.frame(Recall = auprc_t1grs_npod_full_obj$curve_data[,1], Precision = auprc_t1grs_npod_full_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_npod_full <- data.frame(Recall = auprc_grs2_npod_full_obj$curve_data[,1], Precision = auprc_grs2_npod_full_obj$curve_data[,2], Model = "GRS2")
  pr_df_npod_full <- bind_rows(pr_prob_df_npod_full, pr_grs2_df_npod_full)
  baseline_pr_npod_full <- sum(ALL_npod_full$DISEASE == 2) / nrow(ALL_npod_full)
  plot_pr_npod_full <- ggplot(pr_df_npod_full, aes(x = Recall, y = Precision, color = Model)) + geom_path(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_npod_full, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="NPOD: Full Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0, 1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_npod_full_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_npod_full_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_npod_full, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_npod_full)
  # #ggsave('figs_may/fig3/NPOD_auprc_full.pdf', plot = plot_pr_npod_full, width = 6, height = 5)
}


# --- NPOD: MHC Only Model AUPRC ---
cat("\n--- NPOD AUPRC: MHC Only Model ---\n")
all_results[["NPOD_AUPRC"]][["MHCOnly_Model"]] <- list()

ALL_npod_mhc <- ALL_npod_merged
ALL_npod_mhc$DISEASE <- ifelse(ALL_npod_mhc$Condition == 'T1D', 2, 1) # Case (T1D) = 2, Control = 1
ALL_npod_mhc$DISEASE <- as.numeric(as.character(ALL_npod_mhc$DISEASE))
ALL_npod_mhc <- ALL_npod_mhc[complete.cases(ALL_npod_mhc$DISEASE, ALL_npod_mhc$T1GRS_MHC, ALL_npod_mhc$Total_sum_HLA), ]

auprc_t1grs_npod_mhc_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_mhc$T1GRS_MHC[ALL_npod_mhc$DISEASE == 1],
                                               scores_class1 = ALL_npod_mhc$T1GRS_MHC[ALL_npod_mhc$DISEASE == 2],
                                               curve = TRUE, n_boot=100)
auprc_grs2_npod_mhc_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_mhc$Total_sum_HLA[ALL_npod_mhc$DISEASE == 1],
                                              scores_class1 = ALL_npod_mhc$Total_sum_HLA[ALL_npod_mhc$DISEASE == 2],
                                              curve = TRUE, n_boot=100)

all_results[["NPOD_AUPRC"]][["MHCOnly_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_npod_mhc_obj
all_results[["NPOD_AUPRC"]][["MHCOnly_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_npod_mhc_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_npod_mhc_obj$auprc, auprc_t1grs_npod_mhc_obj$ci_low, auprc_t1grs_npod_mhc_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_npod_mhc_obj$auprc, auprc_grs2_npod_mhc_obj$ci_low, auprc_grs2_npod_mhc_obj$ci_high))

# Plotting (ggsave commented out)
if(!is.null(auprc_t1grs_npod_mhc_obj$curve_data) && !is.null(auprc_grs2_npod_mhc_obj$curve_data)){
  pr_prob_df_npod_mhc <- data.frame(Recall = auprc_t1grs_npod_mhc_obj$curve_data[,1], Precision = auprc_t1grs_npod_mhc_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_npod_mhc <- data.frame(Recall = auprc_grs2_npod_mhc_obj$curve_data[,1], Precision = auprc_grs2_npod_mhc_obj$curve_data[,2], Model = "GRS2")
  pr_df_npod_mhc <- bind_rows(pr_prob_df_npod_mhc, pr_grs2_df_npod_mhc)
  baseline_pr_npod_mhc <- sum(ALL_npod_mhc$DISEASE == 2) / nrow(ALL_npod_mhc)
  plot_pr_npod_mhc <- ggplot(pr_df_npod_mhc, aes(x = Recall, y = Precision, color = Model)) + geom_path(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_npod_mhc, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="NPOD: MHC Only Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0, 1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_npod_mhc_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_npod_mhc_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_npod_mhc, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_npod_mhc)
  # #ggsave('figs_may/fig3/NPOD_auprc_MHCOnly.pdf', plot = plot_pr_npod_mhc, width = 6, height = 5)
}

# --- NPOD: Non-MHC Model AUPRC ---
cat("\n--- NPOD AUPRC: Non-MHC Model ---\n")
all_results[["NPOD_AUPRC"]][["NonMHC_Model"]] <- list()

ALL_npod_nonmhc <- ALL_npod_merged
ALL_npod_nonmhc$DISEASE <- ifelse(ALL_npod_nonmhc$Condition == 'T1D', 2, 1) # Case (T1D) = 2, Control = 1
ALL_npod_nonmhc$DISEASE <- as.numeric(as.character(ALL_npod_nonmhc$DISEASE))
ALL_npod_nonmhc <- ALL_npod_nonmhc[complete.cases(ALL_npod_nonmhc$DISEASE, ALL_npod_nonmhc$T1GRS_nonMHC, ALL_npod_nonmhc$sum_All5_nonMHC), ]

auprc_t1grs_npod_nonmhc_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_nonmhc$T1GRS_nonMHC[ALL_npod_nonmhc$DISEASE == 1], # Controls
                                                  scores_class1 = ALL_npod_nonmhc$T1GRS_nonMHC[ALL_npod_nonmhc$DISEASE == 2], # Cases
                                                  curve = TRUE, n_boot=100)
auprc_grs2_npod_nonmhc_obj <- calculate_auprc_ci(scores_class0 = ALL_npod_nonmhc$sum_All5_nonMHC[ALL_npod_nonmhc$DISEASE == 1],   # Controls
                                                 scores_class1 = ALL_npod_nonmhc$sum_All5_nonMHC[ALL_npod_nonmhc$DISEASE == 2],   # Cases
                                                 curve = TRUE, n_boot=100)

all_results[["NPOD_AUPRC"]][["NonMHC_Model"]][["T1GRS_AUPRC"]] <- auprc_t1grs_npod_nonmhc_obj
all_results[["NPOD_AUPRC"]][["NonMHC_Model"]][["GRS2_AUPRC"]] <- auprc_grs2_npod_nonmhc_obj

cat("T1GRS AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_t1grs_npod_nonmhc_obj$auprc, auprc_t1grs_npod_nonmhc_obj$ci_low, auprc_t1grs_npod_nonmhc_obj$ci_high))
cat("GRS2 AUPRC:", sprintf("%.3f (95%% CI: %.3f-%.3f)\n", auprc_grs2_npod_nonmhc_obj$auprc, auprc_grs2_npod_nonmhc_obj$ci_low, auprc_grs2_npod_nonmhc_obj$ci_high))

# Plotting (ggsave commented out)
if(!is.null(auprc_t1grs_npod_nonmhc_obj$curve_data) && !is.null(auprc_grs2_npod_nonmhc_obj$curve_data)){
  pr_prob_df_npod_nonmhc <- data.frame(Recall = auprc_t1grs_npod_nonmhc_obj$curve_data[,1], Precision = auprc_t1grs_npod_nonmhc_obj$curve_data[,2], Model = "T1GRS")
  pr_grs2_df_npod_nonmhc <- data.frame(Recall = auprc_grs2_npod_nonmhc_obj$curve_data[,1], Precision = auprc_grs2_npod_nonmhc_obj$curve_data[,2], Model = "GRS2")
  pr_df_npod_nonmhc <- bind_rows(pr_prob_df_npod_nonmhc, pr_grs2_df_npod_nonmhc)
  baseline_pr_npod_nonmhc <- sum(ALL_npod_nonmhc$DISEASE == 2) / nrow(ALL_npod_nonmhc)
  plot_pr_npod_nonmhc <- ggplot(pr_df_npod_nonmhc, aes(x = Recall, y = Precision, color = Model)) + geom_path(linewidth = 0.8) + geom_hline(yintercept = baseline_pr_npod_nonmhc, linetype="dashed", color="grey50") + theme_minimal(base_size = 16) + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.5), axis.line.y = element_line(color = "black", linewidth = 0.5)) + labs(title="NPOD: Non-MHC Model\nAUPRC Curve", x = "Recall", y = "Precision", color = "Model") + ylim(0, 1) + annotate("text", x = 0.65, y = 0.15, label = paste0("T1GRS AUPRC = ", round(auprc_t1grs_npod_nonmhc_obj$auprc, 3)), color = "slateblue", size = 5) + annotate("text", x = 0.65, y = 0.1, label = paste0("GRS2 AUPRC = ", round(auprc_grs2_npod_nonmhc_obj$auprc, 3)), color = "firebrick", size = 5) + annotate("text", x = 0.65, y = 0.05, label = paste0("Baseline = ", round(baseline_pr_npod_nonmhc, 3)), color = "grey50", size = 5) + scale_color_manual(values = c("T1GRS" = "slateblue", "GRS2" = "firebrick"))
  print(plot_pr_npod_nonmhc)
  # #ggsave('figs_may/fig3/NPOD_auprc_nonMHC.pdf', plot = plot_pr_npod_nonmhc, width = 6, height = 5)
}

# --- Print All Stored Results ---
cat("\n\n\n--- Summary of All Results ---\n")
for (dataset_name in names(all_results)) {
  cat(paste0("\n--- Dataset/Analysis: ", gsub("_", " ", dataset_name), " ---\n"))
  dataset_results <- all_results[[dataset_name]]
  for (model_type_name in names(dataset_results)) {
    cat(paste0("\n  Model Type: ", gsub("_", " ", model_type_name), "\n"))
    model_results <- dataset_results[[model_type_name]]
    
    if (!is.null(model_results$T1GRS_AUC)) {
      cat(sprintf("    T1GRS AUC: %.3f (95%% CI: %.3f-%.3f)\n",
                  model_results$T1GRS_AUC$auc, model_results$T1GRS_AUC$ci_low, model_results$T1GRS_AUC$ci_high))
    }
    if (!is.null(model_results$GRS2_AUC)) {
      cat(sprintf("    GRS2 AUC: %.3f (95%% CI: %.3f-%.3f)\n",
                  model_results$GRS2_AUC$auc, model_results$GRS2_AUC$ci_low, model_results$GRS2_AUC$ci_high))
    }
    if (!is.null(model_results$DeLong_P_value_AUC)) {
      cat(sprintf("    DeLong Test P-value (AUCs): %.4f\n", model_results$DeLong_P_value_AUC))
    }
    if (!is.null(model_results$T1GRS_AUPRC)) {
      cat(sprintf("    T1GRS AUPRC: %.3f (95%% CI: %.3f-%.3f)\n",
                  model_results$T1GRS_AUPRC$auprc, model_results$T1GRS_AUPRC$ci_low, model_results$T1GRS_AUPRC$ci_high))
    }
    if (!is.null(model_results$GRS2_AUPRC)) {
      cat(sprintf("    GRS2 AUPRC: %.3f (95%% CI: %.3f-%.3f)\n",
                  model_results$GRS2_AUPRC$auprc, model_results$GRS2_AUPRC$ci_low, model_results$GRS2_AUPRC$ci_high))
    }
  }
}
cat("\n--- End of Analysis Summary ---\n")