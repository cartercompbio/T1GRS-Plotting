# --- Load Libraries ---
library(ComplexHeatmap)
library(readr)
library(dplyr)
library(circlize)
library(grid)
library(tibble) # Needed for column_to_rownames
library(tidyr)  # For pivot_wider

# --- Parameters ---
shap_data_file <- '~/Code/T1D/data2025/FCshapval.txt'
output_dir     <- '~/Code/T1D/figs_may/fig5' # Directory to save plots

# --- Define Colors ---
cluster_colors <- c("C0" = "magenta", "C1" = "#B38D97", "C2" = "#437F97", "C3" = "#849324", "C4" = "#F58814")
unified_color_scale <- c("#99B2DD", "white", "#F06449")

# --- Load and Prepare Data ---
print("Loading and preparing data...")
shap_df <- read.table(shap_data_file, sep='\t', header=T)
shap_df$MHC <- grepl("^rs|HLA|SNPS|intron|exon", shap_df$Feature, ignore.case = TRUE)
ordered_cluster_names <- names(cluster_colors)[names(cluster_colors) %in% unique(shap_df$Cluster)]
shap_df$Cluster <- factor(shap_df$Cluster, levels = ordered_cluster_names, ordered = TRUE)


# --- Helper Function for Feature Selection ---
select_unique_positive_features <- function(data, clusters, num_features = 5) {
  final_selection_list <- list()
  selected_features_pool <- c()
  
  for (cluster_name in clusters) {
    cluster_candidates <- data %>%
      filter(Cluster == cluster_name, MeanSHAPDifference > 0) %>%
      arrange(desc(MeanSHAPDifference))
    
    cluster_unique_features <- data.frame()
    
    for (i in 1:nrow(cluster_candidates)) {
      feature_name <- cluster_candidates$Feature[i]
      if (!feature_name %in% selected_features_pool) {
        cluster_unique_features <- bind_rows(cluster_unique_features, cluster_candidates[i,])
        selected_features_pool <- c(selected_features_pool, feature_name)
      }
      if (nrow(cluster_unique_features) >= num_features) break
    }
    final_selection_list[[cluster_name]] <- cluster_unique_features
  }
  bind_rows(final_selection_list)
}

# --- P-value Star Function ---
cell_fun_stars = function(j, i, x, y, w, h, fill, pval_matrix_ref) {
  p_val <- pval_matrix_ref[i, j]
  if (!is.na(p_val)) {
    stars <- dplyr::case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "")
    if (stars != "") grid.text(stars, x, y, gp = gpar(fontsize = 12))
  }
}

# ----------------------------------------------------
# --- Generate Heatmap for MHC Features ---
# ----------------------------------------------------
print("Generating heatmap for: Top 5 Unique Positive MHC Features per Cluster")

mhc_data_subset <- shap_df %>% filter(MHC == TRUE)
mhc_data_subset <- mhc_data_subset[mhc_data_subset$Gene %in% c('rs1064173','HLA_DRB1*04:01','HLA_DQA1*03:01',
                                                                   'HLA_B*39:06','HLA_DQB1*03:01', 'HLA_B*39', 'HLA_DPB1*15:01', 'HLA_DQB1*03:02',
                                                                   'HLA_DQA1*05:01','HLA_DRB1*03:01'),]
mhc_top_features <- mhc_data_subset

if (nrow(mhc_top_features) > 0) {
  mhc_top_features <- mhc_top_features %>%
    rename(log2fc = MeanSHAPDifference) %>%
    mutate(TopInCluster = Cluster)
  
  mhc_selected_feature_names <- unique(mhc_top_features$Feature)
  mhc_all_values <- shap_df %>% filter(Feature %in% mhc_selected_feature_names)
  
  mhc_log2fc_mat <- mhc_all_values %>%
    select(Feature, Cluster, MeanSHAPDifference) %>%
    pivot_wider(names_from = Cluster, values_from = MeanSHAPDifference, values_fill = 0) %>%
    column_to_rownames("Feature")
  
  mhc_pval_mat <- mhc_all_values %>%
    select(Feature, Cluster, pval) %>%
    pivot_wider(names_from = Cluster, values_from = pval, values_fill = 1) %>%
    column_to_rownames("Feature")
  
  mhc_anno_df <- mhc_top_features %>% select(Feature, Gene, TopInCluster) #%>% column_to_rownames("Feature")
  mhc_cluster_anno <- HeatmapAnnotation(Cluster = ordered_cluster_names, col = list(Cluster = cluster_colors), show_legend=F, show_annotation_name=F)
  
  mhc_col_fun <- colorRamp2(c(-0.25, 0,0.25), unified_color_scale)
  library(dendextend)
  row_dend <- as.dendrogram(hclust(dist(log2(2^(as.matrix(mhc_log2fc_mat))))))
  reversed_row_dend <- (row_dend)
  
  mhc_ht <- Heatmap(log2(2^(as.matrix(mhc_log2fc_mat))), name = "Log2 Fold Change\nSHAP Val.", col = mhc_col_fun,
                    cell_fun = function(j,i,x,y,w,h,fill) cell_fun_stars(j,i,x,y,w,h,fill,as.matrix(mhc_pval_mat)),
                    #right_annotation = mhc_top_in_cluster_anno,
                    bottom_annotation = mhc_cluster_anno,
                    
                    # Use the pre-computed and reversed dendrogram
                    cluster_rows = reversed_row_dend,
                    
                    #row_split = mhc_anno_df_ordered$TopInCluster,
                    row_title_rot = 0, cluster_columns = FALSE, border = T,
                    column_title = "MHC Features",
                    row_labels = unique(mhc_top_features$Gene),
                    row_names_side = "left",
                    
                    # Set to TRUE to see the reversed dendrogram
                    show_row_dend = F,  column_title_gp = gpar(fontface = "bold"),

                    
                    row_title_side = "right",
                    row_title_gp = gpar(col = "black")
  )
  
  # Save plot as PDF with adjusted width
  mhc_plot_width <- 7.5
  mhc_file_path <- file.path(output_dir, "heatmap_positive_mhc_squat.pdf")
  pdf(mhc_file_path, width = mhc_plot_width, height = max(0.18 * nrow(mhc_log2fc_mat_ordered)))
  draw(mhc_ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
  dev.off()
  print(paste("MHC Heatmap saved to:", mhc_file_path))
} else {
  print("No positive MHC features found. Skipping MHC heatmap.")
}


# ----------------------------------------------------
# --- Generate Heatmap for Non-MHC Features ---
# ----------------------------------------------------
print("Generating heatmap for: Top 5 Unique Positive Non-MHC Features per Cluster")

non_mhc_data_subset <- shap_df %>% filter(MHC == FALSE)
non_mhc_data_subset <- non_mhc_data_subset[non_mhc_data_subset$Gene %in% c('PTPN22','INS','IKZF4','DEXI','BACH2', 'CEL', 'GLIS3', 'RUNX3'),]
non_mhc_top_features <- non_mhc_data_subset

if (nrow(non_mhc_top_features) > 0) {
  non_mhc_top_features <- non_mhc_top_features %>%
    rename(log2fc = MeanSHAPDifference) %>%
    mutate(TopInCluster = Cluster)
  
  non_mhc_selected_feature_names <- unique(non_mhc_top_features$Feature)
  non_mhc_all_values <- shap_df %>% filter(Feature %in% non_mhc_selected_feature_names)
  
  non_mhc_log2fc_mat <- non_mhc_all_values %>%
    select(Feature, Cluster, MeanSHAPDifference) %>%
    pivot_wider(names_from = Cluster, values_from = MeanSHAPDifference, values_fill = 0) %>%
    column_to_rownames("Feature")
  
  non_mhc_pval_mat <- non_mhc_all_values %>%
    select(Feature, Cluster, pval) %>%
    pivot_wider(names_from = Cluster, values_from = pval, values_fill = 1) %>%
    column_to_rownames("Feature")
  
  non_mhc_anno_df <- non_mhc_top_features %>% select(Feature, Gene, TopInCluster) #%>% column_to_rownames("Feature")
  non_mhc_cluster_anno <- HeatmapAnnotation(Cluster = ordered_cluster_names, col = list(Cluster = cluster_colors), show_legend=F, show_annotation_name=F)

  non_mhc_col_fun <- colorRamp2(c(-0.1, 0,0.1), unified_color_scale)
  
  non_mhc_ht <- Heatmap(log2(2^(as.matrix(non_mhc_log2fc_mat))), name = "Log2 Fold Change\nSHAP Val.", col = non_mhc_col_fun,
                        cell_fun = function(j,i,x,y,w,h,fill) cell_fun_stars(j,i,x,y,w,h,fill,as.matrix(non_mhc_pval_mat)),
                        #right_annotation = non_mhc_top_in_cluster_anno,
                        bottom_annotation = non_mhc_cluster_anno,
                        cluster_rows = T, 
                        #row_split = non_mhc_anno_df_ordered$TopInCluster,
                        row_title_rot = 0, cluster_columns = FALSE, border = T,
                        column_title_gp = gpar(fontface = "bold"),
                        column_title = "Non-MHC Features",
                        row_labels = unique(non_mhc_top_features$Gene),
                        row_names_side = "left",show_row_dend = F,
                        row_title_side = "right",
                        row_title_gp = gpar(col = "black")
  )
  
  # Save plot as PDF with adjusted width
  non_mhc_plot_width <- 7
  non_mhc_file_path <- file.path(output_dir, "heatmap_positive_non_mhc_squat.pdf")
  pdf(non_mhc_file_path, width = non_mhc_plot_width, height = max(0.16 * nrow(non_mhc_log2fc_mat_ordered)))
  draw(non_mhc_ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
  dev.off()
  print(paste("Non-MHC Heatmap saved to:", non_mhc_file_path))
} else {
  print("No positive Non-MHC features found. Skipping Non-MHC heatmap.")
}

print("R script finished.")



