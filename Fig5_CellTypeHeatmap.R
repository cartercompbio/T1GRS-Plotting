# --- Load Libraries ---
library(ComplexHeatmap)
library(readr)
library(dplyr)
library(circlize)
library(grid)
library(tibble) # Needed for column_to_rownames
library(tidyr)  # For pivot_wider
library(ggplot2) # For bar plot

# --- Parameters ---
enrichment_data_file <- '~/Code/T1D/data2025/mod_cPPA_perm_input_eg_SuSIE_95_CS_meta_T1D_137_signals.tsv.wcppa.pval.csv'
output_dir             <- '~/Code/T1D/figs_may/fig5' # Directory to save plots

# --- Define Colors ---
cluster_colors <- c("C1" = "#B38D97", "C2" = "#437A91", "C3" = "#849324", "C4" = "#F55F14", "C0" = "magenta")
# MODIFICATION: Changed color scale for -log10(p-value)
unified_color_scale <- c("grey95", "#fca311", "#d00000") 

# --- Load and Prepare Data ---
print("Loading and preparing data...")
enrichment_df_raw <- read.csv(enrichment_data_file) %>%
  rename(
    Value = obs.cppa,
    Feature = Cell.Type,
    pval = P.Value
  )

# Add dummy data for C4
all_cell_types <- unique(enrichment_df_raw$Feature)
dummy_c4_df <- tibble(
  Cluster = "C4",
  Feature = all_cell_types,
  Value = 0,
  pval = 1
)
enrichment_df <- bind_rows(enrichment_df_raw, dummy_c4_df)


# Ensure clusters are factored in a specific order
ordered_cluster_names <- c("C1", "C2", "C3", "C4")
enrichment_df$Cluster <- factor(enrichment_df$Cluster, levels = ordered_cluster_names, ordered = TRUE)


# ----------------------------------------------------
# --- Feature Selection ---
# ----------------------------------------------------
print("Selecting top 5 unique enriched cell types per cluster...")
final_selection_list <- list()
selected_features_pool <- c()
num_features <- 5

for (cluster_name in ordered_cluster_names) {
  cluster_candidates <- enrichment_df %>%
    filter(Cluster == cluster_name, Value > 0) %>%
    arrange(desc(Value))
  
  cluster_unique_features <- data.frame()
  
  if (nrow(cluster_candidates) > 0) {
    for (i in 1:nrow(cluster_candidates)) {
      feature_name <- cluster_candidates$Feature[i]
      if (!feature_name %in% selected_features_pool) {
        cluster_unique_features <- bind_rows(cluster_unique_features, cluster_candidates[i,])
        selected_features_pool <- c(selected_features_pool, feature_name)
      }
      if (nrow(cluster_unique_features) >= num_features) break
    }
  }
  final_selection_list[[cluster_name]] <- cluster_unique_features
}
top_features <- bind_rows(final_selection_list)


# ----------------------------------------------------
# --- Generate Heatmap ---
# ----------------------------------------------------
if (nrow(top_features) > 0) {
  print("Generating heatmap...")
  top_features_heatmap <- top_features %>%
    mutate(TopInCluster = Cluster)
  
  selected_feature_names <- unique(top_features_heatmap$Feature)
  all_values_for_features <- enrichment_df %>% 
    filter(Feature %in% selected_feature_names) %>%
    # MODIFICATION: Add -log10(pval) column for coloring
    mutate(log10_pval = -log10(pval))
  
  # MODIFICATION: Create a matrix of -log10(p-values) for heatmap colors
  log10_pval_mat <- all_values_for_features %>%
    select(Feature, Cluster, log10_pval) %>%
    pivot_wider(names_from = Cluster, values_from = log10_pval, values_fill = 0) %>%
    column_to_rownames("Feature")
  
  # Create a matrix of original p-values for the star annotations
  pval_mat <- all_values_for_features %>%
    select(Feature, Cluster, pval) %>%
    pivot_wider(names_from = Cluster, values_from = pval, values_fill = 1) %>%
    column_to_rownames("Feature")
  
  anno_df <- top_features_heatmap %>% select(Feature, TopInCluster) %>% column_to_rownames("Feature")
  feature_order <- top_features_heatmap %>% arrange(TopInCluster, desc(Value)) %>% pull(Feature)
  
  # Order all matrices and dataframes consistently
  log10_pval_mat_ordered <- as.matrix(log10_pval_mat[feature_order, ordered_cluster_names, drop = FALSE])
  pval_mat_ordered <- as.matrix(pval_mat[feature_order, ordered_cluster_names, drop = FALSE])
  anno_df_ordered <- anno_df[feature_order, , drop = FALSE]
  anno_df_ordered$TopInCluster <- factor(anno_df_ordered$TopInCluster, levels = ordered_cluster_names)
  
  cluster_anno <- HeatmapAnnotation(Cluster = ordered_cluster_names, col = list(Cluster = cluster_colors), show_legend=F, show_annotation_name=F)
  top_in_cluster_anno <- rowAnnotation(`Top In` = anno_df_ordered$TopInCluster, col = list(`Top In` = cluster_colors), show_annotation_name=F)
  
  # MODIFICATION: Define the color scale based on the range of -log10(p-values)
  max_log_pval <- max(log10_pval_mat_ordered[is.finite(log10_pval_mat_ordered)], na.rm = TRUE)
  col_fun <- colorRamp2(c(0, -log10(0.05), max_log_pval), unified_color_scale)
  
  # MODIFICATION: Create a custom legend for p-value significance
  significance_legend <- Legend(
    labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1"),
    title = "Significance",
    type = "points",
    pch = c("***", "**", "*", "o"),
    legend_gp = gpar(col = "black"),
    labels_gp = gpar(col = "black", fontsize = 10),
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  ht <- Heatmap(log10_pval_mat_ordered, name = "-log10(P-value)", col = col_fun,column_title_gp = gpar(fontface = "bold"),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  p_val <- pval_mat_ordered[i, j]
                  if (!is.na(p_val)) {
                    stars <- dplyr::case_when(
                      p_val < 0.001 ~ "***", 
                      p_val < 0.01  ~ "**", 
                      p_val < 0.05  ~ "*", 
                      p_val < 0.1   ~ "o",
                      TRUE ~ ""
                    )
                    if (stars != "") grid.text(stars, x, y, gp = gpar(fontsize = 14, col="black"))
                  }
                },
               # right_annotation = top_in_cluster_anno, 
               bottom_annotation = cluster_anno,
                cluster_rows = FALSE, #row_split = anno_df_ordered$TopInCluster,
                row_title_rot = 0, cluster_columns = FALSE, border = TRUE,
                column_title = "Cell Type Enrichment",
                row_names_side = "left",
                row_title_side = "right",
                row_title_gp = gpar(col = "black")
  )
  
  file_path <- file.path(output_dir, "heatmap_cell_type_enrichment_log10pval.pdf")
  pdf(file_path, width = 6, height =4)
  # MODIFICATION: Draw the heatmap and the custom significance legend
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
       annotation_legend_list = list(significance_legend))
  dev.off()
  print(paste("Cell Type Enrichment Heatmap saved to:", file_path))
} else {
  print("No positive features found for heatmap. Skipping.")
}

