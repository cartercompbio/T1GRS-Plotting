# --- Load Libraries ---
# Ensure these packages are installed: install.packages(c("ComplexHeatmap", "readr", "dplyr", "circlize", "tibble"))
library(ComplexHeatmap)
library(readr)
library(dplyr)
library(circlize)
library(grid)
library(tibble) # Needed for column_to_rownames

# --- Parameters ---
# Define file paths (Ensure these match where Python saved the files)
log2fc_file <- '~/Code/T1D/data2025/heatmap_log2fc_matrix.csv'
pval_file   <- '~/Code/T1D/data2025/heatmap_pval_matrix.csv'
anno_file   <- '~/Code/T1D/data2025/heatmap_annotations.csv'
output_dir  <- '~/Code/T1D/figs_may/fig5' # Directory to save plots

# Define colors
# Cluster colors
cluster_colors <- c("C0" = "magenta", "C1" = "#B38D97", "C2" = "#437F97", "C3" = "#849324", "C4" = "#F58814")

mhc_col='#7061CF'
cd4_col='#FFC20A'
panc_col='#36AB9F'
inflam_col='#D00000'

# Combine MHC and Category colors for the unified legend
mhc_color <- c("MHC" = "#7061CF") # Specific color for MHC type
category_colors <- c("Regulatory" = "#FFC20A", "Initiation" = "#36AB9F", "Amplification" = "#D00000", "Uncategorized" = "grey80")
category_colors <- c("Tcell" = "#FFC20A", "Pancreatic" = "#36AB9F", "Innate Immune" = "#D00000", "Uncategorized" = "grey80")

all_type_colors <- c(mhc_color, category_colors) # Used for the combined legend

# Define order for sorting and legends
all_type_order <- c('MHC', "Tcell", "Pancreatic", "Innate Immune", 'Uncategorized') # For combined legend order

category_order <- c("Tcell", "Pancreatic", "Innate Immune", 'Uncategorized') # For sorting non-MHC data

# Heatmap color scales
mhc_color_scale <- c("#99B2DD", "white", "#F06449")
non_mhc_color_scale <- c("#99B2DD", "white", "#F06449")

# --- Load Data ---
log2fc_mat <- as.matrix(read_csv(log2fc_file, col_types = cols(.default = "d", Feature = "c")) %>% tibble::column_to_rownames("Feature"))
pval_mat   <- as.matrix(read_csv(pval_file, col_types = cols(.default = "d", Feature = "c")) %>% tibble::column_to_rownames("Feature"))
anno_df    <- read_csv(anno_file, col_types = cols(Feature = "c", Gene = "c", MHC = "l", Category = "c")) %>% tibble::column_to_rownames("Feature")
anno_df$MHC[startsWith(rownames(anno_df),'rs')] <- TRUE

print("Data loaded successfully.")

# --- Data Validation and Preparation ---
if (!identical(rownames(log2fc_mat), rownames(pval_mat)) || !identical(rownames(log2fc_mat), rownames(anno_df))) {
  stop("Row names (Features) do not match between loaded files.")
}
if(!all(colnames(log2fc_mat) %in% names(cluster_colors))) {
  warning("Cluster names mismatch. Filtering colors.")
  cluster_colors <- cluster_colors[names(cluster_colors) %in% colnames(log2fc_mat)]
  if(length(cluster_colors) != ncol(log2fc_mat)) stop("Cannot resolve cluster color mismatch.")
} else {
  ordered_cluster_names <- intersect(names(cluster_colors), colnames(log2fc_mat))
  log2fc_mat <- log2fc_mat[, ordered_cluster_names, drop = FALSE]
  pval_mat <- pval_mat[, ordered_cluster_names, drop = FALSE]
}

anno_df$Category <- factor(anno_df$Category, levels = category_order, ordered = TRUE)

# Split data
mhc_features    <- rownames(anno_df[anno_df$MHC, ])
non_mhc_features <- rownames(anno_df[!anno_df$MHC, ])

log2fc_mat_mhc <- 2^(log2fc_mat[mhc_features, , drop = FALSE])  #convert this matrix to FC
pval_mat_mhc   <- pval_mat[mhc_features, , drop = FALSE]
anno_df_mhc    <- anno_df[mhc_features, , drop = FALSE]

log2fc_mat_non_mhc <- 2^(log2fc_mat[non_mhc_features, , drop = FALSE]) # convert this to FC
pval_mat_non_mhc   <- pval_mat[non_mhc_features, , drop = FALSE]
anno_df_non_mhc    <- anno_df[non_mhc_features, , drop = FALSE]

# Sort Non-MHC data correctly
if(nrow(anno_df_non_mhc) > 0) {
  anno_df_non_mhc_sorted <- anno_df_non_mhc %>%
    dplyr::arrange(Category, rownames(.))
  sorted_non_mhc_rownames <- rownames(anno_df_non_mhc_sorted)
  log2fc_mat_non_mhc <- log2fc_mat_non_mhc[sorted_non_mhc_rownames, , drop = FALSE]
  pval_mat_non_mhc <- pval_mat_non_mhc[sorted_non_mhc_rownames, , drop = FALSE]
  print("Non-MHC data sorted by Category, then Feature name.")
}

# Sort MHC data alphabetically by Feature name (rownames)
if(nrow(log2fc_mat_mhc) > 0) {
  mhc_sorted_index <- order(rownames(log2fc_mat_mhc))
  log2fc_mat_mhc <- log2fc_mat_mhc[mhc_sorted_index, , drop = FALSE]
  pval_mat_mhc <- pval_mat_mhc[mhc_sorted_index, , drop = FALSE]
  anno_df_mhc <- anno_df_mhc[mhc_sorted_index, , drop = FALSE]
  print("MHC data sorted alphabetically by Feature name.")
}

# --- Define Annotation Objects ---
cluster_annotation <- HeatmapAnnotation(
  Cluster = colnames(log2fc_mat),
  col = list(Cluster = cluster_colors[colnames(log2fc_mat)]),
  show_legend = FALSE, show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "cm")
)

# Row Annotation for MHC plot
if(nrow(anno_df_mhc) > 0) {
  mhc_anno_data <- factor(rep("MHC", nrow(anno_df_mhc)), levels = "MHC")
  mhc_annotation <- rowAnnotation(
    `Gene Type` = mhc_anno_data,
    col = list(`Gene Type` = all_type_colors["MHC"]),
    # --- FIX: Show annotation name for MHC plot ---
    show_annotation_name = FALSE,
    #annotation_name_side = "top", # Position label above the bar
    #annotation_name_gp = gpar(fontsize = 9, fontface="italic"), # Style the label
    # --- Keep legend hidden for this specific annotation ---
    show_legend = FALSE,
    simple_anno_size = unit(0.5, "cm")
  )
} else { mhc_annotation <- NULL }


# Row Annotation (Categories for Non-MHC)
if(nrow(anno_df_non_mhc) > 0) {
  # Ensure the Category factor in the sorted df includes all levels for the legend
  anno_df_non_mhc_sorted$CategoryForAnno <- factor(anno_df_non_mhc_sorted$Category, levels = all_type_order)
  
  category_annotation <- rowAnnotation(
    # Use the new factor with all levels for the annotation itself
    # The actual colors displayed will depend on the data present (non-MHC categories)
    Category = anno_df_non_mhc_sorted$CategoryForAnno,
    # Provide the full color map mapping levels to colors
    col = list(Category = all_type_colors),
    show_annotation_name = FALSE, # No label above the category bar itself
    simple_anno_size = unit(0.5, "cm"),
    # Define the combined legend appearance using the full color map and order
    annotation_legend_param = list(
      title = "Gene Type",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      # Use the full list of types for the legend labels and colors
      labels = names(all_type_colors[all_type_order]), # Ensure order
      legend_gp = gpar(fill = all_type_colors[all_type_order]) # Ensure order
    )
  )
} else { category_annotation <- NULL }

# --- Define P-value Star Function (same as before) ---
cell_fun_stars = function(j, i, x, y, width, height, fill, pval_matrix) {
  if (i > 0 && i <= nrow(pval_matrix) && j > 0 && j <= ncol(pval_matrix)) {
    p_val <- pval_matrix[i, j]
    if (!is.na(p_val)) {
      stars <- dplyr::case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "")
      if (stars != "") grid.text(stars, x, y, gp = gpar(fontsize = 10))
    }
  }
}

# --- Create MHC Heatmap ---
if (nrow(log2fc_mat_mhc) > 0) {
  print("Generating MHC Heatmap...")
  
  log2fc_mat_mhc['rs1064173',]<-log2fc_mat_mhc['rs1064173',]/2.5 # for better viewing contrast only
  
  mhc_limit <- max(abs(range(log2fc_mat_mhc, na.rm = TRUE)), 0.1)
  col_fun_mhc <- colorRamp2(c(0.5, 1, 1.5), mhc_color_scale)
  mhc_row_labels <- anno_df_mhc[rownames(log2fc_mat_mhc), "Gene"]
  
  ht_mhc <- Heatmap(
    log2fc_mat_mhc, name = "log2(FC) MHC", col = col_fun_mhc,
    cell_fun = function(j, i, x, y, width, h, fill) cell_fun_stars(j, i, x, y, width, h, fill, pval_mat_mhc),
    left_annotation = mhc_annotation, # Includes label now
    bottom_annotation = cluster_annotation,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_labels = mhc_row_labels,
    row_names_side = "left", row_names_gp = gpar(fontsize = 8),
    #row_split = anno_df_mhc_sorted$Category, # Split by original category factor
    row_title_gp = gpar(fontsize = 10, fontface = "bold"), row_title_rot = 0, gap = unit(2, "mm"),
    column_names_gp = gpar(fontsize = 9),
    column_title = "MHC SHAP Value Difference", column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.5),
    heatmap_legend_param = list(title = "SHAP Diff.", title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9))
  )
  
  mhc_plot_height <- 3
  mhc_plot_width <- 4.5
  png_file_mhc <- file.path(output_dir, "complexheatmap_mhc_v5.png") # New filename
  tryCatch({
    png(png_file_mhc, width = mhc_plot_width, height = mhc_plot_height, units = "in", res = 300)
    # Draw MHC plot - annotation legend is hidden via show_legend=FALSE in mhc_annotation
    draw(ht_mhc, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off()
    print(paste("MHC Heatmap saved to:", png_file_mhc))
  }, error = function(e) {
    print(paste("Error plotting/saving MHC heatmap:", e$message))
    if(names(dev.cur()) != "null device") { dev.off() }
  })
  
} else { print("No MHC features to plot.") }


# --- Create Non-MHC Heatmap ---
if (nrow(log2fc_mat_non_mhc) > 0) {
  print("Generating Non-MHC Heatmap...")
  
  #log2fc_mat_non_mhc['chr11:2160994:A:T',]<-log2fc_mat_non_mhc['chr11:2160994:A:T',]/1.66 # for better viewing contrast only
  
  non_mhc_limit <- max(abs(range(log2fc_mat_non_mhc, na.rm = TRUE)), 0.1)
  col_fun_non_mhc <- colorRamp2(c(0.95, 1, 1.05), non_mhc_color_scale)
  non_mhc_row_labels <- anno_df_non_mhc_sorted[rownames(log2fc_mat_non_mhc), "Gene"]
  
  ht_non_mhc <- Heatmap(
    log2fc_mat_non_mhc, name = "log2(FC) Non-MHC", col = col_fun_non_mhc,
    cell_fun = function(j, i, x, y, width, h, fill) cell_fun_stars(j, i, x, y, width, h, fill, pval_mat_non_mhc),
    left_annotation = category_annotation, # Annotation object now defines the combined legend
    bottom_annotation = cluster_annotation,
    cluster_rows = FALSE,
    row_split = anno_df_non_mhc_sorted$Category, # Split by original category factor
    #row_title_gp = gpar(fontsize = 10, fontface = "bold"), row_title_rot = 0, gap = unit(2, "mm"),
    #row_title_gp=gpar(NULL),
    cluster_columns = FALSE,
    row_labels = non_mhc_row_labels,
    row_names_side = "left", row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    column_title = "Non-MHC SHAP Value Difference  ", column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    border = TRUE, rect_gp = gpar(col = "grey85", lwd = 0.5),
    # --- FIX: Simplify heatmap legend parameters ---
    heatmap_legend_param = list(
      title = "SHAP Diff.",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
      # Let ComplexHeatmap handle breaks and labels automatically
      # , at = NULL, labels = NULL # Explicitly set to NULL if needed
    )
  )
  
  non_mhc_plot_height <- max(5, 0.22 * nrow(log2fc_mat_non_mhc))
  non_mhc_plot_width <- 4.5
  png_file_non_mhc <- file.path(output_dir, "complexheatmap_non_mhc_v5.png") # New filename
  tryCatch({
    png(png_file_non_mhc, width = non_mhc_plot_width, height = non_mhc_plot_height, units = "in", res = 300)
    # Draw Non-MHC plot - annotation legend is shown via category_annotation params
    draw(ht_non_mhc, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
    dev.off()
    print(paste("Non-MHC Heatmap saved to:", png_file_non_mhc))
  }, error = function(e) {
    print(paste("Error plotting/saving Non-MHC heatmap:", e$message))
    if(names(dev.cur()) != "null device") { dev.off() }
  })
  
} else { print("No Non-MHC features to plot.") }

print("R script finished.")








