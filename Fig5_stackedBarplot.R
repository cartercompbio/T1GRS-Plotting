# fig 5 stacked barplot
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats) # For factor manipulation
library(ggh4x)   # For advanced facet styling (colored strips)

# --- Define Aesthetics ---
mhc_col_r <- '#7061CF'
cd4_col_r <- '#FFC20A'
panc_col_r <- '#36AB9F'
inflam_col_r <- '#D00000'
other_col_r <- "grey"

# Colors are keyed by the display names in 'newnames'
category_plot_colors_r <- c(
  "MHC" = mhc_col_r,
  "T-Cells" = cd4_col_r,
  "Pancreas" = panc_col_r,
  "Innate Immune" = inflam_col_r,
  "Other" = other_col_r
)

# Cluster titles
cluster_titles_r <- c(
  `1` = "Cluster 1:\nT-Cell Regulatory Failure", 
  `2` = "Cluster 2:\nBroad Trigger",
  `3` = "Cluster 3:\nPancreatic Inflam.", 
  `4` = "Cluster 4:\nEarly-Onset MHC Driven"
)

# Banner colors with alpha (0.7 alpha is B3 in hex)
cluster_banner_colors_r_raw <- c(
  `1` = "#B38D97", `2` = "#437A91", 
  `3` = "#849324", `4` = "#F55F14"
)
cluster_banner_colors_r_alpha <- paste0(cluster_banner_colors_r_raw, "B3") # Appending B3 for 70% alpha
names(cluster_banner_colors_r_alpha) <- names(cluster_banner_colors_r_raw)


# Load data
final_normalized_df_to_save <- read.table('~/Code/T1D/data2025/shap_normalized_cat_vals.txt', sep=',', header=TRUE)

# --- Prepare Data ---
# Original column names from the file for categories
categories_in_file <- c("MHC", "CD4_Autoimmunity", "Beta_Dysfunction", "General_Inflammation", "Other") 

# Desired display names and order for plotting and legend
newnames <- c("MHC", "T-Cells", "Pancreas", "Innate Immune", "Other")

# Check if critical columns exist
critical_cols_check <- c("Cluster", categories_in_file)
if (!all(critical_cols_check %in% names(final_normalized_df_to_save))) {
  missing_cols <- critical_cols_check[!critical_cols_check %in% names(final_normalized_df_to_save)]
  stop(paste("One or more critical columns are missing from 'final_normalized_df_to_save'. Missing:", 
             paste(missing_cols, collapse=", ")))
}

plot_df_r <- final_normalized_df_to_save[, c("Cluster", categories_in_file)]

# Rename columns to the desired display names *before* pivoting
# Ensure the order of categories_in_file correctly maps to newnames
colnames(plot_df_r) <- c('Cluster', newnames) 

# Pivot to long format
plot_df_long_r <- plot_df_r %>%
  pivot_longer(cols = all_of(newnames), 
               names_to = "Category", # 'Category' will now have values from newnames
               values_to = "SHAP_Value") %>%
  mutate(
    SHAP_Value = ifelse(is.na(SHAP_Value), 0, SHAP_Value),
    # Ensure 'Category' is a factor with levels from 'newnames' for correct stacking and legend order
    Category = factor(Category, levels = newnames), 
    Cluster_Factor = factor(Cluster),
    dummy_x = "SHAP Contributions",
    # Label for segments (only if abs(SHAP_Value) > 0.05)
    Segment_Label = ifelse(abs(SHAP_Value) > 0.05, sprintf("%.2f", SHAP_Value), "")
  ) %>%
  # Arrange by Cluster_Factor, then by the factor levels of Category to ensure stacking order
  # This helps position_stack for geom_text to work correctly.
  # For positive values, categories with lower factor levels (earlier in 'newnames') are at the bottom.
  arrange(Cluster_Factor, Category)

# Data for total labels above/below each stacked bar
total_shap_labels <- plot_df_long_r %>%
  group_by(Cluster_Factor, dummy_x) %>%
  summarise(
    Total_SHAP = sum(SHAP_Value),
    # Determine the y-position for the total label based on the extent of the stack
    Label_Y_Position = ifelse(Total_SHAP >= 0, 
                              sum(SHAP_Value[SHAP_Value > 0]), # Top of positive stack
                              sum(SHAP_Value[SHAP_Value < 0])), # Bottom of negative stack
    .groups = 'drop'
  ) %>%
  mutate(
    Vjust = ifelse(Total_SHAP >= 0, -0.3, 1.3), # Nudge slightly above/below the bar
    Total_SHAP_Label_Text = sprintf("Total: %.2f", Total_SHAP)
  )

# Prepare strip styles for ggh4x with alpha
strip_background_elements_alpha <- lapply(levels(plot_df_long_r$Cluster_Factor), function(level) {
  element_rect(fill = cluster_banner_colors_r_alpha[level], color = "black", linewidth = 0.9,margin(r=-5,l=-5,b=-10))
})
names(strip_background_elements_alpha) <- levels(plot_df_long_r$Cluster_Factor)

# --- Create Plot ---
n_distinct_clusters <- n_distinct(plot_df_long_r$Cluster_Factor)

plot_df_long_r$Category<-factor(plot_df_long_r$Category,levels=c("T-Cells","Pancreas","Innate Immune","Other","MHC"))
true_stacked_bar_plot_r_final <- ggplot(plot_df_long_r, 
                                        aes(x = dummy_x, y = SHAP_Value, fill = Category)) +
  geom_col(position = "stack", width = 0.55) +
  # Add text labels within segments
  geom_text(aes(label = Segment_Label), 
            position = position_stack(vjust = 0.5), 
            size = 4.2, 
            color = "white", # Consider if white is always visible
            fontface = "bold") +
  # Add total SHAP value labels
  geom_text(data = total_shap_labels,
            aes(x = dummy_x, y = Label_Y_Position, label = Total_SHAP_Label_Text, fill = NULL), # fill=NULL to not inherit
            vjust = total_shap_labels$Vjust, 
            size = 4.2, 
            fontface = "bold") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5,linetype = 2) +
  scale_fill_manual(
    values = category_plot_colors_r, # Should use 'newnames' as keys if colors are defined that way
    labels = newnames, # Labels in legend
    name = "SNP Loci Category",
    guide = guide_legend(reverse = TRUE) # Reverses order in legend
  ) +
  ggh4x::facet_wrap2(
    ~ Cluster_Factor, 
    nrow = 1, 
    strip.position = "top",
    labeller = labeller(Cluster_Factor = cluster_titles_r),
    strip = strip_themed(background_x = strip_background_elements_alpha) 
  ) +
  labs(
    x = NULL, 
    y = "Net SHAP Value" 
  ) +
  theme_bw(base_size = 14) + 
  theme(
    panel.border = element_blank(), 
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.8), 
    axis.line.y.left = element_line(colour = "black", linewidth = 0.8), 
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    
    strip.text = element_text(size = 12, face = "bold", color = "black", margin = margin(t=6, b=6,r=-5,l=-5)), # Adjusted size
    
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14,color='black'), 
    axis.title.y = element_text(margin = margin(r = 10), size=14, face="bold"),
    
    panel.spacing.x = unit(-0.03, "lines"), 
    legend.position = "bottom",
    legend.title = element_text(size=14, face="bold"), # Adjusted size
    legend.text = element_text(size=14,color='black')
  )

print(true_stacked_bar_plot_r_final)

# Save the plot (optional)
output_dir <- "~/Code/T1D/figs_may/fig5/" # Define your output directory
ggsave(file.path(output_dir, "fig5_shap_true_stacked_bars_R_with_labels_alpha.pdf"), 
        plot = true_stacked_bar_plot_r_final, 
        width = 10, 
        height = 10, 
        units = "in")
# print(paste0("R plot saved to: ", file.path(output_dir, "fig5_shap_true_stacked_bars_R_with_labels_alpha.pdf")))




