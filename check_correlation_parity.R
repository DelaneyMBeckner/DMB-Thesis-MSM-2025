# Check_Negative_Correlations.R
# Reports % of negative, near-zero, and positive pairwise correlations
# Separates by transition type with simple visualization

library(dplyr)
library(ggplot2)
library(tidyr)

# ============ CONFIGURATION ============
# Point this to your Results folder
RESULTS_DIR <- "E:/Data_Processing/R/Results"
OUTPUT_DIR <- "E:/Data_Processing/R/Results"

# Animals and conditions to check
ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm9", "mPFCm4")  # transgenic only
CONDITIONS <- c("BL", "SD", "WO")
TRANSITION_TYPES <- c("Wake2NREM", "NREM2Wake")
WINDOW_SIZE <- 9  # or 3

# Threshold for "near zero" (absolute value below this = near zero)
NEAR_ZERO_THRESHOLD <- 0.05

# File suffix (e.g., "", "_all_simple", "_all_simple_filtered")
# Goes between "9ep" and "_ccf_trajectory"
FILE_SUFFIX <- "_all_simple"

# Save outputs?
SAVE_OUTPUT <- FALSE
# =======================================


#' Analyze correlations for one transition type
analyze_transition_type <- function(transition_type, animals = ANIMALS, 
                                    conditions = CONDITIONS, window_size = WINDOW_SIZE,
                                    results_dir = RESULTS_DIR, near_zero_thresh = NEAR_ZERO_THRESHOLD,
                                    file_suffix = FILE_SUFFIX) {
  
  cat("\n============================================\n")
  cat("  ", transition_type, "\n")
  cat("============================================\n\n")
  
  all_corrs <- data.frame()
  stats_summary <- data.frame()
  
  for (animal in animals) {
    for (condition in conditions) {
      
      # Build filename pattern
      recording_id <- paste0(animal, "_", condition)
      filename <- file.path(results_dir, 
                            paste0(recording_id, "_", transition_type, "_", 
                                   window_size, "ep", file_suffix, "_ccf_trajectory.csv"))
      
      # Check if file exists
      if (!file.exists(filename)) {
        cat("File not found:", basename(filename), "\n")
        next
      }
      
      # Load trajectory data
      traj <- read.csv(filename, stringsAsFactors = FALSE)
      
      # Check for correlation column
      if (!"correlation" %in% colnames(traj)) {
        cat("No 'correlation' column in:", basename(filename), "\n")
        next
      }
      
      # Add metadata
      traj$animal <- animal
      traj$condition <- condition
      traj$transition_type <- transition_type
      
      # Categorize correlations
      traj$corr_category <- case_when(
        abs(traj$correlation) < near_zero_thresh ~ "near_zero",
        traj$correlation < 0 ~ "negative",
        traj$correlation > 0 ~ "positive",
        TRUE ~ "zero"
      )
      
      all_corrs <- rbind(all_corrs, traj)
      
      # Calculate stats for this recording
      n_total <- nrow(traj)
      n_negative <- sum(traj$corr_category == "negative", na.rm = TRUE)
      n_near_zero <- sum(traj$corr_category == "near_zero", na.rm = TRUE)
      n_positive <- sum(traj$corr_category == "positive", na.rm = TRUE)
      
      stats_summary <- rbind(stats_summary, data.frame(
        animal = animal,
        condition = condition,
        transition_type = transition_type,
        n_total = n_total,
        n_negative = n_negative,
        n_near_zero = n_near_zero,
        n_positive = n_positive,
        pct_negative = (n_negative / n_total) * 100,
        pct_near_zero = (n_near_zero / n_total) * 100,
        pct_positive = (n_positive / n_total) * 100,
        min_corr = min(traj$correlation, na.rm = TRUE),
        max_corr = max(traj$correlation, na.rm = TRUE),
        mean_corr = mean(traj$correlation, na.rm = TRUE)
      ))
      
      cat(sprintf("%s: %.1f%% neg | %.1f%% near-zero | %.1f%% pos\n", 
                  recording_id, 
                  (n_negative / n_total) * 100,
                  (n_near_zero / n_total) * 100,
                  (n_positive / n_total) * 100))
    }
  }
  
  return(list(
    correlations = all_corrs,
    stats = stats_summary
  ))
}


#' Create summary plots
create_plots <- function(all_data, near_zero_thresh = NEAR_ZERO_THRESHOLD) {
  
  plots <- list()
  
  # 1. Histogram of all correlations, faceted by transition type
  plots$histogram <- ggplot(all_data, aes(x = correlation, fill = transition_type)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-near_zero_thresh, near_zero_thresh), 
               linetype = "dotted", color = "gray40", linewidth = 0.5) +
    facet_wrap(~transition_type, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Wake2NREM" = "#2E86AB", "NREM2Wake" = "#A23B72")) +
    labs(
      title = "Distribution of Pairwise Correlations",
      subtitle = sprintf("Dotted lines: near-zero threshold (|r| < %.2f)", near_zero_thresh),
      x = "Correlation (r)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12, face = "bold"))
  
  # 2. Histogram faceted by animal and condition
  plots$histogram_by_animal <- ggplot(all_data, aes(x = correlation, fill = condition)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    facet_grid(animal ~ transition_type, scales = "free_y") +
    scale_fill_manual(values = c("BL" = "#4ECDC4", "SD" = "#FF6B6B", "WO" = "#95E1D3")) +
    labs(
      title = "Correlation Distribution by Animal and Transition Type",
      x = "Correlation (r)",
      y = "Count",
      fill = "Condition"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 10, face = "bold"))
  
  # 3. Category breakdown bar chart
  category_summary <- all_data %>%
    group_by(transition_type, corr_category) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(transition_type) %>%
    mutate(pct = n / sum(n) * 100)
  
  # Order categories
  category_summary$corr_category <- factor(category_summary$corr_category,
                                           levels = c("negative", "near_zero", "positive"))
  
  plots$category_bars <- ggplot(category_summary, 
                                aes(x = transition_type, y = pct, fill = corr_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), 
              position = position_stack(vjust = 0.5), size = 4, color = "white") +
    scale_fill_manual(
      values = c("negative" = "#E74C3C", "near_zero" = "#95A5A6", "positive" = "#27AE60"),
      labels = c("Negative", "Near Zero", "Positive")
    ) +
    labs(
      title = "Correlation Categories by Transition Type",
      subtitle = sprintf("Near-zero defined as |r| < %.2f", near_zero_thresh),
      x = "",
      y = "Percentage",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 12, face = "bold"))
  
  # 4. Category breakdown by condition
  category_by_cond <- all_data %>%
    group_by(transition_type, condition, corr_category) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(transition_type, condition) %>%
    mutate(pct = n / sum(n) * 100)
  
  category_by_cond$corr_category <- factor(category_by_cond$corr_category,
                                           levels = c("negative", "near_zero", "positive"))
  
  plots$category_by_condition <- ggplot(category_by_cond, 
                                        aes(x = condition, y = pct, fill = corr_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(~transition_type) +
    scale_fill_manual(
      values = c("negative" = "#E74C3C", "near_zero" = "#95A5A6", "positive" = "#27AE60"),
      labels = c("Negative", "Near Zero", "Positive")
    ) +
    labs(
      title = "Correlation Categories by Condition",
      x = "Condition",
      y = "Percentage",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"))
  
  # 5. Category breakdown by animal
  category_by_animal <- all_data %>%
    group_by(transition_type, animal, corr_category) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(transition_type, animal) %>%
    mutate(pct = n / sum(n) * 100)
  
  category_by_animal$corr_category <- factor(category_by_animal$corr_category,
                                             levels = c("negative", "near_zero", "positive"))
  
  plots$category_by_animal <- ggplot(category_by_animal, 
                                     aes(x = animal, y = pct, fill = corr_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", pct)), 
              position = position_stack(vjust = 0.5), size = 3, color = "white") +
    facet_wrap(~transition_type) +
    scale_fill_manual(
      values = c("negative" = "#E74C3C", "near_zero" = "#95A5A6", "positive" = "#27AE60"),
      labels = c("Negative", "Near Zero", "Positive")
    ) +
    labs(
      title = "Correlation Categories by Animal",
      subtitle = sprintf("Near-zero defined as |r| < %.2f", near_zero_thresh),
      x = "Animal",
      y = "Percentage",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plots)
}


# ==================== MAIN ====================

cat("============================================\n")
cat("  CORRELATION SIGN ANALYSIS\n")
cat("  Near-zero threshold: |r| <", NEAR_ZERO_THRESHOLD, "\n")
cat("============================================\n")
cat("\nLooking for files matching pattern:\n")
cat("  {animal}_{condition}_{transition}_{window}ep", FILE_SUFFIX, "_ccf_trajectory.csv\n", sep = "")
cat("  e.g., mPFCf5_BL_Wake2NREM_", WINDOW_SIZE, "ep", FILE_SUFFIX, "_ccf_trajectory.csv\n\n", sep = "")

# Analyze each transition type
results_w2n <- analyze_transition_type("Wake2NREM")
results_n2w <- analyze_transition_type("NREM2Wake")

# Combine all data
all_corrs <- rbind(results_w2n$correlations, results_n2w$correlations)
all_stats <- rbind(results_w2n$stats, results_n2w$stats)

# Overall summary
cat("\n============================================\n")
cat("  OVERALL SUMMARY\n")
cat("============================================\n\n")

# Check if any data was loaded
if (nrow(all_corrs) == 0) {
  cat("ERROR: No data loaded. Check that:\n")
  cat("  1. RESULTS_DIR path is correct\n")
  cat("  2. FILE_SUFFIX matches your files\n")
  cat("  3. Files exist for specified ANIMALS, CONDITIONS, TRANSITION_TYPES\n")
  cat("\nExpected filename pattern:\n")
  cat("  {animal}_{condition}_{transition}_{window}ep{suffix}_ccf_trajectory.csv\n")
  cat("  e.g., mPFCf5_BL_Wake2NREM_9ep_all_simple_ccf_trajectory.csv\n")
  stop("No data to analyze")
}

for (trans_type in TRANSITION_TYPES) {
  subset_data <- all_corrs %>% filter(.data$transition_type == trans_type)
  n_total <- nrow(subset_data)
  
  if (n_total == 0) next
  
  n_neg <- sum(subset_data$corr_category == "negative")
  n_zero <- sum(subset_data$corr_category == "near_zero")
  n_pos <- sum(subset_data$corr_category == "positive")
  
  cat(sprintf("%s (n = %d pairs):\n", trans_type, n_total))
  cat(sprintf("  Negative (r < -%.2f):     %5d (%5.1f%%)\n", 
              NEAR_ZERO_THRESHOLD, n_neg, n_neg/n_total*100))
  cat(sprintf("  Near-zero (|r| < %.2f):   %5d (%5.1f%%)\n", 
              NEAR_ZERO_THRESHOLD, n_zero, n_zero/n_total*100))
  cat(sprintf("  Positive (r > %.2f):      %5d (%5.1f%%)\n", 
              NEAR_ZERO_THRESHOLD, n_pos, n_pos/n_total*100))
  cat(sprintf("  Range: [%.4f to %.4f]\n", 
              min(subset_data$correlation), max(subset_data$correlation)))
  cat(sprintf("  Mean: %.4f\n\n", mean(subset_data$correlation)))
}

# Combined totals
cat("--- COMBINED (both transition types) ---\n")
n_total <- nrow(all_corrs)
n_neg <- sum(all_corrs$corr_category == "negative")
n_zero <- sum(all_corrs$corr_category == "near_zero")
n_pos <- sum(all_corrs$corr_category == "positive")

cat(sprintf("Total pairs: %d\n", n_total))
cat(sprintf("  Negative:   %5d (%5.1f%%)\n", n_neg, n_neg/n_total*100))
cat(sprintf("  Near-zero:  %5d (%5.1f%%)\n", n_zero, n_zero/n_total*100))
cat(sprintf("  Positive:   %5d (%5.1f%%)\n", n_pos, n_pos/n_total*100))

cat("\n--- FOR METHODS SECTION ---\n")
cat(sprintf("\"Pairwise correlations ranged from %.3f to %.3f (mean = %.3f).\n",
            min(all_corrs$correlation), max(all_corrs$correlation), mean(all_corrs$correlation)))
cat(sprintf("Of all ROI pairs, %.1f%% showed negative correlations, %.1f%% were\n",
            n_neg/n_total*100, n_zero/n_total*100))
cat(sprintf("near-zero (|r| < %.2f), and %.1f%% were positive.\"\n",
            NEAR_ZERO_THRESHOLD, n_pos/n_total*100))

# Create plots
cat("\nGenerating plots...\n")
plots <- create_plots(all_corrs, NEAR_ZERO_THRESHOLD)

# Display plots
print(plots$histogram)
print(plots$category_bars)
print(plots$histogram_by_animal)
print(plots$category_by_condition)
print(plots$category_by_animal)

# Save outputs
if (SAVE_OUTPUT && nrow(all_corrs) > 0) {
  cat("\nSaving outputs...\n")
  
  # Save stats summary
  stats_file <- file.path(OUTPUT_DIR, "correlation_sign_summary.csv")
  write.csv(all_stats, stats_file, row.names = FALSE)
  cat("  Saved:", stats_file, "\n")
  
  # Save plots
  ggsave(file.path(OUTPUT_DIR, "correlation_histogram.png"), 
         plots$histogram, width = 8, height = 8, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "correlation_categories.png"), 
         plots$category_bars, width = 8, height = 6, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "correlation_histogram_by_animal.png"), 
         plots$histogram_by_animal, width = 12, height = 8, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "correlation_categories_by_condition.png"), 
         plots$category_by_condition, width = 10, height = 6, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "correlation_categories_by_animal.png"), 
         plots$category_by_animal, width = 10, height = 6, dpi = 300)
  cat("  Saved plots\n")
}

cat("\nDone!\n")

