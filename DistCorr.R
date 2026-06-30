library(tidyverse)
library(ggplot2)
library(broom)
library(corrr)
library(viridis)

# Function to load and process ROI data
load_and_process_data <- function(distance_file, correlation_file, full_field_roi = "C74") {
  
  cat("Loading distance data from:", distance_file, "\n")
  distance_df <- read.csv(distance_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(distance_df), "distance pairs\n")
  
  cat("Loading correlation data from:", correlation_file, "\n")
  corr_df <- read.csv(correlation_file, row.names = 1, stringsAsFactors = FALSE)
  cat("Loaded correlation matrix:", dim(corr_df), "\n")
  
  # Convert correlation matrix to long format
  corr_matrix <- as.matrix(corr_df)
  
  # Get upper triangle indices to avoid duplicates
  upper_tri_indices <- which(upper.tri(corr_matrix, diag = FALSE), arr.ind = TRUE)
  
  correlation_pairs <- data.frame(
    ROI_1 = rownames(corr_matrix)[upper_tri_indices[, 1]],
    ROI_2 = colnames(corr_matrix)[upper_tri_indices[, 2]],
    Correlation = corr_matrix[upper_tri_indices],
    stringsAsFactors = FALSE
  )
  
  cat("Created", nrow(correlation_pairs), "correlation pairs\n")
  
  # Handle Full_Field naming difference
  distance_df$ROI_1 <- gsub("Full_Field", full_field_roi, distance_df$ROI_1)
  distance_df$ROI_2 <- gsub("Full_Field", full_field_roi, distance_df$ROI_2)
  
  # Create order-independent pair keys for merging
  create_pair_key <- function(roi1, roi2) {
    paste(sort(c(roi1, roi2)), collapse = "_")
  }
  
  distance_df$pair_key <- mapply(create_pair_key, distance_df$ROI_1, distance_df$ROI_2)
  correlation_pairs$pair_key <- mapply(create_pair_key, correlation_pairs$ROI_1, correlation_pairs$ROI_2)
  
  # Merge the datasets
  merged_df <- merge(distance_df, correlation_pairs, by = "pair_key", suffixes = c("_dist", "_corr"))
  
  cat("Successfully merged", nrow(merged_df), "ROI pairs\n")
  
  return(merged_df[, c("Distance", "Correlation")])
}

# Function to plot distance vs correlation with regression
plot_distance_correlation <- function(data, title = "ROI Distance vs Cross-Correlation") {
  
  # Remove NaN values
  clean_data <- data[complete.cases(data), ]
  
  # Fit linear model
  lm_fit <- lm(Correlation ~ Distance, data = clean_data)
  lm_summary <- summary(lm_fit)
  
  # Create the plot
  final_plot <- ggplot(clean_data, aes(x = Distance, y = Correlation)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(
      title = title,
      x = "Distance (μm)",
      y = "Cross-Correlation"
    ) +
    annotate("text", 
             x = max(clean_data$Distance) * 0.8, 
             y = max(clean_data$Correlation) * 0.8, 
             label = paste0("R² = ", round(lm_summary$r.squared, 4),
                            "\np = ", format(lm_summary$coefficients[2,4], scientific = TRUE, digits = 2)),
             size = 4, color = "red") +
    theme_classic()
  
  print(final_plot)
  
  # Print detailed statistics
  cat("\nLinear Regression Results:\n")
  cat("========================\n")
  cat("Slope:", format(lm_summary$coefficients[2,1], scientific = TRUE, digits = 6), "\n")
  cat("Intercept:", format(lm_summary$coefficients[1,1], digits = 6), "\n")
  cat("R-squared:", round(lm_summary$r.squared, 6), "\n")
  cat("P-value:", format(lm_summary$coefficients[2,4], scientific = TRUE, digits = 3), "\n")
  cat("Sample size:", nrow(clean_data), "\n")
  
  return(lm_fit)

}

# Function to analyze multiple animals
analyze_multiple_animals <- function(distance_files, correlation_files, animal_names = NULL) {
  
  if (is.null(animal_names)) {
    animal_names <- paste0("Animal_", 1:length(distance_files))
  }
  
  all_data <- list()
  results <- list()
  
  # Process each animal
  for (i in seq_along(distance_files)) {
    cat("\nProcessing", animal_names[i], "...\n")
    cat("===============================\n")
    
    # Load and process data
    data <- load_and_process_data(distance_files[i], correlation_files[i])
    data$Animal <- animal_names[i]
    all_data[[i]] <- data
    
    # Fit linear model
    clean_data <- data[complete.cases(data), ]
    lm_fit <- lm(Correlation ~ Distance, data = clean_data)
    lm_summary <- summary(lm_fit)
    
    # Store results
    results[[i]] <- data.frame(
      Animal = animal_names[i],
      Slope = lm_summary$coefficients[2, 1],
      Intercept = lm_summary$coefficients[1, 1],
      R_squared = lm_summary$r.squared,
      P_value = lm_summary$coefficients[2, 4],
      N = nrow(clean_data)
    )
  }
  
  # Combine all data
  combined_data <- do.call(rbind, all_data)
  results_df <- do.call(rbind, results)
  
  # Create combined plot
  colors <- viridis_d(n = length(animal_names))
  
  p_combined <- ggplot(combined_data, aes(x = Distance, y = Correlation, color = Animal)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, size = 1.2) +
    scale_color_viridis_d() +
    labs(
      title = "ROI Distance vs Cross-Correlation - All Animals",
      x = "Distance (μm)",
      y = "Cross-Correlation"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    ) +
    guides(color = guide_legend(title = "Animal"))
  
  print(p_combined)
  
  # Print results summary
  cat("\nSummary of Results:\n")
  cat("==================\n")
  print(results_df)
  
  # Statistical comparison between animals (if more than 1)
  if (length(animal_names) > 1) {
    cat("\nANCOVA - Testing for differences between animals:\n")
    cat("================================================\n")
    
    # ANCOVA to test if slopes differ between animals
    ancova_model <- lm(Correlation ~ Distance * Animal, data = combined_data)
    ancova_summary <- anova(ancova_model)
    print(ancova_summary)
    
    # Test if slopes are significantly different from zero overall
    overall_model <- lm(Correlation ~ Distance, data = combined_data)
    overall_summary <- summary(overall_model)
    
    cat("\nOverall regression (all animals combined):\n")
    cat("=========================================\n")
    cat("R-squared:", round(overall_summary$r.squared, 6), "\n")
    cat("P-value:", format(overall_summary$coefficients[2, 4], scientific = TRUE, digits = 3), "\n")
    cat("Slope:", format(overall_summary$coefficients[2, 1], scientific = TRUE, digits = 6), "\n")
  }
  
  return(list(
    results = results_df,
    combined_data = combined_data,
    plots = p_combined
  ))
}

# Function to create additional diagnostic plots
create_diagnostic_plots <- function(model_fit, data, title_prefix = "") {
  
  # Extract fitted values and residuals
  fitted_vals <- fitted(model_fit)
  residuals_vals <- residuals(model_fit)
  
  # Create diagnostic plots
  par(mfrow = c(2, 2))
  
  # 1. Residuals vs Fitted
  plot(fitted_vals, residuals_vals, 
       main = paste(title_prefix, "Residuals vs Fitted"),
       xlab = "Fitted Values", ylab = "Residuals",
       pch = 16, col = "steelblue")
  abline(h = 0, col = "red", lwd = 2)
  
  # 2. Q-Q plot
  qqnorm(residuals_vals, main = paste(title_prefix, "Q-Q Plot"),
         pch = 16, col = "steelblue")
  qqline(residuals_vals, col = "red", lwd = 2)
  
  # 3. Scale-Location plot
  sqrt_abs_residuals <- sqrt(abs(residuals_vals))
  plot(fitted_vals, sqrt_abs_residuals,
       main = paste(title_prefix, "Scale-Location"),
       xlab = "Fitted Values", ylab = "√|Residuals|",
       pch = 16, col = "steelblue")
  
  # 4. Cook's distance
  cooks_d <- cooks.distance(model_fit)
  plot(cooks_d, main = paste(title_prefix, "Cook's Distance"),
       ylab = "Cook's Distance", pch = 16, col = "steelblue")
  abline(h = 4/length(cooks_d), col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
}

# Example usage and main analysis
main_analysis <- function() {
  
  cat("ROI Distance vs Cross-Correlation Analysis\n")
  cat("=========================================\n")
  
  # Single animal analysis
  cat("\nSingle Animal Analysis:\n")
  cat("======================\n")
  
  # Replace with your actual file paths
  distance_file <- "E:\\Data_Processing\\Python\\Results\\mPFCf5_roi_distances.csv"
  correlation_file <- "E:\\Data_Processing\\R\\results\\cross_p_results_mPFCf5_BL\\mPFCf5_BL_correlations_NREM.csv"
  
  # Load and process data
  data <- load_and_process_data(distance_file, correlation_file)
  
  # Create regression plot
  model <- plot_distance_correlation(data, title = "mPFCf5 - Distance vs Cross-Correlation")
  
  # Create diagnostic plots
  create_diagnostic_plots(model, data, "mPFCf5 - ")
  
  # Multiple animals example (uncomment and modify as needed)
  # cat("\nMultiple Animals Analysis:\n")
  # cat("=========================\n")
  # 
  # distance_files <- c(
  #   "animal1_roi_distances.csv",
  #   "animal2_roi_distances.csv",
  #   "animal3_roi_distances.csv"
  # )
  # 
  # correlation_files <- c(
  #   "animal1_BL_correlations_NREM.csv", 
  #   "animal2_BL_correlations_NREM.csv",
  #   "animal3_BL_correlations_NREM.csv"
  # )
  # 
  # animal_names <- c("Animal_1", "Animal_2", "Animal_3")
  # 
  # results <- analyze_multiple_animals(distance_files, correlation_files, animal_names)
  
  return(model)
}

# Run the analysis
if (interactive()) {
  model_result <- main_analysis()
}