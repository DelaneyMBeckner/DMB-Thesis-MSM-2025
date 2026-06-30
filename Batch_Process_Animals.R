# ============================================================
# BATCH PROCESSOR: Cross-Animal Comparison
# ============================================================
# Processes multiple animals and creates comparison plots
#
# Usage:
#   source("Batch_Process_Animals.R")
#   # Edit configuration section below
#   # Then run: results <- batch_process_animals()
# ============================================================

library(dplyr)
library(ggplot2)

# Load pipeline script
source("Pipeline_Transition_Analysis.R")

# ============================================================
# CONFIGURATION
# ============================================================

# Animals to analyze
ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9")

# Conditions (experimental days)
CONDITIONS <- c("BL", "SD", "WO")

# Transition types to include
TRANSITION_TYPES <- c("Wake2NREM", "NREM2Wake")

# Window size (3 or 9 epochs)
WINDOW_SIZE <- 9

# Data directories
DATA_DIR <- "E:/Data_Processing/R/Data CSVs"
OUTPUT_DIR <- "E:/Data_Processing/R/Results"

# Processing options
SAVE_INTERMEDIATE <- TRUE  # Save per-animal CSVs?
CCF_USE_PARALLEL <- TRUE   # Use parallel processing for CCF?
CCF_N_CORES <- 6           # Number of cores for CCF

# CCF analysis parameters
RUN_LAG_SWEEP <- FALSE     # Run lag sweep analysis? (slower)
MAX_LAG <- 5               # Maximum lag for lag sweep (if enabled)

# Distance-correlation analysis parameters
BINNING_METHOD <- "percentile"   # "equal_width", "percentile", or "both"
N_DISTANCE_BINS <- 10      # Number of bins for equal-width binning
N_PERCENTILE_BINS <- 5     # Number of bins for percentile binning
FIT_EXPONENTIAL <- FALSE   # Fit exponential model to distance-correlation?

# Plotting options
N_PERCENTILE_POINTS <- 100  # Number of points for continuous percentile plots

# ============================================================
# MAIN BATCH PROCESSING FUNCTION
# ============================================================

#' Process all animals and create comparison plots
batch_process_animals <- function(animals = ANIMALS,
                                  conditions = CONDITIONS,
                                  transition_types = TRANSITION_TYPES,
                                  window_size = WINDOW_SIZE,
                                  data_dir = DATA_DIR,
                                  output_dir = OUTPUT_DIR,
                                  save_intermediate = SAVE_INTERMEDIATE,
                                  ccf_use_parallel = CCF_USE_PARALLEL,
                                  ccf_n_cores = CCF_N_CORES,
                                  run_lag_sweep = RUN_LAG_SWEEP,
                                  max_lag = MAX_LAG,
                                  binning_method = BINNING_METHOD,
                                  n_distance_bins = N_DISTANCE_BINS,
                                  n_percentile_bins = N_PERCENTILE_BINS,
                                  fit_exponential = FIT_EXPONENTIAL,
                                  n_percentile_points = N_PERCENTILE_POINTS) {
  
  cat("====================================================\n")
  cat("  BATCH PROCESSING: CROSS-ANIMAL COMPARISON\n")
  cat("====================================================\n")
  cat("Animals:", paste(animals, collapse = ", "), "\n")
  cat("Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("Transitions:", paste(transition_types, collapse = ", "), "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("Output directory:", output_dir, "\n")
  cat("====================================================\n\n")
  
  # Storage for all animal results
  all_animal_results <- list()
  
  # Process each animal SEQUENTIALLY (parallel happens within CCF)
  for (animal in animals) {
    cat("\n\n")
    cat("####################################################\n")
    cat("# ANIMAL", match(animal, animals), "of", length(animals), "\n")
    cat("####################################################\n")
    
    tryCatch({
      # Run pipeline for this animal
      animal_results <- run_animal_pipeline(
        animal_id = animal,
        conditions = conditions,
        transition_types = transition_types,
        window_size = window_size,
        data_dir = data_dir,
        output_dir = output_dir,
        save_intermediate = save_intermediate,
        ccf_use_parallel = ccf_use_parallel,
        ccf_n_cores = ccf_n_cores,
        run_lag_sweep = run_lag_sweep,
        max_lag = max_lag,
        binning_method = binning_method,
        n_distance_bins = n_distance_bins,
        n_percentile_bins = n_percentile_bins,
        fit_exponential = fit_exponential
      )
      
      # Store results
      all_animal_results[[animal]] <- animal_results
      
    }, error = function(e) {
      warning("Error processing animal ", animal, ": ", e$message)
      cat("\nSkipping", animal, "due to error\n\n")
    })
  }
  
  # Check if we have results
  if (length(all_animal_results) == 0) {
    stop("No animals were successfully processed!")
  }
  
  cat("\n\n")
  cat("====================================================\n")
  cat("  ALL ANIMALS PROCESSED\n")
  cat("====================================================\n")
  cat("Successfully processed:", length(all_animal_results), "animals\n\n")
  
  # Combine results across animals
  cat("Combining results across animals...\n")
  all_event_rates <- bind_rows(lapply(all_animal_results, function(x) x$event_rates))
  all_correlations <- bind_rows(lapply(all_animal_results, function(x) x$correlations))
  all_dist_corr <- bind_rows(lapply(all_animal_results, function(x) x$dist_corr))
  
  cat("  Total ROI-condition combinations:", nrow(all_event_rates), "\n")
  cat("  Total pair-condition combinations:", nrow(all_correlations), "\n")
  cat("  Total distance-correlation pairs:", nrow(all_dist_corr), "\n\n")
  
  # Create comparison plots
  cat("Creating comparison plots...\n\n")
  
  plots <- list()
  
  # Plot 1: Event Rate by Activity Percentile
  cat("  [1/3] Event Rate by Activity Percentile...\n")
  plots$event_rate_percentile <- create_percentile_plot(
    data = all_event_rates,
    x_var = "activity_percentile",
    y_var = "mean_event_rate",
    title = "Event Rate by Activity Percentile Across Conditions",
    x_lab = "Activity Percentile (within animal)",
    y_lab = "Mean Event Rate (events/transition)",
    n_points = n_percentile_points
  )
  
  # Plot 2: Correlation by Correlation Percentile
  cat("  [2/3] Correlation by Correlation Percentile...\n")
  plots$correlation_percentile <- create_percentile_plot(
    data = all_correlations,
    x_var = "correlation_percentile",
    y_var = "mean_correlation",
    title = "Correlation by Correlation Percentile Across Conditions",
    x_lab = "Correlation Percentile (within animal)",
    y_lab = "Mean Correlation",
    n_points = n_percentile_points
  )
  
  # Plot 3: Correlation by Distance Percentile
  cat("  [3/3] Correlation by Distance Percentile...\n")
  plots$dist_corr_percentile <- create_distance_percentile_plot(
    data = all_dist_corr,
    title = "Correlation by Distance Percentile Across Conditions",
    x_lab = "Distance Percentile (within animal)",
    y_lab = "Mean Correlation",
    n_points = n_percentile_points
  )
  
  cat("  Plots created\n\n")
  
  # Save plots
  cat("Saving plots to", output_dir, "...\n")
  
  ggsave(file.path(output_dir, "CrossAnimal_EventRate_by_Percentile.png"),
         plots$event_rate_percentile,
         width = 12, height = 8, dpi = 300)
  
  ggsave(file.path(output_dir, "CrossAnimal_Correlation_by_Percentile.png"),
         plots$correlation_percentile,
         width = 12, height = 8, dpi = 300)
  
  ggsave(file.path(output_dir, "CrossAnimal_Correlation_by_Distance.png"),
         plots$dist_corr_percentile,
         width = 12, height = 8, dpi = 300)
  
  cat("  Plots saved\n\n")
  
  # Save combined data
  if (save_intermediate) {
    cat("Saving combined data tables...\n")
    
    write.csv(all_event_rates,
              file.path(output_dir, "CrossAnimal_EventRates_All.csv"),
              row.names = FALSE)
    
    write.csv(all_correlations,
              file.path(output_dir, "CrossAnimal_Correlations_All.csv"),
              row.names = FALSE)
    
    write.csv(all_dist_corr,
              file.path(output_dir, "CrossAnimal_DistCorr_All.csv"),
              row.names = FALSE)
    
    cat("  Data tables saved\n\n")
  }
  
  cat("====================================================\n")
  cat("  BATCH PROCESSING COMPLETE!\n")
  cat("====================================================\n\n")
  
  # Return everything
  return(list(
    animals = all_animal_results,
    event_rates = all_event_rates,
    correlations = all_correlations,
    dist_corr = all_dist_corr,
    plots = plots
  ))
}

# ============================================================
# PLOTTING HELPER FUNCTIONS
# ============================================================

#' Create percentile comparison plot
create_percentile_plot <- function(data, x_var, y_var, title, x_lab, y_lab, n_points = 100) {
  
  # Bin data into percentile groups for smooth lines
  percentile_bins <- seq(0, 100, length.out = n_points + 1)
  
  plot_data <- data %>%
    mutate(percentile_bin = cut(.data[[x_var]], 
                                breaks = percentile_bins,
                                include.lowest = TRUE,
                                labels = FALSE)) %>%
    group_by(animal, condition, percentile_bin) %>%
    summarize(
      percentile = mean(.data[[x_var]], na.rm = TRUE),
      metric_value = mean(.data[[y_var]], na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = percentile, y = metric_value,
                             color = animal, linetype = condition)) +
    geom_line(linewidth = 1) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
    ) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = "Animal",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

#' Create distance-percentile plot (special handling for distance-correlation data)
create_distance_percentile_plot <- function(data, title, x_lab, y_lab, n_points = 100) {
  
  # Bin data by distance percentile
  percentile_bins <- seq(0, 100, length.out = n_points + 1)
  
  plot_data <- data %>%
    mutate(percentile_bin = cut(distance_percentile,
                                breaks = percentile_bins,
                                include.lowest = TRUE,
                                labels = FALSE)) %>%
    group_by(animal, condition, percentile_bin) %>%
    summarize(
      distance_percentile = mean(distance_percentile, na.rm = TRUE),
      mean_correlation = mean(Correlation, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = distance_percentile, y = mean_correlation,
                             color = animal, linetype = condition)) +
    geom_line(linewidth = 1) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
    ) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = "Animal",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

# ============================================================
# AUTO-RUN (if sourced in interactive mode)
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  BATCH PROCESSOR LOADED\n")
  cat("========================================\n")
  cat("\nTo run batch processing, execute:\n")
  cat("  results <- batch_process_animals()\n\n")
  cat("To customize settings, edit the CONFIGURATION section above,\n")
  cat("or pass parameters to batch_process_animals()\n\n")
}
