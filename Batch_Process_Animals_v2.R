# ============================================================
# BATCH PROCESSOR: Cross-Animal Comparison (v2)
# ============================================================
# Processes multiple animals and creates comparison plots
# - Separates analyses by transition type
# - Supports trajectory analysis
#
# Usage:
#   source("Batch_Process_Animals_v2.R")
#   # Edit configuration section below
#   # Then run: results <- batch_process_animals()
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Load pipeline script
source("Pipeline_Transition_Analysis_v2.R")

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
WINDOW_SIZE <- 3

# Data directories
DATA_DIR <- "E:/Data_Processing/R/Data CSVs"
OUTPUT_DIR <- "E:/Data_Processing/R/Results"

# Processing options
SAVE_INTERMEDIATE <- FALSE  # Save per-animal CSVs?
CCF_USE_PARALLEL <- TRUE   # Use parallel processing for CCF?
CCF_N_CORES <- 6           # Number of cores for CCF

# Analysis modes
RUN_TRAJECTORY_ANALYSIS <- TRUE   # Run epoch-level trajectory analysis?
RUN_WHOLE_WINDOW <- FALSE          # Run whole-window analysis?

# CCF analysis parameters
RUN_LAG_SWEEP <- FALSE     # Run lag sweep analysis? (slower)
MAX_LAG <- 5               # Maximum lag for lag sweep (if enabled)

# Distance-correlation analysis parameters
BINNING_METHOD <- "percentile"   # "equal_width", "percentile", or "both"
N_DISTANCE_BINS <- 10      # Number of bins for equal-width binning
N_PERCENTILE_BINS <- 5     # Number of bins for percentile binning
N_TRAJECTORY_BINS <- 5     # Number of bins for trajectory analysis
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
                                  run_trajectory_analysis = RUN_TRAJECTORY_ANALYSIS,
                                  run_whole_window = RUN_WHOLE_WINDOW,
                                  run_lag_sweep = RUN_LAG_SWEEP,
                                  max_lag = MAX_LAG,
                                  binning_method = BINNING_METHOD,
                                  n_distance_bins = N_DISTANCE_BINS,
                                  n_percentile_bins = N_PERCENTILE_BINS,
                                  n_trajectory_bins = N_TRAJECTORY_BINS,
                                  fit_exponential = FIT_EXPONENTIAL,
                                  n_percentile_points = N_PERCENTILE_POINTS) {
  
  cat("====================================================\n")
  cat("  BATCH PROCESSING: CROSS-ANIMAL COMPARISON\n")
  cat("====================================================\n")
  cat("Animals:", paste(animals, collapse = ", "), "\n")
  cat("Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("Transitions:", paste(transition_types, collapse = ", "), "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("Trajectory analysis:", run_trajectory_analysis, "\n")
  cat("Whole window analysis:", run_whole_window, "\n")
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
        run_trajectory_analysis = run_trajectory_analysis,
        run_whole_window = run_whole_window,
        run_lag_sweep = run_lag_sweep,
        max_lag = max_lag,
        binning_method = binning_method,
        n_distance_bins = n_distance_bins,
        n_percentile_bins = n_percentile_bins,
        n_trajectory_bins = n_trajectory_bins,
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
  
  # Whole-window results
  all_event_rates <- NULL
  all_correlations <- NULL
  all_dist_corr <- NULL
  
  if (run_whole_window) {
    all_event_rates <- bind_rows(lapply(all_animal_results, function(x) x$event_rates))
    all_correlations <- bind_rows(lapply(all_animal_results, function(x) x$correlations))
    all_dist_corr <- bind_rows(lapply(all_animal_results, function(x) x$dist_corr))
    
    cat("  Total ROI-condition combinations:", nrow(all_event_rates), "\n")
    cat("  Total pair-condition combinations:", nrow(all_correlations), "\n")
    cat("  Total distance-correlation pairs:", nrow(all_dist_corr), "\n")
  }
  
  # Trajectory results
  all_trajectory_data <- NULL
  all_trajectory_dist <- NULL
  all_trajectory_events <- NULL
  
  if (run_trajectory_analysis) {
    all_trajectory_data <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_data))
    all_trajectory_dist <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_dist))
    
    # Collect event rate trajectory data (if available from pipeline)
    all_trajectory_events <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_events))
    
    cat("  Total trajectory epoch observations:", nrow(all_trajectory_data), "\n")
    cat("  Total trajectory distance-binned observations:", nrow(all_trajectory_dist), "\n")
    if (!is.null(all_trajectory_events) && nrow(all_trajectory_events) > 0) {
      cat("  Total trajectory event observations:", nrow(all_trajectory_events), "\n")
    }
  }
  
  cat("\n")
  
  # Create comparison plots
  cat("Creating comparison plots...\n\n")
  
  plots <- list()
  
  # ====== WHOLE-WINDOW PLOTS (separated by transition type) ======
  if (run_whole_window) {
    cat("  [WHOLE-WINDOW PLOTS]\n")
    
    for (trans_type in transition_types) {
      cat("    Processing", trans_type, "...\n")
      
      # Filter data for this transition type
      event_rates_trans <- all_event_rates %>% filter(transition_type == trans_type)
      correlations_trans <- all_correlations %>% filter(transition_type == trans_type)
      dist_corr_trans <- all_dist_corr %>% filter(transition_type == trans_type)
      
      # Plot 1: Event Rate by Activity Percentile
      plots[[paste0(trans_type, "_event_rate_percentile")]] <- create_percentile_plot(
        data = event_rates_trans,
        x_var = "activity_percentile",
        y_var = "mean_event_rate",
        title = paste(trans_type, "- Event Rate by Activity Percentile"),
        x_lab = "Activity Percentile (within animal)",
        y_lab = "Mean Event Rate (events/transition)",
        n_points = n_percentile_points
      )
      
      # Plot 2: Correlation by Correlation Percentile
      plots[[paste0(trans_type, "_correlation_percentile")]] <- create_percentile_plot(
        data = correlations_trans,
        x_var = "correlation_percentile",
        y_var = "mean_correlation",
        title = paste(trans_type, "- Correlation by Correlation Percentile"),
        x_lab = "Correlation Percentile (within animal)",
        y_lab = "Mean Correlation",
        n_points = n_percentile_points
      )
      
      # Plot 3: Correlation by Distance Percentile
      plots[[paste0(trans_type, "_dist_corr_percentile")]] <- create_distance_percentile_plot(
        data = dist_corr_trans,
        title = paste(trans_type, "- Correlation by Distance Percentile"),
        x_lab = "Distance Percentile (within animal)",
        y_lab = "Mean Correlation",
        n_points = n_percentile_points
      )
    }
  }
  
  # ====== TRAJECTORY PLOTS (separate for each transition type) ======
  if (run_trajectory_analysis) {
    cat("\n  [TRAJECTORY PLOTS]\n")
    
    # Overall trajectory comparison (combined + individual per transition type)
    cat("    Creating overall trajectory comparisons...\n")
    trajectory_overall_plots <- create_trajectory_comparison_plots(
      data = all_trajectory_data,
      title_base = "Correlation Trajectory Across Animals",
      window_size = window_size
    )
    for (plot_name in names(trajectory_overall_plots)) {
      plots[[paste0("trajectory_overall_", plot_name)]] <- trajectory_overall_plots[[plot_name]]
    }
    
    # Distance-binned trajectory comparison (combined + individual per animal/transition)
    cat("    Creating distance-binned trajectory comparisons...\n")
    trajectory_distance_plots <- create_trajectory_distance_comparison_plots(
      data = all_trajectory_dist,
      title_base = "Correlation Trajectory by Distance Bin",
      window_size = window_size
    )
    for (plot_name in names(trajectory_distance_plots)) {
      plots[[paste0("trajectory_by_distance_", plot_name)]] <- trajectory_distance_plots[[plot_name]]
    }
    
    # Faceted trajectory comparison (combined + individual per distance bin/transition)
    cat("    Creating faceted trajectory comparisons...\n")
    trajectory_faceted_plots <- create_trajectory_faceted_plots(
      data = all_trajectory_dist,
      title_base = "Correlation Trajectory: Distance Bins Compared",
      window_size = window_size
    )
    for (plot_name in names(trajectory_faceted_plots)) {
      plots[[paste0("trajectory_faceted_", plot_name)]] <- trajectory_faceted_plots[[plot_name]]
    }
    
    # Event rate trajectory plots (if event data available)
    if (!is.null(all_trajectory_events) && nrow(all_trajectory_events) > 0) {
      cat("    Creating event rate trajectory comparisons...\n")
      event_rate_plots <- create_event_rate_trajectory_plots(
        data = all_trajectory_events,
        title_base = "Event Rate Trajectory Across Animals",
        window_size = window_size
      )
      for (trans_type in names(event_rate_plots)) {
        plots[[paste0("event_rate_trajectory_", trans_type)]] <- event_rate_plots[[trans_type]]
      }
    }
  }
  
  cat("  All plots created\n\n")
  
  # Save plots
  cat("Saving plots to", output_dir, "...\n")
  
  if (run_whole_window) {
    for (trans_type in transition_types) {
      ggsave(file.path(output_dir, paste0("CrossAnimal_", trans_type, "_EventRate_by_Percentile_", window_size, "ep.png")),
             plots[[paste0(trans_type, "_event_rate_percentile")]],
             width = 12, height = 8, dpi = 300)
      
      ggsave(file.path(output_dir, paste0("CrossAnimal_", trans_type, "_Correlation_by_Percentile_", window_size, "ep.png")),
             plots[[paste0(trans_type, "_correlation_percentile")]],
             width = 12, height = 8, dpi = 300)
      
      ggsave(file.path(output_dir, paste0("CrossAnimal_", trans_type, "_Correlation_by_Distance_", window_size, "ep.png")),
             plots[[paste0(trans_type, "_dist_corr_percentile")]],
             width = 12, height = 8, dpi = 300)
    }
  }
  
  if (run_trajectory_analysis) {
    # Create output subdirectory for individual plots
    individual_dir <- file.path(output_dir, "Individual_Plots")
    dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Overall trajectory (combined - both transition types)
    if (!is.null(plots[["trajectory_overall_combined"]])) {
      ggsave(file.path(output_dir, paste0("CrossAnimal_Trajectory_Overall_", window_size, "ep.png")),
             plots[["trajectory_overall_combined"]],
             width = 14, height = 8, dpi = 300)
    }
    
    for (trans_type in transition_types) {
      # Overall trajectory (individual per transition type)
      if (!is.null(plots[[paste0("trajectory_overall_", trans_type)]])) {
        ggsave(file.path(individual_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Overall_", window_size, "ep.png")),
               plots[[paste0("trajectory_overall_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
      
      # Distance-binned trajectory (combined faceted by animal)
      if (!is.null(plots[[paste0("trajectory_by_distance_combined_", trans_type)]])) {
        ggsave(file.path(output_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_by_Distance_", window_size, "ep.png")),
               plots[[paste0("trajectory_by_distance_combined_", trans_type)]],
               width = 14, height = 10, dpi = 300)
      }
      
      # Faceted trajectory (combined faceted by distance bin)
      if (!is.null(plots[[paste0("trajectory_faceted_combined_", trans_type)]])) {
        ggsave(file.path(output_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Faceted_", window_size, "ep.png")),
               plots[[paste0("trajectory_faceted_combined_", trans_type)]],
               width = 12, height = 14, dpi = 300)
      }
      
      # Event rate trajectory (if available)
      if (!is.null(plots[[paste0("event_rate_trajectory_", trans_type)]])) {
        ggsave(file.path(individual_dir, paste0("CrossAnimal_", trans_type, "_EventRate_Trajectory_", window_size, "ep.png")),
               plots[[paste0("event_rate_trajectory_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
    }
    
    # Save individual animal plots (from distance-binned)
    cat("  Saving individual animal plots...\n")
    animals <- unique(all_trajectory_dist$animal)
    for (trans_type in transition_types) {
      for (anim in animals) {
        plot_key <- paste0("trajectory_by_distance_", anim, "_", trans_type)
        if (!is.null(plots[[plot_key]])) {
          ggsave(file.path(individual_dir, paste0(anim, "_", trans_type, "_Trajectory_by_Distance_", window_size, "ep.png")),
                 plots[[plot_key]],
                 width = 10, height = 8, dpi = 300)
        }
      }
    }
    
    # Save individual distance bin plots (from faceted)
    cat("  Saving individual distance bin plots...\n")
    distance_bins <- unique(all_trajectory_dist$distance_bin_label)
    for (trans_type in transition_types) {
      for (dist_bin in distance_bins) {
        dist_bin_clean <- gsub("-", "", dist_bin)
        plot_key <- paste0("trajectory_faceted_", trans_type, "_", dist_bin_clean)
        if (!is.null(plots[[plot_key]])) {
          ggsave(file.path(individual_dir, paste0("CrossAnimal_", trans_type, "_", dist_bin_clean, "_Trajectory_", window_size, "ep.png")),
                 plots[[plot_key]],
                 width = 10, height = 8, dpi = 300)
        }
      }
    }
  }
  
  cat("  Plots saved\n\n")
  
  # Save combined data
  if (save_intermediate) {
    cat("Saving combined data tables...\n")
    
    if (run_whole_window) {
      write.csv(all_event_rates,
                file.path(output_dir, "CrossAnimal_EventRates_All.csv"),
                row.names = FALSE)
      
      write.csv(all_correlations,
                file.path(output_dir, "CrossAnimal_Correlations_All.csv"),
                row.names = FALSE)
      
      write.csv(all_dist_corr,
                file.path(output_dir, "CrossAnimal_DistCorr_All.csv"),
                row.names = FALSE)
    }
    
    if (run_trajectory_analysis) {
      write.csv(all_trajectory_data,
                file.path(output_dir, "CrossAnimal_Trajectory_All.csv"),
                row.names = FALSE)
      
      write.csv(all_trajectory_dist,
                file.path(output_dir, "CrossAnimal_Trajectory_Distance_All.csv"),
                row.names = FALSE)
    }
    
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
    trajectory_data = all_trajectory_data,
    trajectory_dist = all_trajectory_dist,
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

#' Create trajectory comparison plots across animals (with error bars)
#' Returns both combined faceted plot and individual plots per transition type
create_trajectory_comparison_plots <- function(data, title_base, window_size) {
  
  # Compute mean across all ROI pairs per animal/condition/transition/epoch
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  # Combined faceted plot (both transition types)
  plots$combined <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                             color = animal, linetype = condition)) +
    geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  width = 0.1, alpha = 0.5, linewidth = 0.5) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    facet_wrap(~transition_type, ncol = 2, scales = "free_y") +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
    ) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = title_base,
      subtitle = paste0("Mean ± SE pairwise correlation (", window_size, " epoch window)"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation",
      color = "Animal",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  # Individual plots per transition type
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    p <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_correlation,
                               color = animal, linetype = condition)) +
      geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation),
                    width = 0.1, alpha = 0.5, linewidth = 0.5) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
      ) +
      scale_color_viridis_d(option = "turbo") +
      labs(
        title = paste(title_base, "-", trans_type),
        subtitle = paste0("Mean ± SE pairwise correlation (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
    
    plots[[trans_type]] <- p
  }
  
  return(plots)
}

#' Create distance-binned trajectory comparison plots (with error bars)
#' Returns combined faceted plot AND individual plots per animal/transition type
create_trajectory_distance_comparison_plots <- function(data, title_base, window_size) {
  
  # Compute mean across transitions per animal/condition/transition_type/distance_bin/epoch
  plot_data <- data %>%
    group_by(animal, condition, transition_type, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  animals <- unique(plot_data$animal)
  
  # Combined faceted plots (one per transition type, faceted by animal)
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[paste0("combined_", trans_type)]] <- ggplot(trans_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = distance_bin_label, linetype = condition)) +
      geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation),
                    width = 0.1, alpha = 0.4, linewidth = 0.4) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.5) +
      facet_wrap(~animal, ncol = 2, scales = "free_y") +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
      ) +
      scale_color_viridis_d(option = "plasma") +
      labs(
        title = paste(title_base, "-", trans_type),
        subtitle = paste0("Mean ± SE by distance bin (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Distance Bin",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  }
  
  # Individual plots per animal/transition type
  for (trans_type in transition_types) {
    for (anim in animals) {
      subset_data <- plot_data %>% 
        filter(transition_type == trans_type, animal == anim)
      
      if (nrow(subset_data) == 0) next
      
      p <- ggplot(subset_data, aes(x = epoch_in_window, y = mean_correlation,
                                  color = distance_bin_label, linetype = condition)) +
        geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                          ymax = mean_correlation + se_correlation),
                      width = 0.1, alpha = 0.5, linewidth = 0.5) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2.5) +
        scale_linetype_manual(
          values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
          labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
        ) +
        scale_color_viridis_d(option = "plasma") +
        labs(
          title = paste(anim, "-", trans_type),
          subtitle = paste0("Mean ± SE by distance bin (", window_size, " epoch window)"),
          x = "Epoch Position Relative to Transition",
          y = "Mean Pairwise Correlation",
          color = "Distance Bin",
          linetype = "Condition"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold")
        )
      
      plots[[paste0(anim, "_", trans_type)]] <- p
    }
  }
  
  return(plots)
}

#' Create faceted trajectory plots by distance bin (with error bars)
#' Returns combined faceted plot AND individual plots per distance bin/transition type
create_trajectory_faceted_plots <- function(data, title_base, window_size) {
  
  # Compute mean per animal/condition/transition_type/distance_bin/epoch
  plot_data <- data %>%
    group_by(animal, condition, transition_type, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  distance_bins <- unique(plot_data$distance_bin_label)
  
  # Combined faceted plots (one per transition type, faceted by distance bin)
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[paste0("combined_", trans_type)]] <- ggplot(trans_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = animal, linetype = condition)) +
      geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation),
                    width = 0.1, alpha = 0.4, linewidth = 0.4) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.5) +
      facet_wrap(~distance_bin_label, ncol = 2, scales = "free_y") +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
      ) +
      scale_color_viridis_d(option = "turbo") +
      labs(
        title = paste(title_base, "-", trans_type),
        subtitle = paste0("Mean ± SE by distance bin (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  }
  
  # Individual plots per distance bin/transition type
  for (trans_type in transition_types) {
    for (dist_bin in distance_bins) {
      subset_data <- plot_data %>% 
        filter(transition_type == trans_type, distance_bin_label == dist_bin)
      
      if (nrow(subset_data) == 0) next
      
      p <- ggplot(subset_data, aes(x = epoch_in_window, y = mean_correlation,
                                  color = animal, linetype = condition)) +
        geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                          ymax = mean_correlation + se_correlation),
                      width = 0.1, alpha = 0.5, linewidth = 0.5) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2.5) +
        scale_linetype_manual(
          values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
          labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
        ) +
        scale_color_viridis_d(option = "turbo") +
        labs(
          title = paste(trans_type, "-", dist_bin),
          subtitle = paste0("Mean ± SE pairwise correlation (", window_size, " epoch window)"),
          x = "Epoch Position Relative to Transition",
          y = "Mean Pairwise Correlation",
          color = "Animal",
          linetype = "Condition"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold")
        )
      
      # Clean up distance bin name for file naming
      dist_bin_clean <- gsub("-", "", dist_bin)
      plots[[paste0(trans_type, "_", dist_bin_clean)]] <- p
    }
  }
  
  return(plots)
}


#' Create event rate trajectory plots (one per transition type)
create_event_rate_trajectory_plots <- function(data, title_base, window_size) {
  
  # Compute mean event rate per animal/condition/transition_type/epoch
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      total_events = n(),
      n_rois = n_distinct(Cell_Name),
      n_transitions = n_distinct(transition_id),
      mean_event_rate = total_events / (n_rois * n_transitions),
      .groups = "drop"
    ) %>%
    # Compute SE across animals for each condition/transition/epoch
    group_by(condition, transition_type, epoch_in_window) %>%
    mutate(
      grand_mean = mean(mean_event_rate, na.rm = TRUE),
      se_rate = sd(mean_event_rate, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  # Create separate plot for each transition type
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    p <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_event_rate,
                               color = animal, linetype = condition)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
      ) +
      scale_color_viridis_d(option = "turbo") +
      labs(
        title = paste(title_base, "-", trans_type),
        subtitle = paste0("Mean events per ROI per transition (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Event Rate",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
    
    plots[[trans_type]] <- p
  }
  
  return(plots)
}

# ============================================================
# AUTO-RUN (if sourced in interactive mode)
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  BATCH PROCESSOR V2 LOADED\n")
  cat("========================================\n")
  cat("\nTo run batch processing, execute:\n")
  cat("  results <- batch_process_animals()\n\n")
  cat("To customize settings, edit the CONFIGURATION section above,\n")
  cat("or pass parameters to batch_process_animals()\n\n")
}
