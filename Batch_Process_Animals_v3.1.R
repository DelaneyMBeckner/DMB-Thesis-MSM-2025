# ============================================================
# BATCH PROCESSOR: Cross-Animal Comparison (v3)
# ============================================================
# Processes multiple animals and creates comparison plots
# v3 adds: Activity-based ROI filtering with clear labeling
#
# All output files include "_filtered" suffix
# All plot titles include "(Active ROIs Only)"
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Load pipeline script
source("Pipeline_Transition_Analysis_v3.R")

# ============================================================
# CONFIGURATION
# ============================================================

# Animals to analyze
ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm9", "mPFCm4")

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
SAVE_INTERMEDIATE <- TRUE
CCF_USE_PARALLEL <- TRUE
CCF_N_CORES <- 6

# Analysis modes
RUN_TRAJECTORY_ANALYSIS <- TRUE
RUN_WHOLE_WINDOW <- FALSE

# CCF parameters
RUN_LAG_SWEEP <- FALSE
MAX_LAG <- 5

# Distance-correlation parameters
BINNING_METHOD <- "percentile"
N_DISTANCE_BINS <- 10
N_PERCENTILE_BINS <- 3
N_TRAJECTORY_BINS <- 3
FIT_EXPONENTIAL <- FALSE

# Activity-based ROI filtering
FILTER_BY_ACTIVITY <- TRUE
MIN_EVENTS_BASELINE <- 1    # Minimum events required
BASELINE_EPOCHS <- 18       # Per this many epochs (18~=1 sleep cycle)

# Correlation-based pair filtering
CORRELATION_FILTER <- "off"        # "off", "alone", or "alongside"
CORRELATION_FILTER_METHOD <- "percentile"  # "percentile" or "outlier"
CORRELATION_PERCENTILE <- 5       # Top X% (only used if method = "percentile")

# File naming
FILE_SUFFIX <- "_all_simple_new"         # Suffix for output files

# Plotting options
N_PERCENTILE_POINTS <- 100

# ============================================================
# PLOT AESTHETICS
# ============================================================

ANIMAL_COLORS <- c(
  "mPFCf5" = "#60BBA0",  # Female 1
  "mPFCf6" = "#E48F4E",  # Female 2
  "mPFCm4" = "#EC5F61",  # Male 1
  "mPFCm9" = "#9F9BCA"   # Male 2
)
ANIMAL_LABELS <- c(
  "mPFCf5" = "Female 1",
  "mPFCf6" = "Female 2",
  "mPFCm4" = "Male 1",
  "mPFCm9" = "Male 2"
)
ANIMAL_LABELLER <- as_labeller(ANIMAL_LABELS)

CONDITION_COLORS <- c(
  "BL" = "#2E86AB",  # Baseline
  "SD" = "#A23B72",  # Sleep Deprivation
  "WO" = "#E69F00"   # Recovery
)
CONDITION_LINETYPES <- c("BL" = "solid", "SD" = "dashed", "WO" = "dotted")
CONDITION_LABELS    <- c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Recovery")

DISTANCE_BIN_COLORS <- c("Close" = "#3A0557", "Medium" = "#185FA5", "Far" = "#0F6E56")

# ============================================================
# MAIN BATCH PROCESSING FUNCTION
# ============================================================

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
                                  n_percentile_points = N_PERCENTILE_POINTS,
                                  filter_by_activity = FILTER_BY_ACTIVITY,
                                  min_events_baseline = MIN_EVENTS_BASELINE,
                                  baseline_epochs = BASELINE_EPOCHS,
                                  correlation_filter = CORRELATION_FILTER,
                                  correlation_filter_method = CORRELATION_FILTER_METHOD,
                                  correlation_percentile = CORRELATION_PERCENTILE,
                                  file_suffix = FILE_SUFFIX) {
  
  cat("====================================================\n")
  cat("  BATCH PROCESSING v3: CROSS-ANIMAL COMPARISON\n")
  cat("  (Activity-Filtered Analysis)\n")
  cat("====================================================\n")
  cat("Animals:", paste(animals, collapse = ", "), "\n")
  cat("Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("Transitions:", paste(transition_types, collapse = ", "), "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("Trajectory analysis:", run_trajectory_analysis, "\n")
  cat("Whole window analysis:", run_whole_window, "\n")
  if (filter_by_activity) {
    cat("\n*** ACTIVITY FILTER ENABLED ***\n")
    cat("Threshold:", min_events_baseline, "event(s) per", baseline_epochs, "epochs at BL\n")
  }
  cat("Output directory:", output_dir, "\n")
  cat("====================================================\n\n")
  
  # Storage for all animal results
  all_animal_results <- list()
  
  # Process each animal
  for (animal in animals) {
    cat("\n\n")
    cat("####################################################\n")
    cat("# ANIMAL", match(animal, animals), "of", length(animals), "\n")
    cat("####################################################\n")
    
    tryCatch({
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
        fit_exponential = fit_exponential,
        filter_by_activity = filter_by_activity,
        min_events_baseline = min_events_baseline,
        baseline_epochs = baseline_epochs,
        correlation_filter = correlation_filter,
        correlation_filter_method = correlation_filter_method,
        correlation_percentile = correlation_percentile,
        file_suffix = file_suffix
      )
      
      all_animal_results[[animal]] <- animal_results
      
    }, error = function(e) {
      cat("\n!!! ERROR processing", animal, ":", conditionMessage(e), "\n")
      cat("Continuing with remaining animals...\n\n")
    })
  }
  
  cat("\n====================================================\n")
  cat("  COMBINING RESULTS ACROSS ANIMALS\n")
  cat("====================================================\n\n")
  
  # Print filtering summary
  if (filter_by_activity) {
    cat("FILTERING SUMMARY:\n")
    cat("------------------\n")
    total_original <- 0
    total_active <- 0
    for (animal in names(all_animal_results)) {
      stats <- all_animal_results[[animal]]$filter_stats
      if (!is.null(stats)) {
        cat(sprintf("  %s: %d/%d ROIs retained (%.1f%%)\n", 
                    animal, stats$active_rois, stats$total_rois, stats$pct_retained))
        total_original <- total_original + stats$total_rois
        total_active <- total_active + stats$active_rois
      }
    }
    cat(sprintf("  TOTAL: %d/%d ROIs retained (%.1f%%)\n\n", 
                total_active, total_original, 100 * total_active / total_original))
  }
  
  # Combine whole-window results
  all_event_rates <- NULL
  all_correlations <- NULL
  all_dist_corr <- NULL
  
  if (run_whole_window) {
    all_event_rates <- bind_rows(lapply(all_animal_results, function(x) x$event_rates))
    all_correlations <- bind_rows(lapply(all_animal_results, function(x) x$correlations))
    all_dist_corr <- bind_rows(lapply(all_animal_results, function(x) x$dist_corr))
    
    cat("Combined whole-window data:\n")
    cat("  Event rates:", nrow(all_event_rates), "observations\n")
    cat("  Correlations:", nrow(all_correlations), "observations\n")
  }
  
  # Combine trajectory results
  all_trajectory_data <- NULL
  all_trajectory_dist <- NULL
  all_trajectory_events <- NULL
  
  if (run_trajectory_analysis) {
    all_trajectory_data <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_data))
    all_trajectory_dist <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_dist))
    all_trajectory_events <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_events))
    
    cat("Combined trajectory data:\n")
    cat("  Epoch observations:", nrow(all_trajectory_data), "\n")
    cat("  Distance-binned observations:", nrow(all_trajectory_dist), "\n")
    if (!is.null(all_trajectory_events)) {
      cat("  Event observations:", nrow(all_trajectory_events), "\n")
    }
  }
  
  cat("\n")
  
  # Create plots
  cat("Creating comparison plots...\n\n")
  
  plots <- list()
  filter_label <- ifelse(filter_by_activity, " (Active ROIs Only)", "")
  
  # ====== TRAJECTORY PLOTS ======
  if (run_trajectory_analysis && !is.null(all_trajectory_data)) {
    cat("  [TRAJECTORY PLOTS]\n")
    
    # Overall trajectory
    cat("    Creating overall trajectory comparisons...\n")
    trajectory_overall_plots <- create_trajectory_comparison_plots(
      data = all_trajectory_data,
      title_base = paste0("Correlation Trajectory Across Animals", filter_label),
      window_size = window_size
    )
    for (plot_name in names(trajectory_overall_plots)) {
      plots[[paste0("trajectory_overall_", plot_name)]] <- trajectory_overall_plots[[plot_name]]
    }
    
    # Animal-faceted trajectory (2x2 layout for thesis Figure 8)
    cat("    Creating animal-faceted trajectory plots (2x2 layout)...\n")
    trajectory_by_animal_plots <- create_trajectory_by_animal_plots(
      data = all_trajectory_data,
      title_base = paste0("Correlation Trajectory", filter_label),
      window_size = window_size
    )
    for (plot_name in names(trajectory_by_animal_plots)) {
      plots[[paste0("trajectory_by_animal_", plot_name)]] <- trajectory_by_animal_plots[[plot_name]]
    }
    
    # Distance-binned trajectory
    cat("    Creating distance-binned trajectory comparisons...\n")
    trajectory_distance_plots <- create_trajectory_distance_comparison_plots(
      data = all_trajectory_dist,
      title_base = paste0("Correlation Trajectory by Distance Bin", filter_label),
      window_size = window_size
    )
    for (plot_name in names(trajectory_distance_plots)) {
      plots[[paste0("trajectory_by_distance_", plot_name)]] <- trajectory_distance_plots[[plot_name]]
    }
    
    # Faceted trajectory
    cat("    Creating faceted trajectory comparisons...\n")
    trajectory_faceted_plots <- create_trajectory_faceted_plots(
      data = all_trajectory_dist,
      title_base = paste0("Correlation Trajectory: Distance Bins Compared", filter_label),
      window_size = window_size
    )
    for (plot_name in names(trajectory_faceted_plots)) {
      plots[[paste0("trajectory_faceted_", plot_name)]] <- trajectory_faceted_plots[[plot_name]]
    }
    
    # Event rate trajectory
    if (!is.null(all_trajectory_events) && nrow(all_trajectory_events) > 0) {
      cat("    Creating event rate trajectory comparisons...\n")
      event_rate_plots <- create_event_rate_trajectory_plots(
        data = all_trajectory_events,
        title_base = paste0("Event Rate Trajectory Across Animals", filter_label),
        window_size = window_size
      )
      for (trans_type in names(event_rate_plots)) {
        plots[[paste0("event_rate_trajectory_", trans_type)]] <- event_rate_plots[[trans_type]]
      }
      
      event_rate_by_animal_plots <- create_event_rate_by_animal_plots(
        data = all_trajectory_events,
        title_base = paste0("Event Rate Trajectory Across Animals", filter_label),
        window_size = window_size
      )
      for (trans_type in names(event_rate_by_animal_plots)) {
        plots[[paste0("event_rate_by_animal_", trans_type)]] <- event_rate_by_animal_plots[[trans_type]]
      }
    }
  }
  
  cat("  All plots created\n\n")
  
  # Save plots
  cat("Saving plots to", output_dir, "...\n")
  
  # Ensure directory structure exists
  create_output_directories(output_dir)
  
  # Get between_subject paths
  between_subj_dir <- get_output_path(output_dir, "between_subject", NULL)
  
  if (run_trajectory_analysis) {
    
    # Combined overall plot (both transition types -> between_subject root)
    if (!is.null(plots[["trajectory_overall_combined"]])) {
      ggsave(file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Overall_", window_size, "ep", file_suffix, ".png")),
             plots[["trajectory_overall_combined"]],
             width = 14, height = 8, dpi = 300)
    }
    
    for (trans_type in transition_types) {
      trans_dir <- get_output_path(output_dir, "between_subject", trans_type)
      
      # Individual overall plots
      if (!is.null(plots[[paste0("trajectory_overall_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Overall_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_overall_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
      
      # Animal-faceted plots (2x2 layout - thesis Figure 8)
      if (!is.null(plots[[paste0("trajectory_by_animal_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_ByAnimal_2x2_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_by_animal_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
      
      # Distance-binned combined
      if (!is.null(plots[[paste0("trajectory_by_distance_combined_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_by_Distance_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_by_distance_combined_", trans_type)]],
               width = 14, height = 10, dpi = 300)
      }
      
      # Faceted combined
      if (!is.null(plots[[paste0("trajectory_faceted_combined_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Faceted_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_faceted_combined_", trans_type)]],
               width = 12, height = 14, dpi = 300)
      }
      
      # Event rate
      if (!is.null(plots[[paste0("event_rate_trajectory_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_EventRate_Trajectory_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("event_rate_trajectory_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
      
      # Event rate faceted by animal (2x2 layout)
      if (!is.null(plots[[paste0("event_rate_by_animal_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_EventRate_ByAnimal_2x2_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("event_rate_by_animal_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
    }
    
    # Individual distance bin plots (cross-animal comparisons)
    if (!is.null(all_trajectory_dist)) {
      cat("  Saving individual distance bin plots...\n")
      distance_bins <- unique(all_trajectory_dist$distance_bin_label)
      for (trans_type in transition_types) {
        trans_dir <- get_output_path(output_dir, "between_subject", trans_type)
        for (dist_bin in distance_bins) {
          dist_bin_clean <- gsub("-", "", dist_bin)
          plot_key <- paste0("trajectory_faceted_", trans_type, "_", dist_bin_clean)
          if (!is.null(plots[[plot_key]])) {
            ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_", dist_bin_clean, "_Trajectory_", window_size, "ep", file_suffix, ".png")),
                   plots[[plot_key]],
                   width = 10, height = 8, dpi = 300)
          }
        }
      }
    }
  }
  
  cat("  Plots saved\n\n")
  
  # Save combined data
  if (save_intermediate) {
    cat("Saving combined data tables...\n")
    
    if (run_whole_window && !is.null(all_event_rates)) {
      write.csv(all_event_rates,
                file.path(between_subj_dir, paste0("CrossAnimal_EventRates_All", file_suffix, ".csv")),
                row.names = FALSE)
      write.csv(all_correlations,
                file.path(between_subj_dir, paste0("CrossAnimal_Correlations_All", file_suffix, ".csv")),
                row.names = FALSE)
      write.csv(all_dist_corr,
                file.path(between_subj_dir, paste0("CrossAnimal_DistCorr_All", file_suffix, ".csv")),
                row.names = FALSE)
    }
    
    if (run_trajectory_analysis && !is.null(all_trajectory_data)) {
      write.csv(all_trajectory_data,
                file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_All", file_suffix, ".csv")),
                row.names = FALSE)
      write.csv(all_trajectory_dist,
                file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Distance_All", file_suffix, ".csv")),
                row.names = FALSE)
    }
    
    # Save filtering summary
    if (filter_by_activity) {
      filter_summary <- do.call(rbind, lapply(names(all_animal_results), function(animal) {
        stats <- all_animal_results[[animal]]$filter_stats
        if (!is.null(stats)) {
          data.frame(
            animal = animal,
            total_rois = stats$total_rois,
            active_rois = stats$active_rois,
            excluded_rois = stats$excluded_rois,
            pct_retained = stats$pct_retained,
            threshold = stats$threshold,
            min_events = stats$min_events,
            baseline_epochs = stats$baseline_epochs
          )
        }
      }))
      write.csv(filter_summary,
                file.path(between_subj_dir, paste0("CrossAnimal_Filtering_Summary", file_suffix, ".csv")),
                row.names = FALSE)
    }
    
    cat("  Data tables saved\n\n")
  }
  
  cat("====================================================\n")
  cat("  BATCH PROCESSING COMPLETE\n")
  cat("====================================================\n")
  cat("Animals processed:", length(all_animal_results), "\n")
  if (filter_by_activity) {
    cat("ROI filtering: ENABLED\n")
    cat("  Threshold:", min_events_baseline, "event(s) per", baseline_epochs, "epochs\n")
  }
  cat("Output location:", output_dir, "\n")
  cat("====================================================\n\n")
  
  return(list(
    animal_results = all_animal_results,
    event_rates = all_event_rates,
    correlations = all_correlations,
    dist_corr = all_dist_corr,
    trajectory_data = all_trajectory_data,
    trajectory_dist = all_trajectory_dist,
    trajectory_events = all_trajectory_events,
    plots = plots,
    config = list(
      animals = animals,
      conditions = conditions,
      transition_types = transition_types,
      window_size = window_size,
      filter_by_activity = filter_by_activity,
      min_events_baseline = min_events_baseline,
      baseline_epochs = baseline_epochs
    )
  ))
}


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

#' Create trajectory comparison plots (combined + individual per transition type)
create_trajectory_comparison_plots <- function(data, title_base, window_size) {
  
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  # Combined faceted plot
  plots$combined <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                             color = animal, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                    ymax = mean_correlation + se_correlation,
                    fill = animal), alpha = 0.4, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    facet_wrap(~transition_type, ncol = 2, scales = "free_y") +
    scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
    scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS) +
    scale_fill_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, guide = "none") +
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
    
    plots[[trans_type]] <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_correlation,
                               color = animal, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation,
                      fill = animal), alpha = 0.4, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
      scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS) +
      scale_fill_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, guide = "none") +
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
  }
  
  return(plots)
}


#' Create trajectory plots faceted by animal (2x2 layout with SEM ribbons)
#' For thesis Figure 8: each animal in its own panel, conditions overlaid
create_trajectory_by_animal_plots <- function(data, title_base, window_size) {
  
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  # Condition colors and linetypes
  condition_colors    <- CONDITION_COLORS
  condition_linetypes <- CONDITION_LINETYPES
  condition_labels    <- CONDITION_LABELS
  
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[trans_type]] <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_correlation,
                               color = condition, fill = condition, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.8) +
      facet_wrap(~animal, nrow = 2, ncol = 2, scales = "free_y",
                 labeller = ANIMAL_LABELLER) +
      scale_color_manual(values = condition_colors, labels = condition_labels) +
      scale_fill_manual(values = condition_colors, labels = condition_labels) +
      scale_linetype_manual(values = condition_linetypes, labels = condition_labels) +
      labs(
        title = paste(title_base, "-", trans_type),
        subtitle = paste0("Mean ± SEM pairwise correlation (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Condition",
        fill = "Condition",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.spacing = unit(1, "lines")
      )
  }
  
  return(plots)
}


#' Create distance-binned trajectory comparison plots
create_trajectory_distance_comparison_plots <- function(data, title_base, window_size) {
  
  # Compute bin ranges for subtitle
  #bin_ranges <- data %>%
    #group_by(distance_bin_label) %>%
    #summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    #arrange(bin_min) %>%
    #mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    #pull(range_str)
  #bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
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
  
  # Combined plots (one per transition type, faceted by animal)
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[paste0("combined_", trans_type)]] <- ggplot(trans_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = distance_bin_label, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation,
                      fill = distance_bin_label), alpha = 0.4, color = NA) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.5) +
      facet_wrap(~animal, ncol = 2, scales = "free_y",
                 labeller = ANIMAL_LABELLER) +
      scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
      scale_color_manual(values = DISTANCE_BIN_COLORS) +
      scale_fill_manual(values = DISTANCE_BIN_COLORS, guide = "none") +
      labs(
        title = paste(title_base, "-", trans_type),
        #subtitle = bin_ranges_subtitle,
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Distance Bin",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        #plot.subtitle = element_text(size = 11),
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
      
      plots[[paste0(anim, "_", trans_type)]] <- ggplot(subset_data, 
          aes(x = epoch_in_window, y = mean_correlation,
              color = distance_bin_label, linetype = condition)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
                 color = "gray40", fontface = "bold", size = 3) +
        geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation,
                        fill = distance_bin_label), alpha = 0.4, color = NA) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2.5) +
        scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
        scale_color_manual(values = DISTANCE_BIN_COLORS) +
        scale_fill_manual(values = DISTANCE_BIN_COLORS, guide = "none") +
        labs(
          title = paste(ANIMAL_LABELS[anim], "-", trans_type),
          #subtitle = bin_ranges_subtitle,
          x = "Epoch Position Relative to Transition",
          y = "Mean Pairwise Correlation",
          color = "Distance Bin",
          linetype = "Condition"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          #plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold")
        )
    }
  }
  
  return(plots)
}


#' Create faceted trajectory plots by distance bin
create_trajectory_faceted_plots <- function(data, title_base, window_size) {
  
  # Compute bin ranges for subtitle
  #bin_ranges <- data %>%
    #group_by(distance_bin_label) %>%
    #summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    #arrange(bin_min) %>%
    #mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    #pull(range_str)
  #bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
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
  
  # Combined plots (one per transition type, faceted by distance bin)
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[paste0("combined_", trans_type)]] <- ggplot(trans_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = animal, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation,
                      fill = animal), alpha = 0.4, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.5) +
      facet_wrap(~distance_bin_label, ncol = 2, scales = "free_y") +
      scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
      scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS) +
      scale_fill_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, guide = "none") +
      labs(
        title = paste(title_base, "-", trans_type),
        #subtitle = bin_ranges_subtitle,
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        #plot.subtitle = element_text(size = 11),
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
      
      dist_bin_clean <- gsub("-", "", dist_bin)
      
      plots[[paste0(trans_type, "_", dist_bin_clean)]] <- ggplot(subset_data, 
          aes(x = epoch_in_window, y = mean_correlation,
              color = animal, linetype = condition)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
                 color = "gray40", fontface = "bold", size = 3) +
        geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation,
                        fill = animal), alpha = 0.4, color = NA) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2.5) +
        scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
        scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS) +
        scale_fill_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, guide = "none") +
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
    }
  }
  
  return(plots)
}


#' Create event rate trajectory plots
create_event_rate_trajectory_plots <- function(data, title_base, window_size) {
  
  # Compute per-transition event rate first, then mean +/- SE across transitions
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window, transition_id) %>%
    summarise(
      transition_event_rate = n() / n_distinct(Cell_Name),
      .groups = "drop"
    ) %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(transition_event_rate, na.rm = TRUE),
      se_event_rate   = sd(transition_event_rate, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[trans_type]] <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_event_rate,
                               color = animal, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_event_rate - se_event_rate,
                      ymax = mean_event_rate + se_event_rate,
                      fill = animal), alpha = 0.4, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
      scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS) +
      scale_fill_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, guide = "none") +
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
  }
  
  return(plots)
}


#' Create event rate trajectory plots faceted by animal (2x2 layout)
create_event_rate_by_animal_plots <- function(data, title_base, window_size) {
  
  plot_data <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window, transition_id) %>%
    summarise(
      transition_event_rate = n() / n_distinct(Cell_Name),
      .groups = "drop"
    ) %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(transition_event_rate, na.rm = TRUE),
      se_event_rate   = sd(transition_event_rate, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  transition_types <- unique(plot_data$transition_type)
  
  for (trans_type in transition_types) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[trans_type]] <- ggplot(trans_data, aes(x = epoch_in_window, y = mean_event_rate,
                               color = condition, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5,
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_event_rate - se_event_rate,
                      ymax = mean_event_rate + se_event_rate,
                      fill = condition), alpha = 0.4, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      facet_wrap(~animal, nrow = 2, ncol = 2, scales = "free_y",
                 labeller = ANIMAL_LABELLER) +
      scale_color_manual(values = CONDITION_COLORS, labels = CONDITION_LABELS) +
      scale_fill_manual(values = CONDITION_COLORS, labels = CONDITION_LABELS, guide = "none") +
      scale_linetype_manual(values = CONDITION_LINETYPES, labels = CONDITION_LABELS) +
      labs(
        title = paste(title_base, "-", trans_type, "- By Animal"),
        subtitle = paste0("Mean ± SE events per ROI per transition (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Event Rate",
        color = "Condition",
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
  }
  
  return(plots)
}


# ============================================================
# AUTO-RUN
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  BATCH PROCESSOR v3 LOADED\n")
  cat("  (Activity-Filtered Analysis)\n")
  cat("========================================\n")
  cat("\nTo run batch processing, execute:\n")
  cat("  results <- batch_process_animals()\n\n")
  cat("Default filter: >=1 event per 18 epochs at BL\n\n")
}
