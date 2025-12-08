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
ANIMALS <- c("mPFCm4")

# Conditions (experimental days)
CONDITIONS <- c("BL", "SD", "WO")

# Transition types to include
TRANSITION_TYPES <- c("Wake2NREM", "NREM2Wake")

# File naming
FILE_SUFFIX <- "_virusHighCorr"

# Window size (3 or 9 epochs)
WINDOW_SIZE <- 9

# Data directories
DATA_DIR <- "E:/Data_Processing/R/Data CSVs"
OUTPUT_DIR <- "E:/Data_Processing/R/Results/Trajectory_Analysis"

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
CORRELATION_FILTER <- "alongside"        # "off", "alone", or "alongside"
CORRELATION_FILTER_METHOD <- "percentile"  # "percentile" or "outlier"
CORRELATION_PERCENTILE <- 5       # Top X% (only used if method = "percentile")

      # Suffix for output files

# Plotting options
N_PERCENTILE_POINTS <- 100

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
  if (correlation_filter != "off") {
    cat("\n*** CORRELATION FILTER:", toupper(correlation_filter), "***\n")
    if (correlation_filter_method == "percentile") {
      cat("Method: Top", correlation_percentile, "% of correlated pairs\n")
    } else {
      cat("Method: Outlier detection\n")
    }
    if (correlation_filter == "alone") {
      cat("Output: Only high-correlation pairs (tagged in filenames)\n")
    } else {
      cat("Output: Both all-pairs AND high-correlation pairs\n")
    }
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
  all_trajectory_data_highcorr <- NULL
  all_trajectory_dist_highcorr <- NULL
  
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
    
    # Combine high-correlation filtered data (if using "alongside" mode)
    if (correlation_filter == "alongside") {
      all_trajectory_data_highcorr <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_data_highcorr))
      all_trajectory_dist_highcorr <- bind_rows(lapply(all_animal_results, function(x) x$trajectory_dist_highcorr))
      
      if (!is.null(all_trajectory_data_highcorr) && nrow(all_trajectory_data_highcorr) > 0) {
        cat("  High-corr epoch observations:", nrow(all_trajectory_data_highcorr), "\n")
        cat("  High-corr distance-binned observations:", nrow(all_trajectory_dist_highcorr), "\n")
      }
    }
  }
  
  cat("\n")
  
  # Create plots
  cat("Creating comparison plots...\n\n")
  
  plots <- list()
  
  # Build title suffix and file tag for high-correlation filtering
  # These apply to main plots when "alone", and to separate highcorr plots when "alongside"
  if (correlation_filter %in% c("alone", "alongside")) {
    if (correlation_filter_method == "percentile") {
      highcorr_title_suffix <- paste0(" (Top ", correlation_percentile, "% Corr Pairs)")
      highcorr_file_tag <- paste0("_top", correlation_percentile, "pct")
    } else {
      highcorr_title_suffix <- " (High-Corr Outlier Pairs)"
      highcorr_file_tag <- "_outliers"
    }
  } else {
    highcorr_title_suffix <- ""
    highcorr_file_tag <- ""
  }
  
  # When "alone", main data IS highcorr data, so add tags to main plots
  # When "alongside", main data is all pairs, highcorr plots are separate
  main_title_suffix <- ifelse(correlation_filter == "alone", highcorr_title_suffix, "")
  main_file_tag <- ifelse(correlation_filter == "alone", highcorr_file_tag, "")
  
  # ====== TRAJECTORY PLOTS ======
  if (run_trajectory_analysis && !is.null(all_trajectory_data)) {
    cat("  [TRAJECTORY PLOTS]\n")
    if (correlation_filter == "alone") {
      cat("    (Using high-correlation filtered pairs)\n")
    }
    
    # Overall trajectory
    cat("    Creating overall trajectory comparisons...\n")
    trajectory_overall_plots <- create_trajectory_comparison_plots(
      data = all_trajectory_data,
      title_base = paste0("Correlation Trajectory Across Animals", main_title_suffix),
      window_size = window_size
    )
    for (plot_name in names(trajectory_overall_plots)) {
      plots[[paste0("trajectory_overall_", plot_name)]] <- trajectory_overall_plots[[plot_name]]
    }
    
    # Distance-binned trajectory
    cat("    Creating distance-binned trajectory comparisons...\n")
    trajectory_distance_plots <- create_trajectory_distance_comparison_plots(
      data = all_trajectory_dist,
      title_base = paste0("Correlation Trajectory by Distance Bin", main_title_suffix),
      window_size = window_size
    )
    for (plot_name in names(trajectory_distance_plots)) {
      plots[[paste0("trajectory_by_distance_", plot_name)]] <- trajectory_distance_plots[[plot_name]]
    }
    
    # Faceted trajectory
    cat("    Creating faceted trajectory comparisons...\n")
    trajectory_faceted_plots <- create_trajectory_faceted_plots(
      data = all_trajectory_dist,
      title_base = paste0("Correlation Trajectory: Distance Bins Compared", main_title_suffix),
      window_size = window_size
    )
    for (plot_name in names(trajectory_faceted_plots)) {
      plots[[paste0("trajectory_faceted_", plot_name)]] <- trajectory_faceted_plots[[plot_name]]
    }
    
    # Event rate trajectory (not affected by correlation filtering)
    if (!is.null(all_trajectory_events) && nrow(all_trajectory_events) > 0) {
      cat("    Creating event rate trajectory comparisons...\n")
      event_rate_plot <- create_event_rate_trajectory_plots(
        data = all_trajectory_events,
        title_base = "Event Rate Trajectory Across Animals",
        window_size = window_size
      )
      plots[["event_rate_trajectory"]] <- event_rate_plot
    }
    
    # ====== ADDITIONAL HIGH-CORRELATION FILTERED PLOTS (alongside mode only) ======
    if (correlation_filter == "alongside" && 
        !is.null(all_trajectory_data_highcorr) && nrow(all_trajectory_data_highcorr) > 0) {
      
      cat("  [HIGH-CORRELATION FILTERED PLOTS (alongside mode)]\n")
      
      # Build highcorr label for titles
      if (correlation_filter_method == "percentile") {
        highcorr_title_suffix <- paste0(" (Top ", correlation_percentile, "% Corr Pairs)")
      } else {
        highcorr_title_suffix <- " (High-Corr Outlier Pairs)"
      }
      
      # Overall trajectory - highcorr
      cat("    Creating highcorr overall trajectory comparisons...\n")
      trajectory_overall_plots_hc <- create_trajectory_comparison_plots(
        data = all_trajectory_data_highcorr,
        title_base = paste0("Correlation Trajectory Across Animals", highcorr_title_suffix),
        window_size = window_size
      )
      for (plot_name in names(trajectory_overall_plots_hc)) {
        plots[[paste0("trajectory_overall_highcorr_", plot_name)]] <- trajectory_overall_plots_hc[[plot_name]]
      }
      
      # Distance-binned trajectory - highcorr
      if (!is.null(all_trajectory_dist_highcorr) && nrow(all_trajectory_dist_highcorr) > 0) {
        cat("    Creating highcorr distance-binned trajectory comparisons...\n")
        trajectory_distance_plots_hc <- create_trajectory_distance_comparison_plots(
          data = all_trajectory_dist_highcorr,
          title_base = paste0("Correlation Trajectory by Distance Bin", highcorr_title_suffix),
          window_size = window_size
        )
        for (plot_name in names(trajectory_distance_plots_hc)) {
          plots[[paste0("trajectory_by_distance_highcorr_", plot_name)]] <- trajectory_distance_plots_hc[[plot_name]]
        }
        
        # Faceted trajectory - highcorr
        cat("    Creating highcorr faceted trajectory comparisons...\n")
        trajectory_faceted_plots_hc <- create_trajectory_faceted_plots(
          data = all_trajectory_dist_highcorr,
          title_base = paste0("Correlation Trajectory: Distance Bins Compared", highcorr_title_suffix),
          window_size = window_size
        )
        for (plot_name in names(trajectory_faceted_plots_hc)) {
          plots[[paste0("trajectory_faceted_highcorr_", plot_name)]] <- trajectory_faceted_plots_hc[[plot_name]]
        }
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
      ggsave(file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Overall", main_file_tag, "_", window_size, "ep", file_suffix, ".png")),
             plots[["trajectory_overall_combined"]],
             width = 14, height = 8, dpi = 300)
    }
    
    for (trans_type in transition_types) {
      trans_dir <- get_output_path(output_dir, "between_subject", trans_type)
      
      # Individual overall plots
      if (!is.null(plots[[paste0("trajectory_overall_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Overall", main_file_tag, "_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_overall_", trans_type)]],
               width = 10, height = 8, dpi = 300)
      }
      
      # Distance-binned combined
      if (!is.null(plots[[paste0("trajectory_by_distance_combined_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_by_Distance", main_file_tag, "_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_by_distance_combined_", trans_type)]],
               width = 14, height = 10, dpi = 300)
      }
      
      # Faceted combined
      if (!is.null(plots[[paste0("trajectory_faceted_combined_", trans_type)]])) {
        ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Faceted", main_file_tag, "_", window_size, "ep", file_suffix, ".png")),
               plots[[paste0("trajectory_faceted_combined_", trans_type)]],
               width = 12, height = 14, dpi = 300)
      }
    }
    
    # Event rate trajectory (single faceted plot saved to between_subject root)
    if (!is.null(plots[["event_rate_trajectory"]])) {
      ggsave(file.path(between_subj_dir, paste0("CrossAnimal_EventRate_Trajectory_", window_size, "ep", file_suffix, ".png")),
             plots[["event_rate_trajectory"]],
             width = 12, height = 6, dpi = 300)
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
            ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_", dist_bin_clean, "_Trajectory", main_file_tag, "_", window_size, "ep", file_suffix, ".png")),
                   plots[[plot_key]],
                   width = 10, height = 8, dpi = 300)
          }
        }
      }
    }
    
    # ====== SAVE HIGH-CORRELATION FILTERED PLOTS ======
    if (correlation_filter == "alongside" && 
        !is.null(all_trajectory_data_highcorr) && nrow(all_trajectory_data_highcorr) > 0) {
      
      cat("  Saving high-correlation filtered plots...\n")
      
      # Build file suffix for highcorr
      if (correlation_filter_method == "percentile") {
        highcorr_file_tag <- paste0("_top", correlation_percentile, "pct")
      } else {
        highcorr_file_tag <- "_outliers"
      }
      
      # Combined overall plot - highcorr
      if (!is.null(plots[["trajectory_overall_highcorr_combined"]])) {
        ggsave(file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Overall_highcorr", highcorr_file_tag, "_", window_size, "ep", file_suffix, ".png")),
               plots[["trajectory_overall_highcorr_combined"]],
               width = 14, height = 8, dpi = 300)
      }
      
      for (trans_type in transition_types) {
        trans_dir <- get_output_path(output_dir, "between_subject", trans_type)
        
        # Individual overall plots - highcorr
        if (!is.null(plots[[paste0("trajectory_overall_highcorr_", trans_type)]])) {
          ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Overall_highcorr", highcorr_file_tag, "_", window_size, "ep", file_suffix, ".png")),
                 plots[[paste0("trajectory_overall_highcorr_", trans_type)]],
                 width = 10, height = 8, dpi = 300)
        }
        
        # Distance-binned combined - highcorr
        if (!is.null(plots[[paste0("trajectory_by_distance_highcorr_combined_", trans_type)]])) {
          ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_by_Distance_highcorr", highcorr_file_tag, "_", window_size, "ep", file_suffix, ".png")),
                 plots[[paste0("trajectory_by_distance_highcorr_combined_", trans_type)]],
                 width = 14, height = 10, dpi = 300)
        }
        
        # Faceted combined - highcorr
        if (!is.null(plots[[paste0("trajectory_faceted_highcorr_combined_", trans_type)]])) {
          ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_Trajectory_Faceted_highcorr", highcorr_file_tag, "_", window_size, "ep", file_suffix, ".png")),
                 plots[[paste0("trajectory_faceted_highcorr_combined_", trans_type)]],
                 width = 12, height = 14, dpi = 300)
        }
      }
      
      # Individual distance bin plots - highcorr
      if (!is.null(all_trajectory_dist_highcorr) && nrow(all_trajectory_dist_highcorr) > 0) {
        distance_bins_hc <- unique(all_trajectory_dist_highcorr$distance_bin_label)
        for (trans_type in transition_types) {
          trans_dir <- get_output_path(output_dir, "between_subject", trans_type)
          for (dist_bin in distance_bins_hc) {
            dist_bin_clean <- gsub("-", "", dist_bin)
            plot_key <- paste0("trajectory_faceted_highcorr_", trans_type, "_", dist_bin_clean)
            if (!is.null(plots[[plot_key]])) {
              ggsave(file.path(trans_dir, paste0("CrossAnimal_", trans_type, "_", dist_bin_clean, "_Trajectory_highcorr", highcorr_file_tag, "_", window_size, "ep", file_suffix, ".png")),
                     plots[[plot_key]],
                     width = 10, height = 8, dpi = 300)
            }
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
                file.path(between_subj_dir, paste0("CrossAnimal_Trajectory", main_file_tag, "_All", file_suffix, ".csv")),
                row.names = FALSE)
      write.csv(all_trajectory_dist,
                file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Distance", main_file_tag, "_All", file_suffix, ".csv")),
                row.names = FALSE)
      
      # Save high-correlation filtered data (alongside mode only - separate from main data)
      if (correlation_filter == "alongside" && 
          !is.null(all_trajectory_data_highcorr) && nrow(all_trajectory_data_highcorr) > 0) {
        
        write.csv(all_trajectory_data_highcorr,
                  file.path(between_subj_dir, paste0("CrossAnimal_Trajectory", highcorr_file_tag, "_All", file_suffix, ".csv")),
                  row.names = FALSE)
        write.csv(all_trajectory_dist_highcorr,
                  file.path(between_subj_dir, paste0("CrossAnimal_Trajectory_Distance", highcorr_file_tag, "_All", file_suffix, ".csv")),
                  row.names = FALSE)
      }
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
  if (correlation_filter != "off") {
    cat("Correlation filtering:", toupper(correlation_filter), "\n")
    if (correlation_filter == "alone") {
      cat("  Main output uses high-correlation pairs only\n")
    } else if (correlation_filter == "alongside") {
      cat("  Generated both all-pairs AND high-correlation outputs\n")
    }
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
    trajectory_data_highcorr = all_trajectory_data_highcorr,
    trajectory_dist_highcorr = all_trajectory_dist_highcorr,
    plots = plots,
    config = list(
      animals = animals,
      conditions = conditions,
      transition_types = transition_types,
      window_size = window_size,
      filter_by_activity = filter_by_activity,
      min_events_baseline = min_events_baseline,
      baseline_epochs = baseline_epochs,
      correlation_filter = correlation_filter,
      correlation_filter_method = correlation_filter_method,
      correlation_percentile = correlation_percentile
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
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation",
      color = "Animal",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
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
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  }
  
  return(plots)
}


#' Create distance-binned trajectory comparison plots
create_trajectory_distance_comparison_plots <- function(data, title_base, window_size) {
  
  # Compute bin ranges for subtitle
  bin_ranges <- data %>%
    group_by(distance_bin_label) %>%
    summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    arrange(bin_min) %>%
    mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    pull(range_str)
  bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
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
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Distance Bin",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
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
          x = "Epoch Position Relative to Transition",
          y = "Mean Pairwise Correlation",
          color = "Distance Bin",
          linetype = "Condition"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
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
  bin_ranges <- data %>%
    group_by(distance_bin_label) %>%
    summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    arrange(bin_min) %>%
    mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    pull(range_str)
  bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
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
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Animal",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
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
          x = "Epoch Position Relative to Transition",
          y = "Mean Pairwise Correlation",
          color = "Animal",
          linetype = "Condition"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
          axis.title = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold")
        )
    }
  }
  
  return(plots)
}


#' Create event rate trajectory plots
#' Returns a single faceted plot (faceted by transition type)
create_event_rate_trajectory_plots <- function(data, title_base, window_size) {
  
  # First compute event rate per ROI per animal/condition/transition/epoch
  roi_rates <- data %>%
    group_by(animal, condition, transition_type, epoch_in_window, Cell_Name) %>%
    summarise(
      n_events = n(),
      n_transitions = n_distinct(transition_id),
      events_per_epoch = n_events / n_transitions,
      .groups = "drop"
    )
  
  # Then aggregate across ROIs to get mean ± SE
  plot_data <- roi_rates %>%
    group_by(animal, condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(events_per_epoch, na.rm = TRUE),
      se_event_rate = sd(events_per_epoch, na.rm = TRUE) / sqrt(n()),
      n_rois = n(),
      .groups = "drop"
    )
  
  # Single faceted plot
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_event_rate,
                              color = animal, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_errorbar(aes(ymin = mean_event_rate - se_event_rate,
                      ymax = mean_event_rate + se_event_rate),
                  width = 0.1, alpha = 0.4, linewidth = 0.4) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.5) +
    facet_wrap(~transition_type, ncol = 2, scales = "free_y") +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation", "WO" = "Washout")
    ) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = title_base,
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)",
      color = "Animal",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
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
