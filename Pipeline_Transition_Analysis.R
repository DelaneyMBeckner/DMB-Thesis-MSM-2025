# ============================================================
# SINGLE ANIMAL PIPELINE: Full Transition Analysis
# ============================================================
# Runs complete workflow for one animal: Extract → CCF → DistCorr
# Computes percentile rankings within animal
#
# Called by: Batch_Process_Animals.R
# ============================================================

library(dplyr)
library(ggplot2)

# Load component scripts (only if not already loaded)
if (!exists("extract_transitions_main")) source("Extract_Transitions.R")
if (!exists("analyze_transition_ccf")) source("Transition_CCF.R")
if (!exists("analyze_transition_distcorr")) source("Transition_DistCorr.R")

#' Run complete analysis pipeline for one animal
#'
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param conditions Vector of conditions (e.g., c("BL", "SD", "WO"))
#' @param transition_types Vector of transition types (e.g., c("Wake2NREM", "NREM2Wake"))
#' @param window_size Window size in epochs (3 or 9)
#' @param data_dir Directory containing input CSV files
#' @param output_dir Directory for saving results
#' @param save_intermediate Save per-animal intermediate results?
#' @param ccf_use_parallel Use parallel processing for CCF?
#' @param ccf_n_cores Number of cores for CCF parallelization
#' @param run_lag_sweep Run CCF lag sweep analysis?
#' @param max_lag Maximum lag for lag sweep
#' @param binning_method Binning method for DistCorr ("equal_width", "percentile", "both")
#' @param n_distance_bins Number of bins for equal-width binning
#' @param n_percentile_bins Number of bins for percentile binning
#' @param fit_exponential Fit exponential model in DistCorr?
#'
#' @return List containing:
#'   - event_rates: ROI-level event rates with percentiles
#'   - correlations: Pair-level correlations with percentiles
#'   - dist_corr: Distance-correlation data with percentiles
#'   - metadata: Analysis metadata
run_animal_pipeline <- function(animal_id,
                                conditions = c("BL", "SD", "WO"),
                                transition_types = c("Wake2NREM", "NREM2Wake"),
                                window_size = 9,
                                data_dir = "E:/Data_Processing/R/Data CSVs",
                                output_dir = "E:/Data_Processing/R/Results",
                                save_intermediate = TRUE,
                                ccf_use_parallel = TRUE,
                                ccf_n_cores = 6,
                                run_lag_sweep = FALSE,
                                max_lag = 5,
                                binning_method = "both",
                                n_distance_bins = 10,
                                n_percentile_bins = 5,
                                fit_exponential = FALSE) {
  
  cat("\n========================================\n")
  cat("PROCESSING ANIMAL:", animal_id, "\n")
  cat("========================================\n")
  cat("Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("Transitions:", paste(transition_types, collapse = ", "), "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("CCF parallel:", ccf_use_parallel, "(", ccf_n_cores, "cores )\n")
  cat("CCF lag sweep:", run_lag_sweep, ifelse(run_lag_sweep, paste0("( max_lag=", max_lag, " )"), ""), "\n")
  cat("DistCorr binning:", binning_method, "\n")
  cat("  Equal-width bins:", n_distance_bins, "\n")
  cat("  Percentile bins:", n_percentile_bins, "\n")
  cat("  Exponential fit:", fit_exponential, "\n\n")
  
  # Load distances once for this animal
  distances_file <- file.path(data_dir, paste0(animal_id, "_roi_distances.csv"))
  if (!file.exists(distances_file)) {
    stop("Distance file not found: ", distances_file)
  }
  distances_df <- read.csv(distances_file, stringsAsFactors = FALSE)
  cat("Loaded distances:", nrow(distances_df), "pairs\n\n")
  
  # Storage for all results
  all_event_rates <- list()
  all_correlations <- list()
  all_dist_corr <- list()
  
  # Process each condition
  for (condition in conditions) {
    recording_id <- paste0(animal_id, "_", condition)
    
    cat("--------------------------------------------------\n")
    cat("CONDITION:", condition, "(", recording_id, ")\n")
    cat("--------------------------------------------------\n")
    
    # Load data files for this condition
    traces_file <- file.path(data_dir, paste0(recording_id, "_Traces.csv"))
    
    # Check traces file exists (extract_transitions_main will check the others)
    if (!file.exists(traces_file)) {
      warning("Traces file not found: ", traces_file, " - Skipping ", condition)
      next
    }
    
    # Load traces
    cat("  Loading traces...\n")
    traces_full <- read.csv(traces_file, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(traces_full)[1] <- "Time"
    if (is.na(suppressWarnings(as.numeric(traces_full[1, 1])))) {
      traces_full <- traces_full[-1, ]
    }
    traces_full[] <- lapply(traces_full, as.numeric)
    
    cat("  Loaded:", nrow(traces_full), "timepoints,", ncol(traces_full) - 1, "ROIs\n\n")
    
    # Process each transition type
    for (trans_type in transition_types) {
      
      cat("  >> Transition Type:", trans_type, "\n")
      
      # Step 1: Extract transitions
      cat("     [1/3] Extracting transitions...\n")
      extraction <- extract_transitions_main(
        states_file = paste0(recording_id, "_states_df.csv"),
        traces_file = paste0(recording_id, "_Traces.csv"),
        events_file = paste0(recording_id, "_Events.csv"),
        recording_id = recording_id,
        input_dir = data_dir,
        window_sizes = window_size,
        save_output = FALSE,
        output_dir = output_dir
      )
      
      # Filter metadata for this transition type and window size
      metadata_df <- extraction$metadata %>% 
        filter(transition_type == trans_type, window_size == !!window_size)
      
      n_trans <- nrow(metadata_df)
      cat("     Found", n_trans, "transitions\n")
      
      if (n_trans == 0) {
        warning("No transitions found for ", recording_id, " - ", trans_type)
        next
      }
      
      # Step 2: Compute CCF
      cat("     [2/3] Computing cross-correlations...\n")
      ccf_results <- analyze_transition_ccf(
        recording_id = recording_id,
        transition_type = trans_type,
        window_size = window_size,
        traces_df = traces_full,
        metadata_df = metadata_df,
        run_lag_sweep = run_lag_sweep,
        max_lag = max_lag,
        use_parallel = ccf_use_parallel,
        n_cores = ccf_n_cores,
        save_output = FALSE
      )
      
      cat("     Computed correlations for", length(ccf_results$ccf_results), "transitions\n")
      
      # Step 3: Distance-correlation analysis
      cat("     [3/3] Analyzing distance-correlation...\n")
      distcorr_results <- analyze_transition_distcorr(
        recording_id = recording_id,
        transition_type = trans_type,
        window_size = window_size,
        distances_df = distances_df,
        metadata_df = metadata_df,
        ccf_results_list = ccf_results$ccf_results,
        binning_method = binning_method,
        n_equal_width_bins = n_distance_bins,
        n_percentile_bins = n_percentile_bins,
        fit_exponential = fit_exponential,
        save_output = save_intermediate,
        output_dir = output_dir
      )
      
      # Extract event rates (from extraction results)
      event_rates <- extraction$events %>%
        group_by(Cell_Name) %>%
        summarize(
          mean_event_rate = n() / n_distinct(transition_id),  # Events per transition
          .groups = 'drop'
        ) %>%
        rename(ROI = Cell_Name) %>%
        mutate(
          animal = animal_id,
          condition = condition,
          transition_type = trans_type,
          window_size = window_size
        )
      
      # Extract correlations (from distcorr results)
      correlations <- distcorr_results$results$dist_corr_all %>%
        group_by(ROI_1, ROI_2, Distance) %>%
        summarize(
          mean_correlation = mean(Correlation, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        mutate(
          animal = animal_id,
          condition = condition,
          transition_type = trans_type,
          window_size = window_size
        )
      
      # Store results
      key <- paste(condition, trans_type, sep = "_")
      all_event_rates[[key]] <- event_rates
      all_correlations[[key]] <- correlations
      all_dist_corr[[key]] <- distcorr_results$results$dist_corr_all %>%
        mutate(
          animal = animal_id,
          condition = condition,
          transition_type = trans_type,
          window_size = window_size
        )
      
      cat("     Complete!\n\n")
    }
    
    # Clear large objects
    rm(traces_full)
    gc()
  }
  
  # Combine all results
  cat("Combining results across conditions...\n")
  event_rates_combined <- bind_rows(all_event_rates)
  correlations_combined <- bind_rows(all_correlations)
  dist_corr_combined <- bind_rows(all_dist_corr)
  
  # Compute percentile rankings WITHIN this animal (across all conditions)
  cat("Computing percentile rankings...\n")
  
  # Event rate percentiles
  event_rates_ranked <- event_rates_combined %>%
    arrange(mean_event_rate) %>%
    mutate(
      activity_rank = row_number(),
      activity_percentile = 100 * (activity_rank - 1) / (n() - 1)
    )
  
  # Correlation percentiles
  correlations_ranked <- correlations_combined %>%
    arrange(mean_correlation) %>%
    mutate(
      correlation_rank = row_number(),
      correlation_percentile = 100 * (correlation_rank - 1) / (n() - 1)
    )
  
  # Distance percentiles (rank by distance, not correlation)
  dist_corr_ranked <- dist_corr_combined %>%
    arrange(Distance) %>%
    mutate(
      distance_rank = row_number(),
      distance_percentile = 100 * (distance_rank - 1) / (n() - 1)
    )
  
  cat("Percentile rankings complete\n")
  
  # Save intermediate results if requested
  if (save_intermediate) {
    cat("\nSaving intermediate results...\n")
    
    write.csv(event_rates_ranked, 
              file.path(output_dir, paste0(animal_id, "_event_rates_ranked.csv")),
              row.names = FALSE)
    write.csv(correlations_ranked,
              file.path(output_dir, paste0(animal_id, "_correlations_ranked.csv")),
              row.names = FALSE)
    write.csv(dist_corr_ranked,
              file.path(output_dir, paste0(animal_id, "_dist_corr_ranked.csv")),
              row.names = FALSE)
    
    cat("Saved intermediate results\n")
  }
  
  # Return results
  cat("\n========================================\n")
  cat("ANIMAL", animal_id, "COMPLETE\n")
  cat("========================================\n\n")
  
  return(list(
    animal_id = animal_id,
    event_rates = event_rates_ranked,
    correlations = correlations_ranked,
    dist_corr = dist_corr_ranked,
    metadata = list(
      conditions = conditions,
      transition_types = transition_types,
      window_size = window_size,
      n_rois = length(unique(event_rates_ranked$ROI)),
      n_pairs = nrow(correlations_ranked)
    )
  ))
}
