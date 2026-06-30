# Transition_CCF.R
# Computes pairwise cross-correlations for transition windows
# Uses metadata to dynamically extract trace windows from full traces file
# Now includes epoch-by-epoch trajectory analysis

library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(viridis)
library(parallel)
library(doParallel)
library(foreach)

# ==================== CONFIGURATION ====================
SAVE_OUTPUT <- FALSE  # Set to TRUE to save results
INPUT_DIR <- "E:/Data_Processing/R/Data CSVs"  # Where raw traces files are
OUTPUT_DIR <- "E:/Data_Processing/R/Results"   # Where to save correlation matrices
METADATA_DIR <- "E:/Data_Processing/R/Results" # Where metadata files are (from Extract_Transitions)

# Parallel processing
USE_PARALLEL <- TRUE  # Set to TRUE to use parallel processing
N_CORES <- 6  # Number of cores for parallel processing

# Analysis mode
RUN_TRAJECTORY_ANALYSIS <- TRUE  # Set to TRUE to analyze epoch-by-epoch trajectories
RUN_WHOLE_WINDOW <- FALSE  # Set to TRUE to also compute whole-window correlations

# Lag analysis parameters
RUN_LAG_SWEEP <- FALSE  # Set to TRUE to perform lag sweep (slower)
MAX_LAG <- 10  # Maximum lag for sweep (in timepoints, 1 = 0.1s at 10Hz)

# Epoch parameters
EPOCH_DURATION <- 10  # seconds per epoch

# Input files
RECORDING_ID <- "mPFCf5_BL"
TRANSITION_TYPE <- "Wake2NREM"  # or "NREM2Wake"
WINDOW_SIZE <- 3  # or 9
TRACES_FILE <- "mPFCf5_BL_Traces.csv"  # Full traces file
# =======================================================


#' Load full traces file
load_full_traces <- function(traces_file, input_dir = INPUT_DIR) {
  traces_path <- file.path(input_dir, traces_file)
  cat("Loading full traces from:", traces_path, "\n")
  
  traces_df <- read.csv(traces_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Fix the first column name if it's empty or has special characters
  if (colnames(traces_df)[1] == "" || colnames(traces_df)[1] == "Time(s)/Cell Status") {
    colnames(traces_df)[1] <- "Time"
  }
  
  # Remove header row if first row contains text instead of numbers
  first_val <- suppressWarnings(as.numeric(traces_df[1, 1]))
  if (is.na(first_val)) {
    cat("  Removing header row from traces\n")
    traces_df <- traces_df[-1, ]
  }
  
  # Convert all to numeric
  traces_df[] <- lapply(traces_df, function(x) as.numeric(as.character(x)))
  
  # Verify Time column
  if (!"Time" %in% colnames(traces_df)) {
    stop("Time column not found after processing. Columns: ", paste(colnames(traces_df), collapse = ", "))
  }
  if (all(is.na(traces_df$Time))) {
    stop("Time column contains only NA values after conversion")
  }
  
  cat("Loaded traces:", nrow(traces_df), "timepoints,", ncol(traces_df) - 1, "ROIs\n")
  
  return(traces_df)
}


#' Load metadata (from file OR from variable)
load_metadata <- function(recording_id = NULL, transition_type = NULL, window_size = NULL,
                          metadata_dir = NULL, metadata_df = NULL) {
  
  # If metadata provided directly, use it (and filter if needed)
  if (!is.null(metadata_df)) {
    cat("Using provided metadata (no file loading)\n")
    
    # Auto-filter if condition parameters are provided
    if (!is.null(transition_type) && !is.null(window_size)) {
      cat("  Filtering for:", transition_type, ",", window_size, "epochs\n")
      
      metadata_df <- metadata_df %>%
        filter(transition_type == !!transition_type, window_size == !!window_size)
      
      cat("  After filtering:", nrow(metadata_df), "transitions\n")
    } else {
      cat("  No filtering applied - using all conditions\n")
      cat("  Metadata:", nrow(metadata_df), "transitions\n")
    }
    
    return(metadata_df)
  }
  
  # Otherwise load from file
  if (is.null(recording_id) || is.null(transition_type) || is.null(window_size) || is.null(metadata_dir)) {
    stop("Must provide either (metadata_df) OR (recording_id + transition_type + window_size + metadata_dir)")
  }
  
  metadata_file <- file.path(metadata_dir,
                             paste0(recording_id, "_", transition_type, "_", 
                                    window_size, "ep_metadata.csv"))
  
  cat("Loading metadata from:", metadata_file, "\n")
  metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(metadata_df), "transitions\n")
  
  return(metadata_df)
}


#' Extract trace window for a single transition
extract_trace_window <- function(traces_full, window_start_time, window_end_time) {
  window_traces <- traces_full[traces_full$Time >= window_start_time & 
                                 traces_full$Time < window_end_time, ]
  return(window_traces)
}


#' Compute pairwise cross-correlation at specified lag
compute_pairwise_ccf <- function(data_matrix, lag = 0, max_lag = 10) {
  # data_matrix should have ROIs as columns, time as rows (Time column excluded)
  
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix)
  
  if (is.null(series_names)) {
    series_names <- paste0("Series_", 1:n_series)
  }
  
  # Initialize the correlation matrix
  ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  rownames(ccf_matrix) <- series_names
  colnames(ccf_matrix) <- series_names
  
  # Loop through each pair of series
  for (i in 1:n_series) {
    # Auto-correlation at lag 0 is 1
    if (lag == 0) {
      ccf_matrix[i, i] <- 1
    }
    
    # Only process upper triangle to save computation
    if (i < n_series) {
      for (j in (i+1):n_series) {
        series_i <- as.numeric(data_matrix[, i])
        series_j <- as.numeric(data_matrix[, j])
        
        # Check for valid data
        if (all(is.na(series_i)) || all(is.na(series_j)) || 
            sd(series_i, na.rm = TRUE) == 0 || sd(series_j, na.rm = TRUE) == 0) {
          ccf_matrix[i, j] <- 0
          ccf_matrix[j, i] <- 0
          next
        }
        
        # Calculate CCF at specified lag
        tryCatch({
          ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE, na.action = na.pass)
          
          # Extract correlation at the desired lag
          lag_index <- which(ccf_result$lag == lag)
          
          if (length(lag_index) > 0) {
            correlation <- ccf_result$acf[lag_index]
            ccf_matrix[i, j] <- correlation
            ccf_matrix[j, i] <- correlation  # CCF is symmetric for our purposes
          } else {
            ccf_matrix[i, j] <- 0
            ccf_matrix[j, i] <- 0
          }
          
        }, error = function(e) {
          ccf_matrix[i, j] <<- 0
          ccf_matrix[j, i] <<- 0
        })
      }
    }
  }
  
  return(ccf_matrix)
}


#' Compute epoch-level correlations for a single transition (trajectory analysis)
compute_epoch_level_ccf <- function(traces_full, transition_row, window_size, 
                                    epoch_duration = EPOCH_DURATION, lag = 0, max_lag = 10) {
  
  # Extract the full window
  window_start <- transition_row$window_start_time
  window_end <- transition_row$window_end_time
  
  # Calculate epoch boundaries within the window
  epoch_starts <- seq(window_start, window_end - epoch_duration, by = epoch_duration)
  
  # Ensure we have the correct number of epochs
  if (length(epoch_starts) != window_size) {
    warning("Expected ", window_size, " epochs but found ", length(epoch_starts), 
            " for transition ", transition_row$transition_id)
  }
  
  # Initialize results list
  epoch_ccf_list <- list()
  
  # Compute CCF for each epoch
  for (i in 1:length(epoch_starts)) {
    epoch_start <- epoch_starts[i]
    epoch_end <- epoch_start + epoch_duration
    
    # Extract traces for this epoch
    epoch_traces <- extract_trace_window(traces_full, epoch_start, epoch_end)
    
    # Remove Time column for correlation computation
    data_matrix <- epoch_traces %>% select(-Time)
    
    # Compute CCF matrix for this epoch
    if (nrow(data_matrix) > 1 && ncol(data_matrix) > 1) {
      ccf_matrix <- compute_pairwise_ccf(data_matrix, lag = lag, max_lag = max_lag)
      
      # Convert to long format for easier analysis
      roi_pairs <- expand.grid(
        ROI_1 = rownames(ccf_matrix),
        ROI_2 = colnames(ccf_matrix),
        stringsAsFactors = FALSE
      )
      
      roi_pairs$correlation <- as.vector(ccf_matrix)
      roi_pairs$epoch_in_window <- i - window_size - 1  # Negative numbering: -3, -2, -1
      roi_pairs$transition_id <- transition_row$transition_id
      
      # Remove diagonal (self-correlations)
      roi_pairs <- roi_pairs %>%
        filter(ROI_1 != ROI_2)
      
      epoch_ccf_list[[i]] <- roi_pairs
    }
  }
  
  # Combine all epochs for this transition
  if (length(epoch_ccf_list) > 0) {
    transition_results <- bind_rows(epoch_ccf_list)
    return(transition_results)
  } else {
    return(NULL)
  }
}


#' Compute trajectory analysis for all transitions (parallelized)
compute_trajectory_ccfs <- function(traces_full, metadata, epoch_duration = EPOCH_DURATION,
                                   lag = 0, max_lag = 10, use_parallel = TRUE, n_cores = 6) {
  
  cat("\nComputing epoch-level correlations for", nrow(metadata), "transitions...\n")
  
  if (use_parallel && nrow(metadata) > 1) {
    cat("Using parallel processing with", n_cores, "cores\n")
    
    # Setup parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Export necessary functions and variables to cluster
    clusterExport(cl, c("compute_epoch_level_ccf", "compute_pairwise_ccf", 
                        "extract_trace_window", "traces_full", "epoch_duration",
                        "lag", "max_lag"),
                  envir = environment())
    
    # Load required packages on each worker
    clusterEvalQ(cl, {
      library(dplyr)
    })
    
    # Process transitions in parallel
    trajectory_results_list <- parLapply(cl, 1:nrow(metadata), function(i) {
      compute_epoch_level_ccf(traces_full, metadata[i, ], metadata$window_size[i],
                             epoch_duration, lag, max_lag)
    })
    
    stopCluster(cl)
    
  } else {
    cat("Using sequential processing\n")
    
    trajectory_results_list <- lapply(1:nrow(metadata), function(i) {
      compute_epoch_level_ccf(traces_full, metadata[i, ], metadata$window_size[i],
                             epoch_duration, lag, max_lag)
    })
  }
  
  # Combine all transitions
  trajectory_results <- bind_rows(trajectory_results_list)
  
  # Add metadata information
  trajectory_results <- trajectory_results %>%
    left_join(metadata %>% select(transition_id, transition_type, window_size),
              by = "transition_id")
  
  cat("Computed", nrow(trajectory_results), "pairwise correlations across epochs\n")
  
  return(trajectory_results)
}


#' Compute whole-window correlations (original approach)
compute_transition_ccfs <- function(traces_full, metadata, lag = 0, max_lag = 10,
                                   run_lag_sweep = FALSE, use_parallel = TRUE, n_cores = 6) {
  
  cat("\nComputing whole-window correlations for", nrow(metadata), "transitions...\n")
  
  if (use_parallel && nrow(metadata) > 1) {
    cat("Using parallel processing with", n_cores, "cores\n")
    
    # Setup parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Export necessary functions and variables
    clusterExport(cl, c("compute_pairwise_ccf", "extract_trace_window", 
                        "traces_full", "lag", "max_lag", "run_lag_sweep"),
                  envir = environment())
    
    # Load required packages
    clusterEvalQ(cl, {
      library(dplyr)
    })
    
    # Process transitions in parallel
    ccf_results <- parLapply(cl, 1:nrow(metadata), function(i) {
      trans_row <- metadata[i, ]
      window_traces <- extract_trace_window(traces_full, 
                                            trans_row$window_start_time,
                                            trans_row$window_end_time)
      data_matrix <- window_traces %>% select(-Time)
      
      result <- list(
        transition_id = trans_row$transition_id,
        ccf_matrix_lag0 = compute_pairwise_ccf(data_matrix, lag = lag, max_lag = max_lag)
      )
      
      return(result)
    })
    
    stopCluster(cl)
    
  } else {
    cat("Using sequential processing\n")
    
    ccf_results <- lapply(1:nrow(metadata), function(i) {
      trans_row <- metadata[i, ]
      window_traces <- extract_trace_window(traces_full, 
                                            trans_row$window_start_time,
                                            trans_row$window_end_time)
      data_matrix <- window_traces %>% select(-Time)
      
      result <- list(
        transition_id = trans_row$transition_id,
        ccf_matrix_lag0 = compute_pairwise_ccf(data_matrix, lag = lag, max_lag = max_lag)
      )
      
      return(result)
    })
  }
  
  # Name the list by transition IDs
  names(ccf_results) <- metadata$transition_id
  
  cat("Computed correlations for", length(ccf_results), "transitions\n")
  
  return(ccf_results)
}


#' Summarize trajectory results
summarize_trajectory <- function(trajectory_data) {
  cat("\nComputing trajectory summary statistics...\n")
  
  # Overall summary by epoch position (mean across all ROI pairs and transitions)
  epoch_summary <- trajectory_data %>%
    group_by(epoch_in_window) %>%
    summarise(
      n_observations = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      median_correlation = median(correlation, na.rm = TRUE),
      q25_correlation = quantile(correlation, 0.25, na.rm = TRUE),
      q75_correlation = quantile(correlation, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Per-ROI-pair summary (averaging across transitions, by epoch)
  roi_pair_summary <- trajectory_data %>%
    group_by(ROI_1, ROI_2, epoch_in_window) %>%
    summarise(
      n_transitions = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Per-transition summary (mean correlation per epoch)
  transition_summary <- trajectory_data %>%
    group_by(transition_id, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nEpoch position summary:\n")
  print(epoch_summary)
  
  return(list(
    epoch_summary = epoch_summary,
    roi_pair_summary = roi_pair_summary,
    transition_summary = transition_summary
  ))
}


#' Create trajectory visualizations
create_trajectory_plots <- function(trajectory_data, summary_stats, 
                                   recording_id, transition_type, window_size) {
  cat("\nCreating trajectory visualizations...\n")
  
  plots <- list()
  
  # 1. Epoch trajectory: Mean correlation by epoch position
  plots$epoch_trajectory <- ggplot(summary_stats$epoch_summary, 
                                   aes(x = epoch_in_window, y = mean_correlation)) +
    geom_line(size = 1.2, color = "blue") +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                     ymax = mean_correlation + se_correlation),
                 width = 0.2, color = "blue") +
    labs(
      title = paste(recording_id, "-", transition_type, "- Correlation Trajectory"),
      subtitle = paste(window_size, "epoch window; mean ± SE across all ROI pairs"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 2. Box plot: Distribution of correlations by epoch position
  plots$epoch_distribution <- ggplot(trajectory_data, 
                                     aes(x = factor(epoch_in_window), 
                                         y = correlation)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Correlation Distribution"),
      subtitle = paste(window_size, "epoch window"),
      x = "Epoch Position Relative to Transition",
      y = "Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 3. Individual transition trajectories by transition ID
  plots$transition_trajectories <- ggplot(summary_stats$transition_summary,
                                          aes(x = epoch_in_window, y = mean_correlation,
                                              group = transition_id, 
                                              color = factor(transition_id))) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Correlation Trajectory by Transition"),
      subtitle = paste(window_size, "epoch window"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation",
      color = "Transition ID"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right"
    )

  
  # 4. Heatmap: Selected ROI pair trajectories (top 20 most variable)
  # Calculate variability for each ROI pair across epochs
  roi_pair_variability <- summary_stats$roi_pair_summary %>%
    group_by(ROI_1, ROI_2) %>%
    summarise(variability = sd(mean_correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(variability)) %>%
    head(20)
  
  top_pairs <- roi_pair_variability %>%
    mutate(pair_id = paste(ROI_1, ROI_2, sep = "-"))
  
  heatmap_data <- summary_stats$roi_pair_summary %>%
    filter(paste(ROI_1, ROI_2, sep = "-") %in% top_pairs$pair_id) %>%
    mutate(pair_id = paste(ROI_1, ROI_2, sep = "-"))
  
  if (nrow(heatmap_data) > 0) {
    plots$roi_pair_heatmap <- ggplot(heatmap_data,
                                     aes(x = epoch_in_window, y = pair_id, fill = mean_correlation)) +
      geom_tile() +
      scale_fill_viridis(option = "D") +
      labs(
        title = paste(recording_id, "-", transition_type, "- ROI Pair Trajectories"),
        subtitle = paste("Top 20 most variable pairs;", window_size, "epoch window"),
        x = "Epoch Position Relative to Transition",
        y = "ROI Pair",
        fill = "Mean\nCorrelation"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 6)
      )
  }
  
  cat("Created", length(plots), "trajectory plots\n")
  
  return(plots)
}


#' Statistical testing for trajectory analysis
test_trajectory_differences <- function(trajectory_data) {
  cat("\nTesting for differences across epoch positions...\n")
  
  # ANOVA to test if correlations differ by epoch position
  aov_model <- aov(correlation ~ factor(epoch_in_window), data = trajectory_data)
  aov_summary <- summary(aov_model)
  
  cat("\nANOVA Results:\n")
  print(aov_summary)
  
  # Pairwise comparisons between consecutive epochs
  epoch_positions <- sort(unique(trajectory_data$epoch_in_window))
  pairwise_results <- data.frame()
  
  for (i in 1:(length(epoch_positions) - 1)) {
    epoch1 <- epoch_positions[i]
    epoch2 <- epoch_positions[i + 1]
    
    # Get data for each epoch (matching ROI pairs)
    paired_data <- trajectory_data %>%
      filter(epoch_in_window %in% c(epoch1, epoch2)) %>%
      select(transition_id, ROI_1, ROI_2, epoch_in_window, correlation) %>%
      pivot_wider(names_from = epoch_in_window, 
                  values_from = correlation,
                  names_prefix = "epoch_")
    
    if (nrow(paired_data) > 0) {
      t_test <- t.test(paired_data[[paste0("epoch_", epoch1)]], 
                       paired_data[[paste0("epoch_", epoch2)]], 
                       paired = TRUE)
      
      pairwise_results <- rbind(pairwise_results, data.frame(
        comparison = paste(epoch1, "vs", epoch2),
        t_statistic = t_test$statistic,
        p_value = t_test$p.value,
        mean_diff = t_test$estimate,
        conf_low = t_test$conf.int[1],
        conf_high = t_test$conf.int[2]
      ))
    }
  }
  
  cat("\nPairwise comparisons (paired t-tests):\n")
  print(pairwise_results)
  
  return(list(
    anova = aov_summary,
    pairwise = pairwise_results
  ))
}


#' Save trajectory outputs
save_trajectory_outputs <- function(trajectory_data, summary_stats, plots, stats_tests,
                                   recording_id, transition_type, window_size, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving trajectory outputs to:", output_dir, "\n")
  
  # Base filename
  base_name <- paste0(recording_id, "_", transition_type, "_", window_size, "ep")
  
  # Save detailed trajectory data
  trajectory_file <- file.path(output_dir, paste0(base_name, "_ccf_trajectory.csv"))
  write.csv(trajectory_data, trajectory_file, row.names = FALSE)
  cat("  Saved:", trajectory_file, "\n")
  
  # Save epoch summary
  epoch_summary_file <- file.path(output_dir, paste0(base_name, "_ccf_epoch_summary.csv"))
  write.csv(summary_stats$epoch_summary, epoch_summary_file, row.names = FALSE)
  cat("  Saved:", epoch_summary_file, "\n")
  
  # Save ROI pair summary
  roi_pair_summary_file <- file.path(output_dir, paste0(base_name, "_ccf_roi_pair_summary.csv"))
  write.csv(summary_stats$roi_pair_summary, roi_pair_summary_file, row.names = FALSE)
  cat("  Saved:", roi_pair_summary_file, "\n")
  
  # Save transition summary
  transition_summary_file <- file.path(output_dir, paste0(base_name, "_ccf_transition_summary.csv"))
  write.csv(summary_stats$transition_summary, transition_summary_file, row.names = FALSE)
  cat("  Saved:", transition_summary_file, "\n")
  
  # Save statistical tests
  stats_file <- file.path(output_dir, paste0(base_name, "_ccf_trajectory_stats.csv"))
  write.csv(stats_tests$pairwise, stats_file, row.names = FALSE)
  cat("  Saved:", stats_file, "\n")
  
  # Save plots
  for (plot_name in names(plots)) {
    plot_file <- file.path(output_dir, paste0(base_name, "_ccf_", plot_name, ".png"))
    ggsave(plot_file, plots[[plot_name]], width = 10, height = 8, dpi = 300)
    cat("  Saved:", plot_file, "\n")
  }
  
  cat("\nAll trajectory outputs saved!\n")
}


#' Save whole-window correlation matrices
save_ccf_matrices <- function(ccf_results, recording_id, transition_type, 
                              window_size, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving whole-window correlation matrices to:", output_dir, "\n")
  
  files_saved <- 0
  
  for (trans_id in names(ccf_results)) {
    ccf_matrix_lag0 <- ccf_results[[trans_id]]$ccf_matrix_lag0
    
    filename_lag0 <- file.path(output_dir,
                               paste0(recording_id, "_", transition_type, "_", 
                                      window_size, "ep_transition_", 
                                      sprintf("%03d", as.numeric(trans_id)), "_ccf_lag0.csv"))
    
    write.csv(ccf_matrix_lag0, filename_lag0, row.names = TRUE)
    files_saved <- files_saved + 1
  }
  
  cat("Saved", files_saved, "whole-window correlation matrices\n")
}


#' Main analysis function
analyze_transition_ccf <- function(recording_id = RECORDING_ID,
                                   transition_type = TRANSITION_TYPE,
                                   window_size = WINDOW_SIZE,
                                   traces_file = TRACES_FILE,
                                   input_dir = INPUT_DIR,
                                   metadata_dir = METADATA_DIR,
                                   traces_df = NULL,
                                   metadata_df = NULL,
                                   epoch_duration = EPOCH_DURATION,
                                   lag = 0,
                                   run_trajectory_analysis = RUN_TRAJECTORY_ANALYSIS,
                                   run_whole_window = RUN_WHOLE_WINDOW,
                                   run_lag_sweep = RUN_LAG_SWEEP,
                                   max_lag = MAX_LAG,
                                   use_parallel = USE_PARALLEL,
                                   n_cores = N_CORES,
                                   save_output = SAVE_OUTPUT,
                                   output_dir = OUTPUT_DIR) {
  
  cat("==================================================\n")
  cat("Transition Cross-Correlation Analysis\n")
  cat("==================================================\n")
  
  # Determine data source
  if (!is.null(traces_df) && !is.null(metadata_df)) {
    cat("Data source: Provided data frames (pipeline mode)\n")
    traces_full <- traces_df
    # Extract info from metadata for labeling
    if (is.null(recording_id)) recording_id <- unique(metadata_df$recording_id)[1]
    if (is.null(transition_type)) transition_type <- unique(metadata_df$transition_type)[1]
    if (is.null(window_size)) window_size <- unique(metadata_df$window_size)[1]
  } else {
    cat("Data source: CSV files (standalone mode)\n")
    cat("Input directory:", input_dir, "\n")
    cat("Metadata directory:", metadata_dir, "\n")
    
    # Load traces
    traces_full <- load_full_traces(traces_file, input_dir)
  }
  
  cat("Recording ID:", recording_id, "\n")
  cat("Transition type:", transition_type, "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("Epoch duration:", epoch_duration, "seconds\n")
  cat("CCF parameters: lag =", lag, "\n")
  cat("Trajectory analysis:", run_trajectory_analysis, "\n")
  cat("Whole-window analysis:", run_whole_window, "\n")
  cat("Parallel processing:", use_parallel, "\n")
  if (use_parallel) {
    cat("  Cores:", n_cores, "\n")
  }
  cat("Output directory:", output_dir, "\n")
  cat("Save output:", save_output, "\n")
  cat("==================================================\n\n")
  
  # Load metadata (from file or use provided, with filtering)
  metadata <- load_metadata(recording_id, transition_type, window_size, 
                            metadata_dir, metadata_df)
  
  # Initialize results list
  results <- list()
  
  # Run trajectory analysis if requested
  if (run_trajectory_analysis) {
    cat("\n========== TRAJECTORY ANALYSIS ==========\n")
    
    trajectory_data <- compute_trajectory_ccfs(traces_full, metadata, epoch_duration,
                                              lag, max_lag, use_parallel, n_cores)
    
    summary_stats <- summarize_trajectory(trajectory_data)
    
    plots <- create_trajectory_plots(trajectory_data, summary_stats,
                                    recording_id, transition_type, window_size)
    
    stats_tests <- test_trajectory_differences(trajectory_data)
    
    if (save_output) {
      save_trajectory_outputs(trajectory_data, summary_stats, plots, stats_tests,
                             recording_id, transition_type, window_size, output_dir)
    }
    
    results$trajectory_data <- trajectory_data
    results$summary_stats <- summary_stats
    results$plots <- plots
    results$stats_tests <- stats_tests
  }
  
  # Run whole-window analysis if requested
  if (run_whole_window) {
    cat("\n========== WHOLE-WINDOW ANALYSIS ==========\n")
    
    ccf_results <- compute_transition_ccfs(traces_full, metadata, lag, max_lag,
                                          run_lag_sweep, use_parallel, n_cores)
    
    if (save_output) {
      save_ccf_matrices(ccf_results, recording_id, transition_type, 
                       window_size, output_dir)
    }
    
    results$ccf_results <- ccf_results
  }
  
  results$metadata <- metadata
  
  cat("\n==================================================\n")
  cat("CCF analysis complete!\n")
  cat("  Transitions analyzed:", nrow(metadata), "\n")
  if (run_trajectory_analysis) {
    cat("  Trajectory analysis: COMPLETE\n")
    cat("    Total correlations:", nrow(trajectory_data), "\n")
  }
  if (run_whole_window) {
    cat("  Whole-window analysis: COMPLETE\n")
  }
  cat("==================================================\n")
  
  return(results)
}


# =========================== MAIN ===========================
# Run the analysis
if (sys.nframe() == 0) {
  results <- analyze_transition_ccf()
  
  # Display key plots if trajectory analysis was run
  if (RUN_TRAJECTORY_ANALYSIS && !is.null(results$plots)) {
    print(results$plots$epoch_trajectory)
    print(results$plots$epoch_distribution)
    if (!is.null(results$plots$transition_trajectories)) {
      print(results$plots$transition_trajectories)
    }
  }
}
