# Cross-P_Iterator.R
# Batch processing script for Cross-P.R cross-correlation analysis

# Source your full analysis pipeline
# source("Cross-P.R")

# Function to run batch cross-correlation analysis with custom input and output paths
run_batch_cross_p_analysis <- function(input_path = ".", output_path = ".", base_names = NULL) {
  
  # Check if input path exists
  if (!dir.exists(input_path)) {
    stop("Input path does not exist: ", input_path)
  }
  
  # Create output path if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    cat("Created output directory:", output_path, "\n")
  }
  
  # Set default base names if not provided
  if (is.null(base_names)) {
    base_names <- c(
      "mPFCm4_BL", "mPFCm4_SD", "mPFCm4_WO",
      "mPFCm9_BL", "mPFCm9_SD", "mPFCm9_WO",
      "mPFCf5_BL", "mPFCf5_SD", "mPFCf5_WO",
      "mPFCf6_BL", "mPFCf6_SD", "mPFCf6_WO"
    )
  }
  
  cat("Input path:", normalizePath(input_path), "\n")
  cat("Output path:", normalizePath(output_path), "\n")
  cat("Processing", length(base_names), "datasets\n")
  
  # Store original working directory
  original_wd <- getwd()
  
  # Track processing results
  processed_count <- 0
  failed_count <- 0
  failed_datasets <- character(0)
  
  # Modified load_data function for custom filenames - FIXED VERSION
  load_data_custom <- function(states_file, traces_file) {
    tryCatch({
      # Use absolute paths to avoid working directory issues
      states_file <- normalizePath(states_file, mustWork = FALSE)
      traces_file <- normalizePath(traces_file, mustWork = FALSE)
      
      # Check if files exist
      if (!file.exists(states_file)) {
        stop(paste("File not found:", states_file))
      }
      if (!file.exists(traces_file)) {
        stop(paste("File not found:", traces_file))
      }
      
      cat("Loading data from custom files:\n")
      cat("  States:", states_file, "\n")
      cat("  Traces:", traces_file, "\n")
      
      # Load data
      stages_df <- read.csv(states_file)
      traces_df <- read.csv(traces_file, check.names = FALSE)
      
      # Clean traces data
      colnames(traces_df)[1] <- "Time"
      traces_df$Time <- as.numeric(sub("Time\\(s\\)/Cell Status", "0", traces_df$Time))
      traces_df <- traces_df[-1, ]  # Remove header row
      traces_df[] <- lapply(traces_df, function(x) as.numeric(as.character(x)))
      
      # Standardize column names for stages
      if ("state" %in% colnames(stages_df) && !"state_label" %in% colnames(stages_df)) {
        state_mapping <- c("1" = "Wake", "2" = "NREM", "3" = "REM")
        stages_df$state_label <- state_mapping[as.character(stages_df$state)]
      }
      
      # Ensure time column exists
      time_cols <- c("Time from Start", "Time")
      time_col <- intersect(time_cols, colnames(stages_df))[1]
      
      if (is.na(time_col)) {
        # Find first numeric column as time
        numeric_cols <- sapply(stages_df, is.numeric)
        time_col <- names(numeric_cols)[numeric_cols][1]
        stages_df$`Time from Start` <- stages_df[[time_col]]
        time_col <- "Time from Start"
      }
      
      cat("Data loaded successfully!\n")
      list(stages = stages_df, traces = traces_df, time_col = time_col)
      
    }, error = function(e) {
      cat("ERROR in load_data_custom():", conditionMessage(e), "\n")
      stop(e)
    })
  }
  
  # Modified save_results function with base name
  save_results_with_basename <- function(ccf_results, output_dir, base_name) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    saved_files <- list()
    
    for (stage_name in names(ccf_results)) {
      ccf_data <- ccf_results[[stage_name]]
      
      # Save matrices with base name prefix
      corr_file <- file.path(output_dir, paste0(base_name, "_correlations_", stage_name, ".csv"))
      pval_file <- file.path(output_dir, paste0(base_name, "_pvalues_", stage_name, ".csv"))
      sig_file <- file.path(output_dir, paste0(base_name, "_significant_", stage_name, ".csv"))
      heatmap_file <- file.path(output_dir, paste0(base_name, "_heatmap_", stage_name, ".png"))
      
      write.csv(ccf_data$correlations, corr_file)
      write.csv(ccf_data$p_values, pval_file)
      write.csv(ccf_data$p_values < 0.05, sig_file)
      
      # Save heatmap
      heatmap <- create_heatmap(ccf_data, paste("Cross-Correlation:", base_name, "-", stage_name))
      ggsave(heatmap_file, heatmap, width = 10, height = 8, dpi = 300)
      
      saved_files[[stage_name]] <- list(
        correlations = corr_file,
        p_values = pval_file,
        significant = sig_file,
        heatmap = heatmap_file
      )
      
      cat("Saved results for", stage_name, "\n")
    }
    
    return(saved_files)
  }
  
  # Modified main analysis function with custom parameters - FIXED VERSION
  run_cross_p_analysis_custom <- function(states_file, traces_file, output_dir = "results", 
                                          base_name = "dataset", lag = 0, max_lag = 5, 
                                          alpha = 0.05, n_cores = parallel::detectCores() - 1) {
    tryCatch({
      cat("=== Starting Cross-P Analysis ===\n")
      cat("Dataset:", base_name, "\n")
      cat("States file:", states_file, "\n")
      cat("Traces file:", traces_file, "\n")
      cat("Parameters: lag =", lag, ", max_lag =", max_lag, ", alpha =", alpha, "\n")
      
      # Load and process data WITHOUT changing working directory
      cat("\nLoading data...\n")
      data <- load_data_custom(states_file, traces_file)
      
      # Display sleep stages summary
      stage_col <- if ("state_label" %in% colnames(data$stages)) "state_label" else "state"
      cat("\nSleep stages summary:\n")
      print(table(data$stages[[stage_col]]))
      
      # Segment traces by sleep stage
      cat("\nSegmenting traces by sleep stage...\n")
      stage_traces <- segment_by_stage(data$stages, data$traces, data$time_col)
      
      if (length(stage_traces) == 0) {
        stop("No stage traces generated!")
      }
      
      # Display summary of segmented data
      cat("\nSegmented traces summary:\n")
      for (stage in names(stage_traces)) {
        cat("  ", stage, "stage:", nrow(stage_traces[[stage]]), "time points\n")
      }
      
      # Compute cross-correlations
      cat("\nComputing cross-correlations...\n")
      cat("Using", n_cores, "cores for parallel processing\n")
      ccf_results <- compute_ccf_by_stage(stage_traces, lag, max_lag, n_cores)
      
      if (length(ccf_results) == 0) {
        stop("No CCF results generated!")
      }
      
      # Create output directory for results
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
      }
      
      # Save results and create visualizations
      cat("\nSaving results...\n")
      saved_files <- save_results_with_basename(ccf_results, output_dir, base_name)
      
      # Print summary statistics
      cat("\nGenerating summary statistics...\n")
      print_summary(ccf_results, alpha)
      
      # Save a detailed summary file
      summary_file <- file.path(output_dir, paste0(base_name, "_analysis_summary.txt"))
      cat("Saving analysis summary to:", summary_file, "\n")
      
      # Create summary content
      summary_content <- paste0(
        "Cross-P Analysis Summary\n",
        "========================\n",
        "Dataset: ", base_name, "\n",
        "Analysis Date: ", Sys.time(), "\n",
        "States file: ", basename(states_file), "\n",
        "Traces file: ", basename(traces_file), "\n",
        "Parameters:\n",
        "  - Lag: ", lag, "\n",
        "  - Max lag: ", max_lag, "\n",
        "  - Alpha level: ", alpha, "\n",
        "  - Cores used: ", n_cores, "\n\n",
        "Sleep Stages Distribution:\n"
      )
      
      # Add stage distribution to summary
      stage_table <- table(data$stages[[stage_col]])
      for (stage_name in names(stage_table)) {
        summary_content <- paste0(summary_content, "  ", stage_name, ": ", stage_table[stage_name], " epochs\n")
      }
      
      summary_content <- paste0(summary_content, "\nSegmented Data Points:\n")
      for (stage in names(stage_traces)) {
        summary_content <- paste0(summary_content, "  ", stage, " stage: ", nrow(stage_traces[[stage]]), " time points\n")
      }
      
      # Add significance summary
      summary_content <- paste0(summary_content, "\nSignificant Correlations (p < ", alpha, "):\n")
      for (stage_name in names(ccf_results)) {
        pval_matrix <- ccf_results[[stage_name]]$p_values
        diag(pval_matrix) <- NA
        significant_count <- sum(pval_matrix < alpha, na.rm = TRUE)
        total_count <- sum(!is.na(pval_matrix))
        percentage <- round(100 * significant_count / total_count, 1)
        summary_content <- paste0(summary_content, "  ", stage_name, ": ", significant_count, "/", total_count, " (", percentage, "%)\n")
      }
      
      summary_content <- paste0(summary_content, "\nOutput Files Generated:\n")
      for (stage in names(saved_files)) {
        files <- saved_files[[stage]]
        summary_content <- paste0(summary_content, "  ", stage, " stage:\n")
        summary_content <- paste0(summary_content, "    - Correlations: ", basename(files$correlations), "\n")
        summary_content <- paste0(summary_content, "    - P-values: ", basename(files$p_values), "\n")
        summary_content <- paste0(summary_content, "    - Significant: ", basename(files$significant), "\n")
        summary_content <- paste0(summary_content, "    - Heatmap: ", basename(files$heatmap), "\n")
      }
      
      # Write summary to file
      writeLines(summary_content, summary_file)
      
      cat("\n=== Analysis Complete ===\n")
      
      # Return results
      return(list(
        stage_traces = stage_traces,
        ccf_results = ccf_results,
        saved_files = saved_files,
        summary_file = summary_file
      ))
      
    }, error = function(e) {
      cat("CRITICAL ERROR in run_cross_p_analysis_custom():", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # Loop through datasets
  for (base in base_names) {
    cat("\n====================================\n")
    cat("Processing:", base, "\n")
    cat("====================================\n")
    
    # Build full file paths using absolute paths
    traces_file <- file.path(normalizePath(input_path, mustWork = FALSE), paste0(base, "_Traces.csv"))
    states_file <- file.path(normalizePath(input_path, mustWork = FALSE), paste0(base, "_states_df.csv"))
    
    # Check if files exist
    if (!file.exists(traces_file) || !file.exists(states_file)) {
      cat("WARNING: Missing files for", base, "\n")
      cat("  Expected:", traces_file, "\n")
      cat("  Expected:", states_file, "\n")
      
      # List actual files in directory to help debug
      cat("  Files in input directory:\n")
      actual_files <- list.files(input_path, pattern = paste0(base, ".*\\.csv$"))
      if (length(actual_files) > 0) {
        cat("    Found:", paste(actual_files, collapse = ", "), "\n")
      } else {
        cat("    No matching files found\n")
      }
      
      failed_count <- failed_count + 1
      failed_datasets <- c(failed_datasets, base)
      next
    }
    
    # Create dataset-specific output directory
    dataset_output_dir <- file.path(normalizePath(output_path, mustWork = FALSE), paste0("cross_p_results_", base))
    
    # Run the analysis
    tryCatch({
      results <- run_cross_p_analysis_custom(
        states_file = states_file,
        traces_file = traces_file,
        output_dir = dataset_output_dir,
        base_name = base,
        lag = 0,  # Default parameters - can be modified
        max_lag = 5,
        alpha = 0.05,
        n_cores = parallel::detectCores() - 1
      )
      
      if (!is.null(results)) {
        cat("✅ Successfully processed:", base, "\n")
        processed_count <- processed_count + 1
        if (dir.exists(dataset_output_dir)) {
          files <- list.files(dataset_output_dir)
          cat("Results saved to:", dataset_output_dir, "(", length(files), "files)\n")
          cat("Files generated:", paste(head(files, 5), collapse = ", "))
          if (length(files) > 5) cat(", ...")
          cat("\n")
        }
      } else {
        cat("❌ Failed to process:", base, "\n")
        failed_count <- failed_count + 1
        failed_datasets <- c(failed_datasets, base)
      }
    }, error = function(e) {
      cat("❌ Error processing", base, ":", conditionMessage(e), "\n")
      failed_count <- failed_count + 1
      failed_datasets <- c(failed_datasets, base)
    })
  }
  
  # Print summary
  cat("\n====================================\n")
  cat("🎉 Batch Cross-P processing complete!\n")
  cat("====================================\n")
  cat("Successfully processed:", processed_count, "datasets\n")
  cat("Failed to process:", failed_count, "datasets\n")
  
  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  # Create overall batch summary file
  batch_summary_file <- file.path(output_path, paste0("cross_p_batch_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  batch_summary_content <- paste0(
    "Cross-P Batch Analysis Summary\n",
    "==============================\n",
    "Processing Date: ", Sys.time(), "\n",
    "Input Path: ", normalizePath(input_path), "\n",
    "Output Path: ", normalizePath(output_path), "\n\n",
    "Processing Results:\n",
    "Successfully processed: ", processed_count, " datasets\n",
    "Failed to process: ", failed_count, " datasets\n\n"
  )
  
  if (length(failed_datasets) > 0) {
    batch_summary_content <- paste0(batch_summary_content, "Failed datasets:\n")
    for (failed in failed_datasets) {
      batch_summary_content <- paste0(batch_summary_content, "  - ", failed, "\n")
    }
    batch_summary_content <- paste0(batch_summary_content, "\n")
  }
  
  batch_summary_content <- paste0(batch_summary_content, "Successfully processed datasets:\n")
  successful_datasets <- setdiff(base_names, failed_datasets)
  for (success in successful_datasets) {
    batch_summary_content <- paste0(batch_summary_content, "  - ", success, "\n")
  }
  
  writeLines(batch_summary_content, batch_summary_file)
  cat("Batch summary saved to:", batch_summary_file, "\n")
  
  # Return summary
  return(list(
    processed_count = processed_count,
    failed_count = failed_count,
    failed_datasets = failed_datasets,
    successful_datasets = successful_datasets,
    batch_summary_file = batch_summary_file
  ))
}

# Function to check what files are actually available
check_available_files <- function(input_path = ".", pattern = "\\.csv$") {
  if (!dir.exists(input_path)) {
    cat("Directory does not exist:", input_path, "\n")
    return(NULL)
  }
  
  files <- list.files(input_path, pattern = pattern, full.names = FALSE)
  cat("Files found in", normalizePath(input_path), ":\n")
  if (length(files) > 0) {
    for (file in files) {
      cat("  -", file, "\n")
    }
  } else {
    cat("  No .csv files found\n")
  }
  
  # Try to identify base names from files
  traces_files <- files[grepl("_Traces\\.csv$", files)]
  states_files <- files[grepl("_states_df\\.csv$", files)]
  
  if (length(traces_files) > 0) {
    cat("\nTraces files found:\n")
    traces_bases <- sub("_Traces\\.csv$", "", traces_files)
    for (i in seq_along(traces_files)) {
      cat("  Base:", traces_bases[i], "-> File:", traces_files[i], "\n")
    }
  }
  
  if (length(states_files) > 0) {
    cat("\nStates files found:\n")
    states_bases <- sub("_states_df\\.csv$", "", states_files)
    for (i in seq_along(states_files)) {
      cat("  Base:", states_bases[i], "-> File:", states_files[i], "\n")
    }
  }
  
  # Identify complete pairs
  if (length(traces_files) > 0 && length(states_files) > 0) {
    traces_bases <- sub("_Traces\\.csv$", "", traces_files)
    states_bases <- sub("_states_df\\.csv$", "", states_files)
    complete_pairs <- intersect(traces_bases, states_bases)
    
    cat("\nComplete pairs (both traces and states files):\n")
    if (length(complete_pairs) > 0) {
      for (pair in complete_pairs) {
        cat("  -", pair, "\n")
      }
      cat("\nYou can use these base names in run_batch_cross_p_analysis()\n")
    } else {
      cat("  No complete pairs found\n")
    }
    
    return(complete_pairs)
  }
  
  return(NULL)
}

# Convenience function with common parameter sets
run_batch_cross_p_with_params <- function(input_path = ".", output_path = ".", 
                                          base_names = NULL, lag = 0, max_lag = 5, 
                                          alpha = 0.05, n_cores = parallel::detectCores() - 1) {
  cat("Running Cross-P batch analysis with custom parameters:\n")
  cat("  Lag:", lag, "\n")
  cat("  Max lag:", max_lag, "\n") 
  cat("  Alpha:", alpha, "\n")
  cat("  Cores:", n_cores, "\n\n")
  
  return(run_batch_cross_p_analysis(input_path, output_path, base_names))
}

# Example usage:
# 
# 1. Check what files are available first:
# check_available_files("./Data CSVs")
#
# 2. Run with default settings (input and output in current directory):
# source("Cross-P_Iterator.R")
# results <- run_batch_cross_p_analysis()
#
# 3. Specify custom input path:
# results <- run_batch_cross_p_analysis(input_path = "/path/to/data/files")
#
# 4. Specify both input and output paths:
# results <- run_batch_cross_p_analysis(
#   input_path = "/path/to/data/files", 
#   output_path = "/path/to/results"
# )
#
# 5. Process only specific datasets:
# results <- run_batch_cross_p_analysis(
#   input_path = "/path/to/data/files",
#   output_path = "/path/to/results",
#   base_names = c("mPFCf5_BL", "mPFCf5_SD")
# )
#
# 6. Run with custom parameters:
# results <- run_batch_cross_p_with_params(
#   input_path = "/path/to/data/files",
#   output_path = "/path/to/results",
#   lag = 1,
#   max_lag = 10,
#   alpha = 0.01
# )

# Uncomment the line below to run with default settings:
# results <- run_batch_cross_p_analysis()