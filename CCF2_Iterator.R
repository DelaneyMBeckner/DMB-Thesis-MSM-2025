# improved_batch_run_ccf_analysis.R

# Source your full analysis pipeline
#source("CCF2.R")

# Function to run batch analysis with custom input and output paths
run_batch_ccf_analysis <- function(input_path = ".", output_path = ".", base_names = NULL) {
  
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
  
  # Helper function to load data with custom filenames
  load_data_custom <- function(states_file, traces_file) {
    tryCatch({
      # Check if files exist
      if (!file.exists(states_file)) {
        stop(paste("File not found:", states_file))
      }
      if (!file.exists(traces_file)) {
        stop(paste("File not found:", traces_file))
      }
      
      cat("Loading sleep stages data from:", states_file, "\n")
      stages_df <- read.csv(states_file)
      
      cat("Loading traces data from:", traces_file, "\n")
      traces_df <- read.csv(traces_file, check.names = FALSE)
      
      # Clean up column names in traces data
      colnames(traces_df)[1] <- "Time"
      traces_df$Time <- as.numeric(sub("Time\\(s\\)/Cell Status", "0", traces_df$Time))
      
      # Keep only numeric columns in traces
      header_row <- traces_df[1, ]
      traces_df <- traces_df[-1, ]  # Remove the header row with "undecided"
      
      # Convert trace data to numeric 
      traces_df[] <- lapply(traces_df, function(x) as.numeric(as.character(x)))
      
      cat("Data loaded successfully!\n")
      
      return(list(stages = stages_df, traces = traces_df))
    }, error = function(e) {
      cat("ERROR in load_data_custom():", conditionMessage(e), "\n")
      stop(e)
    })
  }
  
  # Modified save_ccf_matrices function to include base name
  save_ccf_matrices_with_basename <- function(ccf_results, output_dir, base_name) {
    matrices_saved <- list()
    
    for (stage in names(ccf_results)) {
      # Include base name in filename
      filename <- file.path(output_dir, paste0(base_name, "_", stage, "_ccf_matrix.csv"))
      write.csv(ccf_results[[stage]], filename, row.names = TRUE)
      matrices_saved[[stage]] <- filename
      cat("Saved", stage, "CCF matrix to:", filename, "\n")
    }
    
    return(matrices_saved)
  }
  
  # Modified save_heatmaps function to include base name
  save_heatmaps_with_basename <- function(heatmaps, output_dir, base_name) {
    heatmaps_saved <- list()
    
    for (stage in names(heatmaps)) {
      # Include base name in filename
      filename <- file.path(output_dir, paste0(base_name, "_", stage, "_heatmap.png"))
      ggsave(filename, heatmaps[[stage]], width = 12, height = 10, dpi = 300)
      heatmaps_saved[[stage]] <- filename
      cat("Saved", stage, "heatmap to:", filename, "\n")
    }
    
    return(heatmaps_saved)
  }
  
  # Modified main function that accepts custom filenames and output directory
  run_sleep_stage_ccf_analysis_custom <- function(states_file, traces_file, output_dir = "results", 
                                                  base_name = "dataset",
                                                  lag = 0, max_lag = 5, use_parallel = TRUE, 
                                                  n_cores = parallel::detectCores() - 1) {
    tryCatch({
      cat("=== Starting Sleep Stage CCF Analysis ===\n")
      cat("Dataset:", base_name, "\n")
      cat("States file:", states_file, "\n")
      cat("Traces file:", traces_file, "\n")
      
      # Load data with custom filenames
      cat("\nLoading data...\n")
      data <- load_data_custom(states_file, traces_file)
      
      # Make sure we have the required columns in the stages data
      if ("state" %in% colnames(data$stages) && !"state_label" %in% colnames(data$stages)) {
        cat("Creating 'state_label' column from 'state' column\n")
        state_mapping <- c("1" = "Wake", "2" = "NREM", "3" = "REM")
        data$stages$state_label <- state_mapping[as.character(data$stages$state)]
      }
      
      # Make sure we have index column
      if (!"index" %in% colnames(data$stages)) {
        cat("Creating 'index' column for stages data\n")
        data$stages$index <- 0:(nrow(data$stages) - 1)
      }
      
      # Make sure the time column exists and is properly formatted
      time_col <- NULL
      if ("Time from Start" %in% colnames(data$stages)) {
        time_col <- "Time from Start"
      } else if ("Time" %in% colnames(data$stages)) {
        time_col <- "Time"
        data$stages$`Time from Start` <- data$stages$Time
        time_col <- "Time from Start"
      } else {
        # Look for a likely time column
        for (col in colnames(data$stages)) {
          if (is.numeric(data$stages[[col]]) && all(data$stages[[col]] >= 0)) {
            cat("Using column '", col, "' as time column\n", sep="")
            data$stages$`Time from Start` <- data$stages[[col]]
            time_col <- "Time from Start"
            break
          }
        }
      }
      
      if (is.null(time_col)) {
        stop("Could not find or create a suitable time column in stages data")
      }
      
      # Display sleep stages summary
      cat("\nSleep stages summary:\n")
      print(table(data$stages$state_label))
      
      # Segment traces by sleep stage
      cat("\nSegmenting traces by sleep stage...\n")
      stage_traces <- segment_by_stage(data$stages, data$traces)
      
      # Display summary of segmented data
      cat("\nSegmented traces summary:\n")
      if (length(stage_traces) == 0) {
        cat("ERROR: No stage traces were generated!\n")
        stop("No stage traces were generated. Cannot proceed with analysis.")
      }
      
      for (stage in names(stage_traces)) {
        cat(stage, "stage:", nrow(stage_traces[[stage]]), "time points\n")
      }
      
      # Compute cross-correlations for each stage
      cat("\nComputing cross-correlations for each stage...\n")
      if (use_parallel) {
        cat("Using parallel processing with", n_cores, "cores\n")
      } else {
        cat("Using single-core processing\n")
      }
      ccf_results <- compute_ccf_by_stage(stage_traces, lag, max_lag, use_parallel, n_cores)
      
      # Check if we have any results
      if (length(ccf_results) == 0) {
        cat("WARNING: No CCF results were generated. Cannot create visualizations.\n")
        return(list(stage_traces = stage_traces, ccf_results = list(), heatmaps = list()))
      }
      
      # Create output directory for results
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
      }
      
      # Save CCF matrices as CSV files with base name
      cat("\nSaving CCF matrices as CSV files...\n")
      matrices_saved <- save_ccf_matrices_with_basename(ccf_results, output_dir, base_name)
      
      # Create heatmaps
      cat("\nCreating heatmaps...\n")
      heatmaps <- create_stage_heatmaps(ccf_results)
      
      # Save heatmaps with base name
      cat("\nSaving heatmaps...\n")
      heatmaps_saved <- save_heatmaps_with_basename(heatmaps, output_dir, base_name)
      
      # Save a summary file with base name
      summary_file <- file.path(output_dir, paste0(base_name, "_analysis_summary.txt"))
      cat("Saving analysis summary to:", summary_file, "\n")
      
      # Create summary content
      summary_content <- paste0(
        "Sleep Stage CCF Analysis Summary\n",
        "================================\n",
        "Dataset: ", base_name, "\n",
        "Analysis Date: ", Sys.time(), "\n",
        "States file: ", basename(states_file), "\n",
        "Traces file: ", basename(traces_file), "\n\n",
        "Sleep Stages Distribution:\n"
      )
      
      # Add stage distribution to summary
      stage_table <- table(data$stages$state_label)
      for (stage_name in names(stage_table)) {
        summary_content <- paste0(summary_content, "  ", stage_name, ": ", stage_table[stage_name], " epochs\n")
      }
      
      summary_content <- paste0(summary_content, "\nSegmented Data:\n")
      for (stage in names(stage_traces)) {
        summary_content <- paste0(summary_content, "  ", stage, " stage: ", nrow(stage_traces[[stage]]), " time points\n")
      }
      
      summary_content <- paste0(summary_content, "\nOutput Files Generated:\n")
      for (stage in names(matrices_saved)) {
        summary_content <- paste0(summary_content, "  CCF Matrix: ", basename(matrices_saved[[stage]]), "\n")
      }
      for (stage in names(heatmaps_saved)) {
        summary_content <- paste0(summary_content, "  Heatmap: ", basename(heatmaps_saved[[stage]]), "\n")
      }
      
      # Write summary to file
      writeLines(summary_content, summary_file)
      
      cat("\n=== Analysis Complete ===\n")
      
      # Return results
      return(list(
        stage_traces = stage_traces,
        ccf_results = ccf_results,
        heatmaps = heatmaps,
        summary_file = summary_file
      ))
    }, error = function(e) {
      cat("CRITICAL ERROR in run_sleep_stage_ccf_analysis_custom():", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # Loop through datasets
  for (base in base_names) {
    cat("\n====================================\n")
    cat("Processing:", base, "\n")
    cat("====================================\n")
    
    # Build full file paths
    traces_file <- file.path(input_path, paste0(base, "_Traces.csv"))
    states_file <- file.path(input_path, paste0(base, "_states_df.csv"))
    
    # Check if files exist
    if (!file.exists(traces_file) || !file.exists(states_file)) {
      cat("WARNING: Missing files for", base, "\n")
      cat("  Expected:", traces_file, "\n")
      cat("  Expected:", states_file, "\n")
      failed_count <- failed_count + 1
      failed_datasets <- c(failed_datasets, base)
      next
    }
    
    # Create dataset-specific output directory
    dataset_output_dir <- file.path(output_path, paste0("results_", base))
    
    # Run the analysis with custom filenames and output directory
    tryCatch({
      results <- run_sleep_stage_ccf_analysis_custom(
        states_file = states_file,
        traces_file = traces_file,
        output_dir = dataset_output_dir,
        base_name = base,  # Pass the base name for file naming
        use_parallel = TRUE
      )
      
      if (!is.null(results)) {
        cat("✅ Successfully processed:", base, "\n")
        processed_count <- processed_count + 1
        if (dir.exists(dataset_output_dir)) {
          files <- list.files(dataset_output_dir)
          cat("Results saved to:", dataset_output_dir, "(", length(files), "files )\n")
          cat("Files generated:", paste(files, collapse = ", "), "\n")
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
  cat("🎉 Batch processing complete!\n")
  cat("====================================\n")
  cat("Successfully processed:", processed_count, "datasets\n")
  cat("Failed to process:", failed_count, "datasets\n")
  
  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  # Create overall batch summary file
  batch_summary_file <- file.path(output_path, paste0("batch_processing_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  batch_summary_content <- paste0(
    "Batch CCF Analysis Summary\n",
    "=========================\n",
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
  
  # Restore original working directory
  setwd(original_wd)
  
  # Return summary
  return(list(
    processed_count = processed_count,
    failed_count = failed_count,
    failed_datasets = failed_datasets,
    batch_summary_file = batch_summary_file
  ))
}

# Example usage:
# 
# 1. Run with default settings (input and output in current directory):
# source("improved_batch_run_ccf_analysis.R")
# results <- run_batch_ccf_analysis()
#
# 2. Specify custom input path:
# results <- run_batch_ccf_analysis(input_path = "/path/to/data/files")
#
# 3. Specify both input and output paths:
# results <- run_batch_ccf_analysis(
#   input_path = "/path/to/data/files", 
#   output_path = "/path/to/results"
# )
#
# 4. Process only specific datasets:
# results <- run_batch_ccf_analysis(
#   input_path = "/path/to/data/files",
#   output_path = "/path/to/results",
#   base_names = c("mPFCf5_BL", "mPFCf5_SD")
# )

# Uncomment the line below to run with default settings:
# results <- run_batch_ccf_analysis()