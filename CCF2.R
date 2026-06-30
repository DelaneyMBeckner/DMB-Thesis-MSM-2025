# Cross-Correlation Analysis by Sleep Stage
# This script segments time series data by sleep stage and computes 
# cross-correlations for each stage separately

library(data.table)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)

# Set working directory to the script location (uncomment and modify if needed)
# setwd("/path/to/your/data/directory")

# Print current working directory for debugging
cat("Current working directory:", getwd(), "\n")

# Load data files with error handling
load_data <- function() {
  tryCatch({
    # Check if files exist
    states_file <- "mPFCf5_BL_states_df.csv"
    traces_file <- "mPFCf5_BL_Traces.csv"
    
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
    cat("ERROR in load_data():", conditionMessage(e), "\n")
    stop(e)
  })
}

# Function to segment traces by sleep stage
segment_by_stage <- function(stages_df, traces_df) {
  tryCatch({
    # Initialize list to store traces for each stage
    stage_traces <- list()
    
    # Print column names to help debug
    cat("Column names in stages_df:", paste(colnames(stages_df), collapse=", "), "\n")
    
    # Determine which column to use for sleep stages
    stage_col <- NULL
    if ("state_label" %in% colnames(stages_df)) {
      stage_col <- "state_label"
      cat("Using 'state_label' column for sleep stages\n")
    } else if ("state" %in% colnames(stages_df)) {
      stage_col <- "state"
      cat("Using 'state' column for sleep stages\n")
    } else {
      stop("Could not find 'state_label' or 'state' column in stages_df")
    }
    
    # Determine which column to use for time
    time_col <- NULL
    if ("Time from Start" %in% colnames(stages_df)) {
      time_col <- "Time from Start"
      cat("Using 'Time from Start' column for timing\n")
    } else if ("Time" %in% colnames(stages_df)) {
      time_col <- "Time"
      cat("Using 'Time' column for timing\n")
    } else {
      stop("Could not find 'Time from Start' or 'Time' column in stages_df")
    }
    
    # Get unique stages
    unique_stages <- unique(stages_df[[stage_col]])
    cat("Found", length(unique_stages), "unique sleep stages:", paste(unique_stages, collapse=", "), "\n")
    
    # Make sure stages_df is sorted by time
    stages_df <- stages_df[order(stages_df[[time_col]]), ]
    
    # Loop through each sleep stage
    for (stage in unique_stages) {
      cat("Processing stage:", stage, "\n")
      
      # Get time ranges for this stage
      stage_times <- stages_df[stages_df[[stage_col]] == stage, ]
      cat("  Found", nrow(stage_times), "segments for this stage\n")
      
      # Initialize dataframe to store traces for this stage
      stage_data <- data.frame()
      
      # Loop through each time segment of this stage
      for (i in 1:nrow(stage_times)) {
        # Get the start time for this segment
        start_time <- stage_times[[time_col]][i]
        
        # Find the row number in the full dataset where this segment starts
        full_idx <- which(stages_df[[time_col]] == start_time & stages_df[[stage_col]] == stage)
        if (length(full_idx) == 0) {
          cat("    Could not find matching row for segment", i, "- skipping\n")
          next
        }
        current_row <- full_idx[1]  # Use the first match if multiple
        
        # Determine the end time for this segment
        if (current_row < nrow(stages_df)) {
          # If there's a next row, use its time as the end time
          end_time <- stages_df[[time_col]][current_row + 1]
        } else {
          # If this is the last entry, add a default interval (5 seconds)
          end_time <- start_time + 5
        }
        
        cat("  Segment", i, "from", start_time, "to", end_time, "\n")
        
        # Find traces within this time range
        segment <- traces_df[traces_df$Time >= start_time & traces_df$Time < end_time, ]
        
        # Add to stage data if there are any traces in this segment
        if (nrow(segment) > 0) {
          cat("    Found", nrow(segment), "data points\n")
          if (nrow(stage_data) == 0) {
            stage_data <- segment
          } else {
            stage_data <- rbind(stage_data, segment)
          }
        } else {
          cat("    No data points found in this segment\n")
        }
      }
      
      # Only add to result if we found data for this stage
      if (nrow(stage_data) > 0) {
        cat("  Collected", nrow(stage_data), "total data points for stage", stage, "\n")
        stage_traces[[stage]] <- stage_data
      } else {
        cat("  WARNING: No data points found for stage", stage, "\n")
      }
    }
    
    return(stage_traces)
  }, error = function(e) {
    cat("ERROR in segment_by_stage():", conditionMessage(e), "\n")
    stop(e)
  })
}

# Function to save CCF matrices as CSV files
save_ccf_matrices <- function(ccf_results, output_dir = "results") {
  tryCatch({
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }
    
    # Check if ccf_results is empty
    if (length(ccf_results) == 0) {
      cat("WARNING: No CCF results to save!\n")
      return(FALSE)
    }
    
    # Save each CCF matrix
    for (stage_name in names(ccf_results)) {
      ccf_matrix <- ccf_results[[stage_name]]
      
      # Create filename
      filename <- file.path(output_dir, paste0("ccf_matrix_", stage_name, ".csv"))
      
      # Save matrix as CSV
      write.csv(ccf_matrix, filename)
      cat("Saved CCF matrix for", stage_name, "to", filename, "\n")
    }
    
    cat("All CCF matrices saved successfully!\n")
    return(TRUE)
  }, error = function(e) {
    cat("ERROR in save_ccf_matrices():", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Original CCF function from paste.txt
compute_pairwise_ccf <- function(data_matrix, lag = 0, max_lag = 10) {
  tryCatch({
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
      ccf_matrix[i, i] <- 1
      
      # Only process upper triangle to avoid redundant computations
      if (i < n_series) {
        for (j in (i+1):n_series) {
          # Create ts objects for just the pair we're analyzing
          series_i <- as.numeric(data_matrix[, i])
          series_j <- as.numeric(data_matrix[, j])
          
          # Calculate CCF and extract the correlation at specified lag
          ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
          lag_index <- which(ccf_result$lag == lag)
          
          if (length(lag_index) > 0) {
            correlation <- ccf_result$acf[lag_index]
            ccf_matrix[i, j] <- correlation
            ccf_matrix[j, i] <- correlation  # Matrix is symmetric for lag 0
          }
          
          # Clean up to free memory
          rm(series_i, series_j, ccf_result)
          gc()
        }
      }
      
      # Print progress
      if (i %% 5 == 0) {
        cat("Processed", i, "of", n_series, "series\n")
      }
    }
    
    return(ccf_matrix)
  }, error = function(e) {
    cat("ERROR in compute_pairwise_ccf():", conditionMessage(e), "\n")
    stop(e)
  })
}

# Approach 2: Parallel processing version for multi-core systems
compute_pairwise_ccf_parallel <- function(data_matrix, lag = 0, max_lag = 10, n_cores = 4) {
  tryCatch({
    n_series <- ncol(data_matrix)
    series_names <- colnames(data_matrix)
    
    if (is.null(series_names)) {
      series_names <- paste0("Series_", 1:n_series)
    }
    
    # Initialize the correlation matrix
    ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
    rownames(ccf_matrix) <- series_names
    colnames(ccf_matrix) <- series_names
    
    # Set up parallel backend
    cat("Setting up parallel cluster with", n_cores, "cores\n")
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Create pairs for upper triangle
    pairs <- expand.grid(i = 1:n_series, j = 1:n_series)
    pairs <- pairs[pairs$i < pairs$j, ]
    
    cat("Starting parallel computation of", nrow(pairs), "series pairs\n")
    
    # Compute correlations in parallel
    results <- foreach(idx = 1:nrow(pairs), .combine = rbind, 
                       .packages = c("stats"), .errorhandling = "pass") %dopar% {
                         i <- pairs$i[idx]
                         j <- pairs$j[idx]
                         
                         series_i <- as.numeric(data_matrix[, i])
                         series_j <- as.numeric(data_matrix[, j])
                         
                         ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
                         lag_index <- which(ccf_result$lag == lag)
                         
                         if (length(lag_index) > 0) {
                           correlation <- ccf_result$acf[lag_index]
                           return(c(i, j, correlation))
                         } else {
                           return(c(i, j, NA))
                         }
                       }
    
    # Stop cluster
    stopCluster(cl)
    cat("Parallel computation completed\n")
    
    # Check if we got results
    if (is.null(results) || nrow(results) == 0) {
      cat("WARNING: No results from parallel computation!\n")
      diag(ccf_matrix) <- 1
      return(ccf_matrix)
    }
    
    # Fill in the correlation matrix
    for (k in 1:nrow(results)) {
      i <- results[k, 1]
      j <- results[k, 2]
      correlation <- results[k, 3]
      
      if (!is.na(correlation)) {
        ccf_matrix[i, j] <- correlation
        ccf_matrix[j, i] <- correlation  # Matrix is symmetric for lag 0
      }
    }
    
    # Set diagonal to 1
    diag(ccf_matrix) <- 1
    
    return(ccf_matrix)
  }, error = function(e) {
    cat("ERROR in compute_pairwise_ccf_parallel():", conditionMessage(e), "\n")
    stop(e)
  })
}

# Function to compute cross-correlation for each stage
compute_ccf_by_stage <- function(stage_traces, lag = 0, max_lag = 10, use_parallel = TRUE, n_cores = 4) {
  tryCatch({
    results <- list()
    
    for (stage_name in names(stage_traces)) {
      cat("Processing stage:", stage_name, "\n")
      
      # Get traces for this stage
      stage_data <- stage_traces[[stage_name]]
      
      # Exclude Time column and convert to matrix for CCF calculation
      data_matrix <- as.matrix(stage_data[, -1])
      
      # Check if we have enough data points for CCF calculation
      if (nrow(data_matrix) <= max_lag + 1) {
        cat("WARNING: Not enough data points for stage", stage_name, 
            "to calculate CCF with max_lag =", max_lag, "\n")
        next
      }
      
      # Calculate cross-correlation
      if (use_parallel) {
        cat("Using parallel processing with", n_cores, "cores for stage", stage_name, "\n")
        ccf_result <- compute_pairwise_ccf_parallel(data_matrix, lag, max_lag, n_cores)
      } else {
        cat("Using single-core processing for stage", stage_name, "\n")
        ccf_result <- compute_pairwise_ccf(data_matrix, lag, max_lag)
      }
      
      # Store result
      if (!is.null(ccf_result)) {
        results[[stage_name]] <- ccf_result
        cat("Successfully computed CCF matrix for stage", stage_name, "\n")
      } else {
        cat("WARNING: Failed to compute CCF matrix for stage", stage_name, "\n")
      }
    }
    
    if (length(results) == 0) {
      cat("WARNING: No CCF results were computed for any stage!\n")
    }
    
    return(results)
  }, error = function(e) {
    cat("ERROR in compute_ccf_by_stage():", conditionMessage(e), "\n")
    return(list())
  })
}

# Function to create heatmaps for each stage
create_stage_heatmaps <- function(ccf_results) {
  tryCatch({
    heatmaps <- list()
    
    if (length(ccf_results) == 0) {
      cat("WARNING: No CCF results to create heatmaps from!\n")
      return(heatmaps)
    }
    
    for (stage_name in names(ccf_results)) {
      cat("Creating heatmap for stage:", stage_name, "\n")
      ccf_matrix <- ccf_results[[stage_name]]
      
      # Create title for this stage
      title <- paste("Cross-Correlation Heatmap for", stage_name, "Stage")
      
      # Create heatmap
      heatmap <- create_correlation_heatmap(ccf_matrix, title)
      
      # Store heatmap
      heatmaps[[stage_name]] <- heatmap
    }
    
    return(heatmaps)
  }, error = function(e) {
    cat("ERROR in create_stage_heatmaps():", conditionMessage(e), "\n")
    return(list())
  })
}

# Function to create a heatmap from the cross-correlation matrix
create_correlation_heatmap <- function(ccf_matrix, title = "Cross-Correlation Heatmap", 
                                       triangle = "full", show_diagonal = TRUE) {
  tryCatch({
    # Convert correlation matrix to long format for ggplot
    # Explicitly use reshape2::melt to avoid the data.table melt conflict
    melted_ccf <- reshape2::melt(ccf_matrix)
    names(melted_ccf) <- c("Series1", "Series2", "Correlation")
    
    # Filter data based on triangle parameter
    if (triangle == "lower") {
      # Keep only lower triangle (including or excluding diagonal based on show_diagonal)
      if (show_diagonal) {
        melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) >= as.numeric(melted_ccf$Series2), ]
      } else {
        melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) > as.numeric(melted_ccf$Series2), ]
      }
    } else if (triangle == "upper") {
      # Keep only upper triangle (including or excluding diagonal based on show_diagonal)
      if (show_diagonal) {
        melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) <= as.numeric(melted_ccf$Series2), ]
      } else {
        melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) < as.numeric(melted_ccf$Series2), ]
      }
    }
    # If triangle == "full", keep all data
    
    # Create a heatmap using ggplot2
    heatmap <- ggplot(melted_ccf, aes(x = Series2, y = Series1, fill = Correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, limits = c(-1, 1)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = title, x = "", y = "")
    
    # Only add text labels if there aren't too many variables
    if (nrow(ccf_matrix) <= 25) {
      heatmap <- heatmap + 
        geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3)
    }
    
    return(heatmap)
  }, error = function(e) {
    cat("ERROR in create_correlation_heatmap():", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Function to save heatmaps as PDF files
save_heatmaps <- function(heatmaps, output_dir = "results") {
  tryCatch({
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }
    
    if (length(heatmaps) == 0) {
      cat("WARNING: No heatmaps to save!\n")
      return(FALSE)
    }
    
    # Save each heatmap
    for (stage_name in names(heatmaps)) {
      heatmap <- heatmaps[[stage_name]]
      
      if (is.null(heatmap)) {
        cat("WARNING: Heatmap for stage", stage_name, "is NULL, skipping...\n")
        next
      }
      
      # Create filename (changed extension from .pdf to .png)
      filename <- file.path(output_dir, paste0("ccf_heatmap_", stage_name, ".png"))
      
      # Save heatmap
      tryCatch({
        ggsave(filename, heatmap, width = 10, height = 8, dpi = 300)
        cat("Saved heatmap for", stage_name, "to", filename, "\n")
      }, error = function(e) {
        cat("ERROR saving heatmap for", stage_name, ":", conditionMessage(e), "\n")
      })
    }
    
    cat("All heatmaps saved successfully!\n")
    return(TRUE)
  }, error = function(e) {
    cat("ERROR in save_heatmaps():", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Main function to run the analysis
run_sleep_stage_ccf_analysis <- function(lag = 0, max_lag = 5, use_parallel = TRUE, n_cores = parallel::detectCores() - 1) {
  tryCatch({
    cat("=== Starting Sleep Stage CCF Analysis ===\n")
    cat("Current working directory:", getwd(), "\n")
    
    # Check required files before starting
    required_files <- c("mPFCf5_BL_states_df.csv", "mPFCf5_BL_Traces.csv")
    missing_files <- required_files[!file.exists(required_files)]
    
    if (length(missing_files) > 0) {
      cat("ERROR: The following required files are missing:\n")
      for (file in missing_files) {
        cat("  -", file, "\n")
      }
      stop("Cannot proceed without required input files")
    }
    
    # Load data
    cat("\nLoading data...\n")
    data <- load_data()
    
    # Make sure we have the required columns in the stages data
    if ("state" %in% colnames(data$stages) && !"state_label" %in% colnames(data$stages)) {
      # Create state_label column from state if it doesn't exist
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
      # Rename to expected column name if needed
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
      # Print more diagnostic information to help identify the problem
      cat("First few rows of stages data:\n")
      print(head(data$stages))
      cat("\nFirst few rows of traces data:\n")
      print(head(data$traces))
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
    results_dir <- "results"
    
    # Save CCF matrices as CSV files
    cat("\nSaving CCF matrices as CSV files...\n")
    matrices_saved <- save_ccf_matrices(ccf_results, results_dir)
    
    # Create heatmaps
    cat("\nCreating heatmaps...\n")
    heatmaps <- create_stage_heatmaps(ccf_results)
    
    # Save heatmaps
    cat("\nSaving heatmaps...\n")
    heatmaps_saved <- save_heatmaps(heatmaps, results_dir)
    
    cat("\n=== Analysis Complete ===\n")
    
    # Return results
    return(list(
      stage_traces = stage_traces,
      ccf_results = ccf_results,
      heatmaps = heatmaps
    ))
  }, error = function(e) {
    cat("CRITICAL ERROR in run_sleep_stage_ccf_analysis():", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Run the analysis with parallel processing (default)
cat("========================================\n")
cat("Starting CCF Analysis Script\n")
cat("----------------------------------------\n")