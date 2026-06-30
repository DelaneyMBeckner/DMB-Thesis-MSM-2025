# Debug Script for CCF Analysis Issues
# Add this to the beginning of your analysis to diagnose problems

# Function to check if files exist and are readable
check_input_files <- function(input_path, base_names) {
  cat("=== FILE EXISTENCE CHECK ===\n")
  
  for (base in base_names) {
    cat("\nChecking dataset:", base, "\n")
    
    traces_file <- file.path(input_path, paste0(base, "_Traces.csv"))
    states_file <- file.path(input_path, paste0(base, "_states_df.csv"))
    
    cat("Traces file:", traces_file, "\n")
    cat("  Exists:", file.exists(traces_file), "\n")
    if (file.exists(traces_file)) {
      cat("  Size:", file.size(traces_file), "bytes\n")
    }
    
    cat("States file:", states_file, "\n")
    cat("  Exists:", file.exists(states_file), "\n")
    if (file.exists(states_file)) {
      cat("  Size:", file.size(states_file), "bytes\n")
    }
  }
}

# Function to test data loading for a single dataset
test_data_loading <- function(states_file, traces_file) {
  cat("=== TESTING DATA LOADING ===\n")
  cat("States file:", states_file, "\n")
  cat("Traces file:", traces_file, "\n")
  
  tryCatch({
    # Test states file
    if (!file.exists(states_file)) {
      stop("States file does not exist")
    }
    
    cat("Loading states file...\n")
    stages_df <- read.csv(states_file)
    cat("States loaded - dimensions:", dim(stages_df), "\n")
    cat("States columns:", paste(colnames(stages_df), collapse = ", "), "\n")
    cat("First few rows of states:\n")
    print(head(stages_df, 3))
    
    # Test traces file
    if (!file.exists(traces_file)) {
      stop("Traces file does not exist")
    }
    
    cat("\nLoading traces file...\n")
    traces_df <- read.csv(traces_file, check.names = FALSE)
    cat("Traces loaded - dimensions:", dim(traces_df), "\n")
    cat("Traces columns (first 10):", paste(head(colnames(traces_df), 10), collapse = ", "), "\n")
    
    # Check for time column
    cat("Checking for time column...\n")
    time_candidates <- c("Time", "time", "Time(s)", "Time(s)/Cell Status", "Time.from.Start")
    found_time <- intersect(colnames(traces_df), time_candidates)
    cat("Found time columns:", paste(found_time, collapse = ", "), "\n")
    
    return(list(stages = stages_df, traces = traces_df, success = TRUE))
    
  }, error = function(e) {
    cat("ERROR in data loading:", conditionMessage(e), "\n")
    return(list(success = FALSE, error = conditionMessage(e)))
  })
}

# Function to test the segmentation process
test_segmentation <- function(stages_df, traces_df) {
  cat("=== TESTING SEGMENTATION ===\n")
  
  tryCatch({
    # Check required columns in stages
    required_cols <- c("state_label", "state", "Time.from.Start", "Time")
    available_cols <- intersect(required_cols, colnames(stages_df))
    cat("Available required columns:", paste(available_cols, collapse = ", "), "\n")
    
    # Check state column
    if ("state_label" %in% colnames(stages_df)) {
      cat("State labels found:", paste(unique(stages_df$state_label), collapse = ", "), "\n")
      cat("State label counts:\n")
      print(table(stages_df$state_label))
    } else if ("state" %in% colnames(stages_df)) {
      cat("State values found:", paste(unique(stages_df$state), collapse = ", "), "\n")
      cat("State counts:\n")
      print(table(stages_df$state))
    } else {
      cat("WARNING: No state column found!\n")
    }
    
    # Check time alignment
    if ("Time.from.Start" %in% colnames(stages_df)) {
      time_col <- "Time.from.Start"
    } else if ("Time" %in% colnames(stages_df)) {
      time_col <- "Time"
    } else {
      cat("WARNING: No time column found in stages!\n")
      return(FALSE)
    }
    
    cat("Time range in stages:", range(stages_df[[time_col]], na.rm = TRUE), "\n")
    
    # Check traces time column
    traces_time_col <- NULL
    if ("Time" %in% colnames(traces_df)) {
      traces_time_col <- "Time"
    } else {
      traces_time_col <- colnames(traces_df)[1]  # Assume first column is time
    }
    
    cat("Using traces time column:", traces_time_col, "\n")
    cat("Time range in traces:", range(traces_df[[traces_time_col]], na.rm = TRUE), "\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in segmentation test:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Main debugging function
debug_ccf_analysis <- function(input_path = ".", base_names = NULL) {
  cat("🔍 DEBUGGING CCF ANALYSIS\n")
  cat("================================\n")
  
  # Set default base names if not provided
  if (is.null(base_names)) {
    base_names <- c("mPFCm4_BL", "mPFCm4_SD", "mPFCm4_WO")  # Test with just a few
  }
  
  # Check working directory
  cat("Current working directory:", getwd(), "\n")
  
  # Fix path handling - use forward slashes for R compatibility
  normalized_path <- normalizePath(input_path, winslash = "/", mustWork = FALSE)
  cat("Input path:", normalized_path, "\n")
  cat("Input path exists:", dir.exists(normalized_path), "\n")
  
  # List all CSV files in input directory
  if (dir.exists(normalized_path)) {
    csv_files <- list.files(normalized_path, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
    cat("CSV files found:", length(csv_files), "\n")
    if (length(csv_files) > 0) {
      cat("Files:\n")
      for (f in csv_files) {
        cat("  -", f, "\n")
      }
    }
    
    # Also show all files to help identify naming issues
    all_files <- list.files(normalized_path, full.names = FALSE)
    if (length(all_files) > length(csv_files)) {
      cat("All files in directory (", length(all_files), " total):\n")
      for (f in head(all_files, 20)) {  # Show first 20 files
        cat("  -", f, "\n")
      }
      if (length(all_files) > 20) {
        cat("  ... and", length(all_files) - 20, "more files\n")
      }
    }
  } else {
    cat("ERROR: Input directory does not exist or is not accessible\n")
    cat("Check if the path is correct and you have read permissions\n")
    return(invisible(NULL))
  }
  
  # Check file existence with fixed paths
  check_input_files(normalized_path, base_names)
  
  # Test loading the first available dataset
  cat("\n=== TESTING FIRST AVAILABLE DATASET ===\n")
  for (base in base_names) {
    traces_file <- file.path(normalized_path, paste0(base, "_Traces.csv"))
    states_file <- file.path(normalized_path, paste0(base, "_states_df.csv"))
    
    cat("Looking for:\n")
    cat("  Traces:", traces_file, "\n")
    cat("  States:", states_file, "\n")
    
    if (file.exists(traces_file) && file.exists(states_file)) {
      cat("Testing dataset:", base, "\n")
      
      # Test data loading
      load_result <- test_data_loading(states_file, traces_file)
      
      if (load_result$success) {
        # Test segmentation
        seg_result <- test_segmentation(load_result$stages, load_result$traces)
        
        if (seg_result) {
          cat("✅ Dataset", base, "passed basic tests\n")
          
          # Try running just the clustering part if function exists
          cat("\n=== TESTING CLUSTERING ===\n")
          if (exists("perform_whole_day_clustering")) {
            tryCatch({
              clustering_result <- perform_whole_day_clustering(load_result$traces)
              if (!is.null(clustering_result)) {
                cat("✅ Clustering successful, found", length(clustering_result$trace_order), "traces\n")
                cat("First 5 traces:", paste(head(clustering_result$trace_order, 5), collapse = ", "), "\n")
              } else {
                cat("❌ Clustering failed\n")
              }
            }, error = function(e) {
              cat("❌ Clustering error:", conditionMessage(e), "\n")
            })
          } else {
            cat("ℹ️ perform_whole_day_clustering function not loaded - skipping clustering test\n")
            cat("Make sure to source your main CCF3.R file first\n")
          }
        } else {
          cat("❌ Dataset", base, "failed segmentation test\n")
        }
      } else {
        cat("❌ Dataset", base, "failed loading test:", load_result$error, "\n")
      }
      
      break  # Only test the first available dataset
    }
  }
  
  # If no files found, provide helpful suggestions
  if (length(csv_files) == 0) {
    cat("\n=== TROUBLESHOOTING SUGGESTIONS ===\n")
    cat("No CSV files found. Possible issues:\n")
    cat("1. Files are in a different directory - specify full path:\n")
    cat("   debug_ccf_analysis(input_path = 'C:/path/to/your/data')\n")
    cat("2. File names don't match expected pattern:\n")
    cat("   Expected: basename_Traces.csv and basename_states_df.csv\n")
    cat("3. Files have different extensions or capitalization\n")
    cat("4. Working directory is not where you think it is\n")
    cat("   Current directory: ", getwd(), "\n")
  }
  
  cat("\n=== DEBUGGING COMPLETE ===\n")
}

# Quick test function for a single dataset
quick_test_dataset <- function(base_name, input_path = ".") {
  cat("🚀 QUICK TEST FOR:", base_name, "\n")
  
  # Fix path handling
  normalized_path <- normalizePath(input_path, winslash = "/", mustWork = FALSE)
  
  traces_file <- file.path(normalized_path, paste0(base_name, "_Traces.csv"))
  states_file <- file.path(normalized_path, paste0(base_name, "_states_df.csv"))
  
  cat("Files:\n")
  cat("  Traces:", traces_file, "- exists:", file.exists(traces_file), "\n")
  cat("  States:", states_file, "- exists:", file.exists(states_file), "\n")
  
  if (file.exists(traces_file) && file.exists(states_file)) {
    if (exists("run_sleep_stage_ccf_analysis_consistent")) {
      result <- tryCatch({
        # Quick test of the full pipeline
        run_sleep_stage_ccf_analysis_consistent(
          states_file = states_file,
          traces_file = traces_file,
          output_dir = paste0("test_", base_name),
          base_name = base_name
        )
      }, error = function(e) {
        cat("❌ Full pipeline error:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (!is.null(result)) {
        cat("✅ Pipeline completed successfully!\n")
      } else {
        cat("❌ Pipeline returned NULL\n")
      }
    } else {
      cat("ℹ️ run_sleep_stage_ccf_analysis_consistent function not loaded\n")
      cat("Make sure to source your main CCF3.R file first\n")
    }
  }
}

# Run debugging automatically when script is sourced
cat("Debug functions loaded. To debug your analysis, run:\n")
cat("  debug_ccf_analysis()  # Full debugging\n")
cat("  quick_test_dataset('mPFCm4_BL')  # Test single dataset\n")