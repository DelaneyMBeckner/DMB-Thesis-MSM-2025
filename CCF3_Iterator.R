# Modified CCF3_Iterator.R - Batch processing with consistent trace ordering

# Load required libraries
library(cluster)
library(dendextend)
library(ggdendro)
library(ggplot2)
library(parallel)

# Source the modified analysis pipeline (CCF3.R should be loaded first)
# source("CCF3.R")

# Main batch function with consistent ordering
run_batch_ccf_analysis_consistent <- function(input_path = ".", 
                                              output_path = ".", 
                                              base_names = NULL,
                                              clustering_method = "ward.D2",
                                              distance_metric = "correlation") {
  
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
  
  cat("=== CCF Analysis with Consistent Trace Ordering ===\n")
  cat("Input path:", normalizePath(input_path), "\n")
  cat("Output path:", normalizePath(output_path), "\n")
  cat("Processing", length(base_names), "datasets\n")
  cat("Clustering method:", clustering_method, "\n")
  cat("Distance metric:", distance_metric, "\n")
  
  # Store original working directory
  original_wd <- getwd()
  
  # Track processing results
  processed_count <- 0
  failed_count <- 0
  failed_datasets <- character(0)
  all_results <- list()
  
  # Enhanced save function for consistent ordering results
  save_all_results_consistent <- function(results, output_dir, base_name) {
    saved_files <- list()
    
    if (is.null(results) || length(results) == 0) {
      cat("WARNING: No results to save for", base_name, "\n")
      return(saved_files)
    }
    
    # Files are already saved in the main analysis function, just track them
    if (!is.null(results$saved_files)) {
      saved_files <- results$saved_files
      
      # Count different types of files
      matrix_files <- length(saved_files$matrices)
      heatmap_files <- length(saved_files$heatmaps)
      other_files <- sum(!is.null(saved_files$dendrogram), 
                         !is.null(saved_files$trace_order),
                         !is.null(saved_files$correlation_matrix))
      
      cat("Saved files for", base_name, ":\n")
      cat("  CCF matrices:", matrix_files, "\n")
      cat("  Heatmaps:", heatmap_files, "\n")
      cat("  Other files:", other_files, "\n")
    }
    
    return(saved_files)
  }
  
  # Main processing loop
  for (base in base_names) {
    cat("\n=============================================\n")
    cat("Processing:", base, "\n")
    cat("=============================================\n")
    
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
    
    # Run the consistent analysis
    tryCatch({
      results <- run_sleep_stage_ccf_analysis_consistent(
        states_file = states_file,
        traces_file = traces_file,
        output_dir = dataset_output_dir,
        base_name = base,
        use_parallel = TRUE,
        clustering_method = clustering_method,
        distance_metric = distance_metric
      )
      
      if (!is.null(results)) {
        cat("✅ Successfully processed:", base, "\n")
        processed_count <- processed_count + 1
        all_results[[base]] <- results
        
        # Save and track results
        saved_files_info <- save_all_results_consistent(results, dataset_output_dir, base)
        
        if (dir.exists(dataset_output_dir)) {
          files <- list.files(dataset_output_dir)
          cat("Results saved to:", dataset_output_dir, "(", length(files), "files )\n")
          
          # Show specific file types
          ccf_files <- length(grep("ccf_matrix", files))
          heatmap_files <- length(grep("heatmap", files))
          dendro_files <- length(grep("dendrogram", files))
          other_files <- length(files) - ccf_files - heatmap_files - dendro_files
          
          cat("File breakdown: CCF matrices(", ccf_files, "), Heatmaps(", heatmap_files, 
              "), Dendrograms(", dendro_files, "), Other(", other_files, ")\n")
          
          # Show trace order information
          if (!is.null(results$trace_order)) {
            cat("Trace order established:", length(results$trace_order), "traces\n")
            cat("First 5 traces in order:", paste(head(results$trace_order, 5), collapse = ", "), "\n")
          }
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
    
    # Clean up memory after each dataset
    gc()
  }
  
  # Create comprehensive batch summary
  cat("\n=============================================\n")
  cat("🎉 Batch processing with consistent ordering complete!\n")
  cat("=============================================\n")
  cat("Successfully processed:", processed_count, "datasets\n")
  cat("Failed to process:", failed_count, "datasets\n")
  
  if (length(failed_datasets) > 0) {
    cat("Failed datasets:", paste(failed_datasets, collapse = ", "), "\n")
  }
  
  # Create detailed batch summary file
  batch_summary_file <- file.path(output_path, paste0("batch_consistent_analysis_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  batch_summary_content <- paste0(
    "Batch CCF Analysis with Consistent Trace Ordering Summary\n",
    "========================================================\n",
    "Processing Date: ", Sys.time(), "\n",
    "Input Path: ", normalizePath(input_path), "\n",
    "Output Path: ", normalizePath(output_path), "\n",
    "Clustering Method: ", clustering_method, "\n",
    "Distance Metric: ", distance_metric, "\n\n",
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
  
  # Add successfully processed datasets with details
  batch_summary_content <- paste0(batch_summary_content, "Successfully processed datasets:\n")
  successful_datasets <- setdiff(base_names, failed_datasets)
  for (success in successful_datasets) {
    batch_summary_content <- paste0(batch_summary_content, "  - ", success, "\n")
    
    # Add trace order information if available
    if (success %in% names(all_results) && !is.null(all_results[[success]]$trace_order)) {
      n_traces <- length(all_results[[success]]$trace_order)
      batch_summary_content <- paste0(batch_summary_content, "    * Traces clustered: ", n_traces, "\n")
      batch_summary_content <- paste0(batch_summary_content, "    * Order: ", 
                                      paste(head(all_results[[success]]$trace_order, 10), collapse = ", "))
      if (n_traces > 10) {
        batch_summary_content <- paste0(batch_summary_content, ", ... (", n_traces - 10, " more)")
      }
      batch_summary_content <- paste0(batch_summary_content, "\n")
    }
  }
  
  # Add methodological details
  batch_summary_content <- paste0(batch_summary_content, "\nMethodological Details:\n")
  batch_summary_content <- paste0(batch_summary_content, "- Consistent trace ordering applied across all sleep stages\n")
  batch_summary_content <- paste0(batch_summary_content, "- Whole-day hierarchical clustering performed first\n")
  batch_summary_content <- paste0(batch_summary_content, "- All heatmaps use identical trace ordering for comparability\n")
  batch_summary_content <- paste0(batch_summary_content, "- Cross-correlations computed separately for each sleep stage\n")
  batch_summary_content <- paste0(batch_summary_content, "- Consistent color scales applied across stages\n")
  
  writeLines(batch_summary_content, batch_summary_file)
  cat("Batch summary saved to:", batch_summary_file, "\n")
  
  # Create comparison summary across all datasets
  if (processed_count > 1) {
    comparison_file <- file.path(output_path, paste0("cross_dataset_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    comparison_content <- paste0(
      "Cross-Dataset Comparison Summary\n",
      "===============================\n",
      "Analysis Date: ", Sys.time(), "\n",
      "Total Datasets Processed: ", processed_count, "\n\n"
    )
    
    # Add trace count comparison
    comparison_content <- paste0(comparison_content, "Trace Count Comparison:\n")
    for (success in successful_datasets) {
      if (success %in% names(all_results) && !is.null(all_results[[success]]$trace_order)) {
        n_traces <- length(all_results[[success]]$trace_order)
        comparison_content <- paste0(comparison_content, "  ", success, ": ", n_traces, " traces\n")
      }
    }
    
    comparison_content <- paste0(comparison_content, "\nNote: Each dataset uses its own optimized trace ordering\n")
    comparison_content <- paste0(comparison_content, "based on whole-day clustering of that dataset's traces.\n")
    comparison_content <- paste0(comparison_content, "This ensures optimal visualization while maintaining\n")
    comparison_content <- paste0(comparison_content, "consistency within each dataset across sleep stages.\n")
    
    writeLines(comparison_content, comparison_file)
    cat("Cross-dataset comparison saved to:", comparison_file, "\n")
  }
  
  # Restore original working directory
  setwd(original_wd)
  
  # Return comprehensive summary
  return(list(
    processed_count = processed_count,
    failed_count = failed_count,
    failed_datasets = failed_datasets,
    successful_datasets = successful_datasets,
    all_results = all_results,
    batch_summary_file = batch_summary_file,
    clustering_method = clustering_method,
    distance_metric = distance_metric
  ))
}

# Utility function to compare trace orders across datasets
compare_trace_orders <- function(batch_results, output_path = ".") {
  if (is.null(batch_results$all_results) || length(batch_results$all_results) == 0) {
    cat("No results available for comparison\n")
    return(NULL)
  }
  
  cat("Comparing trace orders across datasets...\n")
  
  # Extract trace orders
  trace_orders <- list()
  for (dataset in names(batch_results$all_results)) {
    if (!is.null(batch_results$all_results[[dataset]]$trace_order)) {
      trace_orders[[dataset]] <- batch_results$all_results[[dataset]]$trace_order
    }
  }
  
  if (length(trace_orders) < 2) {
    cat("Need at least 2 datasets for comparison\n")
    return(NULL)
  }
  
  # Create comparison matrix
  comparison_file <- file.path(output_path, paste0("trace_order_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  
  # Find common traces across datasets
  all_traces <- unique(unlist(trace_orders))
  
  # Create matrix showing position of each trace in each dataset
  comparison_matrix <- matrix(NA, nrow = length(all_traces), ncol = length(trace_orders))
  rownames(comparison_matrix) <- all_traces
  colnames(comparison_matrix) <- names(trace_orders)
  
  for (dataset in names(trace_orders)) {
    for (i in seq_along(trace_orders[[dataset]])) {
      trace <- trace_orders[[dataset]][i]
      comparison_matrix[trace, dataset] <- i
    }
  }
  
  write.csv(comparison_matrix, comparison_file, row.names = TRUE)
  cat("Trace order comparison saved to:", comparison_file, "\n")
  
  return(comparison_matrix)
}

# Example usage functions with documentation
print_usage_examples <- function() {
  cat("\n=== CCF3_Iterator Usage Examples ===\n\n")
  
  cat("1. Basic usage with default settings:\n")
  cat("   source('CCF3.R')\n")
  cat("   source('CCF3_Iterator.R')\n")
  cat("   results <- run_batch_ccf_analysis_consistent()\n\n")
  
  cat("2. Specify custom input and output paths:\n")
  cat("   results <- run_batch_ccf_analysis_consistent(\n")
  cat("     input_path = '/path/to/data/files',\n")
  cat("     output_path = '/path/to/results'\n")
  cat("   )\n\n")
  
  cat("3. Process only specific datasets:\n")
  cat("   results <- run_batch_ccf_analysis_consistent(\n")
  cat("     input_path = '/path/to/data',\n")
  cat("     output_path = '/path/to/results',\n")
  cat("     base_names = c('mPFCf5_BL', 'mPFCf5_SD', 'mPFCf5_WO')\n")
  cat("   )\n\n")
  
  cat("4. Use different clustering method:\n")
  cat("   results <- run_batch_ccf_analysis_consistent(\n")
  cat("     clustering_method = 'complete',\n")
  cat("     distance_metric = 'euclidean'\n")
  cat("   )\n\n")
  
  cat("5. Compare trace orders across datasets:\n")
  cat("   comparison <- compare_trace_orders(results, '/path/to/output')\n\n")
  
  cat("Available clustering methods: 'ward.D2', 'complete', 'average', 'single'\n")
  cat("Available distance metrics: 'correlation', 'euclidean', 'manhattan'\n\n")
}

# Print usage examples when the script is loaded
cat("CCF3_Iterator.R loaded successfully!\n")
cat("This script provides batch processing with consistent trace ordering.\n")
print_usage_examples()

# Uncomment the line below to run with default settings:
# results <- run_batch_ccf_analysis_consistent()