# Modified CCF3.R - Hierarchical Clustering on Whole Day Traces
# This version clusters the entire day's traces once, then applies consistent ordering

library(cluster)
library(dendextend)
library(ggdendro)
library(corrplot)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

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
    if ("Time.from.Start" %in% colnames(stages_df)) {
      time_col <- "Time.from.Start"
      cat("Using 'Time.from.Start' column for timing\n")
    } else if ("Time" %in% colnames(stages_df)) {
      time_col <- "Time"
      cat("Using 'Time' column for timing\n")
    } else {
      stop("Could not find 'Time.from.Start' or 'Time' column in stages_df")
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

# Original CCF function from CCF2.R
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

# Parallel processing version for multi-core systems
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

# Function to perform hierarchical clustering on whole day traces
perform_whole_day_clustering <- function(traces_df, method = "ward.D2", 
                                         distance_metric = "correlation") {
  tryCatch({
    cat("Performing whole-day hierarchical clustering with method:", method, 
        "and distance metric:", distance_metric, "\n")
    
    # Remove time column and convert to numeric matrix
    time_cols <- c("Time", "time", "Time(s)", "Time(s)/Cell Status")
    trace_data <- traces_df[, !names(traces_df) %in% time_cols, drop = FALSE]
    
    # Convert to numeric matrix (skip first row if it contains headers)
    if (any(is.na(suppressWarnings(as.numeric(trace_data[1, ]))))) {
      trace_data <- trace_data[-1, ]
    }
    
    # Convert all columns to numeric
    trace_matrix <- apply(trace_data, 2, function(x) as.numeric(as.character(x)))
    colnames(trace_matrix) <- colnames(trace_data)
    
    # Remove any columns with all NAs or constant values
    valid_cols <- apply(trace_matrix, 2, function(x) {
      !all(is.na(x)) && var(x, na.rm = TRUE) > 0
    })
    
    if (sum(valid_cols) < 2) {
      stop("Insufficient valid traces for clustering (need at least 2)")
    }
    
    trace_matrix <- trace_matrix[, valid_cols, drop = FALSE]
    cat("Using", ncol(trace_matrix), "valid traces for clustering\n")
    
    # Compute correlation matrix between traces (columns)
    cor_matrix <- cor(trace_matrix, use = "pairwise.complete.obs")
    
    # Handle any NAs in correlation matrix
    if (any(is.na(cor_matrix))) {
      warning("NAs found in correlation matrix, using complete observations only")
      cor_matrix[is.na(cor_matrix)] <- 0
    }
    
    # Convert correlation matrix to distance matrix
    if (distance_metric == "correlation") {
      # Use 1 - |correlation| as distance
      dist_matrix <- as.dist(1 - abs(cor_matrix))
    } else if (distance_metric == "euclidean") {
      # Use Euclidean distance on correlation values
      dist_matrix <- dist(cor_matrix, method = "euclidean")
    } else if (distance_metric == "manhattan") {
      # Use Manhattan distance on correlation values
      dist_matrix <- dist(cor_matrix, method = "manhattan")
    } else {
      stop("Unsupported distance metric. Use 'correlation', 'euclidean', or 'manhattan'")
    }
    
    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = method)
    hc$labels <- colnames(cor_matrix)
    
    # Get the ordering from the dendrogram
    trace_order <- hc$labels[hc$order]
    
    cat("Hierarchical clustering completed. Trace order determined.\n")
    
    return(list(
      hclust_result = hc,
      correlation_matrix = cor_matrix,
      trace_order = trace_order,
      valid_traces = colnames(trace_matrix)
    ))
    
  }, error = function(e) {
    cat("ERROR in perform_whole_day_clustering():", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Function to reorder CCF matrix according to clustering order
reorder_ccf_matrix <- function(ccf_matrix, trace_order) {
  tryCatch({
    # Find intersection of available traces in CCF matrix and desired order
    available_traces <- intersect(trace_order, rownames(ccf_matrix))
    available_traces <- intersect(available_traces, colnames(ccf_matrix))
    
    if (length(available_traces) < 2) {
      warning("Insufficient traces available for reordering")
      return(ccf_matrix)
    }
    
    # Reorder both rows and columns
    reordered_matrix <- ccf_matrix[available_traces, available_traces, drop = FALSE]
    
    return(reordered_matrix)
    
  }, error = function(e) {
    cat("ERROR in reorder_ccf_matrix():", conditionMessage(e), "\n")
    return(ccf_matrix)
  })
}

# Function to create dendrogram from whole day clustering
create_whole_day_dendrogram <- function(clustering_result, title = "Whole Day Trace Clustering", 
                                        rotate_labels = TRUE, label_size = 6) {
  tryCatch({
    if (is.null(clustering_result) || is.null(clustering_result$hclust_result)) {
      cat("ERROR: Clustering result is NULL\n")
      return(NULL)
    }
    
    hc_result <- clustering_result$hclust_result
    
    # Convert to dendrogram object
    dend <- as.dendrogram(hc_result)
    
    # Create the plot using ggdendro
    dend_data <- dendro_data(dend)
    
    # Create base plot
    p <- ggplot() +
      geom_segment(data = dend_data$segments, 
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "steelblue", size = 0.5) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = ifelse(rotate_labels, 90, 0), 
                                   hjust = ifelse(rotate_labels, 1, 0.5),
                                   vjust = ifelse(rotate_labels, 0.5, 1),
                                   size = label_size),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "gray30"),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = title,
        x = "Traces (ROIs)",
        y = "Distance"
      )
    
    # Add labels
    p <- p + geom_text(data = dend_data$labels, 
                       aes(x = x, y = y - max(dend_data$segments$y) * 0.02, label = label),
                       hjust = ifelse(rotate_labels, 1.1, 0.5),
                       vjust = ifelse(rotate_labels, 0.5, 1.2),
                       angle = ifelse(rotate_labels, 90, 0),
                       size = label_size/2.5,
                       color = "gray20")
    
    return(p)
    
  }, error = function(e) {
    cat("ERROR in create_whole_day_dendrogram():", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Modified function to create consistently ordered heatmaps
create_stage_heatmaps_consistent <- function(ccf_results, trace_order, 
                                             use_consistent_scale = TRUE) {
  tryCatch({
    if (length(ccf_results) == 0) {
      cat("WARNING: No CCF results to create heatmaps from!\n")
      return(list())
    }
    
    heatmaps <- list()
    
    # Determine consistent color scale across all stages if requested
    global_range <- NULL
    if (use_consistent_scale) {
      all_values <- unlist(lapply(ccf_results, function(x) as.vector(x)))
      all_values <- all_values[!is.na(all_values)]
      if (length(all_values) > 0) {
        global_range <- range(all_values)
        cat("Using consistent color scale across stages:", 
            round(global_range[1], 3), "to", round(global_range[2], 3), "\n")
      }
    }
    
    for (stage_name in names(ccf_results)) {
      cat("Creating consistently ordered heatmap for stage:", stage_name, "\n")
      
      # Reorder the CCF matrix according to clustering
      ccf_matrix <- reorder_ccf_matrix(ccf_results[[stage_name]], trace_order)
      
      if (nrow(ccf_matrix) < 2 || ncol(ccf_matrix) < 2) {
        cat("WARNING: Insufficient data for", stage_name, "heatmap\n")
        next
      }
      
      # Convert matrix to long format for ggplot
      ccf_long <- expand.grid(
        ROI1 = factor(rownames(ccf_matrix), levels = rownames(ccf_matrix)),
        ROI2 = factor(colnames(ccf_matrix), levels = colnames(ccf_matrix))
      )
      ccf_long$CCF <- as.vector(ccf_matrix)
      
      # Create heatmap
      p <- ggplot(ccf_long, aes(x = ROI1, y = ROI2, fill = CCF)) +
        geom_tile(color = "white", size = 0.1) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.title = element_text(size = 10),
          panel.grid = element_blank()
        ) +
        labs(
          title = paste("Cross-Correlation Heatmap -", stage_name, "Stage"),
          subtitle = "Traces ordered by whole-day hierarchical clustering",
          x = "ROI",
          y = "ROI",
          fill = "CCF"
        ) +
        coord_fixed()
      
      # Apply color scale
      if (!is.null(global_range)) {
        p <- p + scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = 0,
          limits = global_range,
          na.value = "gray90"
        )
      } else {
        p <- p + scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = 0,
          na.value = "gray90"
        )
      }
      
      heatmaps[[stage_name]] <- p
    }
    
    return(heatmaps)
    
  }, error = function(e) {
    cat("ERROR in create_stage_heatmaps_consistent():", conditionMessage(e), "\n")
    return(list())
  })
}

# Fixed function to create properly aligned clustered heatmap with dendrogram
# Fixed function to create properly aligned clustered heatmap with dendrogram
create_whole_day_clustered_heatmap_aligned <- function(clustering_result, title = "Whole Day Trace Correlation with Clustering") {
  tryCatch({
    if (is.null(clustering_result) || is.null(clustering_result$hclust_result) || 
        is.null(clustering_result$correlation_matrix)) {
      cat("ERROR: Clustering result or correlation matrix is NULL\n")
      return(NULL)
    }
    
    hc_result <- clustering_result$hclust_result
    cor_matrix <- clustering_result$correlation_matrix
    
    # Reorder correlation matrix according to clustering
    trace_order <- hc_result$labels[hc_result$order]
    reordered_cor_matrix <- cor_matrix[trace_order, trace_order, drop = FALSE]
    
    # Convert correlation matrix to long format for ggplot
    # Keep same order for both axes (no reversal)
    cor_long <- expand.grid(
      ROI1 = factor(rownames(reordered_cor_matrix), levels = rownames(reordered_cor_matrix)),
      ROI2 = factor(colnames(reordered_cor_matrix), levels = rownames(reordered_cor_matrix)) # Same order as ROI1
    )
    cor_long$Correlation <- as.vector(reordered_cor_matrix)
    
    # Create the main heatmap
    heatmap_plot <- ggplot(cor_long, aes(x = ROI1, y = ROI2, fill = Correlation)) +
      geom_tile(color = "white", size = 0.1) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Correlation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 5, "pt")
      ) +
      labs(x = "ROI", y = "ROI") +
      coord_fixed() +
      # KEY FIX: Set exact limits and expansion to match dendrogram
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0))
    
    # Create dendrogram data
    dend <- as.dendrogram(hc_result)
    dend_data <- dendro_data(dend)
    
    # KEY FIX: Map dendrogram positions to heatmap tile centers
    n_labels <- length(trace_order)
    
    # The dendrogram leaf positions (1, 2, 3, ..., n) need to align with 
    # the discrete factor positions in the heatmap
    # ggplot discrete scales place factors at integer positions: 1, 2, 3, etc.
    
    # Extract leaf labels in dendrogram order for positioning
    leaf_labels <- dend_data$labels
    leaf_labels$label <- trace_order[leaf_labels$x]  # Map positions to actual ROI names
    # Remove "C" prefix from labels for dendrogram only
    leaf_labels$label <- gsub("^C", "", leaf_labels$label)
    
    top_dendro_plot <- ggplot() +
      geom_segment(data = dend_data$segments, 
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "steelblue", size = 0.5) +
      # Add leaf labels at the bottom of dendrogram (horizontal)
      geom_text(data = leaf_labels, 
                aes(x = x, y = y - 0.05, label = label),
                angle = 0, hjust = 0.5, vjust = 1, size = 3, color = "black") +
      theme_void() +
      theme(
        plot.margin = margin(5, 5, 0, 5, "pt")
      ) +
      # KEY FIX: Match the heatmap's discrete scale exactly
      scale_x_continuous(
        limits = c(0.5, n_labels + 0.5), 
        expand = c(0, 0),
        breaks = 1:n_labels,
        labels = NULL
      ) +
      scale_y_continuous(expand = c(0.05, 0.05))  # Add slight padding for labels
    
    # Load required libraries
    if (!require(gridExtra, quietly = TRUE)) {
      stop("gridExtra package is required but not installed")
    }
    if (!require(grid, quietly = TRUE)) {
      stop("grid package is required but not installed")
    }
    
    # KEY FIX: Use gtable for precise alignment
    if (!require(gtable, quietly = TRUE)) {
      stop("gtable package is required but not installed")
    }
    
    # Convert plots to grobs for precise control
    dendro_grob <- ggplotGrob(top_dendro_plot)
    heatmap_grob <- ggplotGrob(heatmap_plot)
    
    # Ensure the panel widths match exactly
    # Find the panel columns in both grobs
    dendro_panel_col <- which(dendro_grob$layout$name == "panel")
    heatmap_panel_col <- which(heatmap_grob$layout$name == "panel")
    
    # Match the widths
    max_width <- unit.pmax(dendro_grob$widths, heatmap_grob$widths)
    dendro_grob$widths <- max_width
    heatmap_grob$widths <- max_width
    
    # Create the combined plot
    combined_plot <- gtable_rbind(dendro_grob, heatmap_grob, size = "max")
    
    # Add title
    title_grob <- textGrob(title, gp = gpar(fontsize = 14, fontface = "bold"))
    combined_plot <- gtable_add_rows(combined_plot, unit(1, "line"), 0)
    combined_plot <- gtable_add_grob(combined_plot, title_grob, 1, 1, 1, ncol(combined_plot))
    
    # Draw the plot
    grid.newpage()
    grid.draw(combined_plot)
    
    return(combined_plot)
    
  }, error = function(e) {
    cat("ERROR in create_whole_day_clustered_heatmap_aligned():", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Alternative simpler approach if gtable causes issues
create_whole_day_clustered_heatmap_simple <- function(clustering_result, title = "Whole Day Trace Correlation with Clustering") {
  tryCatch({
    if (is.null(clustering_result) || is.null(clustering_result$hclust_result) || 
        is.null(clustering_result$correlation_matrix)) {
      cat("ERROR: Clustering result or correlation matrix is NULL\n")
      return(NULL)
    }
    
    hc_result <- clustering_result$hclust_result
    cor_matrix <- clustering_result$correlation_matrix
    
    # Reorder correlation matrix according to clustering
    trace_order <- hc_result$labels[hc_result$order]
    reordered_cor_matrix <- cor_matrix[trace_order, trace_order, drop = FALSE]
    
    # Convert correlation matrix to long format for ggplot
    cor_long <- expand.grid(
      ROI1 = factor(rownames(reordered_cor_matrix), levels = rownames(reordered_cor_matrix)),
      ROI2 = factor(colnames(reordered_cor_matrix), levels = rownames(reordered_cor_matrix))
    )
    cor_long$Correlation <- as.vector(reordered_cor_matrix)
    
    # Create the main heatmap with precise positioning
    heatmap_plot <- ggplot(cor_long, aes(x = ROI1, y = ROI2, fill = Correlation)) +
      geom_tile(color = "white", size = 0.1) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Correlation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 5, "pt")
      ) +
      labs(x = "ROI", y = "ROI") +
      coord_fixed() +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0))
    
    # Create dendrogram with exact alignment
    dend <- as.dendrogram(hc_result)
    dend_data <- dendro_data(dend)
    n_labels <- length(trace_order)
    
    # Extract and prepare leaf labels
    leaf_labels <- dend_data$labels
    leaf_labels$label <- trace_order[leaf_labels$x]  # Map positions to actual ROI names
    # Remove "C" prefix from labels for dendrogram only
    leaf_labels$label <- gsub("^C", "", leaf_labels$label)
    
    top_dendro_plot <- ggplot() +
      geom_segment(data = dend_data$segments, 
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "steelblue", size = 0.5) +
      # Add leaf labels at the bottom of dendrogram (horizontal)
      geom_text(data = leaf_labels, 
                aes(x = x, y = y - 0.05, label = label),
                angle = 0, hjust = 0.5, vjust = 1, size = 3, color = "black") +
      theme_void() +
      theme(
        plot.margin = margin(5, 5, 0, 5, "pt"),
        panel.background = element_blank(),
        plot.background = element_blank()
      ) +
      # Match heatmap positioning exactly
      scale_x_continuous(
        limits = c(0.5, n_labels + 0.5), 
        expand = c(0, 0)
      ) +
      scale_y_continuous(expand = c(0.05, 0.05))  # Add padding for labels
    
    # Use cowplot for better alignment control
    if (!require(cowplot, quietly = TRUE)) {
      message("cowplot not available, using gridExtra")
      combined_plot <- grid.arrange(
        top_dendro_plot,
        heatmap_plot,
        ncol = 1, nrow = 2,
        heights = c(1, 4),
        top = textGrob(title, gp = gpar(fontsize = 14, fontface = "bold"))
      )
    } else {
      # Use cowplot for better alignment
      combined_plot <- plot_grid(
        top_dendro_plot,
        heatmap_plot,
        ncol = 1, nrow = 2,
        rel_heights = c(1, 4),
        align = "v",
        axis = "lr"
      )
      
      # Add title
      titled_plot <- plot_grid(
        ggdraw() + draw_label(title, fontface = "bold", size = 14),
        combined_plot,
        ncol = 1, rel_heights = c(0.1, 1)
      )
      
      return(titled_plot)
    }
    
    return(combined_plot)
    
  }, error = function(e) {
    cat("ERROR in create_whole_day_clustered_heatmap_simple():", conditionMessage(e), "\n")
    return(NULL)
  })
}

cat("Improved aligned clustered heatmap functions loaded!\n")
cat("Try create_whole_day_clustered_heatmap_aligned() first.\n") 
cat("If alignment issues persist, try create_whole_day_clustered_heatmap_simple().\n")
cat("These versions ensure dendrogram branches align perfectly with heatmap cells.\n")

# Modified main analysis function with whole-day clustering
run_sleep_stage_ccf_analysis_consistent <- function(states_file = "states_df.csv", 
                                                    traces_file = "Traces.csv",
                                                    output_dir = "results", 
                                                    base_name = "dataset",
                                                    lag = 0, max_lag = 5, 
                                                    use_parallel = TRUE, 
                                                    n_cores = parallel::detectCores() - 1,
                                                    clustering_method = "ward.D2",
                                                    distance_metric = "correlation") {
  tryCatch({
    cat("=== Starting Sleep Stage CCF Analysis with Consistent Ordering ===\n")
    cat("Dataset:", base_name, "\n")
    cat("Clustering method:", clustering_method, "\n")
    cat("Distance metric:", distance_metric, "\n")
    
    # Load data (assuming load_data function exists from original pipeline)
    cat("\nLoading data...\n")
    data <- load_data_custom(states_file, traces_file)
    
    # Data preprocessing (same as original)
    if ("state" %in% colnames(data$stages) && !"state_label" %in% colnames(data$stages)) {
      cat("Creating 'state_label' column from 'state' column\n")
      state_mapping <- c("1" = "Wake", "2" = "NREM", "3" = "REM")
      data$stages$state_label <- state_mapping[as.character(data$stages$state)]
    }
    
    if (!"index" %in% colnames(data$stages)) {
      cat("Creating 'index' column for stages data\n")
      data$stages$index <- 0:(nrow(data$stages) - 1)
    }
    
    # Handle time column
    time_col <- NULL
    if ("Time.from.Start" %in% colnames(data$stages)) {
      time_col <- "Time.from.Start"
    } else if ("Time" %in% colnames(data$stages)) {
      time_col <- "Time"
      data$stages$`Time.from.Start` <- data$stages$Time
      time_col <- "Time.from.Start"
    } else {
      for (col in colnames(data$stages)) {
        if (is.numeric(data$stages[[col]]) && all(data$stages[[col]] >= 0)) {
          cat("Using column '", col, "' as time column\n", sep="")
          data$stages$`Time.from.Start` <- data$stages[[col]]
          time_col <- "Time.from.Start"
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
    
    # STEP 1: Perform whole-day clustering
    cat("\n=== STEP 1: Performing whole-day trace clustering ===\n")
    clustering_result <- perform_whole_day_clustering(
      data$traces, 
      method = clustering_method,
      distance_metric = distance_metric
    )
    
    if (is.null(clustering_result)) {
      stop("Whole-day clustering failed!")
    }
    
    cat("Clustering completed. Trace order established for consistent visualization.\n")
    trace_order <- clustering_result$trace_order
    cat("Trace order (first 10):", paste(head(trace_order, 10), collapse = ", "), "\n")
    
    # STEP 2: Segment traces by sleep stage
    cat("\n=== STEP 2: Segmenting traces by sleep stage ===\n")
    stage_traces <- segment_by_stage(data$stages, data$traces)
    
    if (length(stage_traces) == 0) {
      stop("No stage traces were generated. Cannot proceed with analysis.")
    }
    
    for (stage in names(stage_traces)) {
      cat(stage, "stage:", nrow(stage_traces[[stage]]), "time points\n")
    }
    
    # STEP 3: Compute cross-correlations for each stage (without additional clustering)
    cat("\n=== STEP 3: Computing cross-correlations by stage ===\n")
    if (use_parallel) {
      cat("Using parallel processing with", n_cores, "cores\n")
    } else {
      cat("Using single-core processing\n")
    }
    ccf_results <- compute_ccf_by_stage(stage_traces, lag, max_lag, use_parallel, n_cores)
    
    if (length(ccf_results) == 0) {
      cat("WARNING: No CCF results were generated.\n")
      return(list(stage_traces = stage_traces, ccf_results = list()))
    }
    
    # STEP 4: Create consistently ordered visualizations
    cat("\n=== STEP 4: Creating consistently ordered visualizations ===\n")
    
    # Create output directory
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }
    
    # Create whole-day clustering dendrogram
    whole_day_dendrogram <- create_whole_day_dendrogram(
      clustering_result, 
      title = paste("Whole Day Trace Clustering -", base_name)
    )
    
    # Create consistently ordered heatmaps
    consistent_heatmaps <- create_stage_heatmaps_consistent(ccf_results, trace_order)
    
    # Create whole-day clustered heatmap with dendrogram
    whole_day_clustered_heatmap <- create_whole_day_clustered_heatmap_aligned(
      clustering_result, 
      title = paste("Whole Day Trace Correlation with Clustering -", base_name)
    )
    
    # STEP 5: Save all results
    cat("\n=== STEP 5: Saving results ===\n")
    
    # Save whole-day clustering results
    if (!is.null(clustering_result$correlation_matrix)) {
      cor_file <- file.path(output_dir, paste0(base_name, "_whole_day_correlation_matrix.csv"))
      write.csv(clustering_result$correlation_matrix, cor_file, row.names = TRUE)
      cat("Saved whole-day correlation matrix to:", cor_file, "\n")
    }
    
    # Save trace order
    order_file <- file.path(output_dir, paste0(base_name, "_trace_order.csv"))
    write.csv(data.frame(Position = 1:length(trace_order), Trace = trace_order), 
              order_file, row.names = FALSE)
    cat("Saved trace order to:", order_file, "\n")
    
    # Save whole-day dendrogram
    if (!is.null(whole_day_dendrogram)) {
      dend_file <- file.path(output_dir, paste0(base_name, "_whole_day_dendrogram.png"))
      ggsave(dend_file, whole_day_dendrogram, width = 14, height = 8, dpi = 300)
      cat("Saved whole-day dendrogram to:", dend_file, "\n")
    }
    
    # Save whole-day clustered heatmap
    if (!is.null(whole_day_clustered_heatmap)) {
      clustered_heatmap_file <- file.path(output_dir, paste0(base_name, "_whole_day_clustered_heatmap.png"))
      ggsave(clustered_heatmap_file, whole_day_clustered_heatmap, width = 16, height = 10, dpi = 300)
      cat("Saved whole-day clustered heatmap to:", clustered_heatmap_file, "\n")
    }
    
    # Save consistently ordered CCF matrices
    saved_matrices <- list()
    for (stage in names(ccf_results)) {
      # Reorder matrix according to clustering
      reordered_matrix <- reorder_ccf_matrix(ccf_results[[stage]], trace_order)
      
      matrix_file <- file.path(output_dir, paste0(base_name, "_", stage, "_ccf_matrix_ordered.csv"))
      write.csv(reordered_matrix, matrix_file, row.names = TRUE)
      saved_matrices[[stage]] <- matrix_file
      cat("Saved", stage, "ordered CCF matrix to:", matrix_file, "\n")
    }
    
    # Save consistently ordered heatmaps
    saved_heatmaps <- list()
    for (stage in names(consistent_heatmaps)) {
      heatmap_file <- file.path(output_dir, paste0(base_name, "_", stage, "_heatmap_ordered.png"))
      ggsave(heatmap_file, consistent_heatmaps[[stage]], width = 12, height = 10, dpi = 300)
      saved_heatmaps[[stage]] <- heatmap_file
      cat("Saved", stage, "ordered heatmap to:", heatmap_file, "\n")
    }
    
    # Create comprehensive summary
    summary_file <- file.path(output_dir, paste0(base_name, "_consistent_analysis_summary.txt"))
    summary_content <- paste0(
      "Sleep Stage CCF Analysis with Consistent Ordering\n",
      "================================================\n",
      "Dataset: ", base_name, "\n",
      "Analysis Date: ", Sys.time(), "\n",
      "States file: ", basename(states_file), "\n",
      "Traces file: ", basename(traces_file), "\n\n",
      "Analysis Parameters:\n",
      "- Clustering method: ", clustering_method, "\n",
      "- Distance metric: ", distance_metric, "\n",
      "- Lag: ", lag, "\n",
      "- Max lag: ", max_lag, "\n",
      "- Parallel processing: ", use_parallel, "\n\n",
      "Sleep Stages Distribution:\n"
    )
    
    stage_table <- table(data$stages$state_label)
    for (stage_name in names(stage_table)) {
      summary_content <- paste0(summary_content, "  ", stage_name, ": ", 
                                stage_table[stage_name], " epochs\n")
    }
    
    summary_content <- paste0(summary_content, "\nTrace Ordering:\n")
    summary_content <- paste0(summary_content, "Total traces clustered: ", length(trace_order), "\n")
    summary_content <- paste0(summary_content, "Trace order: ", paste(trace_order, collapse = ", "), "\n\n")
    
    summary_content <- paste0(summary_content, "Key Features:\n")
    summary_content <- paste0(summary_content, "- All heatmaps use consistent trace ordering\n")
    summary_content <- paste0(summary_content, "- Ordering based on whole-day trace similarities\n")
    summary_content <- paste0(summary_content, "- No per-stage clustering performed\n")
    summary_content <- paste0(summary_content, "- Consistent color scales across stages\n")
    
    writeLines(summary_content, summary_file)
    cat("Saved analysis summary to:", summary_file, "\n")
    
    cat("\n=== Analysis with Consistent Ordering Complete ===\n")
    
    return(list(
      stage_traces = stage_traces,
      ccf_results = ccf_results,
      whole_day_clustering = clustering_result,
      trace_order = trace_order,
      whole_day_dendrogram = whole_day_dendrogram,
      whole_day_clustered_heatmap = whole_day_clustered_heatmap,
      consistent_heatmaps = consistent_heatmaps,
      summary_file = summary_file,
      saved_files = list(
        matrices = saved_matrices,
        heatmaps = saved_heatmaps,
        dendrogram = if(!is.null(whole_day_dendrogram)) dend_file else NULL,
        trace_order = order_file,
        correlation_matrix = if(exists("cor_file")) cor_file else NULL
      )
    ))
    
  }, error = function(e) {
    cat("CRITICAL ERROR in run_sleep_stage_ccf_analysis_consistent():", conditionMessage(e), "\n")
    return(NULL)
  })
}

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
    # Read with stringsAsFactors = FALSE to preserve original data
    traces_df <- read.csv(traces_file, check.names = FALSE, stringsAsFactors = FALSE)
    
    cat("Original traces dimensions:", nrow(traces_df), "x", ncol(traces_df), "\n")
    cat("First few column names:", paste(head(colnames(traces_df), 5), collapse = ", "), "\n")
    
    # Fix the empty first column name
    if (colnames(traces_df)[1] == "" || colnames(traces_df)[1] == "X" || 
        grepl("^X\\d*$", colnames(traces_df)[1])) {
      colnames(traces_df)[1] <- "Time"
      cat("Renamed first column to 'Time'\n")
    }
    
    # Print first few values of the first column to debug
    cat("First 5 values in first column:", paste(head(traces_df[,1], 5), collapse = ", "), "\n")
    
    # Remove header rows - this CSV has two header rows
    # First row: column names (already used by read.csv)
    # Second row: "Time(s)/Cell Status, undecided, undecided, ..."
    first_row_first_col <- as.character(traces_df[1, 1])
    cat("First data row, first column:", first_row_first_col, "\n")
    
    # Check if this looks like a header row
    if (grepl("Time.*Cell.*Status|undecided", first_row_first_col, ignore.case = TRUE)) {
      cat("Removing header row with cell status information\n")
      traces_df <- traces_df[-1, ]
    }
    
    # Check if there's still a problematic second row
    if (nrow(traces_df) > 0) {
      second_row_first_col <- as.character(traces_df[1, 1])
      if (grepl("Time.*Cell.*Status|undecided|^[A-Za-z]", second_row_first_col, ignore.case = TRUE) &&
          is.na(suppressWarnings(as.numeric(second_row_first_col)))) {
        cat("Removing second header row:", second_row_first_col, "\n")
        traces_df <- traces_df[-1, ]
      }
    }
    
    # Reset row names after removing header rows
    rownames(traces_df) <- NULL
    
    # Now the first column should contain actual time data starting with "0"
    if (nrow(traces_df) > 0) {
      cat("After header removal, first value:", traces_df[1, 1], "\n")
    }
    
    # Convert time column to numeric with better error handling
    time_col_name <- colnames(traces_df)[1]
    cat("Converting time column to numeric...\n")
    
    # Convert time column to numeric
    original_time_values <- traces_df[[time_col_name]]
    numeric_time_values <- suppressWarnings(as.numeric(as.character(original_time_values)))
    
    # Check conversion success
    conversion_success_rate <- sum(!is.na(numeric_time_values)) / length(numeric_time_values)
    cat("Time conversion success rate:", round(conversion_success_rate * 100, 2), "%\n")
    
    if (conversion_success_rate < 0.9) {
      cat("Low conversion success rate. Showing problematic values:\n")
      problematic_indices <- which(is.na(numeric_time_values))
      if (length(problematic_indices) > 0) {
        cat("Problematic values (first 10):", 
            paste(head(original_time_values[problematic_indices], 10), collapse = ", "), "\n")
      }
    }
    
    traces_df[[time_col_name]] <- numeric_time_values
    
    # Remove rows with invalid time values
    valid_time_rows <- !is.na(traces_df[[time_col_name]])
    if (sum(valid_time_rows) < nrow(traces_df)) {
      invalid_count <- sum(!valid_time_rows)
      cat("Removing", invalid_count, "rows with invalid time values\n")
      traces_df <- traces_df[valid_time_rows, ]
    }
    
    # Rename time column to standard name
    colnames(traces_df)[1] <- "Time"
    
    # Convert all trace columns to numeric (skip Time column)
    trace_cols <- setdiff(colnames(traces_df), "Time")
    cat("Converting", length(trace_cols), "trace columns to numeric...\n")
    
    # Convert trace columns to numeric
    for (col in trace_cols) {
      original_values <- traces_df[[col]]
      numeric_values <- suppressWarnings(as.numeric(as.character(original_values)))
      traces_df[[col]] <- numeric_values
    }
    
    # Remove any remaining rows with all NA trace values
    trace_matrix <- as.matrix(traces_df[, trace_cols])
    valid_trace_rows <- apply(trace_matrix, 1, function(x) !all(is.na(x)))
    if (sum(valid_trace_rows) < nrow(traces_df)) {
      invalid_count <- sum(!valid_trace_rows)
      cat("Removing", invalid_count, "rows with all NA trace values\n")
      traces_df <- traces_df[valid_trace_rows, ]
    }
    
    # Final validation
    if (nrow(traces_df) == 0) {
      stop("No valid data rows remain after processing")
    }
    
    if (all(is.na(traces_df$Time))) {
      stop("All time values are NA after conversion")
    }
    
    # Check for infinite values
    if (any(is.infinite(traces_df$Time))) {
      cat("WARNING: Found infinite values in Time column, removing them\n")
      finite_time_rows <- is.finite(traces_df$Time)
      traces_df <- traces_df[finite_time_rows, ]
    }
    
    # Final check - make sure we have valid data
    if (nrow(traces_df) < 10) {
      stop("Insufficient data after processing (less than 10 rows)")
    }
    
    time_range <- range(traces_df$Time, na.rm = TRUE)
    cat("Final data loaded successfully!\n")
    cat("Traces dimensions:", nrow(traces_df), "x", ncol(traces_df), "\n")
    cat("Time range:", time_range[1], "to", time_range[2], "\n")
    cat("Number of trace columns:", length(trace_cols), "\n")
    
    # Show a sample of the final data
    cat("Sample of processed data:\n")
    print(head(traces_df[, 1:min(5, ncol(traces_df))]))
    
    return(list(stages = stages_df, traces = traces_df))
    
  }, error = function(e) {
    cat("ERROR in load_data_custom():", conditionMessage(e), "\n")
    stop(e)
  })
}