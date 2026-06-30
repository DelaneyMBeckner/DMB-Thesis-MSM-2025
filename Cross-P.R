# Optimized Cross-Correlation Analysis by Sleep Stage
# Modified version with dual color scale for heatmaps

library(data.table)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)

# Load and validate data
load_data <- function() {
  files <- c(states = "mPFCf5_BL_states_df.csv", traces = "mPFCf5_BL_Traces.csv")
  
  if (!all(file.exists(files))) {
    stop("Missing files: ", paste(files[!file.exists(files)], collapse = ", "))
  }
  
  cat("Loading data from:", getwd(), "\n")
  
  # Load data
  stages_df <- read.csv(files["states"])
  traces_df <- read.csv(files["traces"], check.names = FALSE)
  
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
}

# Segment traces by sleep stage
segment_by_stage <- function(stages_df, traces_df, time_col) {
  stage_col <- if ("state_label" %in% colnames(stages_df)) "state_label" else "state"
  stages_df <- stages_df[order(stages_df[[time_col]]), ]
  
  unique_stages <- unique(stages_df[[stage_col]])
  cat("Processing", length(unique_stages), "stages:", paste(unique_stages, collapse = ", "), "\n")
  
  stage_traces <- list()
  
  for (stage in unique_stages) {
    stage_times <- stages_df[stages_df[[stage_col]] == stage, ]
    stage_data <- data.frame()
    
    for (i in 1:nrow(stage_times)) {
      start_time <- stage_times[[time_col]][i]
      
      # Find end time
      full_idx <- which(stages_df[[time_col]] == start_time & stages_df[[stage_col]] == stage)[1]
      end_time <- if (full_idx < nrow(stages_df)) {
        stages_df[[time_col]][full_idx + 1]
      } else {
        start_time + 5  # Default 5-second interval
      }
      
      # Extract segment
      segment <- traces_df[traces_df$Time >= start_time & traces_df$Time < end_time, ]
      stage_data <- rbind(stage_data, segment)
    }
    
    if (nrow(stage_data) > 0) {
      stage_traces[[stage]] <- stage_data
      cat("Stage", stage, ":", nrow(stage_data), "data points\n")
    }
  }
  
  stage_traces
}

# Parallel CCF computation with p-values
compute_ccf_parallel <- function(data_matrix, lag = 0, max_lag = 10, n_cores = parallel::detectCores() - 1) {
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix) %||% paste0("Series_", 1:n_series)
  
  # Initialize matrices
  ccf_matrix <- diag(n_series)
  pval_matrix <- matrix(0, n_series, n_series)
  dimnames(ccf_matrix) <- dimnames(pval_matrix) <- list(series_names, series_names)
  
  # Set up parallel processing
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  # Create pairs (upper triangle only)
  pairs <- expand.grid(i = 1:n_series, j = 1:n_series)
  pairs <- pairs[pairs$i < pairs$j, ]
  
  cat("Computing", nrow(pairs), "correlations with", n_cores, "cores\n")
  
  # Parallel computation
  results <- foreach(idx = 1:nrow(pairs), .combine = rbind, 
                     .packages = "stats", .errorhandling = "pass") %dopar% {
                       i <- pairs$i[idx]
                       j <- pairs$j[idx]
                       
                       series_i <- as.numeric(data_matrix[, i])
                       series_j <- as.numeric(data_matrix[, j])
                       
                       # Remove NAs
                       valid <- !is.na(series_i) & !is.na(series_j)
                       series_i <- series_i[valid]
                       series_j <- series_j[valid]
                       
                       if (length(series_i) < 3) return(c(i, j, NA, NA))
                       
                       if (lag == 0) {
                         cor_test <- cor.test(series_i, series_j, method = "pearson")
                         c(i, j, cor_test$estimate, cor_test$p.value)
                       } else {
                         ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
                         lag_index <- which(ccf_result$lag == lag)
                         
                         if (length(lag_index) > 0) {
                           correlation <- ccf_result$acf[lag_index]
                           n <- length(series_i) - abs(lag)
                           
                           if (n > 3) {
                             t_stat <- correlation * sqrt((n - 2) / (1 - correlation^2))
                             p_value <- 2 * (1 - pt(abs(t_stat), df = n - 2))
                           } else {
                             p_value <- NA
                           }
                           c(i, j, correlation, p_value)
                         } else {
                           c(i, j, NA, NA)
                         }
                       }
                     }
  
  # Fill matrices
  if (!is.null(results) && nrow(results) > 0) {
    for (k in 1:nrow(results)) {
      i <- results[k, 1]
      j <- results[k, 2]
      correlation <- results[k, 3]
      p_value <- results[k, 4]
      
      if (!is.na(correlation)) {
        ccf_matrix[i, j] <- ccf_matrix[j, i] <- correlation
        pval_matrix[i, j] <- pval_matrix[j, i] <- p_value
      }
    }
  }
  
  list(correlations = ccf_matrix, p_values = pval_matrix)
}

# Process all stages
compute_ccf_by_stage <- function(stage_traces, lag = 0, max_lag = 10, n_cores = parallel::detectCores() - 1) {
  results <- list()
  
  for (stage_name in names(stage_traces)) {
    cat("Processing stage:", stage_name, "\n")
    stage_data <- stage_traces[[stage_name]]
    data_matrix <- as.matrix(stage_data[, -1])  # Exclude Time column
    
    if (nrow(data_matrix) <= max_lag + 1) {
      cat("Insufficient data for stage", stage_name, "\n")
      next
    }
    
    results[[stage_name]] <- compute_ccf_parallel(data_matrix, lag, max_lag, n_cores)
  }
  
  results
}

# Create heatmap with dual color scale
create_heatmap <- function(ccf_data, title, alpha = 0.05) {
  ccf_matrix <- ccf_data$correlations
  pval_matrix <- ccf_data$p_values
  
  # Prepare data with separate correlation values for significant/non-significant
  plot_data <- data.frame(
    Series1 = rep(rownames(ccf_matrix), each = ncol(ccf_matrix)),
    Series2 = rep(colnames(ccf_matrix), nrow(ccf_matrix)),
    Correlation = as.vector(ccf_matrix),
    P_Value = as.vector(pval_matrix),
    Significant = as.vector(pval_matrix) < alpha
  )
  
  # Create separate correlation values for each significance level
  plot_data$Corr_Sig <- ifelse(plot_data$Significant, plot_data$Correlation, NA)
  plot_data$Corr_NonSig <- ifelse(!plot_data$Significant, plot_data$Correlation, NA)
  
  # Create the plot with dual color scales
  p <- ggplot(plot_data, aes(x = Series2, y = Series1)) +
    # Non-significant correlations: cyan-black-magenta scale
    geom_tile(aes(fill = Corr_NonSig), data = plot_data[!is.na(plot_data$Corr_NonSig), ]) +
    scale_fill_gradient2(low = "cyan", mid = "black", high = "magenta", 
                         midpoint = 0, limits = c(-1, 1),
                         name = "Non-significant\nCorrelation", na.value = "transparent") +
    # Add second fill scale for significant correlations
    ggnewscale::new_scale_fill() +
    geom_tile(aes(fill = Corr_Sig), data = plot_data[!is.na(plot_data$Corr_Sig), ]) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1),
                         name = "Significant\nCorrelation", na.value = "transparent") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right") +
    labs(title = paste(title, "\n(Blue-Red: p <", alpha, "; Cyan-Magenta: p ≥", alpha, ")"), 
         x = "", y = "")
  
  # Add text labels for small matrices
  if (nrow(ccf_matrix) <= 25) {
    p <- p + 
      geom_text(aes(label = sprintf("%.2f", Correlation)), 
                color = "white", size = 2.5, vjust = -0.2, fontface = "bold") +
      geom_text(aes(label = ifelse(Significant & P_Value < 0.001, "***",
                                   ifelse(Significant & P_Value < 0.01, "**",
                                          ifelse(Significant & P_Value < 0.05, "*", "")))),
                color = "white", size = 3, vjust = 1.2, fontface = "bold")
  }
  
  return(p)
}

# Alternative version without ggnewscale dependency
create_heatmap_alternative <- function(ccf_data, title, alpha = 0.05) {
  ccf_matrix <- ccf_data$correlations
  pval_matrix <- ccf_data$p_values
  
  # Prepare data
  plot_data <- data.frame(
    Series1 = rep(rownames(ccf_matrix), each = ncol(ccf_matrix)),
    Series2 = rep(colnames(ccf_matrix), nrow(ccf_matrix)),
    Correlation = as.vector(ccf_matrix),
    P_Value = as.vector(pval_matrix),
    Significant = as.vector(pval_matrix) < alpha
  )
  
  # Create custom color mapping
  plot_data$Color_Category <- ifelse(plot_data$Significant, "Significant", "Non-significant")
  plot_data$Color_Value <- ifelse(plot_data$Significant, 
                                  plot_data$Correlation + 2,  # Shift to 1-3 range
                                  plot_data$Correlation)       # Keep in -1 to 1 range
  
  # Define color breaks and labels
  breaks <- c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
  colors <- c("cyan", "cyan4", "black", "magenta4", "magenta",  # Non-significant
              "blue", "lightblue", "white", "pink", "red")      # Significant (shifted)
  
  p <- ggplot(plot_data, aes(x = Series2, y = Series1, fill = Color_Value)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors, 
                         breaks = breaks,
                         labels = c("-1", "-0.5", "0", "0.5", "1", "-1", "-0.5", "0", "0.5", "1"),
                         limits = c(-1, 3),
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste(title, "\n(Blue-Red: p <", alpha, "; Cyan-Magenta: p ≥", alpha, ")"), 
         x = "", y = "")
  
  # Add text labels for small matrices
  if (nrow(ccf_matrix) <= 25) {
    p <- p + 
      geom_text(aes(label = sprintf("%.2f", Correlation)), 
                color = "white", size = 2.5, vjust = -0.2, fontface = "bold") +
      geom_text(aes(label = ifelse(Significant & P_Value < 0.001, "***",
                                   ifelse(Significant & P_Value < 0.01, "**",
                                          ifelse(Significant & P_Value < 0.05, "*", "")))),
                color = "white", size = 3, vjust = 1.2, fontface = "bold")
  }
  
  return(p)
}

# Save results
save_results <- function(ccf_results, output_dir = "results") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (stage_name in names(ccf_results)) {
    ccf_data <- ccf_results[[stage_name]]
    
    # Save matrices
    write.csv(ccf_data$correlations, 
              file.path(output_dir, paste0("correlations_", stage_name, ".csv")))
    write.csv(ccf_data$p_values, 
              file.path(output_dir, paste0("pvalues_", stage_name, ".csv")))
    write.csv(ccf_data$p_values < 0.05, 
              file.path(output_dir, paste0("significant_", stage_name, ".csv")))
    
    # Try to create heatmap with ggnewscale, fall back to alternative if not available
    tryCatch({
      library(ggnewscale)
      heatmap <- create_heatmap(ccf_data, paste("Cross-Correlation:", stage_name))
    }, error = function(e) {
      cat("ggnewscale not available, using alternative heatmap function\n")
      heatmap <- create_heatmap_alternative(ccf_data, paste("Cross-Correlation:", stage_name))
    })
    
    ggsave(file.path(output_dir, paste0("heatmap_", stage_name, ".png")), 
           heatmap, width = 12, height = 8, dpi = 300)
    
    cat("Saved results for", stage_name, "\n")
  }
}

# Print summary statistics
print_summary <- function(ccf_results, alpha = 0.05) {
  cat("\nSignificant correlations summary (p <", alpha, "):\n")
  
  for (stage_name in names(ccf_results)) {
    pval_matrix <- ccf_results[[stage_name]]$p_values
    corr_matrix <- ccf_results[[stage_name]]$correlations
    
    diag(pval_matrix) <- NA  # Exclude diagonal
    significant_count <- sum(pval_matrix < alpha, na.rm = TRUE)
    total_count <- sum(!is.na(pval_matrix))
    
    cat(stage_name, ":", significant_count, "/", total_count, 
        "(", round(100 * significant_count / total_count, 1), "%)\n")
    
    # Show top correlations
    if (significant_count > 0) {
      sig_mask <- pval_matrix < alpha & !is.na(pval_matrix)
      sig_corrs <- corr_matrix[sig_mask]
      sig_pvals <- pval_matrix[sig_mask]
      
      top_indices <- order(abs(sig_corrs), decreasing = TRUE)[1:min(3, length(sig_corrs))]
      sig_positions <- which(sig_mask, arr.ind = TRUE)
      
      cat("  Top correlations:\n")
      for (i in top_indices) {
        row_idx <- sig_positions[i, 1]
        col_idx <- sig_positions[i, 2]
        cat("   ", rownames(corr_matrix)[row_idx], "vs", 
            colnames(corr_matrix)[col_idx], ": r =", 
            round(sig_corrs[i], 3), ", p =", 
            format(sig_pvals[i], scientific = TRUE, digits = 3), "\n")
      }
    }
  }
}

# Main analysis function
run_analysis <- function(lag = 0, max_lag = 5, alpha = 0.05, n_cores = parallel::detectCores() - 1) {
  cat("=== Sleep Stage CCF Analysis ===\n")
  
  # Load and process data
  data <- load_data()
  stage_traces <- segment_by_stage(data$stages, data$traces, data$time_col)
  
  if (length(stage_traces) == 0) {
    stop("No stage traces generated!")
  }
  
  # Compute correlations
  cat("\nComputing cross-correlations...\n")
  ccf_results <- compute_ccf_by_stage(stage_traces, lag, max_lag, n_cores)
  
  if (length(ccf_results) == 0) {
    stop("No CCF results generated!")
  }
  
  # Save results and create visualizations
  save_results(ccf_results)
  print_summary(ccf_results, alpha)
  
  cat("\n=== Analysis Complete ===\n")
  list(stage_traces = stage_traces, ccf_results = ccf_results)
}

# Execute analysis
if (file.exists("mPFCf5_BL_states_df.csv") && file.exists("mPFCf5_BL_Traces.csv")) {
  results <- run_analysis(lag = 0, max_lag = 5, alpha = 0.05)
  
  if (dir.exists("results")) {
    files <- list.files("results")
    cat("\nResults saved:", length(files), "files\n")
    cat(paste("  -", files, collapse = "\n"), "\n")
  }
} else {
  cat("ERROR: Required input files not found!\n")
  cat("Need: mPFCf5_BL_states_df.csv, mPFCf5_BL_Traces.csv\n")
}