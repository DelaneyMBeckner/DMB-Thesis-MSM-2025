# Epoch_CCF_Heatmaps.R
# Generate cross-correlation matrix heatmaps for specific epochs and ROIs
# Standalone script leveraging patterns from pipeline v3

library(ggplot2)
library(reshape2)
library(viridis)

# ==================== CONFIGURATION ====================
INPUT_DIR <- "E:/Data_Processing/R/Data CSVs"
OUTPUT_DIR <- "E:/Data_Processing/R/Results"
SAVE_OUTPUT <- TRUE

# Recording specification
TRACES_FILE <- "mPFCf5_BL_Traces.csv"
RECORDING_ID <- "mPFCf5_BL"

# Epoch and ROI selection
EPOCHS_TO_ANALYZE <- c(22)
ROIS_TO_INCLUDE <- c("C60", "C02", "C33")

# Timing parameters (consistent with pipeline v3)
EPOCH_DURATION <- 10  # seconds
SAMPLING_RATE <- 10   # Hz (0.1 sec per sample)
# =======================================================


#' Load traces file (adapted from Transition_CCF_v3)
load_traces <- function(traces_file, input_dir = INPUT_DIR) {
  traces_path <- file.path(input_dir, traces_file)
  cat("Loading traces from:", traces_path, "\n")
  
  traces_df <- read.csv(traces_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Fix first column name
  if (colnames(traces_df)[1] == "" || grepl("Time", colnames(traces_df)[1], ignore.case = TRUE)) {
    colnames(traces_df)[1] <- "Time"
  }
  
  # Remove header row if present
  first_val <- suppressWarnings(as.numeric(traces_df[1, 1]))
  if (is.na(first_val)) {
    cat("  Removing header row\n")
    traces_df <- traces_df[-1, ]
  }
  
  # Convert to numeric
  traces_df[] <- lapply(traces_df, function(x) as.numeric(as.character(x)))
  
  cat("Loaded:", nrow(traces_df), "timepoints,", ncol(traces_df) - 1, "ROIs\n")
  return(traces_df)
}


#' Extract traces for a specific epoch
extract_epoch <- function(traces_df, epoch_num, epoch_duration = EPOCH_DURATION) {
  # Epoch numbering starts at 1
  # Epoch 1: time 0 to 10 sec
  # Epoch N: time (N-1)*10 to N*10 sec
  start_time <- (epoch_num - 1) * epoch_duration
  end_time <- epoch_num * epoch_duration
  
  epoch_data <- traces_df[traces_df$Time >= start_time & traces_df$Time < end_time, ]
  
  cat("  Epoch", epoch_num, ": time", start_time, "-", end_time, "sec,", 
      nrow(epoch_data), "samples\n")
  
  return(epoch_data)
}


#' Compute pairwise cross-correlation matrix (from CCF3.R)
compute_pairwise_ccf <- function(data_matrix, lag = 0, max_lag = 10) {
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix)
  
  ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  rownames(ccf_matrix) <- series_names
  colnames(ccf_matrix) <- series_names
  
  for (i in 1:n_series) {
    ccf_matrix[i, i] <- 1
    
    if (i < n_series) {
      for (j in (i+1):n_series) {
        series_i <- as.numeric(data_matrix[, i])
        series_j <- as.numeric(data_matrix[, j])
        
        ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
        lag_index <- which(ccf_result$lag == lag)
        
        if (length(lag_index) > 0) {
          ccf_matrix[i, j] <- ccf_result$acf[lag_index]
          ccf_matrix[j, i] <- ccf_result$acf[lag_index]
        }
      }
    }
  }
  
  return(ccf_matrix)
}


#' Plot CCF heatmap for a single epoch
plot_ccf_heatmap <- function(ccf_matrix, epoch_num, recording_id) {
  # Convert to long format
  ccf_long <- melt(ccf_matrix)
  colnames(ccf_long) <- c("ROI_1", "ROI_2", "Correlation")
  
  # Create heatmap
  p <- ggplot(ccf_long, aes(x = ROI_2, y = ROI_1, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 5, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.75, 0.75), oob = scales::squish) +
    labs(title = paste0(recording_id, " - Epoch ", epoch_num),
         subtitle = "Cross-correlation at lag 0",
         x = "", y = "") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      aspect.ratio = 1
    ) +
    coord_fixed()
  
  return(p)
}


# ==================== MAIN EXECUTION ====================

cat("\n========================================\n")
cat("Epoch Cross-Correlation Heatmaps\n")
cat("Recording:", RECORDING_ID, "\n")
cat("Epochs:", paste(EPOCHS_TO_ANALYZE, collapse = ", "), "\n")
cat("ROIs:", paste(ROIS_TO_INCLUDE, collapse = ", "), "\n")
cat("========================================\n\n")

# Load traces
traces_full <- load_traces(TRACES_FILE, INPUT_DIR)

# Check that requested ROIs exist
available_rois <- colnames(traces_full)[colnames(traces_full) != "Time"]
missing_rois <- setdiff(ROIS_TO_INCLUDE, available_rois)

if (length(missing_rois) > 0) {
  cat("WARNING: ROIs not found in data:", paste(missing_rois, collapse = ", "), "\n")
  cat("Available ROIs:", paste(head(available_rois, 20), collapse = ", "), "...\n")
  ROIS_TO_INCLUDE <- intersect(ROIS_TO_INCLUDE, available_rois)
}

if (length(ROIS_TO_INCLUDE) < 2) {
  stop("Need at least 2 ROIs for correlation analysis")
}

# Subset to requested ROIs
traces_subset <- traces_full[, c("Time", ROIS_TO_INCLUDE)]
cat("\nAnalyzing", length(ROIS_TO_INCLUDE), "ROIs:", paste(ROIS_TO_INCLUDE, collapse = ", "), "\n\n")

# Initialize storage for plots and matrices
heatmaps <- list()
ccf_matrices <- list()

# Process each epoch
for (epoch in EPOCHS_TO_ANALYZE) {
  cat("Processing epoch", epoch, "...\n")
  
  # Extract epoch data
  epoch_data <- extract_epoch(traces_subset, epoch)
  
  if (nrow(epoch_data) < 10) {
    cat("  WARNING: Insufficient data for epoch", epoch, "- skipping\n")
    next
  }
  
  # Compute CCF matrix (exclude Time column)
  data_matrix <- epoch_data[, ROIS_TO_INCLUDE, drop = FALSE]
  ccf_matrix <- compute_pairwise_ccf(data_matrix)
  
  # Store matrix
  ccf_matrices[[paste0("Epoch_", epoch)]] <- ccf_matrix
  
  # Create heatmap
  heatmaps[[paste0("Epoch_", epoch)]] <- plot_ccf_heatmap(ccf_matrix, epoch, RECORDING_ID)
  
  cat("  CCF matrix:\n")
  print(round(ccf_matrix, 3))
  cat("\n")
}

# Display plots
cat("\n========================================\n")
cat("Displaying heatmaps...\n")
cat("========================================\n")

for (epoch_name in names(heatmaps)) {
  print(heatmaps[[epoch_name]])
  readline(prompt = paste0("Press [Enter] for next epoch (", epoch_name, " shown)..."))
}

# Save outputs if enabled
if (SAVE_OUTPUT && length(heatmaps) > 0) {
  cat("\nSaving outputs...\n")
  
  # Combined plot
  if (length(heatmaps) == 4) {
    library(gridExtra)
    combined_plot <- grid.arrange(
      heatmaps[[1]], heatmaps[[2]], 
      heatmaps[[3]], heatmaps[[4]], 
      ncol = 2,
      top = paste0(RECORDING_ID, " - Epochs ", 
                   paste(EPOCHS_TO_ANALYZE, collapse = ", "),
                   " - ROIs: ", paste(ROIS_TO_INCLUDE, collapse = ", "))
    )
    
    combined_file <- file.path(OUTPUT_DIR, 
                               paste0(RECORDING_ID, "_CCF_Epochs",
                                      min(EPOCHS_TO_ANALYZE), "-", max(EPOCHS_TO_ANALYZE),
                                      "_", paste(ROIS_TO_INCLUDE, collapse = "-"), ".png"))
    ggsave(combined_file, combined_plot, width = 10, height = 10, dpi = 150)
    cat("Saved combined plot:", combined_file, "\n")
  }
  
  # Individual epoch plots
  for (i in seq_along(heatmaps)) {
    epoch <- EPOCHS_TO_ANALYZE[i]
    epoch_file <- file.path(OUTPUT_DIR,
                            paste0(RECORDING_ID, "_CCF_Epoch", epoch,
                                   "_", paste(ROIS_TO_INCLUDE, collapse = "-"), ".png"))
    ggsave(epoch_file, heatmaps[[i]], width = 6, height = 5, dpi = 150)
    cat("Saved:", epoch_file, "\n")
  }
  
  # Save CCF matrices as CSV
  for (epoch_name in names(ccf_matrices)) {
    matrix_file <- file.path(OUTPUT_DIR,
                             paste0(RECORDING_ID, "_", epoch_name, 
                                    "_CCF_", paste(ROIS_TO_INCLUDE, collapse = "-"), ".csv"))
    write.csv(ccf_matrices[[epoch_name]], matrix_file)
    cat("Saved matrix:", matrix_file, "\n")
  }
}

cat("\n========================================\n")
cat("Analysis complete.\n")
cat("========================================\n")
