# Transition_DistCorr_v3.R
# Analyzes relationship between spatial distance and functional correlation
# during state transitions. Uses correlation matrices from Transition_CCF.R
# and pairwise distances between ROI centroids.
# v3: Adds file_suffix parameter for filtered analysis labeling

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(colorspace)
library(ggnewscale)

# ==================== CONFIGURATION ====================
SAVE_OUTPUT <- FALSE  # Set to TRUE to save results
INPUT_DIR <- "E:/Data_Processing/R/Data CSVs"  # Where distance files are
CCF_DIR <- "E:/Data_Processing/R/Results"      # Where CCF matrices are (from Transition_CCF)
METADATA_DIR <- "E:/Data_Processing/R/Results" # Where metadata files are
OUTPUT_DIR <- "E:/Data_Processing/R/Results"   # Where to save distance-correlation results

# Analysis mode
RUN_TRAJECTORY_ANALYSIS <- TRUE  # Set to TRUE for epoch-by-epoch trajectory analysis
RUN_WHOLE_WINDOW <- FALSE  # Set to TRUE for whole-window analysis (original mode)

# Analysis parameters
BINNING_METHOD <- "both"  # "equal_width", "percentile", or "both"
N_DISTANCE_BINS <- 10  # Number of equal-width bins (for equal_width method)
N_PERCENTILE_BINS <- 3  # Number of percentile bins (for percentile method, whole-window)
N_TRAJECTORY_BINS <- 3  # Number of percentile bins for trajectory analysis
FIT_EXPONENTIAL <- FALSE  # Set to TRUE to fit exponential decay (slower)

# File labeling (v3)
FILE_SUFFIX <- ""  # Set to "_filtered" when using activity-filtered ROIs

# Input files
RECORDING_ID <- "mPFCf5_BL"
TRANSITION_TYPE <- "Wake2NREM"  # or "NREM2Wake"
WINDOW_SIZE <- 3  # or 9
DISTANCES_FILE <- "mPFCf5_roi_distances.csv"  # Pairwise ROI distances
# =======================================================


#' Load pairwise distances between ROIs
load_distances <- function(distances_file, input_dir = INPUT_DIR) {
  distances_path <- file.path(input_dir, distances_file)
  cat("Loading pairwise distances from:", distances_path, "\n")
  
  distances_df <- read.csv(distances_path, stringsAsFactors = FALSE)
  
  # Verify required columns
  required_cols <- c("ROI_1", "ROI_2", "Distance")
  if (!all(required_cols %in% colnames(distances_df))) {
    stop("Distance file must have columns: ROI_1, ROI_2, Distance")
  }
  
  cat("Loaded", nrow(distances_df), "distance pairs\n")
  
  return(distances_df)
}


#' Load metadata for transitions
load_metadata <- function(recording_id, transition_type, window_size, 
                         metadata_dir = METADATA_DIR, metadata_df = NULL) {
  
  # If data frame provided, use it
  if (!is.null(metadata_df)) {
    cat("Using provided metadata data frame\n")
    return(metadata_df)
  }
  
  # Otherwise load from file
  metadata_file <- file.path(metadata_dir,
                            paste0(recording_id, "_", transition_type, "_",
                                  window_size, "ep_metadata.csv"))
  
  cat("Loading metadata from:", metadata_file, "\n")
  metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat("Loaded metadata for", nrow(metadata_df), "transitions\n")
  
  return(metadata_df)
}


#' Load CCF correlation matrix for a single transition
load_transition_ccf <- function(recording_id, transition_type, window_size, 
                               transition_id, ccf_dir = CCF_DIR) {
  
  # Construct filename following naming convention
  ccf_file <- file.path(ccf_dir,
                       paste0(recording_id, "_", transition_type, "_",
                             window_size, "ep_transition_",
                             sprintf("%03d", transition_id), "_ccf_lag0.csv"))
  
  if (!file.exists(ccf_file)) {
    warning("CCF file not found: ", ccf_file)
    return(NULL)
  }
  
  # Load correlation matrix
  ccf_matrix <- read.csv(ccf_file, row.names = 1, check.names = FALSE, 
                        stringsAsFactors = FALSE)
  
  return(as.matrix(ccf_matrix))
}


#' Convert correlation matrix to long format (pairwise)
matrix_to_pairs <- function(corr_matrix, include_diagonal = FALSE) {
  
  n_rois <- nrow(corr_matrix)
  roi_names <- rownames(corr_matrix)
  
  # Get indices for upper triangle (avoid duplicates)
  if (include_diagonal) {
    indices <- which(upper.tri(corr_matrix, diag = TRUE), arr.ind = TRUE)
  } else {
    indices <- which(upper.tri(corr_matrix, diag = FALSE), arr.ind = TRUE)
  }
  
  # Create pairs data frame
  pairs_df <- data.frame(
    ROI_1 = roi_names[indices[, 1]],
    ROI_2 = roi_names[indices[, 2]],
    Correlation = corr_matrix[indices],
    stringsAsFactors = FALSE
  )
  
  return(pairs_df)
}


#' Create order-independent pair key for merging
create_pair_key <- function(roi1, roi2) {
  paste(sort(c(roi1, roi2)), collapse = "_")
}


#' Merge correlation pairs with distance data
merge_distance_correlation <- function(corr_pairs, distances_df) {
  
  # Create pair keys for both datasets
  corr_pairs$pair_key <- mapply(create_pair_key, 
                                corr_pairs$ROI_1, corr_pairs$ROI_2)
  distances_df$pair_key <- mapply(create_pair_key, 
                                  distances_df$ROI_1, distances_df$ROI_2)
  
  # Merge on pair_key
  merged_df <- merge(distances_df, corr_pairs, by = "pair_key", 
                    suffixes = c("_dist", "_corr"))
  
  # Keep only essential columns
  result <- data.frame(
    ROI_1 = merged_df$ROI_1_dist,
    ROI_2 = merged_df$ROI_2_dist,
    Distance = merged_df$Distance,
    Correlation = merged_df$Correlation,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


#' Merge trajectory data with distances
merge_trajectory_distances <- function(trajectory_data, distances_df) {
  cat("\nMerging trajectory data with distances...\n")
  
  # Create pair keys for both datasets
  trajectory_data$pair_key <- mapply(create_pair_key, 
                                     trajectory_data$ROI_1, trajectory_data$ROI_2)
  distances_df$pair_key <- mapply(create_pair_key, 
                                  distances_df$ROI_1, distances_df$ROI_2)
  
  # Merge on pair_key
  merged_df <- merge(trajectory_data, distances_df, by = "pair_key")
  
  # Keep essential columns and clean up
  result <- merged_df %>%
    select(transition_id, epoch_in_window, ROI_1 = ROI_1.x, ROI_2 = ROI_2.x,
           correlation, Distance, transition_type, window_size) %>%
    filter(!is.na(Distance), !is.na(correlation))
  
  cat("Merged data: ", nrow(result), "observations\n")
  cat("  Epochs:", length(unique(result$epoch_in_window)), "\n")
  cat("  Transitions:", length(unique(result$transition_id)), "\n")
  
  return(result)
}


#' Bin distances using percentile method and compute trajectory summaries
bin_trajectory_by_distance <- function(trajectory_dist_data, n_bins = 3) {
  cat("\nBinning trajectory data by distance (", n_bins, "percentile bins)...\n", sep = "")
  
  # Compute percentile breaks on unique distances
  unique_distances <- unique(trajectory_dist_data$Distance)
  percentile_breaks <- quantile(unique_distances, 
                                probs = seq(0, 1, length.out = n_bins + 1),
                                na.rm = TRUE)
  percentile_breaks <- unique(percentile_breaks)
  
  cat("Distance bin boundaries:\n")
  print(percentile_breaks)
  
  # Assign bins
  trajectory_dist_data$distance_bin <- cut(trajectory_dist_data$Distance,
                                          breaks = percentile_breaks,
                                          include.lowest = TRUE,
                                          labels = FALSE)  # Numeric first for sorting
  
  # Convert to descriptive factor labels
  if (n_bins == 3) {
    bin_labels <- c("Close", "Medium", "Far")
  } else if (n_bins == 5) {
    bin_labels <- c("Close", "Medium-Close", "Medium", "Medium-Far", "Far")
  } else {
    # Generate generic labels for other bin counts
    bin_labels <- paste("Bin", 1:n_bins)
  }
  
  trajectory_dist_data$distance_bin_label <- factor(trajectory_dist_data$distance_bin,
                                                    levels = 1:n_bins,
                                                    labels = bin_labels[1:n_bins])
  
  # Remove any NA bins
  trajectory_dist_data <- trajectory_dist_data %>%
    filter(!is.na(distance_bin))
  
  # Summary by bin and epoch (mean across all transitions)
  bin_epoch_summary <- trajectory_dist_data %>%
    group_by(distance_bin, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      median_correlation = median(correlation, na.rm = TRUE),
      n_pairs = n(),
      min_distance = min(Distance),
      max_distance = max(Distance),
      mean_distance = mean(Distance),
      .groups = "drop"
    )
  
  # Summary by bin, epoch, AND transition (for individual trajectory lines)
  transition_summary <- trajectory_dist_data %>%
    group_by(distance_bin, distance_bin_label, epoch_in_window, transition_id) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )
  
  # Overall bin characteristics
  bin_characteristics <- trajectory_dist_data %>%
    group_by(distance_bin, distance_bin_label) %>%
    summarise(
      bin_min = min(Distance),
      bin_max = max(Distance),
      bin_center = mean(Distance),
      n_total_pairs = n(),
      .groups = "drop"
    )
  
  cat("Created summaries for", n_bins, "distance bins\n")
  
  return(list(
    binned_data = trajectory_dist_data,
    bin_epoch_summary = bin_epoch_summary,
    transition_summary = transition_summary,
    bin_characteristics = bin_characteristics,
    percentile_breaks = percentile_breaks
  ))
}


#' Create trajectory plots with distance binning
create_trajectory_distance_plots <- function(binned_results, recording_id, 
                                            transition_type, window_size) {
  cat("\nCreating trajectory-distance plots...\n")
  
  plots <- list()
  
  # Define color scheme for all plots
  library(colorspace)
  n_bins <- nlevels(binned_results$binned_data$distance_bin_label)
  base_hues <- seq(0, 300, length.out = n_bins + 1)[1:n_bins]  # Spread across color wheel
  
  # Create colors for mean lines (full saturation)
  bin_colors_mean <- hcl(h = base_hues, c = 80, l = 65)
  names(bin_colors_mean) <- levels(binned_results$binned_data$distance_bin_label)
  
  # 1. Distance distribution with bin boundaries
  distance_data <- binned_results$binned_data %>%
    select(Distance) %>%
    distinct()
  
  plots$distance_distribution <- ggplot(distance_data, aes(x = Distance)) +
    geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
    geom_vline(xintercept = binned_results$percentile_breaks, 
               color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Distance Distribution"),
      x = "Distance (pixels)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12)
    )
  
  # 2. Spaghetti plot: Correlation trajectories binned by distance (simplified)
  # Show only distance bin means with variance ribbons
  
  plots$trajectory_spaghetti <- ggplot(binned_results$bin_epoch_summary) +
    # Transition marker
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    # Variance ribbons for mean trajectories (show ± SE)
    geom_ribbon(aes(x = epoch_in_window, 
                    ymin = mean_correlation - se_correlation,
                    ymax = mean_correlation + se_correlation,
                    fill = distance_bin_label),
                alpha = 0.2) +
    # Mean trajectories
    geom_line(aes(x = epoch_in_window, y = mean_correlation,
                  group = distance_bin_label, color = distance_bin_label),
             size = 1.5) +  # Reduced from 2.5
    geom_point(aes(x = epoch_in_window, y = mean_correlation,
                   color = distance_bin_label),
              size = 3) +  # Reduced from 4
    scale_color_manual(
      values = bin_colors_mean,
      name = "Distance Bin",
      breaks = names(bin_colors_mean)
    ) +
    scale_fill_manual(
      values = bin_colors_mean,
      guide = "none"
    ) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Correlation Trajectory by Distance"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    )

  
  # 3. Faceted view: Separate panel for each distance bin with discrete transition colors
  # Create color mapping for transitions (discrete palette)
  trans_ids <- sort(unique(binned_results$transition_summary$transition_id))
  n_trans_total <- length(trans_ids)
  
  # Use discrete colors from viridis palette
  library(viridis)
  transition_colors_discrete <- viridis(n_trans_total, option = "plasma")
  names(transition_colors_discrete) <- as.character(trans_ids)
  
  # Merge with transition summary
  transition_summary_colored <- binned_results$transition_summary %>%
    mutate(transition_id_factor = factor(transition_id))
  
  plots$trajectory_faceted <- ggplot() +
    # Transition marker
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    # Individual transitions (discrete colors with legend)
    geom_line(data = transition_summary_colored,
             aes(x = epoch_in_window, y = mean_correlation,
                 group = transition_id, color = transition_id_factor),
             size = 1.0, alpha = 0.8) +
    scale_color_manual(
      values = transition_colors_discrete,
      name = "Transition ID",
      labels = trans_ids
    ) +
    # Mean trajectories (gray, bold) - add as second layer
    ggnewscale::new_scale_color() +
    geom_line(data = binned_results$bin_epoch_summary,
             aes(x = epoch_in_window, y = mean_correlation),
             color = "gray30", size = 1.8) +
    geom_point(data = binned_results$bin_epoch_summary,
              aes(x = epoch_in_window, y = mean_correlation),
              color = "gray30", size = 3) +
    facet_wrap(~distance_bin_label, ncol = 2) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Correlation by Distance Bin"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right"
    )
  
  # 4. Distance bin characteristics
  bin_char_plot <- binned_results$bin_characteristics %>%
    mutate(distance_range = paste0(round(bin_min), "-", round(bin_max), " px"))
  
  plots$bin_characteristics <- ggplot(bin_char_plot,
                                     aes(x = distance_bin_label, y = bin_center,
                                         fill = distance_bin_label)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_errorbar(aes(ymin = bin_min, ymax = bin_max), width = 0.3) +
    scale_fill_manual(values = bin_colors_mean, guide = "none") +
    geom_text(aes(label = distance_range), vjust = -0.5, size = 3) +
    labs(
      title = paste(recording_id, "- Distance Bin Characteristics"),
      x = "Distance Bin",
      y = "Distance (pixels)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  cat("Created", length(plots), "trajectory-distance plots\n")
  
  return(plots)
}


#' Compute binned averages of correlation vs distance
compute_binned_averages <- function(dist_corr_df, binning_method = "both",
                                   n_equal_width_bins = N_DISTANCE_BINS,
                                   n_percentile_bins = N_PERCENTILE_BINS,
                                   equal_width_breaks = NULL,
                                   percentile_breaks = NULL) {
  
  # Remove NA values
  clean_df <- dist_corr_df[complete.cases(dist_corr_df), ]
  
  if (nrow(clean_df) == 0) {
    warning("No valid distance-correlation pairs found")
    return(NULL)
  }
  
  results <- list()
  
  # Equal-width binning
  if (binning_method %in% c("equal_width", "both")) {
    # Use provided breaks or compute them (for backward compatibility)
    if (is.null(equal_width_breaks)) {
      equal_width_breaks <- seq(min(clean_df$Distance), max(clean_df$Distance), 
                                length.out = n_equal_width_bins + 1)
    }
    
    clean_df$distance_bin_ew <- cut(clean_df$Distance, 
                                    breaks = equal_width_breaks, 
                                    include.lowest = TRUE,
                                    labels = FALSE)  # Use numeric labels
    
    binned_ew <- clean_df %>%
      group_by(distance_bin_ew) %>%
      summarize(
        bin_id = first(distance_bin_ew),  # Use bin number as ID
        bin_min = min(Distance),
        bin_max = max(Distance),
        bin_center = mean(Distance),
        mean_correlation = mean(Correlation),
        sd_correlation = sd(Correlation),
        se_correlation = sd(Correlation) / sqrt(n()),
        n_pairs = n(),
        .groups = 'drop'
      ) %>%
      mutate(binning_method = "equal_width") %>%
      select(-distance_bin_ew)  # Remove the temporary column
    
    results$equal_width <- binned_ew
  }
  
  # Percentile-based binning (equal-count)
  if (binning_method %in% c("percentile", "both")) {
    # Use provided breaks or compute them (for backward compatibility)
    if (is.null(percentile_breaks)) {
      percentile_breaks <- quantile(clean_df$Distance, 
                                    probs = seq(0, 1, length.out = n_percentile_bins + 1),
                                    na.rm = TRUE)
      percentile_breaks <- unique(percentile_breaks)
    }
    
    if (length(percentile_breaks) < 2) {
      warning("Insufficient unique distance values for percentile binning")
    } else {
      clean_df$distance_bin_pct <- cut(clean_df$Distance,
                                       breaks = percentile_breaks,
                                       include.lowest = TRUE,
                                       labels = FALSE)  # Use numeric labels
      
      binned_pct <- clean_df %>%
        group_by(distance_bin_pct) %>%
        summarize(
          bin_id = first(distance_bin_pct),  # Use bin number as ID
          bin_min = min(Distance),
          bin_max = max(Distance),
          bin_center = mean(Distance),
          mean_correlation = mean(Correlation),
          sd_correlation = sd(Correlation),
          se_correlation = sd(Correlation) / sqrt(n()),
          n_pairs = n(),
          .groups = 'drop'
        ) %>%
        mutate(binning_method = "percentile") %>%
        select(-distance_bin_pct)  # Remove the temporary column
      
      results$percentile <- binned_pct
    }
  }
  
  # Return based on method
  if (binning_method == "equal_width") {
    return(results$equal_width)
  } else if (binning_method == "percentile") {
    return(results$percentile)
  } else {
    return(results)  # Return list with both
  }
}


#' Fit linear model to distance-correlation relationship
fit_linear_model <- function(dist_corr_df) {
  
  # Remove NA values
  clean_df <- dist_corr_df[complete.cases(dist_corr_df), ]
  
  if (nrow(clean_df) < 3) {
    warning("Insufficient data points for regression (n = ", nrow(clean_df), ")")
    return(NULL)
  }
  
  # Fit linear model
  lm_fit <- lm(Correlation ~ Distance, data = clean_df)
  lm_summary <- summary(lm_fit)
  
  # Extract key statistics
  result <- list(
    model = lm_fit,
    slope = coef(lm_fit)[2],
    intercept = coef(lm_fit)[1],
    r_squared = lm_summary$r.squared,
    p_value = coef(summary(lm_fit))[2, 4],
    n = nrow(clean_df)
  )
  
  return(result)
}


#' Fit exponential decay model (optional)
fit_exponential_model <- function(dist_corr_df) {
  
  # Remove NA values and ensure positive correlations for exponential fit
  clean_df <- dist_corr_df[complete.cases(dist_corr_df), ]
  
  if (nrow(clean_df) < 3) {
    warning("Insufficient data points for exponential fit")
    return(NULL)
  }
  
  # Exponential decay model: y = a * exp(-b * x) + c
  # Use nls for nonlinear least squares
  tryCatch({
    # Initial parameter guesses
    a_init <- max(clean_df$Correlation) - min(clean_df$Correlation)
    b_init <- 0.001
    c_init <- min(clean_df$Correlation)
    
    exp_fit <- nls(Correlation ~ a * exp(-b * Distance) + c,
                  data = clean_df,
                  start = list(a = a_init, b = b_init, c = c_init),
                  control = nls.control(maxiter = 100))
    
    # Calculate pseudo R-squared
    ss_res <- sum(residuals(exp_fit)^2)
    ss_tot <- sum((clean_df$Correlation - mean(clean_df$Correlation))^2)
    pseudo_r_squared <- 1 - (ss_res / ss_tot)
    
    result <- list(
      model = exp_fit,
      a = coef(exp_fit)[1],
      b = coef(exp_fit)[2],
      c = coef(exp_fit)[3],
      pseudo_r_squared = pseudo_r_squared,
      n = nrow(clean_df)
    )
    
    return(result)
    
  }, error = function(e) {
    warning("Exponential fit failed: ", e$message)
    return(NULL)
  })
}


#' Analyze distance-correlation for all transitions
analyze_all_transitions <- function(recording_id, transition_type, window_size,
                                   distances_file, 
                                   input_dir = INPUT_DIR,
                                   ccf_dir = CCF_DIR,
                                   metadata_dir = METADATA_DIR,
                                   distances_df = NULL,
                                   metadata_df = NULL,
                                   ccf_results_list = NULL,
                                   binning_method = "both",
                                   n_equal_width_bins = N_DISTANCE_BINS,
                                   n_percentile_bins = N_PERCENTILE_BINS,
                                   fit_exponential = FIT_EXPONENTIAL) {
  
  cat("\n========================================\n")
  cat("DISTANCE-CORRELATION ANALYSIS\n")
  cat("========================================\n")
  cat("Recording:", recording_id, "\n")
  cat("Transition:", transition_type, "\n")
  cat("Window:", window_size, "epochs\n")
  cat("Binning method:", binning_method, "\n")
  if (binning_method %in% c("equal_width", "both")) {
    cat("  Equal-width bins:", n_equal_width_bins, "\n")
  }
  if (binning_method %in% c("percentile", "both")) {
    cat("  Percentile bins:", n_percentile_bins, "\n")
  }
  cat("Exponential fit:", ifelse(fit_exponential, "ENABLED", "DISABLED"), "\n")
  cat("CCF data source:", ifelse(is.null(ccf_results_list), "FILES", "MEMORY"), "\n\n")
  
  # Load distances
  if (is.null(distances_df)) {
    distances_df <- load_distances(distances_file, input_dir)
  } else {
    cat("Using provided distances data frame\n")
  }
  
  # Load metadata
  metadata_df <- load_metadata(recording_id, transition_type, window_size, 
                              metadata_dir, metadata_df)
  
  # Initialize results storage
  all_dist_corr <- list()
  all_binned_ew <- list()
  all_binned_pct <- list()
  linear_fits <- list()
  exp_fits <- list()
  
  # COMPUTE GLOBAL DISTANCE BREAKS ONCE (not per-transition)
  cat("\nComputing global distance breaks...\n")
  distance_range <- range(distances_df$Distance)
  
  # Equal-width breaks: fixed distance intervals
  equal_width_breaks <- NULL
  if (binning_method %in% c("equal_width", "both")) {
    equal_width_breaks <- seq(distance_range[1], distance_range[2], 
                              length.out = n_equal_width_bins + 1)
    cat("  Equal-width bins:", n_equal_width_bins, "bins from", 
        round(distance_range[1], 1), "to", round(distance_range[2], 1), "pixels\n")
  }
  
  # Percentile breaks: based on global distance distribution
  percentile_breaks <- NULL
  if (binning_method %in% c("percentile", "both")) {
    percentile_breaks <- quantile(distances_df$Distance, 
                                  probs = seq(0, 1, length.out = n_percentile_bins + 1),
                                  na.rm = TRUE)
    percentile_breaks <- unique(percentile_breaks)
    cat("  Percentile bins:", length(percentile_breaks) - 1, "bins (equal-count)\n")
  }
  
  # Process each transition
  cat("\nProcessing transitions...\n")
  n_transitions <- nrow(metadata_df)
  
  for (i in 1:n_transitions) {
    trans_id <- metadata_df$transition_id[i]
    
    cat("  Transition", trans_id, "(", i, "/", n_transitions, ")...\n")
    
    # Load CCF matrix - check memory first, then file
    ccf_matrix <- NULL
    
    if (!is.null(ccf_results_list)) {
      # Try to get from memory (ccf_results_list from Transition_CCF.R)
      trans_id_str <- as.character(trans_id)
      if (trans_id_str %in% names(ccf_results_list)) {
        ccf_matrix <- ccf_results_list[[trans_id_str]]$ccf_matrix_lag0
        if (!is.null(ccf_matrix)) {
          cat("    Using CCF matrix from memory\n")
        }
      }
    }
    
    # Fall back to loading from file if not in memory
    if (is.null(ccf_matrix)) {
      ccf_matrix <- load_transition_ccf(recording_id, transition_type, window_size,
                                       trans_id, ccf_dir)
      if (!is.null(ccf_matrix)) {
        cat("    Loaded CCF matrix from file\n")
      }
    }
    
    if (is.null(ccf_matrix)) {
      cat("    Skipping - CCF matrix not found (memory or file)\n")
      next
    }
    
    # Convert to pairwise format
    corr_pairs <- matrix_to_pairs(ccf_matrix, include_diagonal = FALSE)
    
    # Merge with distances
    dist_corr <- merge_distance_correlation(corr_pairs, distances_df)
    
    if (nrow(dist_corr) == 0) {
      cat("    Skipping - no matching distance pairs\n")
      next
    }
    
    # Add transition identifier
    dist_corr$transition_id <- trans_id
    
    # Store results
    all_dist_corr[[as.character(trans_id)]] <- dist_corr
    
    # Compute binned averages using GLOBAL breaks (same for all transitions)
    binned_result <- compute_binned_averages(dist_corr, 
                                            binning_method = binning_method,
                                            n_equal_width_bins = n_equal_width_bins,
                                            n_percentile_bins = n_percentile_bins,
                                            equal_width_breaks = equal_width_breaks,
                                            percentile_breaks = percentile_breaks)
    
    if (!is.null(binned_result)) {
      # Handle different return structures based on method
      if (binning_method == "both") {
        if (!is.null(binned_result$equal_width)) {
          binned_result$equal_width$transition_id <- trans_id
          all_binned_ew[[as.character(trans_id)]] <- binned_result$equal_width
        }
        if (!is.null(binned_result$percentile)) {
          binned_result$percentile$transition_id <- trans_id
          all_binned_pct[[as.character(trans_id)]] <- binned_result$percentile
        }
      } else if (binning_method == "equal_width") {
        binned_result$transition_id <- trans_id
        all_binned_ew[[as.character(trans_id)]] <- binned_result
      } else if (binning_method == "percentile") {
        binned_result$transition_id <- trans_id
        all_binned_pct[[as.character(trans_id)]] <- binned_result
      }
    }
    
    # Fit linear model
    linear_fit <- fit_linear_model(dist_corr)
    if (!is.null(linear_fit)) {
      linear_fits[[as.character(trans_id)]] <- linear_fit
      cat("    Linear fit: slope =", 
          format(linear_fit$slope, scientific = TRUE, digits = 3),
          ", RÂ² =", round(linear_fit$r_squared, 4),
          ", p =", format(linear_fit$p_value, scientific = TRUE, digits = 2), "\n")
    }
    
    # Fit exponential model (if requested)
    if (fit_exponential) {
      exp_fit <- fit_exponential_model(dist_corr)
      if (!is.null(exp_fit)) {
        exp_fits[[as.character(trans_id)]] <- exp_fit
        cat("    Exponential fit: pseudo-RÂ² =", 
            round(exp_fit$pseudo_r_squared, 4), "\n")
      }
    }
  }
  
  cat("\nCompleted analysis for", length(all_dist_corr), "transitions\n")
  
  # Combine all data
  combined_dist_corr <- do.call(rbind, all_dist_corr)
  
  combined_binned_ew <- NULL
  combined_binned_pct <- NULL
  
  if (length(all_binned_ew) > 0) {
    combined_binned_ew <- do.call(rbind, all_binned_ew)
  }
  if (length(all_binned_pct) > 0) {
    combined_binned_pct <- do.call(rbind, all_binned_pct)
  }
  
  # Aggregate statistics across transitions
  cat("\nComputing aggregate statistics...\n")
  
  aggregate_stats <- list(
    mean_slope = mean(sapply(linear_fits, function(x) x$slope)),
    sd_slope = sd(sapply(linear_fits, function(x) x$slope)),
    mean_r_squared = mean(sapply(linear_fits, function(x) x$r_squared)),
    median_r_squared = median(sapply(linear_fits, function(x) x$r_squared)),
    n_transitions = length(linear_fits)
  )
  
  cat("  Mean slope:", format(aggregate_stats$mean_slope, scientific = TRUE, digits = 3), 
      "Â±", format(aggregate_stats$sd_slope, scientific = TRUE, digits = 2), "\n")
  cat("  Mean RÂ²:", round(aggregate_stats$mean_r_squared, 4), "\n")
  
  # Return results
  results <- list(
    dist_corr_all = combined_dist_corr,
    binned_equal_width = combined_binned_ew,
    binned_percentile = combined_binned_pct,
    linear_fits = linear_fits,
    aggregate_stats = aggregate_stats,
    metadata = metadata_df,
    binning_method = binning_method
  )
  
  if (fit_exponential) {
    results$exp_fits <- exp_fits
  }
  
  return(results)
}


#' Create visualizations for distance-correlation analysis
create_distcorr_plots <- function(results, recording_id, transition_type, 
                                 window_size, fit_exponential = FIT_EXPONENTIAL) {
  
  cat("\nCreating visualizations...\n")
  
  plots <- list()
  
  dist_corr_all <- results$dist_corr_all
  binned_ew <- results$binned_equal_width
  binned_pct <- results$binned_percentile
  linear_fits <- results$linear_fits
  aggregate_stats <- results$aggregate_stats
  binning_method <- results$binning_method
  
  # 1. Scatter plot: all transitions combined
  cat("  1. Combined scatter plot with regression\n")
  
  # Fit overall linear model
  overall_fit <- fit_linear_model(dist_corr_all)
  
  if (!is.null(overall_fit)) {
    plots$scatter_combined <- ggplot(dist_corr_all, 
                                     aes(x = Distance, y = Correlation)) +
      geom_point(alpha = 0.2, size = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
      annotate("text", 
               x = max(dist_corr_all$Distance) * 0.7,
               y = max(dist_corr_all$Correlation) * 0.9,
               label = paste0("R² = ", round(overall_fit$r_squared, 4), "\n",
                            "p = ", format(overall_fit$p_value, scientific = TRUE, digits = 2), "\n",
                            "slope = ", format(overall_fit$slope, scientific = TRUE, digits = 3), "\n",
                            "n = ", overall_fit$n, " pairs, ",
                            aggregate_stats$n_transitions, " transitions"),
               size = 3.5, hjust = 0) +
      labs(
        title = paste(recording_id, "-", transition_type, 
                     "- Distance vs Correlation (All Transitions)"),
        x = "Distance (pixels)",
        y = "Correlation (Pearson's r)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title = element_text(size = 11)
      )
  }
  
  # 2. Binned averages - Equal-width
  if (!is.null(binned_ew) && binning_method %in% c("equal_width", "both")) {
    cat("  2. Binned averages plot (equal-width)\n")
    
    # Aggregate bins across transitions using bin_id
    # Now all transitions should have all bins (same global breaks)
    binned_aggregate_ew <- binned_ew %>%
      group_by(bin_id) %>%
      summarize(
        bin_center = mean(bin_center),
        mean_correlation = mean(mean_correlation),
        sd_correlation = sd(mean_correlation),
        se_aggregate = sd(mean_correlation) / sqrt(n()),
        n_transitions = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        # Replace NA (from single-transition bins) with 0
        sd_correlation = ifelse(is.na(sd_correlation), 0, sd_correlation),
        se_aggregate = ifelse(is.na(se_aggregate), 0, se_aggregate)
      )
    
    plots$binned_averages_equal_width <- ggplot(binned_aggregate_ew, 
                                               aes(x = bin_center, y = mean_correlation)) +
      geom_ribbon(aes(ymin = mean_correlation - sd_correlation,
                     ymax = mean_correlation + sd_correlation),
                 fill = "steelblue", alpha = 0.2) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_point(size = 2.5, color = "steelblue") +
      geom_errorbar(aes(ymin = mean_correlation - se_aggregate,
                       ymax = mean_correlation + se_aggregate),
                   width = 50, color = "steelblue", linewidth = 0.8) +
      labs(
        title = paste(recording_id, "-", transition_type, 
                     "- Binned Distance vs Correlation (Equal-Width)"),
        x = "Distance (pixels)",
        y = "Mean Correlation"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title = element_text(size = 11)
      )
  }
  
  # 3. Binned averages - Percentile
  if (!is.null(binned_pct) && binning_method %in% c("percentile", "both")) {
    cat("  3. Binned averages plot (percentile)\n")
    
    # Aggregate bins across transitions using bin_id
    binned_aggregate_pct <- binned_pct %>%
      group_by(bin_id) %>%
      summarize(
        bin_center = mean(bin_center),
        mean_correlation = mean(mean_correlation),
        sd_correlation = sd(mean_correlation),
        se_aggregate = sd(mean_correlation) / sqrt(n()),
        n_transitions = n(),
        n_pairs_per_bin = mean(n_pairs),
        .groups = 'drop'
      ) %>%
      mutate(
        # Replace NA (from single-transition bins) with 0
        sd_correlation = ifelse(is.na(sd_correlation), 0, sd_correlation),
        se_aggregate = ifelse(is.na(se_aggregate), 0, se_aggregate)
      )
    
    plots$binned_averages_percentile <- ggplot(binned_aggregate_pct, 
                                              aes(x = bin_center, y = mean_correlation)) +
      geom_ribbon(aes(ymin = mean_correlation - sd_correlation,
                     ymax = mean_correlation + sd_correlation),
                 fill = "darkgreen", alpha = 0.2) +
      geom_line(color = "darkgreen", linewidth = 1) +
      geom_point(size = 2.5, color = "darkgreen") +
      geom_errorbar(aes(ymin = mean_correlation - se_aggregate,
                       ymax = mean_correlation + se_aggregate),
                   width = 50, color = "darkgreen", linewidth = 0.8) +
      labs(
        title = paste(recording_id, "-", transition_type, 
                     "- Binned Distance vs Correlation (Percentile)"),
        x = "Distance (pixels)",
        y = "Mean Correlation"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title = element_text(size = 11)
      )
  }
  
  # 4. Combined binning comparison (if both methods used)
  if (!is.null(binned_ew) && !is.null(binned_pct) && binning_method == "both") {
    cat("  4. Binning methods comparison plot\n")
    
    # Prepare data for comparison - aggregate by bin_id
    binned_aggregate_ew <- binned_ew %>%
      group_by(bin_id) %>%
      summarize(
        bin_center = mean(bin_center),
        mean_correlation = mean(mean_correlation),
        sd_correlation = sd(mean_correlation),
        se_aggregate = sd(mean_correlation) / sqrt(n()),
        .groups = 'drop'
      ) %>%
      mutate(
        # Replace NA (from single-transition bins) with 0
        sd_correlation = ifelse(is.na(sd_correlation), 0, sd_correlation),
        se_aggregate = ifelse(is.na(se_aggregate), 0, se_aggregate),
        method = "Equal-Width"
      )
    
    binned_aggregate_pct <- binned_pct %>%
      group_by(bin_id) %>%
      summarize(
        bin_center = mean(bin_center),
        mean_correlation = mean(mean_correlation),
        sd_correlation = sd(mean_correlation),
        se_aggregate = sd(mean_correlation) / sqrt(n()),
        .groups = 'drop'
      ) %>%
      mutate(
        # Replace NA (from single-transition bins) with 0
        sd_correlation = ifelse(is.na(sd_correlation), 0, sd_correlation),
        se_aggregate = ifelse(is.na(se_aggregate), 0, se_aggregate),
        method = "Percentile"
      )
    
    combined_binned <- rbind(binned_aggregate_ew, binned_aggregate_pct)
    
    plots$binning_comparison <- ggplot(combined_binned,
                                      aes(x = bin_center, y = mean_correlation,
                                          color = method, fill = method)) +
      geom_ribbon(aes(ymin = mean_correlation - sd_correlation,
                     ymax = mean_correlation + sd_correlation),
                 alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean_correlation - se_aggregate,
                       ymax = mean_correlation + se_aggregate),
                   width = 50, alpha = 0.7) +
      scale_color_manual(values = c("Equal-Width" = "steelblue", 
                                   "Percentile" = "darkgreen")) +
      scale_fill_manual(values = c("Equal-Width" = "steelblue",
                                  "Percentile" = "darkgreen")) +
      labs(
        title = paste(recording_id, "-", transition_type, 
                     "- Binning Methods Comparison"),
        x = "Distance (pixels)",
        y = "Mean Correlation",
        color = "Binning Method",
        fill = "Binning Method"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title = element_text(size = 11),
        legend.position = "bottom"
      )
  }
  
  # 5. Per-transition slopes distribution
  cat("  5. Slope distribution histogram\n")
  
  slopes_df <- data.frame(
    transition_id = as.numeric(names(linear_fits)),
    slope = sapply(linear_fits, function(x) x$slope),
    r_squared = sapply(linear_fits, function(x) x$r_squared)
  )
  
  plots$slope_distribution <- ggplot(slopes_df, aes(x = slope)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    geom_vline(xintercept = aggregate_stats$mean_slope, 
              color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text",
            x = aggregate_stats$mean_slope,
            y = Inf,
            label = paste("Mean =", format(aggregate_stats$mean_slope, 
                                          scientific = TRUE, digits = 3)),
            hjust = -0.1, vjust = 1.5, size = 3.5, color = "red") +
    labs(
      title = paste(recording_id, "-", transition_type, 
                   "- Distribution of Slopes Across Transitions"),
      x = "Slope (correlation / pixel)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot"
    )
  
  # 6. RÂ² distribution
  cat("  6. RÂ² distribution histogram\n")
  
  plots$r_squared_distribution <- ggplot(slopes_df, aes(x = r_squared)) +
    geom_histogram(bins = 20, fill = "darkgreen", color = "white") +
    geom_vline(xintercept = aggregate_stats$mean_r_squared,
              color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text",
            x = aggregate_stats$mean_r_squared,
            y = Inf,
            label = paste("Mean =", round(aggregate_stats$mean_r_squared, 4)),
            hjust = -0.1, vjust = 1.5, size = 3.5, color = "red") +
    labs(
      title = paste(recording_id, "-", transition_type, "- Distribution of R² Across Transitions"),
      x = "R² (goodness of fit)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot"
    )
  
  # 7. Color-coded scatter plot (all transitions on one plot)
  cat("  7. Color-coded scatter plot (all transitions)\n")
  
  # Limit to reasonable number of transitions for color coding
  n_to_show <- min(20, length(unique(dist_corr_all$transition_id)))
  transitions_to_show <- unique(dist_corr_all$transition_id)[1:n_to_show]
  
  dist_corr_subset <- dist_corr_all %>%
    filter(transition_id %in% transitions_to_show) %>%
    mutate(transition_id = factor(transition_id))
  
  plots$scatter_colored <- ggplot(dist_corr_subset, 
                                 aes(x = Distance, y = Correlation, color = transition_id)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = paste(recording_id, "-", transition_type, 
                   "- Distance vs Correlation (Color-Coded by Transition)"),
      x = "Distance (pixels)",
      y = "Correlation",
      color = "Transition ID"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      legend.position = "right"
    )
  
  # 8. Correlation decay curve (if exponential fit requested)
  if (fit_exponential && !is.null(results$exp_fits)) {
    cat("  8. Exponential decay comparison\n")
    
    # Get first successful exponential fit for demonstration
    exp_fit <- results$exp_fits[[1]]
    
    if (!is.null(exp_fit)) {
      # Generate prediction curve
      x_seq <- seq(min(dist_corr_all$Distance), 
                  max(dist_corr_all$Distance), 
                  length.out = 100)
      y_pred <- exp_fit$a * exp(-exp_fit$b * x_seq) + exp_fit$c
      
      pred_df <- data.frame(Distance = x_seq, Correlation = y_pred)
      
      plots$exponential_decay <- ggplot(dist_corr_all, 
                                       aes(x = Distance, y = Correlation)) +
        geom_point(alpha = 0.2, size = 0.5, color = "steelblue") +
        geom_line(data = pred_df, color = "red", linewidth = 1) +
        annotate("text",
                x = max(dist_corr_all$Distance) * 0.7,
                y = max(dist_corr_all$Correlation) * 0.9,
                label = paste0("y = ", round(exp_fit$a, 4), 
                             " × exp(-", format(exp_fit$b, scientific = TRUE, digits = 2),
                             " × x) + ", round(exp_fit$c, 4), "\n",
                             "Pseudo-R² = ", round(exp_fit$pseudo_r_squared, 4)),
                size = 3.5, hjust = 0) +
        labs(
          title = paste(recording_id, "-", transition_type, "- Exponential Decay Fit"),
          x = "Distance (pixels)",
          y = "Correlation"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.title.position = "plot"
        )
    }
  }
  
  cat("Created", length(plots), "plots\n")
  
  return(plots)
}


#' Save results to CSV files
save_distcorr_results <- function(results, recording_id, transition_type, 
                                 window_size, output_dir, file_suffix = "") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving results to:", output_dir, "\n")
  
  files_saved <- 0
  
  # 1. Save combined distance-correlation pairs
  filename_all <- file.path(output_dir,
                           paste0(recording_id, "_", transition_type, "_",
                                 window_size, "ep", file_suffix, "_distcorr_all_pairs.csv"))
  write.csv(results$dist_corr_all, filename_all, row.names = FALSE)
  cat("  Saved:", basename(filename_all), "\n")
  files_saved <- files_saved + 1
  
  # 2. Save binned averages (equal-width)
  if (!is.null(results$binned_equal_width)) {
    filename_binned_ew <- file.path(output_dir,
                                paste0(recording_id, "_", transition_type, "_",
                                      window_size, "ep", file_suffix, "_distcorr_binned_equal_width.csv"))
    write.csv(results$binned_equal_width, filename_binned_ew, row.names = FALSE)
    cat("  Saved:", basename(filename_binned_ew), "\n")
    files_saved <- files_saved + 1
  }
  
  # 3. Save binned averages (percentile)
  if (!is.null(results$binned_percentile)) {
    filename_binned_pct <- file.path(output_dir,
                                 paste0(recording_id, "_", transition_type, "_",
                                       window_size, "ep", file_suffix, "_distcorr_binned_percentile.csv"))
    write.csv(results$binned_percentile, filename_binned_pct, row.names = FALSE)
    cat("  Saved:", basename(filename_binned_pct), "\n")
    files_saved <- files_saved + 1
  }
  
  # 4. Save per-transition regression statistics
  linear_stats <- data.frame(
    transition_id = as.numeric(names(results$linear_fits)),
    slope = sapply(results$linear_fits, function(x) x$slope),
    intercept = sapply(results$linear_fits, function(x) x$intercept),
    r_squared = sapply(results$linear_fits, function(x) x$r_squared),
    p_value = sapply(results$linear_fits, function(x) x$p_value),
    n_pairs = sapply(results$linear_fits, function(x) x$n)
  )
  
  filename_stats <- file.path(output_dir,
                             paste0(recording_id, "_", transition_type, "_",
                                   window_size, "ep", file_suffix, "_distcorr_linear_stats.csv"))
  write.csv(linear_stats, filename_stats, row.names = FALSE)
  cat("  Saved:", basename(filename_stats), "\n")
  files_saved <- files_saved + 1
  
  # 5. Save aggregate statistics
  aggregate_df <- data.frame(
    metric = c("mean_slope", "sd_slope", "mean_r_squared", 
              "median_r_squared", "n_transitions"),
    value = c(results$aggregate_stats$mean_slope,
             results$aggregate_stats$sd_slope,
             results$aggregate_stats$mean_r_squared,
             results$aggregate_stats$median_r_squared,
             results$aggregate_stats$n_transitions)
  )
  
  filename_agg <- file.path(output_dir,
                           paste0(recording_id, "_", transition_type, "_",
                                 window_size, "ep", file_suffix, "_distcorr_aggregate_stats.csv"))
  write.csv(aggregate_df, filename_agg, row.names = FALSE)
  cat("  Saved:", basename(filename_agg), "\n")
  files_saved <- files_saved + 1
  
  # 6. Save exponential fit results (if available)
  if (!is.null(results$exp_fits)) {
    exp_stats <- data.frame(
      transition_id = as.numeric(names(results$exp_fits)),
      a = sapply(results$exp_fits, function(x) x$a),
      b = sapply(results$exp_fits, function(x) x$b),
      c = sapply(results$exp_fits, function(x) x$c),
      pseudo_r_squared = sapply(results$exp_fits, function(x) x$pseudo_r_squared),
      n_pairs = sapply(results$exp_fits, function(x) x$n)
    )
    
    filename_exp <- file.path(output_dir,
                             paste0(recording_id, "_", transition_type, "_",
                                   window_size, "ep", file_suffix, "_distcorr_exp_stats.csv"))
    write.csv(exp_stats, filename_exp, row.names = FALSE)
    cat("  Saved:", basename(filename_exp), "\n")
    files_saved <- files_saved + 1
  }
  
  cat("Total files saved:", files_saved, "\n")
}


#' Save plots to PNG files
save_distcorr_plots <- function(plots, recording_id, transition_type, 
                               window_size, output_dir, file_suffix = "") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving plots to:", output_dir, "\n")
  
  files_saved <- 0
  
  for (plot_name in names(plots)) {
    filename <- file.path(output_dir,
                         paste0(recording_id, "_", transition_type, "_",
                               window_size, "ep", file_suffix, "_distcorr_", plot_name, ".png"))
    
    ggsave(filename, plots[[plot_name]], width = 10, height = 8, dpi = 300)
    cat("  Saved:", basename(filename), "\n")
    files_saved <- files_saved + 1
  }
  
  cat("Total plots saved:", files_saved, "\n")
}


#' Main analysis function (pipeline-ready)
analyze_transition_distcorr <- function(recording_id = RECORDING_ID,
                                       transition_type = TRANSITION_TYPE,
                                       window_size = WINDOW_SIZE,
                                       distances_file = DISTANCES_FILE,
                                       input_dir = INPUT_DIR,
                                       ccf_dir = CCF_DIR,
                                       metadata_dir = METADATA_DIR,
                                       distances_df = NULL,
                                       metadata_df = NULL,
                                       ccf_results_list = NULL,
                                       trajectory_data = NULL,
                                       run_trajectory_analysis = RUN_TRAJECTORY_ANALYSIS,
                                       run_whole_window = RUN_WHOLE_WINDOW,
                                       n_trajectory_bins = N_TRAJECTORY_BINS,
                                       binning_method = BINNING_METHOD,
                                       n_equal_width_bins = N_DISTANCE_BINS,
                                       n_percentile_bins = N_PERCENTILE_BINS,
                                       fit_exponential = FIT_EXPONENTIAL,
                                       save_output = SAVE_OUTPUT,
                                       output_dir = OUTPUT_DIR,
                                       file_suffix = FILE_SUFFIX) {
  
  cat("\n")
  cat("====================================================\n")
  cat("  TRANSITION DISTANCE-CORRELATION ANALYSIS (v3)\n")
  cat("====================================================\n\n")
  
  # Initialize output
  output <- list()
  
  # Load distances (needed for both modes)
  if (is.null(distances_df)) {
    distances_df <- load_distances(distances_file, input_dir)
  } else {
    cat("Using provided distances data frame\n")
  }
  
  # Trajectory analysis mode
  if (run_trajectory_analysis) {
    cat("\n========== TRAJECTORY ANALYSIS ==========\n")
    
    if (is.null(trajectory_data)) {
      stop("Trajectory analysis requires trajectory_data from CCF analysis")
    }
    
    # Merge trajectory data with distances
    trajectory_dist <- merge_trajectory_distances(trajectory_data, distances_df)
    
    # Bin by distance and compute summaries
    binned_results <- bin_trajectory_by_distance(trajectory_dist, n_trajectory_bins)
    
    # Create plots
    trajectory_plots <- create_trajectory_distance_plots(binned_results, 
                                                         recording_id,
                                                         transition_type, 
                                                         window_size)
    
    # Save outputs if requested
    if (save_output) {
      cat("\nSaving trajectory analysis outputs...\n")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Save binned data
      write.csv(binned_results$binned_data,
               file.path(output_dir, paste0(recording_id, "_", transition_type, "_",
                                           window_size, "ep", file_suffix, "_distcorr_trajectory_binned.csv")),
               row.names = FALSE)
      
      # Save summaries
      write.csv(binned_results$bin_epoch_summary,
               file.path(output_dir, paste0(recording_id, "_", transition_type, "_",
                                           window_size, "ep", file_suffix, "_distcorr_trajectory_bin_epoch_summary.csv")),
               row.names = FALSE)
      
      write.csv(binned_results$bin_characteristics,
               file.path(output_dir, paste0(recording_id, "_", transition_type, "_",
                                           window_size, "ep", file_suffix, "_distcorr_trajectory_bin_characteristics.csv")),
               row.names = FALSE)
      
      # Save plots
      for (plot_name in names(trajectory_plots)) {
        filename <- file.path(output_dir,
                             paste0(recording_id, "_", transition_type, "_",
                                   window_size, "ep", file_suffix, "_distcorr_traj_", plot_name, ".png"))
        ggsave(filename, trajectory_plots[[plot_name]], width = 10, height = 8, dpi = 300)
        cat("  Saved:", basename(filename), "\n")
      }
    }
    
    output$trajectory_results <- binned_results
    output$trajectory_plots <- trajectory_plots
  }
  
  # Whole-window analysis mode (original)
  if (run_whole_window) {
    cat("\n========== WHOLE-WINDOW ANALYSIS ==========\n")
    
    # Analyze all transitions
    results <- analyze_all_transitions(
      recording_id = recording_id,
      transition_type = transition_type,
      window_size = window_size,
      distances_file = distances_file,
      input_dir = input_dir,
      ccf_dir = ccf_dir,
      metadata_dir = metadata_dir,
      distances_df = distances_df,
      metadata_df = metadata_df,
      ccf_results_list = ccf_results_list,
      binning_method = binning_method,
      n_equal_width_bins = n_equal_width_bins,
      n_percentile_bins = n_percentile_bins,
      fit_exponential = fit_exponential
    )
    
    # Create visualizations
    plots <- create_distcorr_plots(
      results = results,
      recording_id = recording_id,
      transition_type = transition_type,
      window_size = window_size,
      fit_exponential = fit_exponential
    )
    
    # Save outputs if requested
    if (save_output) {
      save_distcorr_results(results, recording_id, transition_type, 
                           window_size, output_dir, file_suffix)
      save_distcorr_plots(plots, recording_id, transition_type, 
                         window_size, output_dir, file_suffix)
    }
    
    output$results <- results
    output$plots <- plots
  }
  
  cat("\n====================================================\n")
  cat("  ANALYSIS COMPLETE\n")
  if (run_trajectory_analysis) cat("  - Trajectory analysis: COMPLETE\n")
  if (run_whole_window) cat("  - Whole-window analysis: COMPLETE\n")
  cat("====================================================\n")
  
  return(output)
}


# ==================== SCRIPT EXECUTION ====================
if (!interactive()) {
  # Run analysis if script is sourced
  analysis_output <- analyze_transition_distcorr(
    recording_id = RECORDING_ID,
    transition_type = TRANSITION_TYPE,
    window_size = WINDOW_SIZE,
    distances_file = DISTANCES_FILE,
    save_output = SAVE_OUTPUT
  )
}
