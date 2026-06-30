# ============================================================
# SINGLE ANIMAL PIPELINE: Full Transition Analysis (v3)
# ============================================================
# Runs complete workflow for one animal: Extract Ã¢â€ â€™ CCF Ã¢â€ â€™ DistCorr
# v3 adds: Activity-based ROI filtering using baseline data
#
# Called by: Batch_Process_Animals_v3.R
# ============================================================

library(dplyr)
library(ggplot2)

# Load component scripts (only if not already loaded)
if (!exists("extract_transitions_main")) source("Extract_Transitions_v3.R")
if (!exists("analyze_transition_ccf")) source("Transition_CCF_v3.R")
if (!exists("analyze_transition_distcorr")) source("Transition_DistCorr_v3.R")
if (!exists("analyze_transition_event_rates")) source("Transition_EventRate_v3.R")


create_output_directories <- function(output_dir) {
  analysis_types <- c("event", "correlation", "distance_correlation", 
                      "within_subject", "between_subject")
  transition_types <- c("Wake2NREM", "NREM2Wake")
  
  for (analysis in analysis_types) {
    # Create base analysis directory
    dir.create(file.path(output_dir, analysis), recursive = TRUE, showWarnings = FALSE)
    
    # Create transition type subdirectories
    for (trans in transition_types) {
      dir.create(file.path(output_dir, analysis, trans), recursive = TRUE, showWarnings = FALSE)
    }
  }
}


#' Get output directory path for a specific analysis type and transition
#' 
#' @param output_dir Base output directory
#' @param analysis_type One of: "event", "correlation", "distance_correlation", 
#'                      "within_subject", "between_subject"
#' @param transition_type "Wake2NREM", "NREM2Wake", or NULL for combined/both
#' @return Full path to appropriate output directory
get_output_path <- function(output_dir, analysis_type, transition_type = NULL) {
  if (is.null(transition_type)) {
    # Combined files go in the analysis type root
    return(file.path(output_dir, analysis_type))
  } else {
    return(file.path(output_dir, analysis_type, transition_type))
  }
}


# ============================================================
# ACTIVITY-BASED ROI FILTERING
# ============================================================

#' Compute event counts per ROI from events file
count_events_per_roi <- function(events_file, data_dir) {
  events_path <- file.path(data_dir, events_file)
  
  if (!file.exists(events_path)) {
    warning("Events file not found: ", events_path)
    return(NULL)
  }
  
  events_df <- read.csv(events_path, stringsAsFactors = FALSE)
  
  # Handle different column name formats
  cell_col <- NULL
  if ("Cell.Name" %in% colnames(events_df)) {
    cell_col <- "Cell.Name"
  } else if (" Cell Name" %in% colnames(events_df)) {
    cell_col <- " Cell Name"
  } else if ("Cell Name" %in% colnames(events_df)) {
    cell_col <- "Cell Name"
  } else if ("Cell_Name" %in% colnames(events_df)) {
    cell_col <- "Cell_Name"
  }
  
  if (is.null(cell_col)) {
    warning("Could not find Cell Name column in events file")
    return(NULL)
  }
  
  event_counts <- events_df %>%
    group_by(ROI = .data[[cell_col]]) %>%
    summarise(event_count = n(), .groups = "drop")
  
  event_counts$ROI <- trimws(event_counts$ROI)
  return(event_counts)
}


#' Identify active ROIs based on baseline activity
identify_active_rois <- function(animal_id, data_dir, 
                                  min_events = 1, 
                                  baseline_epochs = 18,
                                  total_recording_epochs = 360) {
  
  bl_events_file <- paste0(animal_id, "_BL_Events.csv")
  event_counts <- count_events_per_roi(bl_events_file, data_dir)
  
  if (is.null(event_counts)) {
    warning("Could not compute event counts for ", animal_id)
    return(list(active_rois = NULL, stats = NULL, roi_details = NULL))
  }
  
  # Scale threshold: min_events per baseline_epochs Ã¢â€ â€™ total recording
  threshold <- min_events * (total_recording_epochs / baseline_epochs)
  
  # Create detailed ROI report
  roi_details <- event_counts %>%
    mutate(
      threshold = threshold,
      kept = event_count >= threshold,
      status = ifelse(kept, "ACTIVE", "EXCLUDED"),
      events_per_epoch = round(event_count / total_recording_epochs, 4),
      events_per_minute = round(event_count / (total_recording_epochs * 10 / 60), 4)
    ) %>%
    arrange(desc(event_count))
  
  active_rois <- roi_details %>%
    filter(kept) %>%
    pull(ROI)
  
  stats <- list(
    total_rois = nrow(event_counts),
    active_rois = length(active_rois),
    excluded_rois = nrow(event_counts) - length(active_rois),
    threshold = threshold,
    min_events = min_events,
    baseline_epochs = baseline_epochs,
    pct_retained = round(100 * length(active_rois) / nrow(event_counts), 1)
  )
  
  return(list(active_rois = active_rois, stats = stats, roi_details = roi_details))
}


#' Filter traces to only include active ROIs
filter_traces_by_rois <- function(traces_df, active_rois) {
  keep_cols <- c("Time", intersect(colnames(traces_df), active_rois))
  return(traces_df[, keep_cols, drop = FALSE])
}


#' Filter distances to only include active ROI pairs
filter_distances_by_rois <- function(distances_df, active_rois) {
  distances_df %>%
    filter(ROI_1 %in% active_rois & ROI_2 %in% active_rois)
}


#' Compute distance bin characteristics for an animal
#' 
#' @param distances_df Distance data frame with ROI_1, ROI_2, Distance columns
#' @param n_bins Number of percentile bins (default 3: Close, Medium, Far)
#' @return Data frame with bin characteristics (for inclusion in results table)
compute_distance_bin_characteristics <- function(distances_df, n_bins = 3) {
  
  unique_distances <- unique(distances_df$Distance)
  
  # Compute percentile breaks
  percentile_breaks <- quantile(unique_distances, 
                                probs = seq(0, 1, length.out = n_bins + 1),
                                na.rm = TRUE)
  percentile_breaks <- unique(percentile_breaks)
  
  # Assign bins to each distance
  distances_df$distance_bin <- cut(distances_df$Distance,
                                   breaks = percentile_breaks,
                                   include.lowest = TRUE,
                                   labels = FALSE)
  
  # Create bin labels
  if (n_bins == 3) {
    bin_labels <- c("Close", "Medium", "Far")
  } else if (n_bins == 5) {
    bin_labels <- c("Close", "Medium-Close", "Medium", "Medium-Far", "Far")
  } else {
    bin_labels <- paste("Bin", 1:n_bins)
  }
  
  # Compute bin characteristics
  bin_characteristics <- distances_df %>%
    filter(!is.na(distance_bin)) %>%
    group_by(distance_bin) %>%
    summarise(
      bin_min = min(Distance),
      bin_max = max(Distance),
      bin_center = mean(Distance),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      bin_label = bin_labels[distance_bin],
      bin_range = paste0(round(bin_min, 1), "-", round(bin_max, 1))
    ) %>%
    select(bin = distance_bin, bin_label, bin_min, bin_max, bin_center, bin_range, n_pairs)
  
  return(list(
    bin_characteristics = bin_characteristics,
    percentile_breaks = percentile_breaks
  ))
}


# ============================================================
# INDIVIDUAL ANIMAL PLOTTING FUNCTIONS
# ============================================================

#' Create trajectory plot for a single animal (all conditions overlaid)
create_animal_trajectory_plot <- function(trajectory_data, animal_id, window_size, 
                                          filter_label = "") {
  
  plot_data <- trajectory_data %>%
    group_by(condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      n_pairs = n(),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                              color = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  width = 0.1, alpha = 0.5, linewidth = 0.5) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    facet_wrap(~transition_type, ncol = 2, scales = "free_y") +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72", "WO" = "#F18F01"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep", "WO" = "Washout")
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep", "WO" = "Washout")
    ) +
    labs(
      title = paste0(animal_id, " - Correlation Trajectory"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation",
      color = "Condition",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.title.position = "plot",
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}


#' Create distance-binned trajectory plot for a single animal
create_animal_trajectory_by_distance_plot <- function(trajectory_dist, animal_id, 
                                                       window_size, filter_label = "") {
  
  plot_data <- trajectory_dist %>%
    group_by(condition, transition_type, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plots <- list()
  
  for (trans_type in unique(plot_data$transition_type)) {
    trans_data <- plot_data %>% filter(transition_type == trans_type)
    
    plots[[trans_type]] <- ggplot(trans_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = distance_bin_label, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                        ymax = mean_correlation + se_correlation),
                    width = 0.1, alpha = 0.4, linewidth = 0.4) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2) +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Dep", "WO" = "Washout")
      ) +
      scale_color_viridis_d(option = "plasma") +
      labs(
        title = paste0(animal_id, " - ", trans_type, " - Correlation by Distance"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation",
        color = "Distance Bin",
        linetype = "Condition"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.title.position = "plot",
        axis.title = element_text(size = 11),
        legend.position = "right"
      )
  }
  
  return(plots)
}


#' Create event rate trajectory plot for a single animal (faceted by transition type)
create_animal_event_rate_plot <- function(trajectory_events, animal_id, 
                                          window_size, filter_label = "") {
  
  if (is.null(trajectory_events) || nrow(trajectory_events) == 0) {
    return(NULL)
  }
  
  # Calculate mean and SE per ROI per epoch per condition per transition
  roi_epoch_data <- trajectory_events %>%
    group_by(condition, transition_type, epoch_in_window, Cell_Name) %>%
    summarise(
      events_per_epoch = n() / n_distinct(transition_id),
      .groups = "drop"
    )
  
  # Then summarize across ROIs to get mean ± SE
  plot_data <- roi_epoch_data %>%
    group_by(condition, transition_type, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(events_per_epoch, na.rm = TRUE),
      se_event_rate = sd(events_per_epoch, na.rm = TRUE) / sqrt(n()),
      n_rois = n(),
      .groups = "drop"
    )
  
  # Create single faceted plot
 p <- ggplot(plot_data, 
      aes(x = epoch_in_window, y = mean_event_rate,
          color = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_errorbar(aes(ymin = mean_event_rate - se_event_rate,
                      ymax = mean_event_rate + se_event_rate),
                  width = 0.1, alpha = 0.5, linewidth = 0.5) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    facet_wrap(~transition_type, ncol = 2, scales = "free_y") +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72", "WO" = "#F18F01"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep", "WO" = "Washout")
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed", "WO" = "dotted"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep", "WO" = "Washout")
    ) +
    labs(
      title = paste0(animal_id, " - Event Rate Trajectory"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)",
      color = "Condition",
      linetype = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.title.position = "plot",
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}
# ============================================================
# SIMPLIFIED PLOTTING FUNCTIONS - WAKE→NREM FOCUS
# ============================================================
# Add these to Pipeline_Transition_Analysis_v3.R
# After the existing plotting functions (around line 355)

#' Create SIMPLIFIED trajectory plot (BL vs SD only, Wake→NREM default)
create_animal_trajectory_plot_simplified <- function(trajectory_data, animal_id, window_size, 
                                                     filter_label = "", 
                                                     transition_filter = "Wake2NREM") {
  
  # Filter to BL and SD only, and one transition type
  plot_data <- trajectory_data %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      n_pairs = n(),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                             color = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  width = 0.1, alpha = 0.5, linewidth = 0.5) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation")
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation")
    ) +
    labs(
      title = paste0(animal_id, " - Sleep Onset Transitions", filter_label),
      subtitle = paste0("Mean ± SE pairwise correlation (", window_size, " epoch window)"),
      x = "Epoch Position Relative to Sleep Onset",
      y = "Mean Pairwise Correlation",
      color = "Condition",
      linetype = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  return(p)
}


#' Create SIMPLIFIED distance-binned plot (BL only, Wake→NREM default)
create_animal_trajectory_by_distance_plot_simplified <- function(trajectory_dist, animal_id, 
                                                                 window_size, filter_label = "",
                                                                 transition_filter = "Wake2NREM") {
  
  # Compute bin ranges for subtitle
  bin_ranges <- trajectory_dist %>%
    group_by(distance_bin_label) %>%
    summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    arrange(bin_min) %>%
    mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    pull(range_str)
  bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
  # Filter to BL only and one transition type
  plot_data <- trajectory_dist %>%
    filter(condition == "BL") %>%
    filter(transition_type == transition_filter) %>%
    group_by(distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, 
              aes(x = epoch_in_window, y = mean_correlation,
                  color = distance_bin_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_errorbar(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  width = 0.1, alpha = 0.4, linewidth = 0.4) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = paste0(animal_id, " - Sleep Onset - Baseline", filter_label),
      subtitle = bin_ranges_subtitle,
      x = "Epoch Position Relative to Sleep Onset",
      y = "Mean Pairwise Correlation",
      color = "Distance Bin"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(size = 12)
    )
  
  return(p)
}


#' Create SIMPLIFIED event rate plot (BL vs SD only, Wake→NREM default)
create_animal_event_rate_plot_simplified <- function(trajectory_events, animal_id, 
                                                     window_size, filter_label = "",
                                                     transition_filter = "Wake2NREM") {
  
  if (is.null(trajectory_events) || nrow(trajectory_events) == 0) {
    return(NULL)
  }
  
  # Filter to BL and SD only, and one transition type
  plot_data <- trajectory_events %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      total_events = n(),
      n_rois = n_distinct(Cell_Name),
      n_transitions = n_distinct(transition_id),
      mean_event_rate = total_events / (n_rois * n_transitions),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, 
              aes(x = epoch_in_window, y = mean_event_rate,
                  color = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation")
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Deprivation")
    ) +
    labs(
      title = paste0(animal_id, " - Sleep Onset Transitions", filter_label),
      subtitle = paste0("Mean events per ROI per transition (", window_size, " epoch window)"),
      x = "Epoch Position Relative to Sleep Onset",
      y = "Mean Event Rate",
      color = "Condition",
      linetype = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(size = 12)
    )
  
  return(p)
}

# ============================================================
# CORRELATION-BASED PAIR FILTERING
# ============================================================

#' Identify highly correlated pairs based on pre-transition correlation
#' 
#' @param trajectory_data Combined trajectory data with correlation per epoch
#' @param window_size Window size (determines which epochs are "pre-transition")
#' @param method "percentile" for top X%, "outlier" for boxplot upper fence
#' @param percentile Top percentile to keep (e.g., 10 = top 10%)
#' @return List with high_corr_pairs, threshold, and filter statistics
identify_high_correlation_pairs <- function(trajectory_data, 
                                            window_size,
                                            method = "percentile",
                                            percentile = 10) {
  

  cat("\n[CORRELATION FILTER] Identifying highly correlated pairs...\n")
  
 # Pre-transition epochs are all negative epochs: -window_size to -1
  # (Epoch +1 is post-transition, epoch 0 doesn't exist)
  pre_transition_epochs <- -window_size:-1
  
  cat("  Pre-transition epochs:", paste(pre_transition_epochs, collapse = ", "), "\n")
  cat("  Available epochs in data:", paste(sort(unique(trajectory_data$epoch_in_window)), collapse = ", "), "\n")
  
  # Create pair key for grouping
  pair_data <- trajectory_data %>%
    mutate(pair_key = paste(pmin(ROI_1, ROI_2), pmax(ROI_1, ROI_2), sep = "_"))
  
  # Calculate mean pre-transition correlation per pair (across all conditions)
  pair_pre_corr <- pair_data %>%
    filter(epoch_in_window %in% pre_transition_epochs) %>%
    group_by(pair_key, ROI_1, ROI_2) %>%
    summarise(
      mean_pre_corr = mean(correlation, na.rm = TRUE),
      n_observations = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_pre_corr))
  
  cat("  Total pairs:", nrow(pair_pre_corr), "\n")
  cat("  Mean pre-transition correlation:", round(mean(pair_pre_corr$mean_pre_corr), 4), "\n")
  cat("  SD:", round(sd(pair_pre_corr$mean_pre_corr), 4), "\n")
  
  # Calculate threshold based on method
  if (method == "outlier") {
    q1 <- quantile(pair_pre_corr$mean_pre_corr, 0.25, na.rm = TRUE)
    q3 <- quantile(pair_pre_corr$mean_pre_corr, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    threshold <- q3 + 1.5 * iqr
    cat("  Method: Outlier (Q3 + 1.5*IQR)\n")
    cat("  Q1:", round(q1, 4), " Q3:", round(q3, 4), " IQR:", round(iqr, 4), "\n")
  } else {
    # Percentile method
    threshold <- quantile(pair_pre_corr$mean_pre_corr, 1 - percentile/100, na.rm = TRUE)
    cat("  Method: Top", percentile, "percentile\n")
  }
  
  cat("  Threshold:", round(threshold, 4), "\n")
  
  # Identify pairs above threshold
  high_corr_pairs <- pair_pre_corr %>%
    filter(mean_pre_corr >= threshold)
  
  cat("  Pairs above threshold:", nrow(high_corr_pairs), 
      "(", round(100 * nrow(high_corr_pairs) / nrow(pair_pre_corr), 1), "%)\n")
  
  # Return results
  stats <- list(
    total_pairs = nrow(pair_pre_corr),
    high_corr_pairs = nrow(high_corr_pairs),
    pct_retained = round(100 * nrow(high_corr_pairs) / nrow(pair_pre_corr), 1),
    threshold = threshold,
    method = method,
    percentile = if(method == "percentile") percentile else NA,
    pre_transition_epochs = pre_transition_epochs
  )
  
  return(list(
    high_corr_pairs = high_corr_pairs$pair_key,
    pair_details = high_corr_pairs,
    all_pair_corr = pair_pre_corr,
    stats = stats
  ))
}


#' Filter trajectory data to only include high-correlation pairs
filter_trajectory_by_pairs <- function(trajectory_data, high_corr_pairs) {
  trajectory_data %>%
    mutate(pair_key = paste(pmin(ROI_1, ROI_2), pmax(ROI_1, ROI_2), sep = "_")) %>%
    filter(pair_key %in% high_corr_pairs) %>%
    select(-pair_key)
}


# ============================================================
# MAIN PIPELINE FUNCTION
# ============================================================

run_animal_pipeline <- function(animal_id,
                                conditions = c("BL", "SD", "WO"),
                                transition_types = c("Wake2NREM", "NREM2Wake"),
                                window_size = 3,
                                data_dir = "E:/Data_Processing/R/Data CSVs",
                                output_dir = "E:/Data_Processing/R/Results",
                                save_intermediate = TRUE,
                                force_save_intermediate = FALSE,  # Override auto-suppression when filtering
                                ccf_use_parallel = TRUE,
                                ccf_n_cores = 6,
                                run_trajectory_analysis = TRUE,
                                run_whole_window = FALSE,
                                run_lag_sweep = FALSE,
                                max_lag = 5,
                                binning_method = "percentile",
                                n_distance_bins = 10,
                                n_percentile_bins = 3,
                                n_trajectory_bins = 3,
                                fit_exponential = FALSE,
                                filter_by_activity = TRUE,
                                min_events_baseline = 1,
                                baseline_epochs = 18,
                                correlation_filter = "off",
                                correlation_filter_method = "percentile",
                                correlation_percentile = 10,
                                file_suffix = "_filtered") {
  
  cat("\n========================================\n")
  cat("PROCESSING ANIMAL:", animal_id, "\n")
  cat("========================================\n")
  cat("Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("Transitions:", paste(transition_types, collapse = ", "), "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("CCF parallel:", ccf_use_parallel, "(", ccf_n_cores, "cores )\n")
  
  # Auto-suppress intermediate saving when correlation filtering is enabled
  # (Intermediate files show all pairs; final output shows filtered pairs - confusing)
  if (correlation_filter != "off" && save_intermediate && !force_save_intermediate) {
    cat("\n[NOTE] Suppressing intermediate sub-script outputs (correlation_filter =", 
        correlation_filter, ")\n")
    cat("       Final filtered outputs will be saved to within_subject/\n")
    cat("       (Use force_save_intermediate = TRUE to override)\n")
    save_intermediate <- FALSE
  }
  
  # Create output directory structure
  create_output_directories(output_dir)
  
  # ============================================================
  # STEP 0: ACTIVITY-BASED ROI FILTERING
  # ============================================================
  
  active_rois <- NULL
  filter_stats <- NULL
  roi_details <- NULL
  
  if (filter_by_activity) {
    cat("\n[ACTIVITY FILTER] Identifying active ROIs from baseline...\n")
    
    filter_result <- identify_active_rois(
      animal_id = animal_id,
      data_dir = data_dir,
      min_events = min_events_baseline,
      baseline_epochs = baseline_epochs
    )
    
    active_rois <- filter_result$active_rois
    filter_stats <- filter_result$stats
    roi_details <- filter_result$roi_details
    
    if (is.null(active_rois) || length(active_rois) == 0) {
      warning("No active ROIs found for ", animal_id, " - skipping")
      return(NULL)
    }
    
    cat("  Threshold:", filter_stats$threshold, "events in BL recording\n")
    cat("  (Based on", min_events_baseline, "event(s) per", baseline_epochs, "epochs)\n")
    cat("  Active ROIs:", filter_stats$active_rois, "of", filter_stats$total_rois, 
        "(", filter_stats$pct_retained, "% retained)\n")
    cat("  Excluded:", filter_stats$excluded_rois, "ROIs\n\n")
  }
  
  # Load and filter distances
  distances_file <- file.path(data_dir, paste0(animal_id, "_roi_distances.csv"))
  if (!file.exists(distances_file)) {
    stop("Distance file not found: ", distances_file)
  }
  distances_df <- read.csv(distances_file, stringsAsFactors = FALSE)
  cat("Loaded distances:", nrow(distances_df), "pairs\n")
  
  if (filter_by_activity) {
    distances_df <- filter_distances_by_rois(distances_df, active_rois)
    cat("After filtering:", nrow(distances_df), "pairs\n")
  }
  
  # Compute distance bin characteristics (same for all conditions since distances don't change)
  distance_bin_info <- compute_distance_bin_characteristics(distances_df, n_trajectory_bins)
  cat("\nDistance bin boundaries (", n_trajectory_bins, " bins):\n", sep = "")
  print(distance_bin_info$bin_characteristics)
  cat("\n")
  
  # Storage for results
  all_event_rates <- list()
  all_correlations <- list()
  all_dist_corr <- list()
  all_trajectory_data <- list()
  all_trajectory_dist <- list()
  all_trajectory_events <- list()
  all_trajectory_event_rates <- list()
  
  # Process each condition
  for (condition in conditions) {
    recording_id <- paste0(animal_id, "_", condition)
    
    cat("--------------------------------------------------\n")
    cat("CONDITION:", condition, "(", recording_id, ")\n")
    cat("--------------------------------------------------\n")
    
    traces_file <- file.path(data_dir, paste0(recording_id, "_Traces.csv"))
    
    if (!file.exists(traces_file)) {
      warning("Traces file not found: ", traces_file, " - Skipping ", condition)
      next
    }
    
    # Load traces
    cat("  Loading traces...\n")
    traces_full <- read.csv(traces_file, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(traces_full)[1] <- "Time"
    if (is.na(suppressWarnings(as.numeric(traces_full[1, 1])))) {
      traces_full <- traces_full[-1, ]
    }
    traces_full[] <- lapply(traces_full, as.numeric)
    
    original_n_rois <- ncol(traces_full) - 1
    cat("  Loaded:", nrow(traces_full), "timepoints,", original_n_rois, "ROIs\n")
    
    # Apply ROI filter
    if (filter_by_activity) {
      traces_full <- filter_traces_by_rois(traces_full, active_rois)
      cat("  After filtering:", ncol(traces_full) - 1, "ROIs\n\n")
    }
    
    # Process each transition type
    for (trans_type in transition_types) {
      
      cat("  >> Transition Type:", trans_type, "\n")
      
      # Extract transitions
      cat("     [1/3] Extracting transitions...\n")
      extraction <- extract_transitions_main(
        states_file = paste0(recording_id, "_states_df.csv"),
        traces_file = paste0(recording_id, "_Traces.csv"),
        events_file = paste0(recording_id, "_Events.csv"),
        recording_id = recording_id,
        input_dir = data_dir,
        window_sizes = window_size,
        save_output = save_intermediate,
        output_dir = output_dir,
        file_suffix = file_suffix
      )
      
      metadata_df <- extraction$metadata %>% 
        filter(transition_type == trans_type, window_size == !!window_size)
      
      n_trans <- nrow(metadata_df)
      cat("     Found", n_trans, "transitions\n")
      
      if (n_trans == 0) {
        warning("No transitions found for ", recording_id, " - ", trans_type)
        next
      }
      
      # ==== WHOLE-WINDOW ANALYSIS ====
      if (run_whole_window) {
        cat("     [WHOLE-WINDOW] Computing correlations...\n")
        
        ccf_results <- analyze_transition_ccf(
          recording_id = recording_id,
          transition_type = trans_type,
          window_size = window_size,
          traces_df = traces_full,
          metadata_df = metadata_df,
          run_lag_sweep = run_lag_sweep,
          max_lag = max_lag,
          use_parallel = ccf_use_parallel,
          n_cores = ccf_n_cores,
          save_output = save_intermediate,
          output_dir = get_output_path(output_dir, "correlation", trans_type),
          file_suffix = file_suffix
        )
        
        distcorr_results <- analyze_transition_distcorr(
          recording_id = recording_id,
          transition_type = trans_type,
          window_size = window_size,
          distances_df = distances_df,
          metadata_df = metadata_df,
          ccf_results_list = ccf_results$ccf_results,
          binning_method = binning_method,
          n_equal_width_bins = n_distance_bins,
          n_percentile_bins = n_percentile_bins,
          fit_exponential = fit_exponential,
          save_output = save_intermediate,
          output_dir = get_output_path(output_dir, "distance_correlation", trans_type),
          file_suffix = file_suffix
        )
        
        # Filter events to active ROIs
        if (filter_by_activity) {
          filtered_events <- extraction$events %>% filter(Cell_Name %in% active_rois)
        } else {
          filtered_events <- extraction$events
        }
        
        event_rates <- filtered_events %>%
          group_by(Cell_Name) %>%
          summarize(mean_event_rate = n() / n_distinct(transition_id), .groups = 'drop') %>%
          rename(ROI = Cell_Name) %>%
          mutate(animal = animal_id, condition = condition, 
                 transition_type = trans_type, window_size = window_size)
        
        correlations <- distcorr_results$results$dist_corr_all %>%
          group_by(ROI_1, ROI_2, Distance) %>%
          summarize(mean_correlation = mean(Correlation, na.rm = TRUE), .groups = 'drop') %>%
          mutate(animal = animal_id, condition = condition,
                 transition_type = trans_type, window_size = window_size)
        
        key <- paste(condition, trans_type, sep = "_")
        all_event_rates[[key]] <- event_rates
        all_correlations[[key]] <- correlations
        all_dist_corr[[key]] <- distcorr_results$results$dist_corr_all %>%
          mutate(animal = animal_id, condition = condition,
                 transition_type = trans_type, window_size = window_size)
      }
      
      # ==== TRAJECTORY ANALYSIS ====
      if (run_trajectory_analysis) {
        cat("     [TRAJECTORY] Running epoch-level analysis...\n")
        
        ccf_traj <- analyze_transition_ccf(
          recording_id = recording_id,
          transition_type = trans_type,
          window_size = window_size,
          traces_df = traces_full,
          metadata_df = metadata_df,
          run_trajectory_analysis = TRUE,
          run_whole_window = FALSE,
          use_parallel = ccf_use_parallel,
          n_cores = ccf_n_cores,
          save_output = save_intermediate,
          output_dir = get_output_path(output_dir, "correlation", trans_type),
          file_suffix = file_suffix
        )
        
        distcorr_traj <- analyze_transition_distcorr(
          recording_id = recording_id,
          transition_type = trans_type,
          window_size = window_size,
          distances_df = distances_df,
          trajectory_data = ccf_traj$trajectory_data,
          run_trajectory_analysis = TRUE,
          run_whole_window = FALSE,
          n_trajectory_bins = n_trajectory_bins,
          save_output = save_intermediate,
          output_dir = get_output_path(output_dir, "distance_correlation", trans_type),
          file_suffix = file_suffix
        )
        
        # Event rate analysis
        cat("     [EVENT RATES] Analyzing single-cell event rates...\n")
        
        # Filter events to active ROIs first
        if (filter_by_activity) {
          filtered_events <- extraction$events %>% filter(Cell_Name %in% active_rois)
        } else {
          filtered_events <- extraction$events
        }
        
        eventrate_results <- analyze_transition_event_rates(
          recording_id = recording_id,
          transition_type = trans_type,
          window_size = window_size,
          events_df = filtered_events,
          metadata_df = metadata_df,
          save_output = save_intermediate,
          output_dir = get_output_path(output_dir, "event", trans_type),
          file_suffix = file_suffix
        )
        
        key <- paste(condition, trans_type, sep = "_")
        all_trajectory_data[[key]] <- ccf_traj$trajectory_data %>%
          mutate(animal = animal_id, condition = condition,
                 transition_type = trans_type, window_size = window_size)
        
        all_trajectory_dist[[key]] <- distcorr_traj$trajectory_results$binned_data %>%
          mutate(animal = animal_id, condition = condition,
                 transition_type = trans_type, window_size = window_size)
        
        all_trajectory_events[[key]] <- filtered_events %>%
          mutate(animal = animal_id, condition = condition)
        
        all_trajectory_event_rates[[key]] <- eventrate_results$event_rates %>%
          mutate(animal = animal_id, condition = condition, transition_type = trans_type)
      }
      
      cat("     Complete!\n\n")
    }
    
    rm(traces_full)
    gc()
  }
  
  # Combine results
  cat("Combining results across conditions...\n")
  
  event_rates_ranked <- NULL
  correlations_ranked <- NULL
  dist_corr_ranked <- NULL
  
  if (run_whole_window && length(all_event_rates) > 0) {
    event_rates_combined <- bind_rows(all_event_rates)
    correlations_combined <- bind_rows(all_correlations)
    dist_corr_combined <- bind_rows(all_dist_corr)
    
    event_rates_ranked <- event_rates_combined %>%
      arrange(mean_event_rate) %>%
      mutate(activity_rank = row_number(),
             activity_percentile = 100 * (activity_rank - 1) / (n() - 1))
    
    correlations_ranked <- correlations_combined %>%
      arrange(mean_correlation) %>%
      mutate(correlation_rank = row_number(),
             correlation_percentile = 100 * (correlation_rank - 1) / (n() - 1))
    
    dist_corr_ranked <- dist_corr_combined %>%
      arrange(Distance) %>%
      mutate(distance_rank = row_number(),
             distance_percentile = 100 * (distance_rank - 1) / (n() - 1))
  }
  
  trajectory_data_combined <- NULL
  trajectory_dist_combined <- NULL
  trajectory_events_combined <- NULL
  trajectory_event_rates_combined <- NULL
  
  if (run_trajectory_analysis && length(all_trajectory_data) > 0) {
    trajectory_data_combined <- bind_rows(all_trajectory_data)
    trajectory_dist_combined <- bind_rows(all_trajectory_dist)
    trajectory_events_combined <- bind_rows(all_trajectory_events)
    trajectory_event_rates_combined <- bind_rows(all_trajectory_event_rates)
    cat("  Trajectory data combined\n")
  }
  
  # ============================================================
  # CORRELATION-BASED PAIR FILTERING
  # ============================================================
  
  trajectory_data_highcorr <- NULL
  trajectory_dist_highcorr <- NULL
  corr_filter_results <- NULL
  highcorr_tag <- ""
  main_file_tag <- ""  # Tag for main plots when correlation_filter == "alone"
  
  if (run_trajectory_analysis && correlation_filter != "off" && !is.null(trajectory_data_combined)) {
    
    # Identify high-correlation pairs
    corr_filter_results <- identify_high_correlation_pairs(
      trajectory_data = trajectory_data_combined,
      window_size = window_size,
      method = correlation_filter_method,
      percentile = correlation_percentile
    )
    
    # Create tag for file naming
    if (correlation_filter_method == "percentile") {
      highcorr_tag <- paste0("_top", correlation_percentile, "pct")
    } else {
      highcorr_tag <- "_outlier"
    }
    
    # Set main_file_tag for "alone" mode (main plots ARE the highcorr plots)
    if (correlation_filter == "alone") {
      main_file_tag <- highcorr_tag
    }
    
    # Filter trajectory data to high-correlation pairs only
    trajectory_data_highcorr <- filter_trajectory_by_pairs(
      trajectory_data_combined, 
      corr_filter_results$high_corr_pairs
    )
    
    trajectory_dist_highcorr <- filter_trajectory_by_pairs(
      trajectory_dist_combined,
      corr_filter_results$high_corr_pairs
    )
    
    cat("  High-correlation trajectory data:", nrow(trajectory_data_highcorr), "rows\n")
    cat("  High-correlation trajectory dist:", nrow(trajectory_dist_highcorr), "rows\n")
    
    # Warn if no pairs found
    if (nrow(trajectory_data_highcorr) == 0) {
      warning("No high-correlation pairs found - highcorr plots will be skipped")
    }
    
    # If "alone" mode, replace the main data with filtered data
    if (correlation_filter == "alone") {
      trajectory_data_combined <- trajectory_data_highcorr
      trajectory_dist_combined <- trajectory_dist_highcorr
      cat("  Mode: 'alone' - using only high-correlation pairs\n")
    } else {
      cat("  Mode: 'alongside' - will generate both all-pairs and high-correlation plots\n")
    }
  }
  
  # Save intermediate results
  if (save_intermediate) {
    cat("\nSaving intermediate results...\n")
    suffix <- file_suffix
    filter_label <- ifelse(filter_by_activity, " (Active ROIs Only)", "")
    
    # Get within_subject output path
    within_subj_dir <- get_output_path(output_dir, "within_subject", NULL)
    
    # ---- FILTERING REPORT ----
    if (filter_by_activity && !is.null(roi_details)) {
      cat("  Saving filtering report...\n")
      
      # Add animal ID to report
      roi_report <- roi_details %>%
        mutate(animal = animal_id) %>%
        select(animal, ROI, event_count, events_per_minute, threshold, status, kept)
      
      write.csv(roi_report, 
                file.path(within_subj_dir, paste0(animal_id, "_ROI_filtering_report.csv")),
                row.names = FALSE)
      cat("    Saved:", paste0(animal_id, "_ROI_filtering_report.csv"), "\n")
      
      # Summary stats file
      summary_df <- data.frame(
        animal = animal_id,
        total_rois = filter_stats$total_rois,
        active_rois = filter_stats$active_rois,
        excluded_rois = filter_stats$excluded_rois,
        pct_retained = filter_stats$pct_retained,
        threshold_events = filter_stats$threshold,
        min_events_setting = filter_stats$min_events,
        baseline_epochs_setting = filter_stats$baseline_epochs
      )
      write.csv(summary_df,
                file.path(within_subj_dir, paste0(animal_id, "_filtering_summary.csv")),
                row.names = FALSE)
      cat("    Saved:", paste0(animal_id, "_filtering_summary.csv"), "\n")
    }
    
    # ---- DATA CSVs ----
    if (run_whole_window && !is.null(event_rates_ranked)) {
      write.csv(event_rates_ranked, 
                file.path(within_subj_dir, paste0(animal_id, "_event_rates_ranked_", window_size, "ep", suffix, ".csv")),
                row.names = FALSE)
      write.csv(correlations_ranked,
                file.path(within_subj_dir, paste0(animal_id, "_correlations_ranked_", window_size, "ep", suffix, ".csv")),
                row.names = FALSE)
      write.csv(dist_corr_ranked,
                file.path(within_subj_dir, paste0(animal_id, "_dist_corr_ranked_", window_size, "ep", suffix, ".csv")),
                row.names = FALSE)
    }
    
    if (run_trajectory_analysis && !is.null(trajectory_data_combined)) {
      write.csv(trajectory_data_combined,
                file.path(within_subj_dir, paste0(animal_id, "_trajectory_data", main_file_tag, "_", window_size, "ep", suffix, ".csv")),
                row.names = FALSE)
      write.csv(trajectory_dist_combined,
                file.path(within_subj_dir, paste0(animal_id, "_trajectory_dist", main_file_tag, "_", window_size, "ep", suffix, ".csv")),
                row.names = FALSE)
      
      # Save event rate trajectory data (not affected by correlation filtering)
      if (!is.null(trajectory_event_rates_combined) && nrow(trajectory_event_rates_combined) > 0) {
        write.csv(trajectory_event_rates_combined,
                  file.path(within_subj_dir, paste0(animal_id, "_trajectory_eventrates_", window_size, "ep", suffix, ".csv")),
                  row.names = FALSE)
      }
      cat("  Saved trajectory data CSVs\n")
      
      # Save high-correlation pair data if filtering was applied
      if (correlation_filter != "off" && !is.null(corr_filter_results)) {
        # Save pair correlation report
        write.csv(corr_filter_results$all_pair_corr,
                  file.path(within_subj_dir, paste0(animal_id, "_pair_correlations_", window_size, "ep", suffix, ".csv")),
                  row.names = FALSE)
        write.csv(corr_filter_results$pair_details,
                  file.path(within_subj_dir, paste0(animal_id, "_high_corr_pairs_", window_size, "ep", suffix, ".csv")),
                  row.names = FALSE)
        
        # Save filter summary
        filter_summary <- data.frame(
          animal = animal_id,
          method = corr_filter_results$stats$method,
          percentile = corr_filter_results$stats$percentile,
          threshold = corr_filter_results$stats$threshold,
          total_pairs = corr_filter_results$stats$total_pairs,
          high_corr_pairs = corr_filter_results$stats$high_corr_pairs,
          pct_retained = corr_filter_results$stats$pct_retained,
          pre_transition_epochs = paste(corr_filter_results$stats$pre_transition_epochs, collapse = ",")
        )
        write.csv(filter_summary,
                  file.path(within_subj_dir, paste0(animal_id, "_corr_filter_summary_", window_size, "ep", suffix, ".csv")),
                  row.names = FALSE)
        cat("  Saved correlation filter reports\n")
        
        # Save high-correlation filtered data if in alongside mode and data exists
        if (correlation_filter == "alongside" && !is.null(trajectory_data_highcorr) && nrow(trajectory_data_highcorr) > 0) {
          write.csv(trajectory_data_highcorr,
                    file.path(within_subj_dir, paste0(animal_id, "_trajectory_data_highcorr", highcorr_tag, "_", window_size, "ep", suffix, ".csv")),
                    row.names = FALSE)
          write.csv(trajectory_dist_highcorr,
                    file.path(within_subj_dir, paste0(animal_id, "_trajectory_dist_highcorr", highcorr_tag, "_", window_size, "ep", suffix, ".csv")),
                    row.names = FALSE)
          cat("  Saved high-correlation trajectory data CSVs\n")
        }
      }
    }
    
    # ---- INDIVIDUAL ANIMAL PLOTS ----
    if (run_trajectory_analysis && !is.null(trajectory_data_combined)) {
      cat("  Generating individual animal plots...\n")
      
      # Determine label for plots
      plot_label <- filter_label
      if (correlation_filter == "alone") {
        if (correlation_filter_method == "percentile") {
          method_desc <- paste0("Top ", correlation_percentile, "%")
        } else {
          method_desc <- "Outliers"
        }
        plot_label <- paste0(filter_label, " (High-Corr: ", method_desc, ")")
      }
      
      # Overall trajectory plot (combined transition types -> within_subject root)
      traj_plot <- create_animal_trajectory_plot(
        trajectory_data_combined, animal_id, window_size, plot_label
      )
      ggsave(file.path(within_subj_dir, paste0(animal_id, "_Trajectory_Overall", main_file_tag, "_", window_size, "ep", suffix, ".png")),
             traj_plot, width = 12, height = 6, dpi = 300)
      cat("    Saved:", paste0(animal_id, "_Trajectory_Overall", main_file_tag, "_", window_size, "ep", suffix, ".png"), "\n")
      
      # Distance-binned trajectory plots (one per transition type)
      if (!is.null(trajectory_dist_combined)) {
        dist_plots <- create_animal_trajectory_by_distance_plot(
          trajectory_dist_combined, animal_id, window_size, plot_label
        )
        for (trans_type in names(dist_plots)) {
          trans_dir <- get_output_path(output_dir, "within_subject", trans_type)
          ggsave(file.path(trans_dir, paste0(animal_id, "_", trans_type, "_Trajectory_by_Distance", main_file_tag, "_", 
                                              window_size, "ep", suffix, ".png")),
                 dist_plots[[trans_type]], width = 10, height = 8, dpi = 300)
          cat("    Saved:", paste0(animal_id, "_", trans_type, "_Trajectory_by_Distance", main_file_tag, "_", 
                                   window_size, "ep", suffix, ".png"), "\n")
        }
      }
      
      # Event rate trajectory plots (now faceted, single plot)
      if (!is.null(trajectory_events_combined) && nrow(trajectory_events_combined) > 0) {
        event_plot <- create_animal_event_rate_plot(
          trajectory_events_combined, animal_id, window_size, filter_label
        )
        if (!is.null(event_plot)) {
          ggsave(file.path(within_subj_dir, paste0(animal_id, "_EventRate_Trajectory_", 
                                                    window_size, "ep", suffix, ".png")),
                 event_plot, width = 12, height = 6, dpi = 300)
          cat("    Saved:", paste0(animal_id, "_EventRate_Trajectory_", 
                                   window_size, "ep", suffix, ".png"), "\n")
        }
      }
      
      # ---- HIGH-CORRELATION PLOTS (alongside mode) ----
      if (correlation_filter == "alongside" && !is.null(trajectory_data_highcorr) && nrow(trajectory_data_highcorr) > 0) {
        cat("  Generating high-correlation pair plots...\n")
        
        # Create descriptive label for plots
        if (correlation_filter_method == "percentile") {
          method_desc <- paste0("Top ", correlation_percentile, "%")
        } else {
          method_desc <- "Outliers"
        }
        highcorr_label <- paste0(filter_label, " (High-Corr: ", method_desc, ")")
        
        # Overall trajectory plot - high corr
        traj_plot_hc <- create_animal_trajectory_plot(
          trajectory_data_highcorr, animal_id, window_size, highcorr_label
        )
        ggsave(file.path(within_subj_dir, paste0(animal_id, "_Trajectory_Overall_HighCorr", highcorr_tag, "_", window_size, "ep", suffix, ".png")),
               traj_plot_hc, width = 12, height = 6, dpi = 300)
        cat("    Saved:", paste0(animal_id, "_Trajectory_Overall_HighCorr", highcorr_tag, "_", window_size, "ep", suffix, ".png"), "\n")
        
        # Distance-binned trajectory plots - high corr
        if (!is.null(trajectory_dist_highcorr) && nrow(trajectory_dist_highcorr) > 0) {
          dist_plots_hc <- create_animal_trajectory_by_distance_plot(
            trajectory_dist_highcorr, animal_id, window_size, highcorr_label
          )
          for (trans_type in names(dist_plots_hc)) {
            trans_dir <- get_output_path(output_dir, "within_subject", trans_type)
            ggsave(file.path(trans_dir, paste0(animal_id, "_", trans_type, "_Trajectory_by_Distance_HighCorr", highcorr_tag, "_", 
                                                window_size, "ep", suffix, ".png")),
                   dist_plots_hc[[trans_type]], width = 10, height = 8, dpi = 300)
            cat("    Saved:", paste0(animal_id, "_", trans_type, "_Trajectory_by_Distance_HighCorr", highcorr_tag, "_", 
                                     window_size, "ep", suffix, ".png"), "\n")
          }
        }
      }
    }
    
    cat("  All outputs saved!\n")
  }
  
  cat("\n========================================\n")
  cat("ANIMAL", animal_id, "COMPLETE\n")
  if (filter_by_activity) {
    cat("(Filtered:", filter_stats$active_rois, "of", filter_stats$total_rois, "ROIs)\n")
  }
  if (correlation_filter != "off" && !is.null(corr_filter_results)) {
    cat("(High-corr pairs:", corr_filter_results$stats$high_corr_pairs, "of", 
        corr_filter_results$stats$total_pairs, ")\n")
  }
  cat("========================================\n\n")
  
  return(list(
    animal_id = animal_id,
    event_rates = event_rates_ranked,
    correlations = correlations_ranked,
    dist_corr = dist_corr_ranked,
    trajectory_data = trajectory_data_combined,
    trajectory_dist = trajectory_dist_combined,
    trajectory_events = trajectory_events_combined,
    trajectory_event_rates = trajectory_event_rates_combined,
    trajectory_data_highcorr = trajectory_data_highcorr,
    trajectory_dist_highcorr = trajectory_dist_highcorr,
    corr_filter_results = corr_filter_results,
    filter_stats = filter_stats,
    roi_details = roi_details,
    distance_bin_characteristics = distance_bin_info$bin_characteristics,
    metadata = list(
      conditions = conditions,
      transition_types = transition_types,
      window_size = window_size,
      n_rois_original = if(!is.null(filter_stats)) filter_stats$total_rois else NA,
      n_rois_active = if(!is.null(filter_stats)) filter_stats$active_rois else NA,
      n_pairs = nrow(distances_df),
      filtered = filter_by_activity,
      correlation_filter = correlation_filter,
      correlation_filter_method = correlation_filter_method,
      correlation_percentile = correlation_percentile,
      distance_percentile_breaks = distance_bin_info$percentile_breaks
    )
  ))
}


# ============================================================
# CONVENIENCE WRAPPER: Run single animal
# ============================================================

#' Quick wrapper to run pipeline for a single animal
#' 
#' @param animal_id Animal ID (e.g., "mPFCm4")
#' @param ... Additional arguments passed to run_animal_pipeline
#' @return Pipeline results
run_single_animal <- function(animal_id, 
                              data_dir = "E:/Data_Processing/R/Data CSVs",
                              output_dir = "E:/Data_Processing/R/Results",
                              window_size = 3,
                              filter_by_activity = TRUE,
                              file_suffix = "_filtered",
                              ...) {
  
  cat("\n")
  cat("############################################\n")
  cat("# SINGLE ANIMAL ANALYSIS: ", animal_id, "\n")
  cat("############################################\n")
  
  results <- run_animal_pipeline(
    animal_id = animal_id,
    data_dir = data_dir,
    output_dir = output_dir,
    window_size = window_size,
    filter_by_activity = filter_by_activity,
    file_suffix = file_suffix,
    save_intermediate = TRUE,
    ...
  )
  
  cat("\n############################################\n")
  cat("# OUTPUT FILES:\n")
  cat("############################################\n")
  cat("  ", animal_id, "_ROI_filtering_report.csv\n")
  cat("  ", animal_id, "_filtering_summary.csv\n")
  cat("  ", animal_id, "_Trajectory_Overall_", window_size, "ep", file_suffix, ".png\n")
  cat("  ", animal_id, "_Wake2NREM_Trajectory_by_Distance_", window_size, "ep", file_suffix, ".png\n")
  cat("  ", animal_id, "_NREM2Wake_Trajectory_by_Distance_", window_size, "ep", file_suffix, ".png\n")
  cat("  ", animal_id, "_Wake2NREM_EventRate_Trajectory_", window_size, "ep", file_suffix, ".png\n")
  cat("  ", animal_id, "_NREM2Wake_EventRate_Trajectory_", window_size, "ep", file_suffix, ".png\n")
  cat("############################################\n\n")
  
  return(results)
}


# ============================================================
# SCRIPT INFO
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  Pipeline_Transition_Analysis_v3.R\n")
  cat("  (With Activity Filtering & Individual Plots)\n")
  cat("========================================\n")
  cat("\nUsage:\n")
  cat("  # Single animal:\n")
  cat("  results <- run_single_animal('mPFCm4')\n")
  cat("\n")
  cat("  # Full pipeline call:\n")
  cat("  results <- run_animal_pipeline(\n")
  cat("    animal_id = 'mPFCf5',\n")
  cat("    filter_by_activity = TRUE\n")
  cat("  )\n")
  cat("========================================\n\n")
}

