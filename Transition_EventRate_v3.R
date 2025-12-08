# Transition_EventRate_v3.R
# Analyzes event rates in transition windows at the single-cell level
# Compares event rates across epoch positions approaching transitions
# v3: Adds file_suffix parameter for filtered analysis labeling

library(dplyr)
library(ggplot2)
library(tidyr)

# ==================== CONFIGURATION ====================
SAVE_OUTPUT <- TRUE  # Set to TRUE to save results
INPUT_DIR <- "E:/Data_Processing/R/Results"  # Where Extract_Transitions output is
OUTPUT_DIR <- "E:/Data_Processing/R/Results"  # Where to save event rate results
EPOCH_DURATION <- 10  # seconds per epoch

# File labeling (v3)
FILE_SUFFIX <- ""  # Set to "_filtered" when using activity-filtered ROIs

# Input files (from Extract_Transitions.R output)
RECORDING_ID <- "mPFCf5_BL"
TRANSITION_TYPE <- "Wake2NREM"  # or "NREM2Wake"
WINDOW_SIZE <- 3  # or 9
# =======================================================


#' Load events and metadata (from files OR from variables)
load_event_data <- function(recording_id = NULL, transition_type = NULL, window_size = NULL, 
                           input_dir = NULL, events_df = NULL, metadata_df = NULL) {
  
  # If data frames provided directly, use them (and filter if needed)
  if (!is.null(events_df) && !is.null(metadata_df)) {
    cat("Using provided data frames (no file loading)\n")
    
    # Auto-filter if condition parameters are provided
    if (!is.null(transition_type) && !is.null(window_size)) {
      cat("  Filtering for:", transition_type, ",", window_size, "epochs\n")
      
      events_df <- events_df %>%
        filter(transition_type == !!transition_type, window_size == !!window_size)
      
      metadata_df <- metadata_df %>%
        filter(transition_type == !!transition_type, window_size == !!window_size)
      
      cat("  After filtering: ", nrow(events_df), "events,", nrow(metadata_df), "transitions\n")
    } else {
      cat("  No filtering applied - using all conditions\n")
      cat("  Events:", nrow(events_df), "rows\n")
      cat("  Metadata:", nrow(metadata_df), "transitions\n")
    }
    
    return(list(events = events_df, metadata = metadata_df))
  }
  
  # Otherwise load from files (already filtered by filename)
  if (is.null(recording_id) || is.null(transition_type) || is.null(window_size) || is.null(input_dir)) {
    stop("Must provide either (events_df + metadata_df) OR (recording_id + transition_type + window_size + input_dir)")
  }
  
  # Construct filenames
  events_file <- file.path(input_dir, 
                           paste0(recording_id, "_", transition_type, "_", 
                                 window_size, "ep_events.csv"))
  metadata_file <- file.path(input_dir,
                            paste0(recording_id, "_", transition_type, "_", 
                                  window_size, "ep_metadata.csv"))
  
  cat("Loading events from:", events_file, "\n")
  events_df <- read.csv(events_file, stringsAsFactors = FALSE)
  
  cat("Loading metadata from:", metadata_file, "\n")
  metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  cat("Loaded", nrow(events_df), "events across", nrow(metadata_df), "transitions\n")
  
  return(list(events = events_df, metadata = metadata_df))
}


#' Calculate event rate per ROI per epoch position per transition
calculate_event_rates <- function(events_df, metadata_df, epoch_duration = EPOCH_DURATION) {
  cat("\nCalculating event rates...\n")
  
  # Count events per ROI per epoch_in_window per transition
  event_counts <- events_df %>%
    group_by(transition_id, epoch_in_window, Cell_Name) %>%
    summarise(event_count = n(), .groups = "drop")
  
  # Create a complete grid of all combinations
  # (some ROIs might have 0 events in some epochs)
  all_rois <- unique(events_df$Cell_Name)
  all_epochs <- unique(events_df$epoch_in_window)
  all_transitions <- unique(metadata_df$transition_id)
  
  complete_grid <- expand.grid(
    transition_id = all_transitions,
    epoch_in_window = all_epochs,
    Cell_Name = all_rois,
    stringsAsFactors = FALSE
  )
  
  # Merge with actual counts, filling 0 for missing combinations
  event_rates <- complete_grid %>%
    left_join(event_counts, by = c("transition_id", "epoch_in_window", "Cell_Name")) %>%
    mutate(event_count = ifelse(is.na(event_count), 0, event_count)) %>%
    mutate(
      events_per_epoch = event_count,  # Direct count per epoch
      events_per_second = event_count / epoch_duration,
      events_per_minute = event_count / epoch_duration * 60
    )
  
  # Add metadata info
  event_rates <- event_rates %>%
    left_join(metadata_df %>% select(transition_id, transition_type, window_size),
              by = "transition_id")
  
  cat("Calculated event rates for", length(all_rois), "ROIs across", 
      length(all_transitions), "transitions\n")
  
  return(event_rates)
}


#' Compute summary statistics across transitions
summarize_event_rates <- function(event_rates) {
  cat("\nComputing summary statistics...\n")
  
  # Summarize by ROI and epoch position (averaging across transitions)
  roi_summary <- event_rates %>%
    group_by(Cell_Name, epoch_in_window) %>%
    summarise(
      n_transitions = n(),
      mean_events_per_epoch = mean(events_per_epoch),
      sd_events_per_epoch = sd(events_per_epoch),
      se_events_per_epoch = sd(events_per_epoch) / sqrt(n()),
      median_events_per_epoch = median(events_per_epoch),
      .groups = "drop"
    )
  
  # Overall summary by epoch position (averaging across ROIs and transitions)
  epoch_summary <- event_rates %>%
    group_by(epoch_in_window) %>%
    summarise(
      n_observations = n(),
      mean_events_per_epoch = mean(events_per_epoch),
      sd_events_per_epoch = sd(events_per_epoch),
      se_events_per_epoch = sd(events_per_epoch) / sqrt(n()),
      median_events_per_epoch = median(events_per_epoch),
      .groups = "drop"
    )
  
  cat("\nEpoch position summary (events per epoch):\n")
  print(epoch_summary)
  
  return(list(
    roi_summary = roi_summary,
    epoch_summary = epoch_summary
  ))
}


#' Create visualizations of event rates
create_event_rate_plots <- function(event_rates, summary_stats, recording_id, 
                                    transition_type, window_size) {
  cat("\nCreating visualizations...\n")
  
  plots <- list()
  
  # 1. Line plot: Mean event rate by epoch position (averaged across ROIs)
  plots$epoch_trajectory <- ggplot(summary_stats$epoch_summary, 
                                   aes(x = epoch_in_window, y = mean_events_per_epoch)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_line(size = 1.2, color = "blue") +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean_events_per_epoch - se_events_per_epoch,
                     ymax = mean_events_per_epoch + se_events_per_epoch),
                 width = 0.2, color = "blue") +
    labs(
      title = paste(recording_id, "-", transition_type, "- Event Rate Trajectory"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 2. Box plot: Distribution of event rates by epoch position
  plots$epoch_distribution <- ggplot(event_rates, 
                                     aes(x = factor(epoch_in_window), 
                                         y = events_per_epoch)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 0.5) +  # height = 0: no vertical jitter
    labs(
      title = paste(recording_id, "-", transition_type, "- Event Rate Distribution"),
      x = "Epoch Position Relative to Transition",
      y = "Event Rate (events/epoch)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # 3. Heatmap: ROI-specific trajectories
  # Select top variable ROIs for visualization (or all if < 50)
  roi_variability <- event_rates %>%
    group_by(Cell_Name) %>%
    summarise(cv = sd(events_per_epoch) / (mean(events_per_epoch) + 0.001)) %>%
    arrange(desc(cv))
  
  n_rois_to_plot <- min(50, nrow(roi_variability))
  top_rois <- roi_variability$Cell_Name[1:n_rois_to_plot]
  
  event_rates_subset <- event_rates %>%
    filter(Cell_Name %in% top_rois)
  
  # Average across transitions for each ROI and epoch
  roi_epoch_avg <- event_rates_subset %>%
    group_by(Cell_Name, epoch_in_window) %>%
    summarise(mean_events_per_epoch = mean(events_per_epoch), .groups = "drop")
  
  plots$roi_heatmap <- ggplot(roi_epoch_avg, 
                              aes(x = factor(epoch_in_window), 
                                  y = Cell_Name, 
                                  fill = mean_events_per_epoch)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = paste(recording_id, "-", transition_type, "- ROI Event Rate Trajectories"),
      x = "Epoch Position Relative to Transition",
      y = "ROI",
      fill = "Events/epoch"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 10)
    )
  
  # 4. Individual ROI trajectories (spaghetti plots) - split into high-activity and typical
  
  # Identify high-activity ROIs (those with highest mean or peak activity)
  roi_activity_stats <- event_rates %>%
    group_by(Cell_Name) %>%
    summarise(
      mean_activity = mean(events_per_epoch),
      max_activity = max(events_per_epoch),
      .groups = "drop"
    ) %>%
    arrange(desc(max_activity))
  
  # Select up to 10 high-activity ROIs (those with highest peak values)
  n_high_activity <- min(10, nrow(roi_activity_stats))
  high_activity_rois <- roi_activity_stats$Cell_Name[1:n_high_activity]
  
  # Select random sample of remaining ROIs for "typical" plot
  remaining_rois <- setdiff(unique(event_rates$Cell_Name), high_activity_rois)
  n_typical <- min(10, length(remaining_rois))
  typical_rois <- sample(remaining_rois, n_typical)
  
  # Calculate trajectories with variability for high-activity ROIs
  high_activity_trajectories <- event_rates %>%
    filter(Cell_Name %in% high_activity_rois) %>%
    group_by(Cell_Name, epoch_in_window) %>%
    summarise(
      mean_events = mean(events_per_epoch),
      se_events = sd(events_per_epoch) / sqrt(n()),
      .groups = "drop"
    )
  
  plots$roi_trajectories_high <- ggplot(high_activity_trajectories, 
                                        aes(x = epoch_in_window, 
                                            y = mean_events, 
                                            group = Cell_Name, 
                                            color = Cell_Name,
                                            fill = Cell_Name)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_events - se_events, 
                    ymax = mean_events + se_events),
                alpha = 0.2, color = NA) +
    geom_line(alpha = 0.8, size = 0.8) +
    geom_point(alpha = 0.8, size = 1.5) +
    labs(
      title = paste(recording_id, "-", transition_type, "- High-Activity ROI Trajectories"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  # Calculate trajectories with variability for typical ROIs
  typical_trajectories <- event_rates %>%
    filter(Cell_Name %in% typical_rois) %>%
    group_by(Cell_Name, epoch_in_window) %>%
    summarise(
      mean_events = mean(events_per_epoch),
      se_events = sd(events_per_epoch) / sqrt(n()),
      .groups = "drop"
    )
  
  plots$roi_trajectories_typical <- ggplot(typical_trajectories, 
                                           aes(x = epoch_in_window, 
                                               y = mean_events, 
                                               group = Cell_Name, 
                                               color = Cell_Name,
                                               fill = Cell_Name)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "T0", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_events - se_events, 
                    ymax = mean_events + se_events),
                alpha = 0.2, color = NA) +
    geom_line(alpha = 0.8, size = 0.8) +
    geom_point(alpha = 0.8, size = 1.5) +
    labs(
      title = paste(recording_id, "-", transition_type, "- Typical ROI Trajectories"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  cat("Created", length(plots), "plots\n")
  return(plots)
}


#' Statistical testing for changes across epoch positions
test_epoch_differences <- function(event_rates) {
  cat("\nTesting for differences across epoch positions...\n")
  
  # Check for sufficient data
  if (is.null(event_rates) || nrow(event_rates) == 0) {
    cat("  No event data available - skipping statistical tests\n")
    return(list(anova = NULL, pairwise = data.frame()))
  }
  
  n_epochs <- length(unique(event_rates$epoch_in_window))
  if (n_epochs < 2) {
    cat("  Fewer than 2 epochs with data - skipping statistical tests\n")
    return(list(anova = NULL, pairwise = data.frame()))
  }
  
  # ANOVA to test if event rates differ by epoch position
  # Using mixed model approach since we have repeated measures (same ROIs across epochs)
  
  # Simple ANOVA first (not accounting for repeated measures)
  aov_model <- aov(events_per_epoch ~ factor(epoch_in_window), data = event_rates)
  aov_summary <- summary(aov_model)
  
  cat("\nANOVA Results:\n")
  print(aov_summary)
  
  # Pairwise comparisons between consecutive epochs
  epoch_pairs <- combn(sort(unique(event_rates$epoch_in_window)), 2, simplify = FALSE)
  
  pairwise_results <- data.frame()
  
  for (pair in epoch_pairs) {
    epoch1 <- pair[1]
    epoch2 <- pair[2]
    
    data1 <- event_rates %>% filter(epoch_in_window == epoch1) %>% pull(events_per_epoch)
    data2 <- event_rates %>% filter(epoch_in_window == epoch2) %>% pull(events_per_epoch)
    
    # Paired t-test (pairing by ROI and transition)
    # Need to ensure proper pairing
    paired_data <- event_rates %>%
      filter(epoch_in_window %in% c(epoch1, epoch2)) %>%
      select(transition_id, Cell_Name, epoch_in_window, events_per_epoch) %>%
      pivot_wider(names_from = epoch_in_window, 
                  values_from = events_per_epoch,
                  names_prefix = "epoch_")
    
    if (nrow(paired_data) > 0) {
      t_test <- t.test(paired_data[[paste0("epoch_", epoch1)]], 
                       paired_data[[paste0("epoch_", epoch2)]], 
                       paired = TRUE)
      
      pairwise_results <- rbind(pairwise_results, data.frame(
        comparison = paste(epoch1, "vs", epoch2),
        t_statistic = t_test$statistic,
        p_value = t_test$p.value,
        mean_diff = t_test$estimate,
        conf_low = t_test$conf.int[1],
        conf_high = t_test$conf.int[2]
      ))
    }
  }
  
  cat("\nPairwise comparisons (paired t-tests):\n")
  print(pairwise_results)
  
  return(list(
    anova = aov_summary,
    pairwise = pairwise_results
  ))
}


#' Save outputs
save_event_rate_outputs <- function(event_rates, summary_stats, plots, stats_tests,
                                    recording_id, transition_type, window_size, output_dir,
                                    file_suffix = "") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving outputs to:", output_dir, "\n")
  
  # Base filename with suffix
  base_name <- paste0(recording_id, "_", transition_type, "_", window_size, "ep", file_suffix)
  
  # Save detailed event rates
  rates_file <- file.path(output_dir, paste0(base_name, "_event_rates.csv"))
  write.csv(event_rates, rates_file, row.names = FALSE)
  cat("  Saved:", rates_file, "\n")
  
  # Save ROI summary
  roi_summary_file <- file.path(output_dir, paste0(base_name, "_event_rates_roi_summary.csv"))
  write.csv(summary_stats$roi_summary, roi_summary_file, row.names = FALSE)
  cat("  Saved:", roi_summary_file, "\n")
  
  # Save epoch summary
  epoch_summary_file <- file.path(output_dir, paste0(base_name, "_event_rates_epoch_summary.csv"))
  write.csv(summary_stats$epoch_summary, epoch_summary_file, row.names = FALSE)
  cat("  Saved:", epoch_summary_file, "\n")
  
  # Save statistical tests
  stats_file <- file.path(output_dir, paste0(base_name, "_event_rates_stats.csv"))
  write.csv(stats_tests$pairwise, stats_file, row.names = FALSE)
  cat("  Saved:", stats_file, "\n")
  
  # Save plots
  for (plot_name in names(plots)) {
    plot_file <- file.path(output_dir, paste0(base_name, "_event_rates_", plot_name, ".png"))
    ggsave(plot_file, plots[[plot_name]], width = 10, height = 8, dpi = 300)
    cat("  Saved:", plot_file, "\n")
  }
  
  cat("\nAll event rate outputs saved!\n")
}


#' Main analysis function
analyze_transition_event_rates <- function(recording_id = RECORDING_ID,
                                          transition_type = TRANSITION_TYPE,
                                          window_size = WINDOW_SIZE,
                                          input_dir = INPUT_DIR,
                                          events_df = NULL,
                                          metadata_df = NULL,
                                          save_output = SAVE_OUTPUT,
                                          output_dir = OUTPUT_DIR,
                                          file_suffix = FILE_SUFFIX) {
  
  cat("==================================================\n")
  cat("Transition Event Rate Analysis (v3)\n")
  cat("==================================================\n")
  
  # Determine data source
  if (!is.null(events_df) && !is.null(metadata_df)) {
    cat("Data source: Provided data frames (pipeline mode)\n")
    # Extract info from data for labeling
    if (is.null(recording_id)) recording_id <- unique(metadata_df$recording_id)[1]
    if (is.null(transition_type)) transition_type <- unique(metadata_df$transition_type)[1]
    if (is.null(window_size)) window_size <- unique(metadata_df$window_size)[1]
  } else {
    cat("Data source: CSV files (standalone mode)\n")
    cat("Input directory:", input_dir, "\n")
  }
  
  cat("Recording ID:", recording_id, "\n")
  cat("Transition type:", transition_type, "\n")
  cat("Window size:", window_size, "epochs\n")
  cat("Output directory:", output_dir, "\n")
  cat("Save output:", save_output, "\n")
  cat("==================================================\n\n")
  
  # Step 1: Load data (from files or use provided)
  data <- load_event_data(recording_id, transition_type, window_size, input_dir,
                         events_df, metadata_df)
  
  # Step 2: Calculate event rates
  event_rates <- calculate_event_rates(data$events, data$metadata)
  
  # Step 3: Summarize
  summary_stats <- summarize_event_rates(event_rates)
  
  # Step 4: Visualize
  plots <- create_event_rate_plots(event_rates, summary_stats, 
                                   recording_id, transition_type, window_size)
  
  # Step 5: Statistical tests
  stats_tests <- test_epoch_differences(event_rates)
  
  # Step 6: Save if requested
  if (save_output) {
    save_event_rate_outputs(event_rates, summary_stats, plots, stats_tests,
                           recording_id, transition_type, window_size, output_dir,
                           file_suffix)
  } else {
    cat("\nSAVE_OUTPUT = FALSE: Results returned but not saved to disk\n")
  }
  
  # Return all results
  results <- list(
    event_rates = event_rates,
    summary_stats = summary_stats,
    plots = plots,
    stats_tests = stats_tests
  )
  
  cat("\n==================================================\n")
  cat("Event rate analysis complete!\n")
  cat("  Total observations:", nrow(event_rates), "\n")
  cat("  ROIs analyzed:", length(unique(event_rates$Cell_Name)), "\n")
  cat("  Transitions:", length(unique(event_rates$transition_id)), "\n")
  cat("==================================================\n")
  
  return(results)
}


# =========================== MAIN ===========================
# Run the analysis
if (sys.nframe() == 0) {
  results <- analyze_transition_event_rates()
  
  # Display key plots
  print(results$plots$epoch_trajectory)
  print(results$plots$roi_trajectories_high)
}
