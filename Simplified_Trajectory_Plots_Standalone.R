# ============================================================
# SIMPLIFIED TRAJECTORY PLOTS - STANDALONE VERSION
# ============================================================
# Loads pre-computed trajectory CSVs and generates thesis figures
# 
# Figure 3: Correlation Trajectories (BL vs SD, Wake→NREM)
#   - Panel A: Transgenic (mPFCf5, f6, m9)
#   - Panel B: Viral (mPFCm4)
#
# Figure 4: Distance-Binned Correlations (BL vs SD, Wake→NREM)
#   - Four panels, one per animal
#   - Close & Far bins only (Medium removed)
#   - Color = distance bin, Linetype = condition
#
# Figure S1: Event Rate Trajectories (BL vs SD, Wake→NREM)
#   - Panel A: Transgenic (mPFCf5, f6, m9)
#   - Panel B: Viral (mPFCm4)
#   - Includes error bars (SD)
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# CONFIGURATION
# ============================================================

# Animals by expression type
TRANSGENIC_ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm9")
VIRAL_ANIMALS <- c("mPFCm4")
ALL_ANIMALS <- c(TRANSGENIC_ANIMALS, VIRAL_ANIMALS)

# Data paths - EDIT THESE
DATA_DIR <- "E:/Data_Processing/R/Results/within_subject"
WINDOW_SIZE <- 9  # epochs
FILE_SUFFIX <- "_all_simple"  # e.g., "_filtered" or ""

# ============================================================
# DATA LOADING FUNCTIONS
# ============================================================

#' Load trajectory correlation data for one animal
load_trajectory_data <- function(animal_id, data_dir = DATA_DIR, 
                                  window_size = WINDOW_SIZE, 
                                  suffix = FILE_SUFFIX) {
  filename <- paste0(animal_id, "_trajectory_data_", window_size, "ep", suffix, ".csv")
  filepath <- file.path(data_dir, filename)
  
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(NULL)
  }
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  cat("Loaded", nrow(df), "rows from", filename, "\n")
  return(df)
}

#' Load trajectory distance-binned data for one animal
load_trajectory_dist <- function(animal_id, data_dir = DATA_DIR,
                                  window_size = WINDOW_SIZE,
                                  suffix = FILE_SUFFIX) {
  filename <- paste0(animal_id, "_trajectory_dist_", window_size, "ep", suffix, ".csv")
  filepath <- file.path(data_dir, filename)
  
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(NULL)
  }
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  cat("Loaded", nrow(df), "rows from", filename, "\n")
  return(df)
}

#' Load trajectory event rate data for one animal
load_trajectory_events <- function(animal_id, data_dir = DATA_DIR,
                                    window_size = WINDOW_SIZE,
                                    suffix = FILE_SUFFIX) {
  filename <- paste0(animal_id, "_trajectory_eventrates_", window_size, "ep", suffix, ".csv")
  filepath <- file.path(data_dir, filename)
  
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(NULL)
  }
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  cat("Loaded", nrow(df), "rows from", filename, "\n")
  return(df)
}

#' Load all animals' data of a given type
load_all_animals <- function(animals, loader_func, ...) {
  all_data <- list()
  for (animal in animals) {
    df <- loader_func(animal, ...)
    if (!is.null(df)) {
      all_data[[animal]] <- df
    }
  }
  return(all_data)
}

# ============================================================
# SIMPLIFIED PLOTTING FUNCTIONS
# ============================================================

#' Create SIMPLIFIED correlation trajectory plot (BL vs SD only, Wake→NREM)
#' 
#' @param trajectory_data Data frame with correlation trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return ggplot object
create_correlation_plot_simplified <- function(trajectory_data, animal_id, 
                                                window_size = WINDOW_SIZE,
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
                              color = condition, fill = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                    ymax = mean_correlation + se_correlation),
                alpha = 0.4, linewidth = 0, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
    ) +
    scale_fill_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
    ) +
    labs(
      title = paste0(ANIMAL_LABELS[animal_id], " - Pairwise Correlation", filter_label),
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
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    )
  
  return(p)
}


#' Create correlation trajectory plot SPLIT BY CONDITION
#' 
#' Returns two plots: one for BL only, one for SD only. Same design as
#' create_correlation_plot_simplified. Y-axis is forced to span the full BL+SD
#' range so the two split plots share the same scale as each other and as the
#' combined plot.
#'
#' @param trajectory_data Data frame with correlation trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return Named list of ggplot objects: list(BL = p_bl, SD = p_sd)
create_correlation_plot_by_condition_split <- function(trajectory_data, animal_id,
                                                        window_size = WINDOW_SIZE,
                                                        filter_label = "",
                                                        transition_filter = "Wake2NREM") {
  
  # Summarise full combined dataset (BL+SD) to derive shared y-axis limits
  summary_all <- trajectory_data %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      n_pairs = n(),
      .groups = "drop"
    )
  y_min <- min(summary_all$mean_correlation - summary_all$se_correlation, na.rm = TRUE)
  y_max <- max(summary_all$mean_correlation + summary_all$se_correlation, na.rm = TRUE)
  y_pad <- 0.05 * (y_max - y_min)
  ylim_shared <- c(y_min - y_pad, y_max + y_pad)
  
  cond_labels <- c("BL" = "Baseline", "SD" = "Sleep Dep")
  
  plots <- list()
  for (cond in c("BL", "SD")) {
    plot_data <- summary_all %>% filter(condition == cond)
    
    p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                                color = condition, fill = condition, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  alpha = 0.4, linewidth = 0, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_color_manual(
        values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
        guide = "none"
      ) +
      scale_fill_manual(
        values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
        guide = "none"
      ) +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed"),
        guide = "none"
      ) +
      coord_cartesian(ylim = ylim_shared) +
      labs(
        title = paste0(ANIMAL_LABELS[animal_id],
                       " - Pairwise Correlation (", cond_labels[[cond]], ")",
                       filter_label),
        subtitle = paste0("Mean \u00B1 SE pairwise correlation (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Sleep Onset",
        y = "Mean Pairwise Correlation"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.spacing.y = unit(0.1, "cm")
      )
    
    plots[[cond]] <- p
  }
  
  return(plots)
}


#' Create SIMPLIFIED distance-binned plot (BL vs SD, Wake→NREM, Close & Far only)
#' 
#' @param trajectory_dist Data frame with distance-binned trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return ggplot object
create_distance_plot_simplified <- function(trajectory_dist, animal_id, 
                                             window_size = WINDOW_SIZE,
                                             filter_label = "",
                                             transition_filter = "Wake2NREM") {
  
  # Filter to Close and Far only (remove Medium)
  trajectory_dist <- trajectory_dist %>%
    filter(distance_bin_label %in% c("Close", "Far"))
  
  # Compute bin ranges for subtitle
  bin_ranges <- trajectory_dist %>%
    group_by(distance_bin_label) %>%
    summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop") %>%
    arrange(bin_min) %>%
    mutate(range_str = paste0(distance_bin_label, ": ", round(bin_min), "-", round(bin_max), "px")) %>%
    pull(range_str)
  bin_ranges_subtitle <- paste(bin_ranges, collapse = " | ")
  
  # Filter to BL and SD only, and one transition type
  plot_data <- trajectory_dist %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, 
      aes(x = epoch_in_window, y = mean_correlation,
          color = distance_bin_label, fill = distance_bin_label,
          linetype = condition,
          group = interaction(distance_bin_label, condition))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                    ymax = mean_correlation + se_correlation),
                alpha = 0.4, linewidth = 0, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("Close" = "#3A0557", "Far" = "#0F6E56")) +
    scale_fill_manual(values = c("Close" = "#3A0557", "Far" = "#0F6E56"), guide = "none") +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
    ) +
    labs(
      title = paste0(ANIMAL_LABELS[animal_id], " - Distance-Binned Correlation", filter_label),
      subtitle = bin_ranges_subtitle,
      x = "Epoch Position Relative to Sleep Onset",
      y = "Mean Pairwise Correlation",
      color = "Distance Bin",
      linetype = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    )
  
  return(p)
}


#' Create SIMPLIFIED event rate plot (BL vs SD only, Wake→NREM)
#' 
#' @param trajectory_events Data frame with event rate trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return ggplot object
create_event_rate_plot_simplified <- function(trajectory_events, animal_id, 
                                               window_size = WINDOW_SIZE,
                                               filter_label = "",
                                               transition_filter = "Wake2NREM") {
  
  if (is.null(trajectory_events) || nrow(trajectory_events) == 0) {
    return(NULL)
  }
  
  # Detect event rate column name
  event_col <- if ("events_per_epoch" %in% colnames(trajectory_events)) {
    "events_per_epoch"
  } else if ("event_rate" %in% colnames(trajectory_events)) {
    "event_rate"
  } else if ("mean_events" %in% colnames(trajectory_events)) {
    "mean_events"
  } else {
    stop("Could not find event rate column. Available: ", 
         paste(colnames(trajectory_events), collapse = ", "))
  }
  
  # Filter to BL and SD only, and one transition type
  plot_data <- trajectory_events %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(.data[[event_col]], na.rm = TRUE),
      se_event_rate = sd(.data[[event_col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, 
      aes(x = epoch_in_window, y = mean_event_rate,
          color = condition, fill = condition, linetype = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
             color = "gray40", fontface = "bold", size = 3) +
    geom_ribbon(aes(ymin = mean_event_rate - se_event_rate,
                    ymax = mean_event_rate + se_event_rate),
                alpha = 0.4, linewidth = 0, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
    ) +
    scale_fill_manual(
      values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("BL" = "solid", "SD" = "dashed"),
      labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
    ) +
    labs(
      title = paste0(ANIMAL_LABELS[animal_id], " - Event Rate", filter_label),
      subtitle = paste0("Mean ± SE event rate (", window_size, " epoch window)"),
      x = "Epoch Position Relative to Sleep Onset",
      y = "Mean Event Rate (events/epoch)",
      color = "Condition",
      linetype = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    )
  
  return(p)
}


#' Create distance-binned correlation plot SPLIT BY DISTANCE BIN
#' 
#' Returns two plots: one for Close pairs only, one for Far pairs only. Each shows
#' BL and SD overlaid (BL solid, SD dashed) for a single bin. Same design as
#' create_distance_plot_simplified. Y-axis is forced to span the full BL+SD ×
#' Close+Far range so the two split plots share the same scale as each other
#' and as the combined plot.
#'
#' @param trajectory_dist Data frame with distance-binned trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return Named list of ggplot objects: list(Close = p_close, Far = p_far)
create_distance_plot_by_bin_split <- function(trajectory_dist, animal_id,
                                               window_size = WINDOW_SIZE,
                                               filter_label = "",
                                               transition_filter = "Wake2NREM") {
  
  # Filter to Close and Far only (matches combined plot)
  trajectory_dist <- trajectory_dist %>%
    filter(distance_bin_label %in% c("Close", "Far"))
  
  # Per-bin ranges for subtitles
  bin_ranges <- trajectory_dist %>%
    group_by(distance_bin_label) %>%
    summarise(bin_min = min(Distance), bin_max = max(Distance), .groups = "drop")
  
  # Summarise full combined dataset (BL+SD × Close+Far) to derive shared y-axis limits
  summary_all <- trajectory_dist %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, distance_bin_label, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  y_min <- min(summary_all$mean_correlation - summary_all$se_correlation, na.rm = TRUE)
  y_max <- max(summary_all$mean_correlation + summary_all$se_correlation, na.rm = TRUE)
  y_pad <- 0.05 * (y_max - y_min)
  ylim_shared <- c(y_min - y_pad, y_max + y_pad)
  
  plots <- list()
  for (bin in c("Close", "Far")) {
    br <- bin_ranges %>% filter(distance_bin_label == bin)
    bin_subtitle <- paste0(bin, ": ", round(br$bin_min), "-", round(br$bin_max), "px")
    
    plot_data <- summary_all %>% filter(distance_bin_label == bin)
    
    p <- ggplot(plot_data, 
        aes(x = epoch_in_window, y = mean_correlation,
            color = distance_bin_label, fill = distance_bin_label,
            linetype = condition,
            group = interaction(distance_bin_label, condition))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  alpha = 0.4, linewidth = 0, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c("Close" = "#3A0557", "Far" = "#0F6E56"),
                         guide = "none") +
      scale_fill_manual(values = c("Close" = "#3A0557", "Far" = "#0F6E56"),
                        guide = "none") +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed"),
        labels = c("BL" = "Baseline", "SD" = "Sleep Dep")
      ) +
      coord_cartesian(ylim = ylim_shared) +
      labs(
        title = paste0(ANIMAL_LABELS[animal_id],
                       " - Distance-Binned Correlation (", bin, " Pairs)",
                       filter_label),
        subtitle = bin_subtitle,
        x = "Epoch Position Relative to Sleep Onset",
        y = "Mean Pairwise Correlation",
        linetype = "Condition"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.spacing.y = unit(0.1, "cm")
      )
    
    plots[[bin]] <- p
  }
  
  return(plots)
}


#' Create event rate trajectory plot SPLIT BY CONDITION
#' 
#' Returns two plots: one for BL only, one for SD only. Same design as
#' create_event_rate_plot_simplified. Y-axis is forced to span the full BL+SD
#' range so the two split plots share the same scale as each other and as the
#' combined plot.
#'
#' @param trajectory_events Data frame with event rate trajectory data
#' @param animal_id Animal identifier for title
#' @param window_size Window size in epochs
#' @param filter_label Optional label suffix for title
#' @param transition_filter Transition type to show (default: "Wake2NREM")
#' @return Named list of ggplot objects: list(BL = p_bl, SD = p_sd)
create_event_rate_plot_by_condition_split <- function(trajectory_events, animal_id,
                                                       window_size = WINDOW_SIZE,
                                                       filter_label = "",
                                                       transition_filter = "Wake2NREM") {
  
  if (is.null(trajectory_events) || nrow(trajectory_events) == 0) {
    return(NULL)
  }
  
  # Detect event rate column name
  event_col <- if ("events_per_epoch" %in% colnames(trajectory_events)) {
    "events_per_epoch"
  } else if ("event_rate" %in% colnames(trajectory_events)) {
    "event_rate"
  } else if ("mean_events" %in% colnames(trajectory_events)) {
    "mean_events"
  } else {
    stop("Could not find event rate column. Available: ", 
         paste(colnames(trajectory_events), collapse = ", "))
  }
  
  # Summarise full combined dataset (BL+SD) to derive shared y-axis limits
  summary_all <- trajectory_events %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(.data[[event_col]], na.rm = TRUE),
      se_event_rate = sd(.data[[event_col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  y_min <- min(summary_all$mean_event_rate - summary_all$se_event_rate, na.rm = TRUE)
  y_max <- max(summary_all$mean_event_rate + summary_all$se_event_rate, na.rm = TRUE)
  y_pad <- 0.05 * (y_max - y_min)
  ylim_shared <- c(y_min - y_pad, y_max + y_pad)
  
  cond_labels <- c("BL" = "Baseline", "SD" = "Sleep Dep")
  
  plots <- list()
  for (cond in c("BL", "SD")) {
    plot_data <- summary_all %>% filter(condition == cond)
    
    p <- ggplot(plot_data, 
        aes(x = epoch_in_window, y = mean_event_rate,
            color = condition, fill = condition, linetype = condition)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      annotate("text", x = 0, y = Inf, label = "Sleep Onset", vjust = 2, hjust = 0.5, 
               color = "gray40", fontface = "bold", size = 3) +
      geom_ribbon(aes(ymin = mean_event_rate - se_event_rate,
                      ymax = mean_event_rate + se_event_rate),
                  alpha = 0.4, linewidth = 0, color = NA) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_color_manual(
        values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
        guide = "none"
      ) +
      scale_fill_manual(
        values = c("BL" = "#2E86AB", "SD" = "#A23B72"),
        guide = "none"
      ) +
      scale_linetype_manual(
        values = c("BL" = "solid", "SD" = "dashed"),
        guide = "none"
      ) +
      coord_cartesian(ylim = ylim_shared) +
      labs(
        title = paste0(ANIMAL_LABELS[animal_id],
                       " - Event Rate (", cond_labels[[cond]], ")",
                       filter_label),
        subtitle = paste0("Mean \u00B1 SE event rate (", window_size, " epoch window)"),
        x = "Epoch Position Relative to Sleep Onset",
        y = "Mean Event Rate (events/epoch)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.spacing.y = unit(0.1, "cm")
      )
    
    plots[[cond]] <- p
  }
  
  return(plots)
}


# ============================================================
# STANDALONE WRAPPER FUNCTIONS
# ============================================================

#' Generate correlation trajectory plot from CSV
#' 
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return ggplot object
plot_correlation_standalone <- function(animal_id, 
                                         data_dir = DATA_DIR,
                                         window_size = WINDOW_SIZE,
                                         suffix = FILE_SUFFIX) {
  trajectory_data <- load_trajectory_data(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_data)) return(NULL)
  create_correlation_plot_simplified(trajectory_data, animal_id, window_size)
}

#' Generate correlation plot SPLIT BY CONDITION (BL and SD as separate plots)
#'
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return Named list: list(BL = ggplot, SD = ggplot)
plot_correlation_by_condition_split_standalone <- function(animal_id,
                                                            data_dir = DATA_DIR,
                                                            window_size = WINDOW_SIZE,
                                                            suffix = FILE_SUFFIX) {
  trajectory_data <- load_trajectory_data(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_data)) return(NULL)
  create_correlation_plot_by_condition_split(trajectory_data, animal_id, window_size)
}

#' Generate distance-binned plot from CSV
#' 
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return ggplot object
plot_distance_standalone <- function(animal_id,
                                      data_dir = DATA_DIR,
                                      window_size = WINDOW_SIZE,
                                      suffix = FILE_SUFFIX) {
  trajectory_dist <- load_trajectory_dist(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_dist)) return(NULL)
  create_distance_plot_simplified(trajectory_dist, animal_id, window_size)
}

#' Generate distance-binned plot SPLIT BY DISTANCE BIN (Close and Far as separate plots)
#'
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return Named list: list(Close = ggplot, Far = ggplot)
plot_distance_by_bin_split_standalone <- function(animal_id,
                                                   data_dir = DATA_DIR,
                                                   window_size = WINDOW_SIZE,
                                                   suffix = FILE_SUFFIX) {
  trajectory_dist <- load_trajectory_dist(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_dist)) return(NULL)
  create_distance_plot_by_bin_split(trajectory_dist, animal_id, window_size)
}

#' Generate event rate plot from CSV
#' 
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return ggplot object
plot_event_rate_standalone <- function(animal_id,
                                        data_dir = DATA_DIR,
                                        window_size = WINDOW_SIZE,
                                        suffix = FILE_SUFFIX) {
  trajectory_events <- load_trajectory_events(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_events)) return(NULL)
  create_event_rate_plot_simplified(trajectory_events, animal_id, window_size)
}

#' Generate event rate plot SPLIT BY CONDITION (BL and SD as separate plots)
#'
#' @param animal_id Animal identifier (e.g., "mPFCf5")
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix (e.g., "_filtered")
#' @return Named list: list(BL = ggplot, SD = ggplot)
plot_event_rate_by_condition_split_standalone <- function(animal_id,
                                                           data_dir = DATA_DIR,
                                                           window_size = WINDOW_SIZE,
                                                           suffix = FILE_SUFFIX) {
  trajectory_events <- load_trajectory_events(animal_id, data_dir, window_size, suffix)
  if (is.null(trajectory_events)) return(NULL)
  create_event_rate_plot_by_condition_split(trajectory_events, animal_id, window_size)
}


# ============================================================
# BATCH PLOTTING - GENERATE ALL THESIS FIGURES
# ============================================================

#' Generate all Figure 3 panels (Correlation Trajectories)
#' 
#' Returns a list of plots:
#'   - One per transgenic animal (for Panel A assembly)
#'   - mPFCm4 (for Panel B)
generate_figure3_panels <- function(data_dir = DATA_DIR,
                                     window_size = WINDOW_SIZE,
                                     suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE 3: Correlation Trajectories ===\n")
  cat("Transgenic (Panel A):", paste(TRANSGENIC_ANIMALS, collapse = ", "), "\n")
  cat("Viral (Panel B):", paste(VIRAL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p <- plot_correlation_standalone(animal, data_dir, window_size, suffix)
    if (!is.null(p)) {
      plots[[animal]] <- p
    }
  }
  
  cat("\nGenerated", length(plots), "plots\n")
  return(plots)
}

#' Generate Figure 3 panels SPLIT BY CONDITION (BL and SD as separate plots)
#'
#' Returns a nested list: plots[[animal_id]][[condition]] where condition is
#' "BL" or "SD". Y-axis is forced to the combined BL+SD range for direct
#' comparability across the split panels.
generate_figure3_by_condition_panels <- function(data_dir = DATA_DIR,
                                                  window_size = WINDOW_SIZE,
                                                  suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE 3 (split by condition): Correlation Trajectories ===\n")
  cat("Each animal yields separate Baseline and Sleep Dep plots\n")
  cat("Animals:", paste(ALL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p_list <- plot_correlation_by_condition_split_standalone(animal, data_dir,
                                                              window_size, suffix)
    if (!is.null(p_list)) {
      plots[[animal]] <- p_list
    }
  }
  
  cat("\nGenerated", length(plots), "animals x 2 conditions\n")
  return(plots)
}

#' Generate all Figure 4 panels (Distance-Binned Correlations)
#' 
#' Returns a list with one plot per animal
generate_figure4_panels <- function(data_dir = DATA_DIR,
                                     window_size = WINDOW_SIZE,
                                     suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE 4: Distance-Binned Correlations ===\n")
  cat("Animals:", paste(ALL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p <- plot_distance_standalone(animal, data_dir, window_size, suffix)
    if (!is.null(p)) {
      plots[[animal]] <- p
    }
  }
  
  cat("\nGenerated", length(plots), "plots\n")
  return(plots)
}

#' Generate Figure 4 panels SPLIT BY DISTANCE BIN (Close and Far as separate plots)
#'
#' Returns a nested list: plots[[animal_id]][[bin]] where bin is "Close" or "Far".
#' Each plot shows BL + SD overlaid (linetype distinction) for a single bin.
#' Y-axis is forced to the combined BL+SD × Close+Far range for direct
#' comparability across the split panels.
generate_figure4_by_bin_panels <- function(data_dir = DATA_DIR,
                                            window_size = WINDOW_SIZE,
                                            suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE 4 (split by distance bin): Distance-Binned Correlations ===\n")
  cat("Each animal yields separate Close and Far plots\n")
  cat("Animals:", paste(ALL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p_list <- plot_distance_by_bin_split_standalone(animal, data_dir,
                                                     window_size, suffix)
    if (!is.null(p_list)) {
      plots[[animal]] <- p_list
    }
  }
  
  cat("\nGenerated", length(plots), "animals x 2 bins\n")
  return(plots)
}

#' Generate all Figure S1 panels (Event Rate Trajectories)
#' 
#' Returns a list of plots:
#'   - One per transgenic animal (for Panel A assembly)
#'   - mPFCm4 (for Panel B)
generate_figureS1_panels <- function(data_dir = DATA_DIR,
                                      window_size = WINDOW_SIZE,
                                      suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE S1: Event Rate Trajectories ===\n")
  cat("Transgenic (Panel A):", paste(TRANSGENIC_ANIMALS, collapse = ", "), "\n")
  cat("Viral (Panel B):", paste(VIRAL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p <- plot_event_rate_standalone(animal, data_dir, window_size, suffix)
    if (!is.null(p)) {
      plots[[animal]] <- p
    }
  }
  
  cat("\nGenerated", length(plots), "plots\n")
  return(plots)
}


#' Generate Figure S1 panels SPLIT BY CONDITION (BL and SD as separate plots)
#'
#' Returns a nested list: plots[[animal_id]][[condition]] where condition is
#' "BL" or "SD". Y-axis is forced to the combined BL+SD range for direct
#' comparability across the split panels.
generate_figureS1_by_condition_panels <- function(data_dir = DATA_DIR,
                                                   window_size = WINDOW_SIZE,
                                                   suffix = FILE_SUFFIX) {
  
  cat("\n=== FIGURE S1 (split by condition): Event Rate Trajectories ===\n")
  cat("Each animal yields separate Baseline and Sleep Dep plots\n")
  cat("Animals:", paste(ALL_ANIMALS, collapse = ", "), "\n\n")
  
  plots <- list()
  for (animal in ALL_ANIMALS) {
    cat("Processing", animal, "...\n")
    p_list <- plot_event_rate_by_condition_split_standalone(animal, data_dir,
                                                             window_size, suffix)
    if (!is.null(p_list)) {
      plots[[animal]] <- p_list
    }
  }
  
  cat("\nGenerated", length(plots), "animals x 2 conditions\n")
  return(plots)
}


# ============================================================
# CONDITION-FACETED PLOTS (ALL ANIMALS, FACET BY CONDITION)
# ============================================================
# These plots show all animals on the same axes, faceted by
# condition (BL, SD). This layout makes between-animal variance
# visible from vertical spread, within-animal noise from error
# bars, and condition effects from panel-to-panel comparison.
# ============================================================

# Shared animal color palette (matches cross-animal plots)
ANIMAL_COLORS <- c(
  "mPFCf5" = "#7B2D8E",   # purple
  "mPFCf6" = "#1B9E77",   # teal
  "mPFCm4" = "#D92B2B",   # red (viral)
  "mPFCm9" = "#D4A017"    # gold
)

ANIMAL_LABELS <- c(
  "mPFCf5" = "Female 1",
  "mPFCf6" = "Female 2",
  "mPFCm4" = "Male 1",
  "mPFCm9" = "Male 2"
)

#' Combine multiple animals' data into one data frame with animal_id column
#'
#' @param animals Character vector of animal IDs
#' @param loader_func Function to load one animal's data
#' @param ... Additional arguments passed to loader_func
#' @return Combined data frame with animal_id column, or NULL if no data loaded
combine_animal_data <- function(animals, loader_func, ...) {
  combined <- list()
  for (animal in animals) {
    df <- loader_func(animal, ...)
    if (!is.null(df)) {
      df$animal_id <- animal
      combined[[animal]] <- df
    }
  }
  if (length(combined) == 0) return(NULL)
  bind_rows(combined)
}


#' Correlation trajectory faceted by condition, colored by animal
#'
#' @param combined_data Combined data frame from combine_animal_data()
#' @param window_size Window size in epochs
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @param title_suffix Optional title suffix
#' @return ggplot object
create_correlation_plot_by_condition <- function(combined_data,
                                                  window_size = WINDOW_SIZE,
                                                  transition_filter = "Wake2NREM",
                                                  title_suffix = "") {
  
  plot_data <- combined_data %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(animal_id, condition, epoch_in_window) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(condition = factor(condition, levels = c("BL", "SD"),
                              labels = c("Baseline", "Sleep Deprivation")))
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                              color = animal_id, fill = animal_id,
                              linetype = animal_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                    ymax = mean_correlation + se_correlation),
                alpha = 0.4, linewidth = 0, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ condition, nrow = 1) +
    scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, name = "Animal") +
    scale_fill_manual(values = ANIMAL_COLORS, guide = "none") +
    scale_linetype_manual(
      values = c("mPFCf5" = "solid", "mPFCf6" = "solid",
                 "mPFCm4" = "dashed", "mPFCm9" = "solid"),
      labels = ANIMAL_LABELS,
      name = "Animal"
    ) +
    labs(
      title = paste0("Pairwise Correlation by Condition", title_suffix),
      subtitle = paste0("Shaded regions = \u00B1 SE (", window_size, " epoch window, ",
                        gsub("2", " \u2192 ", transition_filter), ")"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}


#' Event rate trajectory faceted by condition, colored by animal
#'
#' @param combined_data Combined data frame from combine_animal_data()
#' @param window_size Window size in epochs
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @param title_suffix Optional title suffix
#' @return ggplot object
create_event_rate_plot_by_condition <- function(combined_data,
                                                 window_size = WINDOW_SIZE,
                                                 transition_filter = "Wake2NREM",
                                                 title_suffix = "") {
  
  # Detect event rate column name
  event_col <- if ("events_per_epoch" %in% colnames(combined_data)) {
    "events_per_epoch"
  } else if ("event_rate" %in% colnames(combined_data)) {
    "event_rate"
  } else if ("mean_events" %in% colnames(combined_data)) {
    "mean_events"
  } else {
    stop("Could not find event rate column. Available: ",
         paste(colnames(combined_data), collapse = ", "))
  }
  
  plot_data <- combined_data %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(animal_id, condition, epoch_in_window) %>%
    summarise(
      mean_event_rate = mean(.data[[event_col]], na.rm = TRUE),
      se_event_rate = sd(.data[[event_col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(condition = factor(condition, levels = c("BL", "SD"),
                              labels = c("Baseline", "Sleep Deprivation")))
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_event_rate,
                              color = animal_id, fill = animal_id,
                              linetype = animal_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(aes(ymin = mean_event_rate - se_event_rate,
                    ymax = mean_event_rate + se_event_rate),
                alpha = 0.4, linewidth = 0, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ condition, nrow = 1) +
    scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, name = "Animal") +
    scale_fill_manual(values = ANIMAL_COLORS, guide = "none") +
    scale_linetype_manual(
      values = c("mPFCf5" = "solid", "mPFCf6" = "solid",
                 "mPFCm4" = "dashed", "mPFCm9" = "solid"),
      labels = ANIMAL_LABELS,
      name = "Animal"
    ) +
    labs(
      title = paste0("Event Rate by Condition", title_suffix),
      subtitle = paste0("Shaded regions = \u00B1 SE (", window_size, " epoch window, ",
                        gsub("2", " \u2192 ", transition_filter), ")"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Event Rate (events/epoch)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}


#' Distance-binned correlation trajectory faceted by condition, colored by animal
#' Returns one plot per distance bin (Close, Far)
#'
#' @param combined_data Combined data frame from combine_animal_data()
#' @param window_size Window size in epochs
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @param distance_bins Which bins to plot (default: c("Close", "Medium", "Far"))
#' @return Named list of ggplot objects, one per distance bin
create_distance_plot_by_condition <- function(combined_data,
                                               window_size = WINDOW_SIZE,
                                               transition_filter = "Wake2NREM",
                                               distance_bins = c("Close", "Medium", "Far")) {
  
  plots <- list()
  
  for (bin in distance_bins) {
    
    # Compute bin range for subtitle
    bin_range <- combined_data %>%
      filter(distance_bin_label == bin) %>%
      summarise(bin_min = round(min(Distance, na.rm = TRUE)),
                bin_max = round(max(Distance, na.rm = TRUE))) %>%
      mutate(range_str = paste0(bin, " pairs: ", bin_min, "-", bin_max, " px"))
    
    plot_data <- combined_data %>%
      filter(distance_bin_label == bin) %>%
      filter(condition %in% c("BL", "SD")) %>%
      filter(transition_type == transition_filter) %>%
      group_by(animal_id, condition, epoch_in_window) %>%
      summarise(
        mean_correlation = mean(correlation, na.rm = TRUE),
        se_correlation = sd(correlation, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      mutate(condition = factor(condition, levels = c("BL", "SD"),
                                labels = c("Baseline", "Sleep Deprivation")))
    
    p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_correlation,
                                color = animal_id, fill = animal_id,
                                linetype = animal_id)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_ribbon(aes(ymin = mean_correlation - se_correlation,
                      ymax = mean_correlation + se_correlation),
                  alpha = 0.4, linewidth = 0, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~ condition, nrow = 1) +
      scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, name = "Animal") +
      scale_fill_manual(values = ANIMAL_COLORS, guide = "none") +
      scale_linetype_manual(
        values = c("mPFCf5" = "solid", "mPFCf6" = "solid",
                   "mPFCm4" = "dashed", "mPFCm9" = "solid"),
        labels = ANIMAL_LABELS,
        name = "Animal"
      ) +
      labs(
        title = paste0("Distance-Binned Correlation by Condition \u2014 ", bin, " Pairs"),
        subtitle = paste0(bin_range$range_str, " | Shaded regions = \u00B1 SE (", window_size,
                          " epoch window, ", gsub("2", " \u2192 ", transition_filter), ")"),
        x = "Epoch Position Relative to Transition",
        y = "Mean Pairwise Correlation"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    plots[[bin]] <- p
  }
  
  return(plots)
}


# ============================================================
# VARIANCE MAGNITUDE PLOTS (CV as Y-axis)
# ============================================================
# These plots show the coefficient of variation (SD/|mean|) as
# the outcome variable, making variance differences across
# distance bins directly visible — medium bins should sit
# lowest, confirming they carry the cleanest signal.
# ============================================================

#' Plot CV by epoch, colored by animal, faceted by condition
#' Returns one plot per distance bin (transgenic animals only by default)
#'
#' @param combined_data Combined data frame from combine_animal_data()
#' @param window_size Window size in epochs
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @param distance_bins Which bins to plot (default: c("Close", "Medium", "Far"))
#' @return Named list of ggplot objects, one per distance bin
create_variance_plot_by_condition <- function(combined_data,
                                               window_size = WINDOW_SIZE,
                                               transition_filter = "Wake2NREM",
                                               distance_bins = c("Close", "Medium", "Far"),
                                               animal_filter = TRANSGENIC_ANIMALS) {
  
  plots <- list()
  
  for (bin in distance_bins) {
    
    plot_data <- combined_data %>%
      filter(animal_id %in% animal_filter) %>%
      filter(distance_bin_label == bin) %>%
      filter(condition %in% c("BL", "SD")) %>%
      filter(transition_type == transition_filter) %>%
      group_by(animal_id, condition, epoch_in_window) %>%
      summarise(
        mean_corr = mean(correlation, na.rm = TRUE),
        sd_corr = sd(correlation, na.rm = TRUE),
        cv_correlation = sd_corr / abs(mean_corr),
        .groups = "drop"
      ) %>%
      mutate(condition = factor(condition, levels = c("BL", "SD"),
                                labels = c("Baseline", "Sleep Deprivation")))
    
    p <- ggplot(plot_data, aes(x = epoch_in_window, y = cv_correlation,
                                color = animal_id, linetype = animal_id)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~ condition, nrow = 1) +
      scale_color_manual(values = ANIMAL_COLORS, labels = ANIMAL_LABELS, name = "Animal") +
      scale_linetype_manual(
        values = c("mPFCf5" = "solid", "mPFCf6" = "solid",
                   "mPFCm4" = "dashed", "mPFCm9" = "solid"),
        labels = ANIMAL_LABELS,
        name = "Animal"
      ) +
      labs(
        title = paste0("Epoch-to-Epoch CV \u2014 ", bin, " Pairs"),
        subtitle = paste0("Coefficient of variation of pairwise correlation (", window_size,
                          " epoch window, ", gsub("2", " \u2192 ", transition_filter), ")"),
        x = "Epoch Position Relative to Transition",
        y = "CV of Pairwise Correlation"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    plots[[bin]] <- p
  }
  
  return(plots)
}


#' Plot CV with all distance bins on one plot, faceted by condition
#' Shows bin-to-bin variance differences directly
#'
#' @param combined_data Combined data frame from combine_animal_data()
#' @param animal_filter Character vector of animals to include
#' @param window_size Window size in epochs
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @return ggplot object
create_variance_plot_by_bin <- function(combined_data,
                                         animal_filter = TRANSGENIC_ANIMALS,
                                         window_size = WINDOW_SIZE,
                                         transition_filter = "Wake2NREM") {
  
  plot_data <- combined_data %>%
    filter(animal_id %in% animal_filter) %>%
    filter(distance_bin_label %in% c("Close", "Medium", "Far")) %>%
    filter(condition %in% c("BL", "SD")) %>%
    filter(transition_type == transition_filter) %>%
    group_by(distance_bin_label, condition, epoch_in_window) %>%
    summarise(
      mean_corr = mean(correlation, na.rm = TRUE),
      sd_corr = sd(correlation, na.rm = TRUE),
      mean_cv = sd_corr / abs(mean_corr),
      .groups = "drop"
    ) %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD"),
                          labels = c("Baseline", "Sleep Deprivation")),
      distance_bin_label = factor(distance_bin_label,
                                   levels = c("Close", "Medium", "Far"))
    )
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_cv,
                              color = distance_bin_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0) +
    facet_wrap(~ condition, nrow = 1) +
    scale_color_manual(
      values = c("Close" = "#3A0557", "Medium" = "#185FA5", "Far" = "#0F6E56"),
      name = "Distance Bin"
    ) +
    labs(
      title = "Epoch-to-Epoch CV by Distance Bin",
      subtitle = paste0("CV across transgenic animals (", window_size,
                        " epoch window, ", gsub("2", " \u2192 ", transition_filter), ")"),
      x = "Epoch Position Relative to Transition",
      y = "CV of Pairwise Correlation"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}


# ============================================================
# BATCH GENERATION - CONDITION-FACETED FIGURES
# ============================================================

#' Generate all condition-faceted presentation figures
#'
#' @param data_dir Directory containing trajectory CSV files
#' @param window_size Window size in epochs
#' @param suffix File suffix
#' @param transition_filter Transition type (default: "Wake2NREM")
#' @return Named list: $correlation, $event_rate, $distance_close/medium/far, $variance_close/medium/far, $variance_by_bin
generate_condition_faceted_figures <- function(data_dir = DATA_DIR,
                                                window_size = WINDOW_SIZE,
                                                suffix = FILE_SUFFIX,
                                                transition_filter = "Wake2NREM") {
  
  cat("\n=== CONDITION-FACETED PRESENTATION FIGURES ===\n")
  cat("Transition:", transition_filter, "\n")
  cat("Animals:", paste(ALL_ANIMALS, collapse = ", "), "\n\n")
  
  results <- list()
  
  # --- Correlation ---
  cat("Loading correlation data...\n")
  corr_data <- combine_animal_data(ALL_ANIMALS, load_trajectory_data,
                                    data_dir = data_dir, window_size = window_size,
                                    suffix = suffix)
  if (!is.null(corr_data)) {
    results$correlation <- create_correlation_plot_by_condition(
      corr_data, window_size, transition_filter)
    cat("  Correlation plot: done\n")
  }
  
  # --- Event Rate ---
  cat("Loading event rate data...\n")
  event_data <- combine_animal_data(ALL_ANIMALS, load_trajectory_events,
                                     data_dir = data_dir, window_size = window_size,
                                     suffix = suffix)
  if (!is.null(event_data)) {
    results$event_rate <- create_event_rate_plot_by_condition(
      event_data, window_size, transition_filter)
    cat("  Event rate plot: done\n")
  }
  
  # --- Distance-Binned ---
  cat("Loading distance-binned data...\n")
  dist_data <- combine_animal_data(ALL_ANIMALS, load_trajectory_dist,
                                    data_dir = data_dir, window_size = window_size,
                                    suffix = suffix)
  if (!is.null(dist_data)) {
    dist_plots <- create_distance_plot_by_condition(
      dist_data, window_size, transition_filter)
    results$distance_close <- dist_plots$Close
    results$distance_medium <- dist_plots$Medium
    results$distance_far <- dist_plots$Far
    cat("  Distance plots (Close, Medium, Far): done\n")
    
    # --- Variance magnitude ---
    var_plots <- create_variance_plot_by_condition(
      dist_data, window_size, transition_filter)
    results$variance_close <- var_plots$Close
    results$variance_medium <- var_plots$Medium
    results$variance_far <- var_plots$Far
    results$variance_by_bin <- create_variance_plot_by_bin(
      dist_data, TRANSGENIC_ANIMALS, window_size, transition_filter)
    cat("  Variance plots (Close, Medium, Far, combined): done\n")
  }
  
  cat("\nGenerated", length(results), "condition-faceted plots\n")
  return(results)
}


# ============================================================
# USAGE EXAMPLES
# ============================================================

# --- Single animal, view in RStudio ---
# p <- plot_correlation_standalone("mPFCf5")
# print(p)  # View in Plots pane, export manually

# --- All Figure 3 panels ---
# fig3_plots <- generate_figure3_panels()
# print(fig3_plots$mPFCf5)  # View transgenic
# print(fig3_plots$mPFCm4)  # View viral

# --- All Figure 4 panels ---
# fig4_plots <- generate_figure4_panels()
# print(fig4_plots$mPFCf5)

# --- All Figure S1 panels ---
# figS1_plots <- generate_figureS1_panels()
# print(figS1_plots$mPFCf5)

# --- With custom paths ---
# p <- plot_correlation_standalone("mPFCf5", 
#                                   data_dir = "D:/MyData/Results/within_subject",
#                                   window_size = 9,
#                                   suffix = "")

# --- Condition-faceted plots (for presentation) ---
# All animals on same axes, faceted by BL vs SD
# pres_plots <- generate_condition_faceted_figures()
# print(pres_plots$correlation)      # All animals, correlation
# print(pres_plots$event_rate)       # All animals, event rate
# print(pres_plots$distance_close)   # All animals, close pairs
# print(pres_plots$distance_medium)  # All animals, medium pairs
# print(pres_plots$distance_far)     # All animals, far pairs
# print(pres_plots$variance_close)   # CV, close pairs
# print(pres_plots$variance_medium)  # CV, medium pairs
# print(pres_plots$variance_far)     # CV, far pairs
# print(pres_plots$variance_by_bin)  # All bins on one plot (transgenic only)

# --- Single condition-faceted plot with custom transition ---
# corr_data <- combine_animal_data(ALL_ANIMALS, load_trajectory_data)
# p <- create_correlation_plot_by_condition(corr_data, 
#        transition_filter = "NREM2Wake")

cat("\n=== Simplified Trajectory Plots - Standalone ===\n")
cat("Edit DATA_DIR, WINDOW_SIZE, and FILE_SUFFIX at top of script\n\n")

# ============================================================
# AUTO-RUN: Generate and print all figures
# ============================================================

cat("--- Generating per-animal figures ---\n")
fig3_plots <- generate_figure3_panels()
fig3_by_cond_plots <- generate_figure3_by_condition_panels()
fig4_plots <- generate_figure4_panels()
fig4_by_bin_plots <- generate_figure4_by_bin_panels()
figS1_plots <- generate_figureS1_panels()
figS1_by_cond_plots <- generate_figureS1_by_condition_panels()

cat("\n--- Printing per-animal figures ---\n")
for (name in names(fig3_plots)) {
  cat("  Figure 3:", name, "\n")
  print(fig3_plots[[name]])
}
for (animal in names(fig3_by_cond_plots)) {
  for (cond in names(fig3_by_cond_plots[[animal]])) {
    cat("  Figure 3 [", cond, "]:", animal, "\n")
    print(fig3_by_cond_plots[[animal]][[cond]])
  }
}
for (name in names(fig4_plots)) {
  cat("  Figure 4:", name, "\n")
  print(fig4_plots[[name]])
}
for (animal in names(fig4_by_bin_plots)) {
  for (bin in names(fig4_by_bin_plots[[animal]])) {
    cat("  Figure 4 [", bin, "]:", animal, "\n")
    print(fig4_by_bin_plots[[animal]][[bin]])
  }
}
for (name in names(figS1_plots)) {
  cat("  Figure S1:", name, "\n")
  print(figS1_plots[[name]])
}
for (animal in names(figS1_by_cond_plots)) {
  for (cond in names(figS1_by_cond_plots[[animal]])) {
    cat("  Figure S1 [", cond, "]:", animal, "\n")
    print(figS1_by_cond_plots[[animal]][[cond]])
  }
}

cat("\n--- Generating condition-faceted presentation figures ---\n")
pres_plots <- generate_condition_faceted_figures()

cat("\n--- Printing condition-faceted figures ---\n")
for (name in names(pres_plots)) {
  cat("  Presentation:", name, "\n")
  print(pres_plots[[name]])
}

cat("\n=== All figures printed. Scroll through RStudio Plots pane. ===\n")
cat("Objects available: fig3_plots, fig3_by_cond_plots, fig4_plots, fig4_by_bin_plots, figS1_plots, figS1_by_cond_plots, pres_plots\n")
