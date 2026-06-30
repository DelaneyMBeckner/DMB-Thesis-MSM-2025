# ============================================================
# TRAJECTORY STATISTICAL ANALYSIS
# ============================================================
# Mixed-effects models for transition trajectory data
# 2-way: Condition × Epoch (for event rates and correlations)
# 3-way: Condition × Epoch × Distance (for correlations)
#
# Requires output from Pipeline_Transition_Analysis_v3.R
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)  # For p-values in lmer
library(emmeans)   # For post-hoc tests

# ============================================================
# DATA LOADING FUNCTIONS
# ============================================================

#' Load trajectory data for a single animal
#' 
#' @param animal_id Animal identifier (e.g., "mPFCm4")
#' @param data_dir Directory containing trajectory CSVs
#' @param window_size Window size used in pipeline (3 or 9)
#' @param file_suffix Suffix used in pipeline (e.g., "_filtered")
#' @return List with trajectory_data and trajectory_dist
load_animal_trajectory_data <- function(animal_id, 
                                        data_dir = "E:/Data_Processing/R/Results/within_subject",
                                        window_size = 3,
                                        file_suffix = "_filtered") {
  
  traj_file <- file.path(data_dir, paste0(animal_id, "_trajectory_data_", window_size, "ep", file_suffix, ".csv"))
  dist_file <- file.path(data_dir, paste0(animal_id, "_trajectory_dist_", window_size, "ep", file_suffix, ".csv"))
  
  if (!file.exists(traj_file)) {
    stop("Trajectory data file not found: ", traj_file)
  }
  
  trajectory_data <- read.csv(traj_file, stringsAsFactors = FALSE)
  
  trajectory_dist <- NULL
  if (file.exists(dist_file)) {
    trajectory_dist <- read.csv(dist_file, stringsAsFactors = FALSE)
  } else {
    warning("Distance-binned trajectory file not found: ", dist_file)
  }
  
  cat("Loaded trajectory data for", animal_id, "\n")
  cat("  Rows:", nrow(trajectory_data), "\n")
  cat("  Conditions:", paste(unique(trajectory_data$condition), collapse = ", "), "\n")
  cat("  Transition types:", paste(unique(trajectory_data$transition_type), collapse = ", "), "\n")
  
  return(list(
    trajectory_data = trajectory_data,
    trajectory_dist = trajectory_dist
  ))
}


#' Load event rate trajectory data for a single animal
#' 
#' @param animal_id Animal identifier
#' @param data_dir Base results directory
#' @param window_size Window size used (3 or 9)
#' @param file_suffix Suffix used in pipeline
#' @return Event rate trajectory dataframe
load_animal_event_data <- function(animal_id,
                                   data_dir = "E:/Data_Processing/R/Results",
                                   window_size = 3,
                                   file_suffix = "_filtered") {
  
 # Event rate data is saved as combined file in within_subject/
  event_file <- file.path(data_dir, "within_subject",
                          paste0(animal_id, "_trajectory_eventrates_", window_size, "ep", file_suffix, ".csv"))
  
  if (!file.exists(event_file)) {
    warning("Event rate trajectory file not found: ", event_file,
            "\nRe-run pipeline to generate this file.")
    return(NULL)
  }
  
  event_data <- read.csv(event_file, stringsAsFactors = FALSE)
  
  # Add animal column if not present
  if (!"animal" %in% colnames(event_data)) {
    event_data$animal <- animal_id
  }
  
  cat("Loaded event rate data for", animal_id, "\n")
  cat("  Rows:", nrow(event_data), "\n")
  cat("  Conditions:", paste(unique(event_data$condition), collapse = ", "), "\n")
  cat("  Transition types:", paste(unique(event_data$transition_type), collapse = ", "), "\n")
  
  return(event_data)
}


# ============================================================
# 2-WAY MIXED MODEL: CONDITION × EPOCH
# ============================================================

#' Run 2-way mixed-effects ANOVA for correlation trajectory
#' 
#' Model: correlation ~ condition * epoch + (1|ROI_pair) + (1|transition_id)
#' 
#' @param trajectory_data Trajectory dataframe with columns: correlation, condition, 
#'                        epoch_in_window, ROI_1, ROI_2, transition_id
#' @param transition_type Filter to specific transition type (NULL for both)
#' @return List with model, anova_table, and emmeans results
run_2way_correlation_model <- function(trajectory_data, transition_type = NULL) {
  
  cat("\n========================================\n")
  cat("2-WAY MIXED MODEL: Correlation ~ Condition × Epoch\n")
  cat("========================================\n")
  
  # Filter by transition type if specified
  if (!is.null(transition_type)) {
    trajectory_data <- trajectory_data %>% filter(transition_type == !!transition_type)
    cat("Transition type:", transition_type, "\n")
  }
  
  # Create ROI pair identifier
  trajectory_data <- trajectory_data %>%
    mutate(
      ROI_pair = paste(ROI_1, ROI_2, sep = "_"),
      condition = factor(condition, levels = c("BL", "SD", "WO")),
      epoch_factor = factor(epoch_in_window)
    )
  
  cat("Observations:", nrow(trajectory_data), "\n")
  cat("ROI pairs:", n_distinct(trajectory_data$ROI_pair), "\n")
  cat("Transitions:", n_distinct(trajectory_data$transition_id), "\n")
  cat("Conditions:", paste(levels(trajectory_data$condition), collapse = ", "), "\n")
  cat("Epochs:", paste(sort(unique(trajectory_data$epoch_in_window)), collapse = ", "), "\n")
  
  # Fit mixed model
  cat("\nFitting model...\n")
  
  model <- lmer(
    correlation ~ condition * epoch_factor + (1|ROI_pair) + (1|transition_id),
    data = trajectory_data,
    REML = TRUE
  )
  
  # Get ANOVA table with Satterthwaite degrees of freedom
  cat("\nANOVA Table (Type III, Satterthwaite):\n")
  anova_table <- anova(model, type = 3, ddf = "Satterthwaite")
  print(anova_table)
  
  # Convert to dataframe for saving
  anova_df <- as.data.frame(anova_table)
  anova_df$Effect <- rownames(anova_df)
  anova_df <- anova_df %>% select(Effect, everything())
  rownames(anova_df) <- NULL
  
  # Random effects summary
  cat("\nRandom Effects:\n")
  print(VarCorr(model))
  
  # Estimated marginal means
  cat("\nEstimated Marginal Means by Condition:\n")
  emm_condition <- emmeans(model, ~ condition)
  print(summary(emm_condition))
  
  cat("\nEstimated Marginal Means by Epoch:\n
")
  emm_epoch <- emmeans(model, ~ epoch_factor)
  print(summary(emm_epoch))
  
  cat("\nEstimated Marginal Means: Condition × Epoch:\n")
  emm_interaction <- emmeans(model, ~ condition | epoch_factor)
  
  # Post-hoc pairwise comparisons
  cat("\nPost-hoc: Condition comparisons at each epoch:\n")
  posthoc_condition <- pairs(emm_interaction, adjust = "tukey")
  posthoc_df <- as.data.frame(posthoc_condition)
  print(posthoc_df)
  
  # Condition effect across all epochs

  cat("\nPost-hoc: Overall condition comparisons:\n")
  posthoc_overall <- pairs(emm_condition, adjust = "tukey")
  print(posthoc_overall)
  
  return(list(
    model = model,
    anova_table = anova_df,
    emm_condition = emm_condition,
    emm_epoch = emm_epoch,
    emm_interaction = emm_interaction,
    posthoc_by_epoch = posthoc_df,
    posthoc_overall = as.data.frame(posthoc_overall),
    data = trajectory_data
  ))
}


#' Run 2-way mixed-effects ANOVA for event rate trajectory
#' 
#' Model: events_per_epoch ~ condition * epoch + (1|Cell_Name) + (1|transition_id)
#' 
#' @param event_data Event rate dataframe
#' @param transition_type Filter to specific transition type (NULL for both)
#' @return List with model, anova_table, and emmeans results
run_2way_eventrate_model <- function(event_data, transition_type = NULL) {
  
  cat("\n========================================\n")
  cat("2-WAY MIXED MODEL: Event Rate ~ Condition × Epoch\n")
  cat("========================================\n")
  
  # Filter by transition type if specified
  if (!is.null(transition_type)) {
    event_data <- event_data %>% filter(transition_type == !!transition_type)
    cat("Transition type:", transition_type, "\n")
  }
  
  # Prepare factors
  event_data <- event_data %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO")),
      epoch_factor = factor(epoch_in_window)
    )
  
  cat("Observations:", nrow(event_data), "\n")
  cat("ROIs:", n_distinct(event_data$Cell_Name), "\n")
  cat("Transitions:", n_distinct(event_data$transition_id), "\n")
  
  # Fit mixed model
  cat("\nFitting model...\n")
  
  model <- lmer(
    events_per_epoch ~ condition * epoch_factor + (1|Cell_Name) + (1|transition_id),
    data = event_data,
    REML = TRUE
  )
  
  # Get ANOVA table
  cat("\nANOVA Table (Type III, Satterthwaite):\n")
  anova_table <- anova(model, type = 3, ddf = "Satterthwaite")
  print(anova_table)
  
  anova_df <- as.data.frame(anova_table)
  anova_df$Effect <- rownames(anova_df)
  anova_df <- anova_df %>% select(Effect, everything())
  rownames(anova_df) <- NULL
  
  # Random effects
  cat("\nRandom Effects:\n")
  print(VarCorr(model))
  
  # EMMs
  emm_condition <- emmeans(model, ~ condition)
  emm_epoch <- emmeans(model, ~ epoch_factor)
  emm_interaction <- emmeans(model, ~ condition | epoch_factor)
  
  # Post-hoc
  cat("\nPost-hoc: Condition comparisons at each epoch:\n")
  posthoc_condition <- pairs(emm_interaction, adjust = "tukey")
  posthoc_df <- as.data.frame(posthoc_condition)
  print(posthoc_df)
  
  return(list(
    model = model,
    anova_table = anova_df,
    emm_condition = emm_condition,
    emm_epoch = emm_epoch,
    emm_interaction = emm_interaction,
    posthoc_by_epoch = posthoc_df,
    posthoc_overall = as.data.frame(pairs(emm_condition, adjust = "tukey")),
    data = event_data
  ))
}


# ============================================================
# 3-WAY MIXED MODEL: CONDITION × EPOCH × DISTANCE
# ============================================================

#' Run 3-way mixed-effects ANOVA for correlation with distance
#' 
#' Model: correlation ~ condition * epoch * distance_bin + (1|ROI_pair) + (1|transition_id)
#' 
#' @param trajectory_dist Distance-binned trajectory dataframe
#' @param transition_type Filter to specific transition type (NULL for both)
#' @return List with model, anova_table, and emmeans results
run_3way_correlation_model <- function(trajectory_dist, transition_type = NULL) {
  
  cat("\n========================================\n")
  cat("3-WAY MIXED MODEL: Correlation ~ Condition × Epoch × Distance\n")
  cat("========================================\n")
  
  # Filter by transition type if specified
  if (!is.null(transition_type)) {
    trajectory_dist <- trajectory_dist %>% filter(transition_type == !!transition_type)
    cat("Transition type:", transition_type, "\n")
  }
  
  # Create factors
  trajectory_dist <- trajectory_dist %>%
    mutate(
      ROI_pair = paste(ROI_1, ROI_2, sep = "_"),
      condition = factor(condition, levels = c("BL", "SD", "WO")),
      epoch_factor = factor(epoch_in_window),
      distance_bin = factor(distance_bin)
    )
  
  cat("Observations:", nrow(trajectory_dist), "\n")
  cat("ROI pairs:", n_distinct(trajectory_dist$ROI_pair), "\n")
  cat("Distance bins:", paste(levels(trajectory_dist$distance_bin), collapse = ", "), "\n")
  
  # Fit model - may need to simplify random structure if convergence issues
  cat("\nFitting model (this may take a moment)...\n")
  
  # Try full model first
  model <- tryCatch({
    lmer(
      correlation ~ condition * epoch_factor * distance_bin + 
        (1|ROI_pair) + (1|transition_id),
      data = trajectory_dist,
      REML = TRUE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    )
  }, error = function(e) {
    cat("Full model failed, trying simplified random structure...\n")
    lmer(
      correlation ~ condition * epoch_factor * distance_bin + (1|ROI_pair),
      data = trajectory_dist,
      REML = TRUE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    )
  })
  
  # ANOVA table
  cat("\nANOVA Table (Type III, Satterthwaite):\n")
  anova_table <- anova(model, type = 3, ddf = "Satterthwaite")
  print(anova_table)
  
  anova_df <- as.data.frame(anova_table)
  anova_df$Effect <- rownames(anova_df)
  anova_df <- anova_df %>% select(Effect, everything())
  rownames(anova_df) <- NULL
  
  # Random effects
  cat("\nRandom Effects:\n")
  print(VarCorr(model))
  
  # EMMs for key comparisons
  cat("\nEstimated Marginal Means: Condition × Distance:\n")
  emm_cond_dist <- emmeans(model, ~ condition | distance_bin)
  print(summary(emm_cond_dist))
  
  cat("\nEstimated Marginal Means: Epoch × Distance:\n")
  emm_epoch_dist <- emmeans(model, ~ epoch_factor | distance_bin)
  
  # Post-hoc: condition at each distance bin
  cat("\nPost-hoc: Condition comparisons by distance bin:\n")
  posthoc_cond_dist <- pairs(emm_cond_dist, adjust = "tukey")
  print(as.data.frame(posthoc_cond_dist))
  
  return(list(
    model = model,
    anova_table = anova_df,
    emm_cond_dist = emm_cond_dist,
    emm_epoch_dist = emm_epoch_dist,
    emm_full = emmeans(model, ~ condition * epoch_factor * distance_bin),
    posthoc_cond_by_dist = as.data.frame(posthoc_cond_dist),
    data = trajectory_dist
  ))
}


# ============================================================
# VISUALIZATION FUNCTIONS
# ============================================================

#' Create 2-way trajectory plot with significance annotations
#' 
#' @param model_results Output from run_2way_correlation_model or run_2way_eventrate_model
#' @param y_var Variable name for y-axis ("correlation" or "events_per_epoch")
#' @param title Plot title
#' @param transition_type Transition type for subtitle
#' @return ggplot object
plot_2way_trajectory <- function(model_results, y_var = "correlation", 

                                  title = "Trajectory Analysis",
                                  transition_type = NULL) {
  
  data <- model_results$data
  
  # Get EMM summary for plotting
  emm_df <- as.data.frame(model_results$emm_interaction)
  emm_df$epoch_in_window <- as.numeric(as.character(emm_df$epoch_factor))
  
  # Get posthoc results for annotations
  posthoc <- model_results$posthoc_by_epoch
  
  # Identify significant comparisons (p < 0.05)
  sig_comparisons <- posthoc %>%
    filter(p.value < 0.05) %>%
    mutate(
      epoch_in_window = as.numeric(as.character(epoch_factor)),
      label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Base plot
  p <- ggplot(emm_df, aes(x = epoch_in_window, y = emmean, color = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = condition), 
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                       labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    scale_fill_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                      labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    labs(
      title = title,
      subtitle = ifelse(!is.null(transition_type), 
                        paste("Transition:", transition_type), 
                        "All transitions"),
      x = "Epoch Position Relative to Transition",
      y = ifelse(y_var == "correlation", "Estimated Marginal Mean Correlation", 
                 "Estimated Marginal Mean Event Rate"),
      color = "Condition",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Add significance markers if any exist
  if (nrow(sig_comparisons) > 0) {
    # Find y position for annotations (above the data)
    y_max <- max(emm_df$upper.CL, na.rm = TRUE)
    y_range <- y_max - min(emm_df$lower.CL, na.rm = TRUE)
    
    # Add asterisks at top of each epoch with significant differences
    sig_epochs <- sig_comparisons %>%
      group_by(epoch_in_window) %>%
      summarise(label = paste(unique(label), collapse = ""), .groups = "drop")
    
    p <- p + 
      geom_text(data = sig_epochs,
                aes(x = epoch_in_window, y = y_max + 0.02 * y_range, label = label),
                inherit.aes = FALSE, size = 5, color = "black")
  }
  
  return(p)
}


#' Create ANOVA summary table plot
#' 
#' @param anova_df ANOVA table as dataframe
#' @param title Title for the table
#' @return ggplot object
plot_anova_table <- function(anova_df, title = "ANOVA Results") {
  
  # Format p-values
  anova_display <- anova_df %>%
    mutate(
      `F value` = round(`F value`, 2),
      `Pr(>F)` = case_when(
        `Pr(>F)` < 0.001 ~ "< 0.001 ***",
        `Pr(>F)` < 0.01 ~ paste0(round(`Pr(>F)`, 3), " **"),
        `Pr(>F)` < 0.05 ~ paste0(round(`Pr(>F)`, 3), " *"),
        TRUE ~ as.character(round(`Pr(>F)`, 3))
      )
    ) %>%
    select(Effect, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`, `Pr(>F)`)
  
  # Create table using ggplot
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, label = title, size = 5, fontface = "bold", hjust = 0.5) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)
  
  # Add table as text (simple approach)
  table_text <- paste(capture.output(print(anova_display, row.names = FALSE)), collapse = "\n")
  
  p <- p + annotate("text", x = 0.5, y = 0.5, label = table_text, 
                    size = 3, family = "mono", hjust = 0.5, vjust = 0.5)
  
  return(p)
}


#' Create 3-way visualization: faceted by distance bin
#' 
#' @param model_results Output from run_3way_correlation_model
#' @param transition_type Transition type for title
#' @return ggplot object
plot_3way_faceted_distance <- function(model_results, transition_type = NULL) {
  
  data <- model_results$data
  
  # Summarize data for plotting
  plot_data <- data %>%
    group_by(condition, epoch_in_window, distance_bin) %>%
    summarise(
      mean_corr = mean(correlation, na.rm = TRUE),
      se_corr = sd(correlation, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, aes(x = epoch_in_window, y = mean_corr, color = condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = mean_corr - se_corr, ymax = mean_corr + se_corr, fill = condition),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ distance_bin, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                       labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    scale_fill_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                      labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    labs(
      title = "Correlation Trajectory by Distance Bin",
      subtitle = ifelse(!is.null(transition_type),
                        paste("Transition:", transition_type),
                        "All transitions"),
      x = "Epoch Position Relative to Transition",
      y = "Mean Correlation ± SE",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


#' Create 3-way heatmap visualization
#' 
#' Rows: Condition × Distance, Columns: Epoch, Fill: Mean correlation
#' 
#' @param model_results Output from run_3way_correlation_model
#' @param transition_type Transition type for title
#' @return ggplot object
plot_3way_heatmap <- function(model_results, transition_type = NULL) {
  
  data <- model_results$data
  
  # Summarize for heatmap
  heatmap_data <- data %>%
    group_by(condition, epoch_in_window, distance_bin) %>%
    summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    mutate(cond_dist = paste(condition, distance_bin, sep = "\n"))
  
  # Order y-axis: BL bins, then SD bins, then WO bins
  cond_order <- expand.grid(
    distance_bin = unique(heatmap_data$distance_bin),
    condition = c("BL", "SD", "WO")
  ) %>%
    mutate(cond_dist = paste(condition, distance_bin, sep = "\n"))
  
  heatmap_data$cond_dist <- factor(heatmap_data$cond_dist, levels = cond_order$cond_dist)
  
  p <- ggplot(heatmap_data, aes(x = factor(epoch_in_window), y = cond_dist, fill = mean_corr)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_vline(xintercept = which(sort(unique(heatmap_data$epoch_in_window)) == 0) + 0.5,
               linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_viridis_c(option = "D", name = "Mean\nCorrelation") +
    labs(
      title = "Condition × Distance × Epoch Heatmap",
      subtitle = ifelse(!is.null(transition_type),
                        paste("Transition:", transition_type),
                        "All transitions"),
      x = "Epoch Position",
      y = "Condition / Distance Bin"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  return(p)
}


#' Create distance-slope comparison across epochs
#' 
#' Shows how the correlation-distance relationship changes across epoch positions
#' 
#' @param model_results Output from run_3way_correlation_model
#' @param transition_type Transition type for title
#' @return ggplot object
plot_distance_slope_by_epoch <- function(model_results, transition_type = NULL) {
  
  data <- model_results$data
  
  # For each condition and epoch, compute correlation-distance slope
  slope_data <- data %>%
    group_by(condition, epoch_in_window) %>%
    summarise(
      slope = coef(lm(correlation ~ mean_distance))[2],
      slope_se = summary(lm(correlation ~ mean_distance))$coefficients[2, 2],
      r_squared = summary(lm(correlation ~ mean_distance))$r.squared,
      n_pairs = n(),
      .groups = "drop"
    )
  
  p <- ggplot(slope_data, aes(x = epoch_in_window, y = slope, color = condition)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(aes(ymin = slope - slope_se, ymax = slope + slope_se, fill = condition),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                       labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    scale_fill_manual(values = c("BL" = "#2166AC", "SD" = "#B2182B", "WO" = "#1B7837"),
                      labels = c("BL" = "Baseline", "SD" = "Sleep Deprived", "WO" = "Recovery")) +
    labs(
      title = "Distance-Correlation Slope Across Transition",
      subtitle = ifelse(!is.null(transition_type),
                        paste("Transition:", transition_type, "| More negative = stronger distance-dependence"),
                        "More negative slope = stronger distance-dependence"),
      x = "Epoch Position Relative to Transition",
      y = "Correlation-Distance Slope",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}


# ============================================================
# MAIN ANALYSIS FUNCTION
# ============================================================

#' Run complete trajectory statistical analysis for one animal
#' 
#' @param animal_id Animal identifier
#' @param data_dir Base results directory
#' @param output_dir Output directory for stats results
#' @param file_suffix Suffix from pipeline
#' @param run_event_rates Include event rate analysis?
#' @param run_3way Include 3-way distance analysis?
#' @param save_outputs Save CSVs and plots?
#' @return List with all results
run_trajectory_statistics <- function(animal_id,
                                      data_dir = "E:/Data_Processing/R/Results",
                                      output_dir = NULL,
                                      file_suffix = "_filtered",
                                      window_size = 3,
                                      run_event_rates = TRUE,
                                      run_3way = TRUE,
                                      save_outputs = TRUE) {
  
  cat("\n")
  cat("############################################################\n")
  cat("# TRAJECTORY STATISTICAL ANALYSIS: ", animal_id, "\n")
  cat("############################################################\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(data_dir, "statistics")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list(animal_id = animal_id)
  
  # Load trajectory data
  traj_data <- load_animal_trajectory_data(
    animal_id, 
    file.path(data_dir, "within_subject"),
    window_size,
    file_suffix
  )
  
  # ============================================================
  # 2-WAY CORRELATION ANALYSIS
  # ============================================================
  
  cat("\n--- 2-Way Correlation Analysis ---\n")
  
  results$corr_2way <- list()
  results$corr_2way_plots <- list()
  
  for (trans_type in c("Wake2NREM", "NREM2Wake")) {
    cat("\nAnalyzing:", trans_type, "\n")
    
    model_result <- run_2way_correlation_model(traj_data$trajectory_data, trans_type)
    results$corr_2way[[trans_type]] <- model_result
    
    # Create plot
    p <- plot_2way_trajectory(model_result, y_var = "correlation",
                               title = paste(animal_id, "- Correlation Trajectory"),
                               transition_type = trans_type)
    results$corr_2way_plots[[trans_type]] <- p
    
    if (save_outputs) {
      # Save ANOVA table
      write.csv(model_result$anova_table,
                file.path(output_dir, paste0(animal_id, "_", trans_type, 
                                             "_correlation_2way_anova_", window_size, "ep", file_suffix, ".csv")),
                row.names = FALSE)
      
      # Save post-hoc results
      write.csv(model_result$posthoc_by_epoch,
                file.path(output_dir, paste0(animal_id, "_", trans_type,
                                             "_correlation_2way_posthoc_", window_size, "ep", file_suffix, ".csv")),
                row.names = FALSE)
      
      # Save plot
      ggsave(file.path(output_dir, paste0(animal_id, "_", trans_type,
                                          "_correlation_2way_EMM_", window_size, "ep", file_suffix, ".png")),
             p, width = 10, height = 6, dpi = 300)
      
      cat("  Saved 2-way correlation results for", trans_type, "\n")
    }
  }
  
  # ============================================================
  # 2-WAY EVENT RATE ANALYSIS
  # ============================================================
  
  if (run_event_rates) {
    cat("\n--- 2-Way Event Rate Analysis ---\n")
    
    event_data <- load_animal_event_data(animal_id, data_dir, window_size, file_suffix)
    
    if (!is.null(event_data)) {
      results$event_2way <- list()
      results$event_2way_plots <- list()
      
      for (trans_type in c("Wake2NREM", "NREM2Wake")) {
        cat("\nAnalyzing:", trans_type, "\n")
        
        model_result <- run_2way_eventrate_model(event_data, trans_type)
        results$event_2way[[trans_type]] <- model_result
        
        p <- plot_2way_trajectory(model_result, y_var = "events_per_epoch",
                                   title = paste(animal_id, "- Event Rate Trajectory"),
                                   transition_type = trans_type)
        results$event_2way_plots[[trans_type]] <- p
        
        if (save_outputs) {
          write.csv(model_result$anova_table,
                    file.path(output_dir, paste0(animal_id, "_", trans_type,
                                                 "_eventrate_2way_anova_", window_size, "ep", file_suffix, ".csv")),
                    row.names = FALSE)
          
          write.csv(model_result$posthoc_by_epoch,
                    file.path(output_dir, paste0(animal_id, "_", trans_type,
                                                 "_eventrate_2way_posthoc_", window_size, "ep", file_suffix, ".csv")),
                    row.names = FALSE)
          
          ggsave(file.path(output_dir, paste0(animal_id, "_", trans_type,
                                              "_eventrate_2way_EMM_", window_size, "ep", file_suffix, ".png")),
                 p, width = 10, height = 6, dpi = 300)
          
          cat("  Saved 2-way event rate results for", trans_type, "\n")
        }
      }
    }
  }
  
  # ============================================================
  # 3-WAY CORRELATION × DISTANCE ANALYSIS
  # ============================================================
  
  if (run_3way && !is.null(traj_data$trajectory_dist)) {
    cat("\n--- 3-Way Correlation × Distance Analysis ---\n")
    
    results$corr_3way <- list()
    results$corr_3way_plots <- list()
    
    for (trans_type in c("Wake2NREM", "NREM2Wake")) {
      cat("\nAnalyzing:", trans_type, "\n")
      
      model_result <- run_3way_correlation_model(traj_data$trajectory_dist, trans_type)
      results$corr_3way[[trans_type]] <- model_result
      
      # Create plots
      p_facet <- plot_3way_faceted_distance(model_result, trans_type)
      p_heat <- plot_3way_heatmap(model_result, trans_type)
      p_slope <- plot_distance_slope_by_epoch(model_result, trans_type)
      
      results$corr_3way_plots[[trans_type]] <- list(
        faceted = p_facet,
        heatmap = p_heat,
        slope = p_slope
      )
      
      if (save_outputs) {
        write.csv(model_result$anova_table,
                  file.path(output_dir, paste0(animal_id, "_", trans_type,
                                               "_correlation_3way_anova_", window_size, "ep", file_suffix, ".csv")),
                  row.names = FALSE)
        
        write.csv(model_result$posthoc_cond_by_dist,
                  file.path(output_dir, paste0(animal_id, "_", trans_type,
                                               "_correlation_3way_posthoc_", window_size, "ep", file_suffix, ".csv")),
                  row.names = FALSE)
        
        ggsave(file.path(output_dir, paste0(animal_id, "_", trans_type,
                                            "_correlation_3way_faceted_", window_size, "ep", file_suffix, ".png")),
               p_facet, width = 10, height = 8, dpi = 300)
        
        ggsave(file.path(output_dir, paste0(animal_id, "_", trans_type,
                                            "_correlation_3way_heatmap_", window_size, "ep", file_suffix, ".png")),
               p_heat, width = 10, height = 8, dpi = 300)
        
        ggsave(file.path(output_dir, paste0(animal_id, "_", trans_type,
                                            "_correlation_3way_slope_", window_size, "ep", file_suffix, ".png")),
               p_slope, width = 10, height = 6, dpi = 300)
        
        cat("  Saved 3-way correlation results for", trans_type, "\n")
      }
    }
  }
  
  # ============================================================
  # SUMMARY
  # ============================================================
  
  cat("\n############################################################\n")
  cat("# ANALYSIS COMPLETE: ", animal_id, "\n")
  cat("############################################################\n")
  
  if (save_outputs) {
    cat("\nOutput files saved to:", output_dir, "\n")
    cat("\nGenerated files:\n")
    cat("  2-way correlation: *_correlation_2way_anova.csv, *_2way_posthoc.csv, *_2way_EMM.png\n")
    if (run_event_rates) {
      cat("  2-way event rate:  *_eventrate_2way_anova.csv, *_2way_posthoc.csv, *_2way_EMM.png\n")
    }
    if (run_3way) {
      cat("  3-way correlation: *_correlation_3way_anova.csv, *_3way_posthoc.csv\n")
      cat("                     *_3way_faceted.png, *_3way_heatmap.png, *_3way_slope.png\n")
    }
  }
  
  return(results)
}


# ============================================================
# CONVENIENCE WRAPPERS
# ============================================================

#' Run statistics for mPFCm4 (viral injection animal)
run_mPFCm4_stats <- function(...) {
  run_trajectory_statistics("mPFCm4", ...)
}

#' Run statistics for transgenic animals
run_transgenic_stats <- function(data_dir = "E:/Data_Processing/R/Results",
                                  output_dir = NULL,
                                  ...) {
  
  animals <- c("mPFCf5", "mPFCf6", "mPFCm9")
  
  all_results <- list()
  
  for (animal in animals) {
    cat("\n\n========================================\n")
    cat("Processing:", animal, "\n")
    cat("========================================\n\n")
    
    all_results[[animal]] <- run_trajectory_statistics(animal, data_dir, output_dir, ...)
  }
  
  return(all_results)
}


# ============================================================
# SCRIPT INFO
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  Trajectory_Statistics.R\n")
  cat("  Mixed-Effects Models for Trajectories\n")
  cat("========================================\n")
  cat("\nUsage:\n")
  cat("  # Single animal:\n")
  cat("  results <- run_trajectory_statistics('mPFCm4')\n")
  cat("\n")
  cat("  # Transgenic animals only:\
n")
  cat("  results <- run_transgenic_stats()\n")
  cat("\n")
  cat("  # Custom paths:\n")
  cat("  results <- run_trajectory_statistics(\n")
  cat("    'mPFCf5',\n")
  cat("    data_dir = 'path/to/Results',\n")
  cat("    output_dir = 'path/to/stats_output'\n")
  cat("  )\n")
  cat("========================================\n\n")
}
