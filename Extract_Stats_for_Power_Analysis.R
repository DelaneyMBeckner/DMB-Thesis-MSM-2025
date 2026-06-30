# ============================================================
# EXTRACT STATISTICS FOR POWER ANALYSIS
# ============================================================
# Reads pipeline output files and computes power analysis for:
#   1. Event rates
#   2. Pairwise correlations (CCF)
#   3. Distance-binned correlations
#
# For each measure:
#   - Extracts animal-level means per condition
#   - Computes within-animal BL-SD differences
#   - Calculates SD of differences (for paired power analysis)
#   - Runs power analysis and generates figures
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(pwr)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

# Results directory (where pipeline outputs are)
RESULTS_DIR <- "E:/Data_Processing/R/Results"

# Animals to include (transgenic only - exclude mPFCm4)
ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm9")

# Conditions to compare
CONDITIONS <- c("BL", "SD", "WO")

# Analysis parameters
WINDOW_SIZE <- 9
TRANSITION_TYPES <- c("Wake2NREM", "NREM2Wake")
FILE_SUFFIX <- "all_simple"  # suffix used in batch run (no leading underscore)

# Output directory for power analysis results
OUTPUT_DIR <- "E:/Data_Processing/R/Results/Power_Analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Power analysis parameters
EFFECT_SIZES <- c(0.20, 0.30, 0.50)  # 20%, 30%, 50% changes
ALPHA <- 0.05
POWER_TARGET <- 0.80

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Generic file finder for epoch summaries
find_summary_files <- function(results_dir, animals, conditions, 
                                transition_types, window_size, suffix,
                                analysis_type, filename_pattern) {
  
  all_data <- data.frame()
  files_found <- c()
  files_missing <- c()
  
  for (animal in animals) {
    for (cond in conditions) {
      for (trans in transition_types) {
        # Build expected filename using provided pattern
        filename <- sprintf(filename_pattern, animal, cond, trans, window_size, suffix)
        
        # Build filepath based on analysis type
        filepath <- file.path(results_dir, analysis_type, trans, filename)
        
        if (file.exists(filepath)) {
          files_found <- c(files_found, filepath)
          
          df <- read.csv(filepath, stringsAsFactors = FALSE)
          df$animal <- animal
          df$condition <- cond
          df$transition_type <- trans
          
          all_data <- bind_rows(all_data, df)
        } else {
          files_missing <- c(files_missing, filepath)
        }
      }
    }
  }
  
  cat("  Found:", length(files_found), "files\n")
  if (length(files_missing) > 0) {
    cat("  Missing:", length(files_missing), "files\n")
  }
  
  return(list(data = all_data, found = files_found, missing = files_missing))
}


#' Compute animal-level summary from epoch data
compute_animal_summaries <- function(epoch_data, value_col = "mean_correlation") {
  
  # Check if value column exists
  if (!value_col %in% colnames(epoch_data)) {
    warning("Column '", value_col, "' not found. Available: ", 
            paste(colnames(epoch_data), collapse = ", "))
    return(NULL)
  }
  
  # Mean per animal x condition x transition type
  # Include CV (coefficient of variation) for epoch-to-epoch variability
  animal_summary <- epoch_data %>%
    group_by(animal, condition, transition_type) %>%
    summarise(
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      sd_value = sd(.data[[value_col]], na.rm = TRUE),
      cv_pct = 100 * sd(.data[[value_col]], na.rm = TRUE) / abs(mean(.data[[value_col]], na.rm = TRUE)),
      n_epochs = n(),
      .groups = "drop"
    )
  
  # Combined across transition types
  animal_summary_combined <- epoch_data %>%
    group_by(animal, condition) %>%
    summarise(
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      sd_value = sd(.data[[value_col]], na.rm = TRUE),
      cv_pct = 100 * sd(.data[[value_col]], na.rm = TRUE) / abs(mean(.data[[value_col]], na.rm = TRUE)),
      n_epochs = n(),
      .groups = "drop"
    ) %>%
    mutate(transition_type = "Combined")
  
  return(bind_rows(animal_summary, animal_summary_combined))
}


#' Compute CV summary statistics across animals (within-animal epoch-to-epoch)
compute_cv_summary <- function(animal_summary) {
  cv_summary <- animal_summary %>%
    filter(transition_type == "Combined") %>%
    group_by(condition) %>%
    summarise(
      mean_cv = mean(cv_pct, na.rm = TRUE),
      sd_cv = sd(cv_pct, na.rm = TRUE),
      min_cv = min(cv_pct, na.rm = TRUE),
      max_cv = max(cv_pct, na.rm = TRUE),
      n_animals = n(),
      .groups = "drop"
    )
  return(cv_summary)
}


#' Compute between-animal CV (variance of animal means)
#' This is the variability that repeated-measures design controls for
compute_between_animal_cv <- function(animal_summary) {
  between_cv <- animal_summary %>%
    filter(transition_type == "Combined") %>%
    group_by(condition) %>%
    summarise(
      mean_of_means = mean(mean_value, na.rm = TRUE),
      sd_of_means = sd(mean_value, na.rm = TRUE),
      between_animal_cv = 100 * sd(mean_value, na.rm = TRUE) / abs(mean(mean_value, na.rm = TRUE)),
      min_animal_mean = min(mean_value, na.rm = TRUE),
      max_animal_mean = max(mean_value, na.rm = TRUE),
      n_animals = n(),
      .groups = "drop"
    )
  return(between_cv)
}


#' Compute within-animal differences for paired analysis
compute_paired_differences <- function(animal_summary, 
                                        baseline_cond = "BL", 
                                        comparison_cond = "SD") {
  
  # Pivot to wide format
  wide_data <- animal_summary %>%
    filter(condition %in% c(baseline_cond, comparison_cond)) %>%
    select(animal, condition, transition_type, mean_value) %>%
    pivot_wider(
      names_from = condition,
      values_from = mean_value
    )
  
  # Compute differences
  wide_data <- wide_data %>%
    mutate(
      difference = .data[[baseline_cond]] - .data[[comparison_cond]],
      pct_change = 100 * difference / abs(.data[[baseline_cond]])
    )
  
  # Summary statistics of differences
  diff_summary <- wide_data %>%
    group_by(transition_type) %>%
    summarise(
      n_animals = n(),
      mean_BL = mean(.data[[baseline_cond]], na.rm = TRUE),
      sd_BL = sd(.data[[baseline_cond]], na.rm = TRUE),
      mean_diff = mean(difference, na.rm = TRUE),
      sd_diff = sd(difference, na.rm = TRUE),
      mean_pct_change = mean(pct_change, na.rm = TRUE),
      sd_pct_change = sd(pct_change, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    individual = wide_data,
    summary = diff_summary
  ))
}


#' Run power analysis for a single measure
run_power_analysis <- function(baseline_mean, sd_difference, n_observed,
                                effect_sizes = EFFECT_SIZES, 
                                alpha = ALPHA, 
                                power_target = POWER_TARGET) {
  
  # Handle edge cases
  if (is.na(sd_difference) || sd_difference == 0 || is.na(baseline_mean)) {
    return(list(
      results = NULL,
      achieved_power = NULL,
      error = "Insufficient data for power analysis"
    ))
  }
  
  # Calculate required sample sizes
  power_results <- data.frame()
  
  for (effect in effect_sizes) {
    expected_diff <- abs(baseline_mean) * effect
    cohens_d <- expected_diff / sd_difference
    
    pwr_result <- pwr.t.test(
      d = cohens_d,
      power = power_target,
      sig.level = alpha,
      type = "paired"
    )
    
    power_results <- bind_rows(power_results, data.frame(
      Effect_Pct = effect * 100,
      Expected_Diff = round(expected_diff, 4),
      Cohens_d = round(cohens_d, 3),
      N_Animals = ceiling(pwr_result$n)
    ))
  }
  
  # Calculate achieved power at current N
  achieved_power <- data.frame()
  effect_range <- seq(0.1, 1.0, by = 0.05)
  
  for (eff in effect_range) {
    expected_diff <- abs(baseline_mean) * eff
    d <- expected_diff / sd_difference
    pwr <- pwr.t.test(d = d, n = n_observed, sig.level = alpha, type = "paired")$power
    achieved_power <- bind_rows(achieved_power, data.frame(
      Effect_Pct = eff * 100, 
      Power = pwr
    ))
  }
  
  return(list(
    results = power_results,
    achieved_power = achieved_power,
    error = NULL
  ))
}


#' Generate power curve plot
create_power_curve_plot <- function(power_results, n_observed, 
                                     baseline_mean, sd_difference,
                                     title, alpha = ALPHA) {
  
  # Generate power curve data
  generate_power_curve <- function(cohens_d, n_range = 2:40) {
    power_values <- sapply(n_range, function(n) {
      pwr.t.test(d = cohens_d, n = n, sig.level = alpha, type = "paired")$power
    })
    data.frame(n = n_range, power = power_values)
  }
  
  power_curves <- data.frame()
  for (i in 1:nrow(power_results)) {
    curve_data <- generate_power_curve(power_results$Cohens_d[i])
    curve_data$effect <- paste0(power_results$Effect_Pct[i], "%")
    power_curves <- bind_rows(power_curves, curve_data)
  }
  
  p <- ggplot(power_curves, aes(x = n, y = power, color = effect)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
    annotate("text", x = 35, y = 0.82, label = "80% Power", size = 3.5) +
    geom_vline(xintercept = n_observed, linetype = "dotted", 
               color = "red", linewidth = 1) +
    annotate("text", x = n_observed + 1.5, y = 0.15, 
             label = paste0("N=", n_observed, "\n(this study)"), 
             size = 3, color = "red", hjust = 0) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 40, 5)) +
    labs(
      title = title,
      subtitle = paste0("Baseline: ", round(baseline_mean, 4),
                        ", SD of differences: ", round(sd_difference, 4)),
      x = "Number of Animals",
      y = "Statistical Power",
      color = "Effect Size"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}


# ============================================================
# MAIN ANALYSIS
# ============================================================

cat("============================================================\n")
cat("  COMPREHENSIVE POWER ANALYSIS\n")
cat("  Event Rates | Correlations | Distance-Binned Correlations\n")
cat("============================================================\n\n")

# Storage for all results
all_results <- list()

# ============================================================
# 1. PAIRWISE CORRELATIONS (CCF)
# ============================================================

cat("\n=== 1. PAIRWISE CORRELATIONS (CCF) ===\n")

# CCF files are loose in the base data directory, not in subfolders
CCF_DIR <- "E:/Data_Processing/R/Results"

ccf_all_data <- data.frame()
ccf_files_found <- c()

for (animal in ANIMALS) {
  for (cond in CONDITIONS) {
    for (trans in TRANSITION_TYPES) {
      filename <- sprintf("%s_%s_%s_%dep_%s_ccf_epoch_summary.csv", 
                          animal, cond, trans, WINDOW_SIZE, FILE_SUFFIX)
      filepath <- file.path(CCF_DIR, filename)
      
      if (file.exists(filepath)) {
        ccf_files_found <- c(ccf_files_found, filepath)
        df <- read.csv(filepath, stringsAsFactors = FALSE)
        df$animal <- animal
        df$condition <- cond
        df$transition_type <- trans
        ccf_all_data <- bind_rows(ccf_all_data, df)
      }
    }
  }
}

cat("  Found:", length(ccf_files_found), "files\n")

if (nrow(ccf_all_data) > 0) {
  ccf_animal_summary <- compute_animal_summaries(ccf_all_data, "mean_correlation")
  ccf_diff <- compute_paired_differences(ccf_animal_summary, "BL", "SD")
  ccf_cv <- compute_cv_summary(ccf_animal_summary)
  ccf_between_cv <- compute_between_animal_cv(ccf_animal_summary)
  
  # Use Combined or Wake2NREM
  ccf_stats <- ccf_diff$summary %>% filter(transition_type == "Combined")
  if (nrow(ccf_stats) == 0 || is.na(ccf_stats$sd_diff)) {
    ccf_stats <- ccf_diff$summary %>% filter(transition_type == "Wake2NREM")
  }
  
  ccf_power <- run_power_analysis(
    ccf_stats$mean_BL, 
    ccf_stats$sd_diff, 
    ccf_stats$n_animals
  )
  
  all_results$correlation <- list(
    animal_summary = ccf_animal_summary,
    differences = ccf_diff,
    cv_summary = ccf_cv,
    between_animal_cv = ccf_between_cv,
    stats = ccf_stats,
    power = ccf_power
  )
  
  cat("\n  Baseline mean:", round(ccf_stats$mean_BL, 4), "\n")
  cat("  SD of differences:", round(ccf_stats$sd_diff, 4), "\n")
  cat("  N animals:", ccf_stats$n_animals, "\n")
  
  # Report epoch-to-epoch CV
  bl_cv <- ccf_cv %>% filter(condition == "BL") %>% pull(mean_cv)
  if (length(bl_cv) > 0 && !is.na(bl_cv)) {
    cat("  Epoch-to-epoch CV (BL):", round(bl_cv, 1), "%\n")
  }
  
  # Report between-animal CV
  bl_between_cv <- ccf_between_cv %>% filter(condition == "BL") %>% pull(between_animal_cv)
  if (length(bl_between_cv) > 0 && !is.na(bl_between_cv)) {
    cat("  Between-animal CV (BL):", round(bl_between_cv, 1), "%\n")
  }
  
  if (!is.null(ccf_power$results)) {
    cat("\n  Required N for 30% effect:", 
        ccf_power$results$N_Animals[ccf_power$results$Effect_Pct == 30], "\n")
  }
} else {
  cat("  No CCF data found\n")
}


# ============================================================
# 2. EVENT RATES
# ============================================================

cat("\n=== 2. EVENT RATES ===\n")

event_pattern <- "%s_%s_%s_%dep_%s_event_rates_epoch_summary.csv"
event_files <- find_summary_files(
  RESULTS_DIR, ANIMALS, CONDITIONS, TRANSITION_TYPES, 
  WINDOW_SIZE, FILE_SUFFIX, "event", event_pattern
)

if (nrow(event_files$data) > 0) {
  # Event rate files might have different column names
  value_col <- if ("mean_event_rate" %in% colnames(event_files$data)) {
    "mean_event_rate"
  } else if ("mean_events" %in% colnames(event_files$data)) {
    "mean_events"
  } else if ("total_events" %in% colnames(event_files$data)) {
    "total_events"
  } else {
    # Find first numeric column that looks like a mean
    num_cols <- colnames(event_files$data)[sapply(event_files$data, is.numeric)]
    mean_cols <- num_cols[grepl("mean|rate|event", num_cols, ignore.case = TRUE)]
    if (length(mean_cols) > 0) mean_cols[1] else num_cols[1]
  }
  
  cat("  Using column:", value_col, "\n")
  
  event_animal_summary <- compute_animal_summaries(event_files$data, value_col)
  
  if (!is.null(event_animal_summary)) {
    event_diff <- compute_paired_differences(event_animal_summary, "BL", "SD")
    event_cv <- compute_cv_summary(event_animal_summary)
    event_between_cv <- compute_between_animal_cv(event_animal_summary)
    
    event_stats <- event_diff$summary %>% filter(transition_type == "Combined")
    if (nrow(event_stats) == 0 || is.na(event_stats$sd_diff)) {
      event_stats <- event_diff$summary %>% filter(transition_type == "Wake2NREM")
    }
    
    event_power <- run_power_analysis(
      event_stats$mean_BL, 
      event_stats$sd_diff, 
      event_stats$n_animals
    )
    
    all_results$event_rate <- list(
      animal_summary = event_animal_summary,
      differences = event_diff,
      cv_summary = event_cv,
      between_animal_cv = event_between_cv,
      stats = event_stats,
      power = event_power
    )
    
    cat("\n  Baseline mean:", round(event_stats$mean_BL, 4), "\n")
    cat("  SD of differences:", round(event_stats$sd_diff, 4), "\n")
    cat("  N animals:", event_stats$n_animals, "\n")
    
    # Report epoch-to-epoch CV
    bl_cv <- event_cv %>% filter(condition == "BL") %>% pull(mean_cv)
    if (length(bl_cv) > 0 && !is.na(bl_cv)) {
      cat("  Epoch-to-epoch CV (BL):", round(bl_cv, 1), "%\n")
    }
    
    # Report between-animal CV
    bl_between_cv <- event_between_cv %>% filter(condition == "BL") %>% pull(between_animal_cv)
    if (length(bl_between_cv) > 0 && !is.na(bl_between_cv)) {
      cat("  Between-animal CV (BL):", round(bl_between_cv, 1), "%\n")
    }
    
    if (!is.null(event_power$results)) {
      cat("\n  Required N for 30% effect:", 
          event_power$results$N_Animals[event_power$results$Effect_Pct == 30], "\n")
    }
  }
} else {
  cat("  No event rate data found\n")
}


# ============================================================
# 3. DISTANCE-BINNED CORRELATIONS
# ============================================================

cat("\n=== 3. DISTANCE-BINNED CORRELATIONS ===\n")

dist_pattern <- "%s_%s_%s_%dep_%s_distcorr_trajectory_bin_epoch_summary.csv"
dist_files <- find_summary_files(
  RESULTS_DIR, ANIMALS, CONDITIONS, TRANSITION_TYPES, 
  WINDOW_SIZE, FILE_SUFFIX, "distance_correlation", dist_pattern
)

if (nrow(dist_files$data) > 0) {
  # Distance files have multiple bins - analyze each separately
  cat("  Columns found:", paste(colnames(dist_files$data), collapse = ", "), "\n")
  
  # Check for distance bin column
  bin_col <- if ("distance_bin" %in% colnames(dist_files$data)) {
    "distance_bin"
  } else if ("distance_bin_label" %in% colnames(dist_files$data)) {
    "distance_bin_label"
  } else if ("bin" %in% colnames(dist_files$data)) {
    "bin"
  } else {
    NULL
  }
  
  if (!is.null(bin_col)) {
    distance_bins <- unique(dist_files$data[[bin_col]])
    cat("  Distance bins:", paste(distance_bins, collapse = ", "), "\n")
    
    all_results$distance_correlation <- list()
    
    for (dbin in distance_bins) {
      bin_data <- dist_files$data %>% filter(.data[[bin_col]] == dbin)
      
      bin_animal_summary <- compute_animal_summaries(bin_data, "mean_correlation")
      
      if (!is.null(bin_animal_summary) && nrow(bin_animal_summary) > 0) {
        bin_diff <- compute_paired_differences(bin_animal_summary, "BL", "SD")
        
        bin_stats <- bin_diff$summary %>% filter(transition_type == "Combined")
        if (nrow(bin_stats) == 0 || is.na(bin_stats$sd_diff)) {
          bin_stats <- bin_diff$summary %>% filter(transition_type == "Wake2NREM")
        }
        
        if (nrow(bin_stats) > 0) {
          bin_power <- run_power_analysis(
            bin_stats$mean_BL, 
            bin_stats$sd_diff, 
            bin_stats$n_animals
          )
          
          # Compute CV for this distance bin
          bin_cv_summary <- compute_cv_summary(bin_animal_summary)
          bin_between_cv <- compute_between_animal_cv(bin_animal_summary)
          
          all_results$distance_correlation[[as.character(dbin)]] <- list(
            animal_summary = bin_animal_summary,
            differences = bin_diff,
            stats = bin_stats,
            power = bin_power,
            cv_summary = bin_cv_summary,
            between_animal_cv = bin_between_cv
          )
          
          cat("\n  ", dbin, ":\n", sep = "")
          cat("    Baseline mean:", round(bin_stats$mean_BL, 4), "\n")
          cat("    SD of differences:", round(bin_stats$sd_diff, 4), "\n")
          if (!is.null(bin_power$results)) {
            cat("    Required N for 30% effect:", 
                bin_power$results$N_Animals[bin_power$results$Effect_Pct == 30], "\n")
          }
          # Report CV
          bl_cv <- bin_cv_summary$mean_cv[bin_cv_summary$condition == "BL"]
          if (length(bl_cv) > 0 && !is.na(bl_cv)) {
            cat("    Within-animal CV (BL):", round(bl_cv, 1), "%\n")
          }
          bl_between <- bin_between_cv$between_animal_cv[bin_between_cv$condition == "BL"]
          if (length(bl_between) > 0 && !is.na(bl_between)) {
            cat("    Between-animal CV (BL):", round(bl_between, 1), "%\n")
          }
        }
      }
    }
  } else {
    # No bin column - treat as single measure
    dist_animal_summary <- compute_animal_summaries(dist_files$data, "mean_correlation")
    if (!is.null(dist_animal_summary)) {
      dist_diff <- compute_paired_differences(dist_animal_summary, "BL", "SD")
      dist_stats <- dist_diff$summary %>% filter(transition_type == "Combined")
      
      dist_power <- run_power_analysis(
        dist_stats$mean_BL, 
        dist_stats$sd_diff, 
        dist_stats$n_animals
      )
      
      all_results$distance_correlation <- list(
        overall = list(
          animal_summary = dist_animal_summary,
          differences = dist_diff,
          stats = dist_stats,
          power = dist_power
        )
      )
    }
  }
} else {
  cat("  No distance correlation data found\n")
}


# ============================================================
# SAVE RESULTS
# ============================================================

cat("\n\n=== SAVING RESULTS ===\n")

# Save animal summaries
if (!is.null(all_results$correlation)) {
  write.csv(all_results$correlation$animal_summary, 
            file.path(OUTPUT_DIR, "correlation_animal_summaries.csv"), row.names = FALSE)
  write.csv(all_results$correlation$differences$individual, 
            file.path(OUTPUT_DIR, "correlation_BL_SD_differences.csv"), row.names = FALSE)
  write.csv(all_results$correlation$power$results, 
            file.path(OUTPUT_DIR, "correlation_power_results.csv"), row.names = FALSE)
  cat("  Correlation results saved\n")
}

if (!is.null(all_results$event_rate)) {
  write.csv(all_results$event_rate$animal_summary, 
            file.path(OUTPUT_DIR, "event_rate_animal_summaries.csv"), row.names = FALSE)
  write.csv(all_results$event_rate$differences$individual, 
            file.path(OUTPUT_DIR, "event_rate_BL_SD_differences.csv"), row.names = FALSE)
  if (!is.null(all_results$event_rate$power$results)) {
    write.csv(all_results$event_rate$power$results, 
              file.path(OUTPUT_DIR, "event_rate_power_results.csv"), row.names = FALSE)
  }
  cat("  Event rate results saved\n")
}

if (length(all_results$distance_correlation) > 0) {
  # Combine all distance bin results
  dist_power_combined <- data.frame()
  for (bin_name in names(all_results$distance_correlation)) {
    bin_results <- all_results$distance_correlation[[bin_name]]
    if (!is.null(bin_results$power$results)) {
      tmp <- bin_results$power$results
      tmp$Distance_Bin <- bin_name
      dist_power_combined <- bind_rows(dist_power_combined, tmp)
    }
  }
  if (nrow(dist_power_combined) > 0) {
    write.csv(dist_power_combined, 
              file.path(OUTPUT_DIR, "distance_correlation_power_results.csv"), row.names = FALSE)
  }
  cat("  Distance correlation results saved\n")
}


# ============================================================
# GENERATE FIGURES
# ============================================================

cat("\n=== GENERATING FIGURES ===\n")

plots <- list()

# Correlation power curve
if (!is.null(all_results$correlation$power$results)) {
  plots$correlation <- create_power_curve_plot(
    all_results$correlation$power$results,
    all_results$correlation$stats$n_animals,
    all_results$correlation$stats$mean_BL,
    all_results$correlation$stats$sd_diff,
    "Power Analysis: Pairwise Correlation (BL vs SD)"
  )
  ggsave(file.path(OUTPUT_DIR, "Correlation_Power_Curve.png"),
         plots$correlation, width = 10, height = 6, dpi = 300)
  cat("  Correlation power curve saved\n")
}

# Event rate power curve
if (!is.null(all_results$event_rate$power$results)) {
  plots$event_rate <- create_power_curve_plot(
    all_results$event_rate$power$results,
    all_results$event_rate$stats$n_animals,
    all_results$event_rate$stats$mean_BL,
    all_results$event_rate$stats$sd_diff,
    "Power Analysis: Event Rate (BL vs SD)"
  )
  ggsave(file.path(OUTPUT_DIR, "Event_Rate_Power_Curve.png"),
         plots$event_rate, width = 10, height = 6, dpi = 300)
  cat("  Event rate power curve saved\n")
}

# Distance bin comparison (if multiple bins)
if (length(all_results$distance_correlation) > 1) {
  dist_comparison <- data.frame()
  for (bin_name in names(all_results$distance_correlation)) {
    bin_results <- all_results$distance_correlation[[bin_name]]
    if (!is.null(bin_results$power$results)) {
      n_30 <- bin_results$power$results$N_Animals[bin_results$power$results$Effect_Pct == 30]
      dist_comparison <- bind_rows(dist_comparison, data.frame(
        Distance_Bin = bin_name,
        N_for_30pct = n_30,
        Baseline_Mean = bin_results$stats$mean_BL,
        SD_Diff = bin_results$stats$sd_diff
      ))
    }
  }
  
  if (nrow(dist_comparison) > 0) {
    # Map numeric bin labels to meaningful names
    bin_labels <- c("1" = "Close", "2" = "Medium", "3" = "Far",
                    "close" = "Close", "medium" = "Medium", "far" = "Far")
    dist_comparison <- dist_comparison %>%
      mutate(Distance_Bin_Label = ifelse(Distance_Bin %in% names(bin_labels),
                                          bin_labels[Distance_Bin],
                                          Distance_Bin),
             Distance_Bin_Label = factor(Distance_Bin_Label, 
                                          levels = c("Close", "Medium", "Far")))
    
    plots$distance_comparison <- ggplot(dist_comparison, 
                                         aes(x = Distance_Bin_Label, y = N_for_30pct, fill = Distance_Bin_Label)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = N_for_30pct), vjust = -0.5, size = 5, fontface = "bold") +
      geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
      annotate("text", x = 0.6, y = 4.5, label = "This study (N=3)", 
               color = "red", size = 4, hjust = 0) +
      scale_fill_manual(values = c("Close" = "#3A0557", "Medium" = "#185FA5", "Far" = "#0F6E56")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Required Sample Size by Distance Bin (30% Effect)",
        subtitle = "Paired design, 80% power, alpha = 0.05",
        x = "Distance Bin",
        y = "Number of Animals Required"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
      )
    
    ggsave(file.path(OUTPUT_DIR, "Distance_Bin_Power_Comparison.png"),
           plots$distance_comparison, width = 10, height = 6, dpi = 300)
    cat("  Distance bin comparison saved\n")
  }
}

# Combined summary bar chart
summary_data <- data.frame()
if (!is.null(all_results$correlation$power$results)) {
  summary_data <- bind_rows(summary_data, data.frame(
    Measure = "Correlation",
    N_for_30pct = all_results$correlation$power$results$N_Animals[
      all_results$correlation$power$results$Effect_Pct == 30]
  ))
}
if (!is.null(all_results$event_rate$power$results)) {
  summary_data <- bind_rows(summary_data, data.frame(
    Measure = "Event Rate",
    N_for_30pct = all_results$event_rate$power$results$N_Animals[
      all_results$event_rate$power$results$Effect_Pct == 30]
  ))
}

if (nrow(summary_data) > 0) {
  plots$summary <- ggplot(summary_data, aes(x = Measure, y = N_for_30pct, fill = Measure)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = N_for_30pct), vjust = -0.5, size = 5, fontface = "bold") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
    annotate("text", x = 0.6, y = 4.5, label = "This study (N=3)", 
             color = "red", size = 4, hjust = 0) +
    scale_fill_manual(values = c("Correlation" = "#A23B72", "Event Rate" = "#2E86AB")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Required Sample Size for 30% Effect",
      subtitle = "Paired design, 80% power, alpha = 0.05",
      x = NULL,
      y = "Number of Animals Required"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(file.path(OUTPUT_DIR, "Sample_Size_Summary.png"),
         plots$summary, width = 8, height = 6, dpi = 300)
  cat("  Summary comparison saved\n")
}

# Combined figure if we have multiple plots
available_plots <- plots[!sapply(plots, is.null)]
if (length(available_plots) >= 2) {
  main_plots <- available_plots[names(available_plots) %in% c("correlation", "event_rate", "summary")]
  
  if (length(main_plots) >= 2) {
    if (length(main_plots) == 2) {
      p_combined <- main_plots[[1]] + main_plots[[2]] +
        plot_annotation(title = "Power Analysis Summary",
                        theme = theme(plot.title = element_text(size = 18, face = "bold")))
    } else {
      p_combined <- (main_plots[[1]] + main_plots[[2]]) / main_plots[[3]] +
        plot_annotation(title = "Power Analysis Summary",
                        theme = theme(plot.title = element_text(size = 18, face = "bold")))
    }
    
    ggsave(file.path(OUTPUT_DIR, "Power_Analysis_Combined.png"),
           p_combined, width = 14, height = 10, dpi = 300)
    cat("  Combined figure saved\n")
  }
}


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("  POWER ANALYSIS SUMMARY\n")
cat("============================================================\n\n")

cat("MEASURES ANALYZED:\n")

if (!is.null(all_results$correlation$power$results)) {
  cat("\n1. PAIRWISE CORRELATION:\n")
  cat("   Baseline mean:", round(all_results$correlation$stats$mean_BL, 4), "\n")
  cat("   SD of differences:", round(all_results$correlation$stats$sd_diff, 4), "\n")
  # Add CV
  if (!is.null(all_results$correlation$cv_summary)) {
    bl_cv <- all_results$correlation$cv_summary %>% filter(condition == "BL") %>% pull(mean_cv)
    if (length(bl_cv) > 0 && !is.na(bl_cv)) {
      cat("   Epoch-to-epoch CV:", round(bl_cv, 1), "%\n")
    }
  }
  cat("   N for 30% effect:", 
      all_results$correlation$power$results$N_Animals[
        all_results$correlation$power$results$Effect_Pct == 30], "\n")
}

if (!is.null(all_results$event_rate$power$results)) {
  cat("\n2. EVENT RATE:\n")
  cat("   Baseline mean:", round(all_results$event_rate$stats$mean_BL, 4), "\n")
  cat("   SD of differences:", round(all_results$event_rate$stats$sd_diff, 4), "\n")
  # Add CV
  if (!is.null(all_results$event_rate$cv_summary)) {
    bl_cv <- all_results$event_rate$cv_summary %>% filter(condition == "BL") %>% pull(mean_cv)
    if (length(bl_cv) > 0 && !is.na(bl_cv)) {
      cat("   Epoch-to-epoch CV:", round(bl_cv, 1), "%\n")
    }
  }
  cat("   N for 30% effect:", 
      all_results$event_rate$power$results$N_Animals[
        all_results$event_rate$power$results$Effect_Pct == 30], "\n")
}

if (length(all_results$distance_correlation) > 0) {
  cat("\n3. DISTANCE-BINNED CORRELATIONS:\n")
  bin_labels <- c("1" = "Close", "2" = "Medium", "3" = "Far",
                  "close" = "Close", "medium" = "Medium", "far" = "Far")
  for (bin_name in names(all_results$distance_correlation)) {
    bin_results <- all_results$distance_correlation[[bin_name]]
    if (!is.null(bin_results$power$results)) {
      display_name <- ifelse(bin_name %in% names(bin_labels), 
                             bin_labels[bin_name], bin_name)
      cat("   ", display_name, ": N =", 
          bin_results$power$results$N_Animals[bin_results$power$results$Effect_Pct == 30],
          "for 30% effect\n")
    }
  }
}

# Determine maximum N recommendation
max_n <- 0
if (!is.null(all_results$correlation$power$results)) {
  max_n <- max(max_n, all_results$correlation$power$results$N_Animals[
    all_results$correlation$power$results$Effect_Pct == 30])
}
if (!is.null(all_results$event_rate$power$results)) {
  max_n <- max(max_n, all_results$event_rate$power$results$N_Animals[
    all_results$event_rate$power$results$Effect_Pct == 30])
}

cat("\n============================================================\n")
cat("  RECOMMENDATION: N =", max_n, "-", max_n + 3, "animals per group\n")
cat("============================================================\n")

# ============================================================
# CV SUMMARY - For Committee Question on Variability
# ============================================================

cat("\n============================================================\n")
cat("  COEFFICIENT OF VARIATION (CV) SUMMARY\n")
cat("============================================================\n\n")

cv_export <- data.frame()
between_cv_export <- data.frame()

# Within-animal (epoch-to-epoch) CV
cat("--- WITHIN-ANIMAL CV (epoch-to-epoch variability) ---\n\n")

if (!is.null(all_results$correlation$cv_summary)) {
  cat("PAIRWISE CORRELATION:\n")
  print(all_results$correlation$cv_summary %>% 
          mutate(across(where(is.numeric), ~round(., 1))))
  cv_export <- bind_rows(cv_export, 
                          all_results$correlation$cv_summary %>% 
                            mutate(measure = "correlation",
                                   cv_type = "within_animal"))
  cat("\n")
}

if (!is.null(all_results$event_rate$cv_summary)) {
  cat("EVENT RATE:\n")
  print(all_results$event_rate$cv_summary %>% 
          mutate(across(where(is.numeric), ~round(., 1))))
  cv_export <- bind_rows(cv_export, 
                          all_results$event_rate$cv_summary %>% 
                            mutate(measure = "event_rate",
                                   cv_type = "within_animal"))
  cat("\n")
}

# Between-animal CV
cat("\n--- BETWEEN-ANIMAL CV (variance of animal means) ---\n\n")

if (!is.null(all_results$correlation$between_animal_cv)) {
  cat("PAIRWISE CORRELATION:\n")
  print(all_results$correlation$between_animal_cv %>% 
          mutate(across(where(is.numeric), ~round(., 1))))
  between_cv_export <- bind_rows(between_cv_export, 
                                  all_results$correlation$between_animal_cv %>% 
                                    mutate(measure = "correlation"))
  cat("\n")
}

if (!is.null(all_results$event_rate$between_animal_cv)) {
  cat("EVENT RATE:\n")
  print(all_results$event_rate$between_animal_cv %>% 
          mutate(across(where(is.numeric), ~round(., 1))))
  between_cv_export <- bind_rows(between_cv_export, 
                                  all_results$event_rate$between_animal_cv %>% 
                                    mutate(measure = "event_rate"))
  cat("\n")
}

# Distance-binned CV
if (length(all_results$distance_correlation) > 0) {
  bin_labels <- c("1" = "Close", "2" = "Medium", "3" = "Far",
                  "close" = "Close", "medium" = "Medium", "far" = "Far")
  
  cat("\n--- DISTANCE-BINNED CORRELATION CV ---\n\n")
  
  for (bin_name in names(all_results$distance_correlation)) {
    bin_results <- all_results$distance_correlation[[bin_name]]
    display_name <- ifelse(bin_name %in% names(bin_labels), 
                           bin_labels[bin_name], bin_name)
    
    if (!is.null(bin_results$cv_summary)) {
      cat(display_name, " (Within-animal CV):\n", sep = "")
      print(bin_results$cv_summary %>% 
              mutate(across(where(is.numeric), ~round(., 1))))
      cv_export <- bind_rows(cv_export, 
                              bin_results$cv_summary %>% 
                                mutate(measure = paste0("dist_", display_name),
                                       cv_type = "within_animal"))
      cat("\n")
    }
    
    if (!is.null(bin_results$between_animal_cv)) {
      cat(display_name, " (Between-animal CV):\n", sep = "")
      print(bin_results$between_animal_cv %>% 
              mutate(across(where(is.numeric), ~round(., 1))))
      between_cv_export <- bind_rows(between_cv_export, 
                                      bin_results$between_animal_cv %>% 
                                        mutate(measure = paste0("dist_", display_name)))
      cat("\n")
    }
  }
}

# Save CV summaries to CSV
if (nrow(cv_export) > 0) {
  write.csv(cv_export, file.path(OUTPUT_DIR, "CV_Summary_Within_Animal.csv"), row.names = FALSE)
  cat("Within-animal CV saved to:", file.path(OUTPUT_DIR, "CV_Summary_Within_Animal.csv"), "\n")
}

if (nrow(between_cv_export) > 0) {
  write.csv(between_cv_export, file.path(OUTPUT_DIR, "CV_Summary_Between_Animal.csv"), row.names = FALSE)
  cat("Between-animal CV saved to:", file.path(OUTPUT_DIR, "CV_Summary_Between_Animal.csv"), "\n")
}

cat("\n--- INTERPRETATION ---\n")
cat("Within-animal CV:  Epoch-to-epoch measurement noise (~14-28%)\n")
cat("Between-animal CV: Baseline differences across animals (controlled by paired design)\n")
cat("SD of differences: Consistency of BL→SD change (what paired t-test uses)\n")
cat("\nRepeated measures eliminates between-animal variance, focusing on consistent changes.\n")

cat("\nFiles saved to:", OUTPUT_DIR, "\n")
cat(list.files(OUTPUT_DIR), sep = "\n")

cat("\n\n=== SUGGESTED THESIS TEXT ===\n\n")
cat('"Post-hoc power analysis based on observed within-animal variability\n')
cat('indicates that detecting a 30% change with 80% power would require\n')
cat('N = ', max_n, ' animals (paired design). The current study (N = 3) was\n', sep = "")
cat('therefore underpowered for moderate effect sizes. Future studies\n')
cat('should include N = ', max_n, '-', max_n + 3, ' animals for robust detection of\n', sep = "")
cat('sleep deprivation effects on mPFC network dynamics."\n')
