# Trace_SNR_Comparison.R
# Computes per-ROI SNR and compares across animals
# SNR = mean event amplitude / baseline noise (SD of quiescent periods)
# v2: Added activity-based ROI filtering (consistent with Pipeline_Transition_Analysis_v3)

library(tidyverse)

# =============================================================================
# CONFIGURATION
# =============================================================================

animals <- c("mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9")
condition <- "BL"

event_window <- 2  # Seconds around each event to exclude from baseline

base_path <- "E:/Data_Processing/R/Data CSVs"
output_path <- "E:/Data_Processing/R/Results"

# Activity-based filtering parameters (match Pipeline_Transition_Analysis_v3)
FILTER_BY_ACTIVITY <- TRUE
MIN_EVENTS_BASELINE <- 1      # Minimum events per baseline_epochs
BASELINE_EPOCHS <- 18         # Number of epochs for threshold calculation
TOTAL_RECORDING_EPOCHS <- 360 # Total epochs in recording

# =============================================================================
# SNR FUNCTIONS
# =============================================================================

compute_roi_snr <- function(trace, event_times, time_vec, event_window) {
  n_frames <- length(trace)
  
  # Create mask of quiescent periods (not near events)
  quiescent <- rep(TRUE, n_frames)
  
  for (evt in event_times) {
    exclude_idx <- which(abs(time_vec - evt) <= event_window)
    quiescent[exclude_idx] <- FALSE
  }
  
  # Baseline noise = SD during quiescent periods
  if (sum(quiescent) > 10) {
    baseline_sd <- sd(trace[quiescent], na.rm = TRUE)
  } else {
    baseline_sd <- NA
  }
  
  # Mean event amplitude (from trace values at event times)
  if (length(event_times) > 0) {
    event_frames <- sapply(event_times, function(evt) {
      which.min(abs(time_vec - evt))
    })
    event_frames <- event_frames[event_frames >= 1 & event_frames <= n_frames]
    mean_amplitude <- mean(trace[event_frames], na.rm = TRUE)
  } else {
    mean_amplitude <- NA
  }
  
  snr <- mean_amplitude / baseline_sd
  
  return(list(
    baseline_sd = baseline_sd,
    mean_amplitude = mean_amplitude,
    snr = snr,
    n_events = length(event_times)
  ))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

all_results <- list()
all_filter_stats <- list()

for (animal in animals) {
  cat("\n========================================\n")
  cat("Processing", animal, "...\n")
  cat("========================================\n")
  
  traces_file <- file.path(base_path, paste0(animal, "_", condition, "_Traces.csv"))
  events_file <- file.path(base_path, paste0(animal, "_", condition, "_Events.csv"))
  
  if (!file.exists(traces_file) | !file.exists(events_file)) {
    cat("  Files not found, skipping\n")
    next
  }
  
  # Load traces - IDPS format
  traces_raw <- read.csv(traces_file, header = TRUE, check.names = FALSE)
  traces <- traces_raw[-1, ]
  colnames(traces)[1] <- "Time_s"
  traces <- as.data.frame(lapply(traces, as.numeric))
  
  time_vec <- traces$Time_s
  roi_cols <- grep("^C[0-9]+", colnames(traces), value = TRUE)
  
  # Load events
  events <- read.csv(events_file, header = TRUE, check.names = FALSE)
  colnames(events) <- c("Time_s", "Cell_Name", "Value")
  events$Cell_Name <- trimws(events$Cell_Name)
  
  # Apply activity-based filtering (computed directly from loaded events CSV)
  if (FILTER_BY_ACTIVITY) {
    # Compute threshold: min_events per baseline_epochs scaled to full recording
    threshold <- MIN_EVENTS_BASELINE * (TOTAL_RECORDING_EPOCHS / BASELINE_EPOCHS)
    
    # Count events per ROI
    event_counts <- events %>%
      group_by(ROI = Cell_Name) %>%
      summarise(event_count = n(), .groups = "drop")
    
    # Identify active ROIs (those meeting threshold)
    active_rois <- event_counts %>%
      filter(event_count >= threshold) %>%
      pull(ROI)
    
    n_total <- nrow(event_counts)
    n_active <- length(active_rois)
    n_excluded <- n_total - n_active
    pct_retained <- round(100 * n_active / n_total, 1)
    
    cat("  Activity filtering:\n")
    cat("    Threshold:", threshold, "events\n")
    cat("    Active ROIs:", n_active, "of", n_total, 
        "(", pct_retained, "% retained)\n")
    cat("    Excluded:", n_excluded, "ROIs\n")
    
    # Filter roi_cols to only active ROIs
    roi_cols_filtered <- roi_cols[roi_cols %in% active_rois]
    
    all_filter_stats[[animal]] <- data.frame(
      animal = animal,
      total_rois = n_total,
      active_rois = n_active,
      excluded_rois = n_excluded,
      pct_retained = pct_retained,
      threshold = threshold
    )
  } else {
    roi_cols_filtered <- roi_cols
    cat("  Activity filtering: DISABLED\n")
  }
  
  cat("  Processing", length(roi_cols_filtered), "ROIs...\n")
  
  # Process each ROI
  animal_results <- data.frame()
  
  for (roi in roi_cols_filtered) {
    trace <- traces[[roi]]
    roi_events <- events %>% filter(Cell_Name == roi) %>% pull(Time_s)
    snr_result <- compute_roi_snr(trace, roi_events, time_vec, event_window)
    event_amplitudes <- events %>% filter(Cell_Name == roi) %>% pull(Value)
    if (length(event_amplitudes) > 0) {
      snr_result$mean_amplitude <- mean(event_amplitudes, na.rm = TRUE)
      snr_result$snr <- snr_result$mean_amplitude / snr_result$baseline_sd
    }
    animal_results <- rbind(animal_results, data.frame(
      animal = animal,
      roi = roi,
      baseline_sd = snr_result$baseline_sd,
      mean_amplitude = snr_result$mean_amplitude,
      snr = snr_result$snr,
      n_events = snr_result$n_events,
      filtered = TRUE  # All ROIs in results passed the filter
    ))
  }
  
  all_results[[animal]] <- animal_results
  cat("  Processed", nrow(animal_results), "active ROIs\n")
}

results_df <- bind_rows(all_results)
filter_stats_df <- bind_rows(all_filter_stats)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

summary_stats <- results_df %>%
  filter(!is.na(snr) & is.finite(snr) & n_events >= 1) %>%
  group_by(animal) %>%
  summarise(
    n_rois = n(),
    mean_snr = mean(snr, na.rm = TRUE),
    sd_snr = sd(snr, na.rm = TRUE),
    median_snr = median(snr, na.rm = TRUE),
    mean_baseline_sd = mean(baseline_sd, na.rm = TRUE),
    mean_amplitude = mean(mean_amplitude, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== FILTERING SUMMARY ===\n")
print(filter_stats_df)

cat("\n=== SNR SUMMARY STATISTICS (Active ROIs Only) ===\n")
print(summary_stats)

# =============================================================================
# STATISTICAL COMPARISON
# =============================================================================

valid_results <- results_df %>%
  filter(!is.na(snr) & is.finite(snr) & n_events >= 1)

kw_test <- kruskal.test(snr ~ animal, data = valid_results)
cat("\n=== Kruskal-Wallis test for SNR across animals ===\n")
print(kw_test)

if (kw_test$p.value < 0.05) {
  cat("\n=== Pairwise Wilcoxon tests (Holm correction) ===\n")
  pairwise <- pairwise.wilcox.test(valid_results$snr, valid_results$animal, 
                                   p.adjust.method = "holm")
  print(pairwise)
}

# mPFCm4 vs transgenic animals
cat("\n=== mPFCm4 vs transgenic animals ===\n")
valid_results <- valid_results %>%
  mutate(expression_type = ifelse(animal == "mPFCm4", "viral", "transgenic"))

wilcox_viral_vs_transgenic <- wilcox.test(snr ~ expression_type, data = valid_results)
print(wilcox_viral_vs_transgenic)

viral_snr <- valid_results %>% filter(expression_type == "viral") %>% pull(snr)
transgenic_snr <- valid_results %>% filter(expression_type == "transgenic") %>% pull(snr)
cat("Viral median SNR:", median(viral_snr), "\n")
cat("Transgenic median SNR:", median(transgenic_snr), "\n")
cat("Ratio:", median(viral_snr) / median(transgenic_snr), "\n")

# =============================================================================
# VISUALIZATION
# =============================================================================

animal_labels <- c("mPFCm4" = "Male 1", "mPFCm9" = "Male 2",
                   "mPFCf5" = "Female 1", "mPFCf6" = "Female 2")

animal_colors <- c("Male 1"   = "#D92B2B",
                   "Male 2"   = "#D4A017",
                   "Female 1" = "#7B2D8E",
                   "Female 2" = "#1B9E77")

animal_levels <- c("Male 1", "Male 2", "Female 1", "Female 2")

valid_results <- valid_results %>%
  mutate(animal_label = factor(animal_labels[animal], levels = animal_levels))

theme_snr <- theme_minimal() +
  theme(
    plot.title          = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.title.position = "plot",
    axis.title          = element_text(size = 13, face = "bold"),
    axis.text           = element_text(size = 12, face = "bold", color = "#111111"),
    legend.title        = element_text(size = 12, face = "bold"),
    legend.text         = element_text(size = 14),
    panel.grid.major.y = element_line(color = "#444444", linewidth = 0.4, linetype = "dashed"),
    panel.grid.major.x = element_line(color = "#bbbbbb", linewidth = 0.4, linetype = "solid"),
    panel.grid.minor.y = element_line(color = "#2B2B2B", linewidth = 0.2, linetype = "dashed"),
    panel.grid.minor.x = element_line(color = "#cccccc", linewidth = 0.2, linetype = "solid"),
    legend.position     = "none"
  )

p1 <- ggplot(valid_results, aes(x = animal_label, y = snr, fill = animal_label)) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = animal_colors) +
  labs(title = "Signal-to-Noise Ratio by Animal",
       x = "Animal", y = "SNR") +
  theme_snr

p2 <- ggplot(valid_results, aes(x = animal_label, y = baseline_sd, fill = animal_label)) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = animal_colors) +
  labs(title = "Baseline Noise by Animal",
       x = "Animal", y = "Std Dev Outside Events (ΔF/F)") +
  theme_snr

p3 <- ggplot(valid_results, aes(x = animal_label, y = mean_amplitude, fill = animal_label)) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = animal_colors) +
  labs(title = "Mean Event Amplitude by Animal",
       x = "Animal", y = "Mean ΔF/F at Events") +
  theme_snr

p4 <- ggplot(valid_results, aes(x = expression_type, y = snr, fill = expression_type)) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c("viral" = "#e41a1c", "transgenic" = "#377eb8")) +
  labs(title = "SNR: Viral vs Transgenic Expression",
       x = "Expression Method", y = "SNR") +
  theme_snr

print(p1)
print(p2)
print(p3)
print(p4)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Add suffix based on filtering
suffix <- ifelse(FILTER_BY_ACTIVITY, "_filtered", "")

write.csv(results_df, 
          file.path(output_path, paste0("trace_snr_by_roi", suffix, ".csv")), 
          row.names = FALSE)
write.csv(summary_stats, 
          file.path(output_path, paste0("trace_snr_summary", suffix, ".csv")), 
          row.names = FALSE)

if (FILTER_BY_ACTIVITY) {
  write.csv(filter_stats_df,
            file.path(output_path, "trace_snr_filter_summary.csv"),
            row.names = FALSE)
}

# Save plots
ggsave(file.path(output_path, paste0("snr_by_animal", suffix, ".png")), 
       p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_path, paste0("baseline_noise_by_animal", suffix, ".png")), 
       p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_path, paste0("event_amplitude_by_animal", suffix, ".png")), 
       p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_path, paste0("snr_viral_vs_transgenic", suffix, ".png")), 
       p4, width = 8, height = 6, dpi = 300)

cat("\n=== Results saved to:", output_path, "===\n")
cat("Files:\n")
cat("  - trace_snr_by_roi", suffix, ".csv\n", sep = "")
cat("  - trace_snr_summary", suffix, ".csv\n", sep = "")
if (FILTER_BY_ACTIVITY) {
  cat("  - trace_snr_filter_summary.csv\n")
}
cat("  - snr_by_animal", suffix, ".png\n", sep = "")
cat("  - baseline_noise_by_animal", suffix, ".png\n", sep = "")
cat("  - event_amplitude_by_animal", suffix, ".png\n", sep = "")
cat("  - snr_viral_vs_transgenic", suffix, ".png\n", sep = "")
