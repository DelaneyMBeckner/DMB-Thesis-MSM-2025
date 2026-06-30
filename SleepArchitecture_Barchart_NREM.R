#!/usr/bin/env Rscript
# SleepArchitecture_Barchart.R
# Complete sleep architecture analysis:
#   1. Stacked bar chart (% time per state, all animals)
#   2. Variance plots (individual points + mean +/- SEM)
#   3. Bout count and duration analysis
#   4. Delta power analysis from EDFs (optional)
#   5. Statistical tests (RM-ANOVA and Friedman)

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# CONFIGURATION
# ============================================================

INPUT_PATH <- "E:/Data_Processing/R/Data CSVs"
INPUT_PATH_EDF <- "E:/Medial_PreFrontal_Cortex/EDFs"
OUTPUT_PATH <- "E:/Data_Processing/R/sleep_architecture_figures"

ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9")
CONDITIONS <- c("BL", "SD", "WO")

STATE_LABELS <- c("1" = "Wake", "2" = "NREM", "3" = "REM")
STATE_COLORS <- c("Wake" = "#F4D03F", "NREM" = "#3498DB", "REM" = "#E74C3C")

# Toggle: TRUE = facet by condition, FALSE = facet by animal
DAY_FIRST <- TRUE

# Delta power settings
ANALYZE_DELTA <- TRUE  # Set FALSE to skip delta power analysis
CHANNEL <- "EEG1"
DELTA_BAND <- c(0.5, 4.0)  # Hz
EPOCH_DURATION <- 10  # seconds

# Load edfReader only if needed
if (ANALYZE_DELTA) {
  library(edfReader)
}

# ============================================================
# BOUT IDENTIFICATION FUNCTION
# ============================================================

identify_bouts <- function(states_df) {
  states_df <- states_df[order(states_df$Time.from.Start), ]
  
  state_vec <- states_df$state
  changes <- c(TRUE, state_vec[-1] != state_vec[-length(state_vec)])
  bout_id <- cumsum(changes)
  
  bout_df <- data.frame(
    bout_id = bout_id,
    state = state_vec,
    state_label = states_df$state_label
  ) %>%
    group_by(bout_id, state, state_label) %>%
    summarize(length_epochs = n(), .groups = 'drop') %>%
    mutate(length_sec = length_epochs * 10)
  
  return(bout_df)
}

# ============================================================
# DELTA POWER FUNCTION
# ============================================================

compute_delta_power <- function(signal, srate, delta_band = DELTA_BAND) {
  # Compute PSD using spectrum()
  spec <- spectrum(signal, plot = FALSE)
  
  # Convert normalized frequency to Hz
  freq_hz <- spec$freq * srate
  
  # Find delta band indices
  delta_idx <- which(freq_hz >= delta_band[1] & freq_hz <= delta_band[2])
  
  # Mean power in delta band
  delta_power <- mean(spec$spec[delta_idx])
  
  return(delta_power)
}

# ============================================================
# LOAD DATA AND CALCULATE PERCENTAGES
# ============================================================

cat("Loading data...\n")

all_percentages <- data.frame()
bout_summary <- data.frame()
delta_results <- data.frame()

for (animal in ANIMALS) {
  for (condition in CONDITIONS) {
    
    filename <- paste0(animal, "_", condition, "_states_df.csv")
    filepath <- file.path(INPUT_PATH, filename)
    
    if (file.exists(filepath)) {
      states_df <- read.csv(filepath)
      states_df$state_label <- trimws(states_df$state_label)
      
      total_epochs <- nrow(states_df)
      
      # State percentages
      pct <- states_df %>%
        group_by(state_label) %>%
        summarize(n_epochs = n(), .groups = 'drop') %>%
        mutate(
          percent = 100 * n_epochs / total_epochs,
          animal = animal,
          condition = condition
        )
      
      all_percentages <- rbind(all_percentages, pct)
      
      # Bout analysis
      bouts <- identify_bouts(states_df)
      
      for (state in c("Wake", "NREM")) {
        state_bouts <- bouts[bouts$state_label == state, ]
        
        n <- nrow(state_bouts)
        mean_dur <- mean(state_bouts$length_sec)
        sem <- sd(state_bouts$length_sec) / sqrt(n)
        
        bout_summary <- rbind(bout_summary, data.frame(
          animal = animal,
          condition = condition,
          state = state,
          n_bouts = n,
          mean_sec = round(mean_dur, 1),
          sem = round(sem, 1)
        ))
      }
      
      cat("Loaded:", filename, "-", total_epochs, "epochs\n")
      
      # Delta power analysis
      if (ANALYZE_DELTA) {
        edf_file <- file.path(INPUT_PATH_EDF, paste0(animal, "_", condition, ".edf"))
        
        if (file.exists(edf_file)) {
          cat("  Analyzing delta power...\n")
          
          tryCatch({
            # Load EDF
            hdr <- readEdfHeader(edf_file)
            srate <- hdr$sHeaders$sRate[1]
            signals <- readEdfSignals(hdr)
            
            # Find EEG channel
            signal_names <- names(signals)
            eeg_idx <- grep(CHANNEL, signal_names, ignore.case = TRUE)[1]
            
            if (!is.na(eeg_idx)) {
              eeg_signal <- signals[[eeg_idx]]$signal
              
              for (i in 1:nrow(states_df)) {
                epoch_start <- states_df$Time.from.Start[i]
                state_label <- states_df$state_label[i]
                
                # Skip REM
                if (!state_label %in% c("Wake", "NREM")) next
                
                # Get epoch samples
                start_samp <- epoch_start * srate + 1
                end_samp <- (epoch_start + EPOCH_DURATION) * srate
                
                # Check bounds
                if (end_samp > length(eeg_signal)) next
                
                # Extract epoch and compute delta power
                epoch_signal <- eeg_signal[start_samp:end_samp]
                
                delta_power <- tryCatch({
                  compute_delta_power(epoch_signal, srate)
                }, error = function(e) NA)
                
                if (!is.na(delta_power)) {
                  delta_results <- rbind(delta_results, data.frame(
                    animal = animal,
                    condition = condition,
                    epoch = i,
                    time = epoch_start,
                    state = state_label,
                    delta_power = delta_power
                  ))
                }
              }
              cat("    ", sum(delta_results$animal == animal & delta_results$condition == condition), "epochs\n")
            } else {
              cat("  Channel", CHANNEL, "not found\n")
            }
          }, error = function(e) {
            cat("  EDF error:", e$message, "\n")
          })
        } else {
          cat("  Missing EDF:", edf_file, "\n")
        }
      }
      
    } else {
      cat("Missing:", filename, "\n")
    }
  }
}

# Create output directory
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)

# ============================================================
# STACKED BAR CHART
# ============================================================

cat("\nCreating stacked bar chart...\n")

all_percentages$condition <- factor(all_percentages$condition, 
                                    levels = c("BL", "SD", "WO"),
                                    labels = c("Baseline", "Sleep Dep", "Recovery"))
all_percentages$state_label <- factor(all_percentages$state_label,
                                      levels = c("REM", "NREM", "Wake"))

if (DAY_FIRST == TRUE) {
  p <- ggplot(all_percentages, aes(x = animal, y = percent, fill = state_label)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    facet_wrap(~ condition, nrow = 1) +
    scale_fill_manual(values = STATE_COLORS, name = "State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Sleep Architecture by Animal and Condition",
      x = NULL,
      y = "Time (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    )
}

if (DAY_FIRST == FALSE) {
  p <- ggplot(all_percentages, aes(x = condition, y = percent, fill = state_label)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    facet_wrap(~ animal, nrow = 1) +
    scale_fill_manual(values = STATE_COLORS, name = "State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Sleep Architecture by Animal and Condition",
      x = NULL,
      y = "Time (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    )
}

ggsave(file.path(OUTPUT_PATH, "sleep_architecture_stacked.png"), p,
       width = 10, height = 5, dpi = 300)
print(p)

# ============================================================
# VARIANCE PLOT - TIME IN STATE
# ============================================================

cat("Creating variance plots...\n")

# Reset factor levels for variance plot (Wake first for display)
all_percentages$state_label <- factor(all_percentages$state_label,
                                       levels = c("Wake", "NREM", "REM"))

p_variance <- ggplot(all_percentages, aes(x = condition, y = percent)) +
  geom_point(aes(color = animal), size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ state_label, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "Animal") +
  labs(title = "Sleep Architecture by State", x = NULL, y = "Time (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_PATH, "sleep_architecture_variance.png"), p_variance,
       width = 10, height = 4, dpi = 300)
print(p_variance)

# ============================================================
# VARIANCE PLOT - BOUT COUNT
# ============================================================

bout_summary$condition <- factor(bout_summary$condition, 
                                  levels = c("BL", "SD", "WO"),
                                  labels = c("Baseline", "Sleep Dep", "Recovery"))
bout_summary$state <- factor(bout_summary$state, levels = c("Wake", "NREM"))

p_bout_count <- ggplot(bout_summary, aes(x = condition, y = n_bouts)) +
  geom_point(aes(color = animal), size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ state, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "Animal") +
  labs(title = "Bout Count by State", x = NULL, y = "Number of Bouts") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_PATH, "bout_count_variance.png"), p_bout_count,
       width = 8, height = 4, dpi = 300)
print(p_bout_count)

# ============================================================
# VARIANCE PLOT - BOUT DURATION
# ============================================================

p_bout_duration <- ggplot(bout_summary, aes(x = condition, y = mean_sec)) +
  geom_point(aes(color = animal), size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ state, scales = "free_y") +
  scale_color_viridis_d(option = "turbo", name = "Animal") +
  labs(title = "Bout Duration by State", x = NULL, y = "Mean Duration (sec)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_PATH, "bout_duration_variance.png"), p_bout_duration,
       width = 8, height = 4, dpi = 300)
print(p_bout_duration)

# ============================================================
# DELTA POWER PLOTS
# ============================================================

if (ANALYZE_DELTA && nrow(delta_results) > 0) {
  cat("\nCreating delta power plots...\n")
  
  # Log transform
  delta_results$log_delta_power <- log10(delta_results$delta_power)
  
  # Summary by recording
  delta_summary <- delta_results %>%
    group_by(animal, condition, state) %>%
    summarize(
      delta_power = mean(delta_power),
      log_delta_power = mean(log_delta_power),
      n_epochs = n(),
      .groups = 'drop'
    )
  
  # Save summaries
  write.csv(delta_summary, file.path(OUTPUT_PATH, "delta_power_summary.csv"), row.names = FALSE)
  write.csv(delta_results, file.path(OUTPUT_PATH, "delta_power_all_epochs.csv"), row.names = FALSE)
  
  # Format for plotting
  delta_summary$condition <- factor(delta_summary$condition,
                                     levels = c("BL", "SD", "WO"),
                                     labels = c("Baseline", "Sleep Dep", "Recovery"))
  # Filter to NREM only for delta power analysis (Wake not informative for sleep pressure)
  delta_summary <- delta_summary %>% filter(state == "NREM")
  delta_summary$state <- factor(delta_summary$state, levels = c("NREM"))
  
  # Plot: Log delta power
  p_delta <- ggplot(delta_summary, aes(x = condition, y = log_delta_power)) +
    geom_point(aes(color = animal), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    facet_wrap(~ state, scales = "free_y") +
    scale_color_viridis_d(option = "turbo", name = "Animal") +
    labs(
      title = "Delta Power by State (0.5-4 Hz)",
      x = NULL,
      y = expression(Log[10]~Delta~Power)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(OUTPUT_PATH, "delta_power_variance.png"), p_delta,
         width = 10, height = 5, dpi = 300)
  print(p_delta)
  
  # Plot: Raw delta power
  p_delta_raw <- ggplot(delta_summary, aes(x = condition, y = delta_power)) +
    geom_point(aes(color = animal), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    facet_wrap(~ state, scales = "free_y") +
    scale_color_viridis_d(option = "turbo", name = "Animal") +
    labs(
      title = "Delta Power by State (Raw)",
      x = NULL,
      y = "Delta Power"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(OUTPUT_PATH, "delta_power_variance_raw.png"), p_delta_raw,
         width = 10, height = 5, dpi = 300)
  print(p_delta_raw)
  
  # Boxplot: All epochs
  delta_results$condition <- factor(delta_results$condition,
                                     levels = c("BL", "SD", "WO"),
                                     labels = c("Baseline", "Sleep Dep", "Recovery"))
  delta_results <- delta_results %>% filter(state == "NREM")
  delta_results$state <- factor(delta_results$state, levels = c("NREM"))
  
  p_delta_box <- ggplot(delta_results, aes(x = condition, y = log_delta_power, fill = animal)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~ state) +
    scale_fill_viridis_d(option = "turbo", name = "Animal") +
    labs(
      title = "Delta Power Distribution (All Epochs)",
      x = NULL,
      y = expression(Log[10]~Delta~Power)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(OUTPUT_PATH, "delta_power_boxplot.png"), p_delta_box,
         width = 10, height = 5, dpi = 300)
  print(p_delta_box)
}

# ============================================================
# SUMMARY TABLE - BOUT STATISTICS
# ============================================================

cat("\nSaving bout summary table...\n")

bout_table <- bout_summary %>%
  select(animal, condition, state, mean_sec, sem) %>%
  pivot_wider(names_from = state, values_from = c(mean_sec, sem))

write.csv(bout_table, file.path(OUTPUT_PATH, "bout_length_summary.csv"), row.names = FALSE)
print(bout_table)

# ============================================================
# STATISTICAL TESTS
# ============================================================

cat("\n===================================================\n")
cat("STATISTICAL ANALYSIS - SLEEP ARCHITECTURE\n")
cat("===================================================\n")
cat("n =", length(ANIMALS), "animals, 3 conditions (repeated measures)\n")
cat("Parametric: Repeated measures ANOVA\n")
cat("Non-parametric: Friedman test\n")
cat("===================================================\n")

# Need raw factor levels for stats
all_percentages_stats <- all_percentages
all_percentages_stats$condition <- factor(all_percentages_stats$condition,
                                           levels = c("Baseline", "Sleep Dep", "Recovery"))

bout_summary_stats <- bout_summary
bout_summary_stats$condition <- factor(bout_summary_stats$condition,
                                        levels = c("Baseline", "Sleep Dep", "Recovery"))

# 1. TIME IN STATE
cat("\n=== TIME IN STATE ===\n")

for (state in c("Wake", "NREM", "REM")) {
  state_data <- all_percentages_stats %>% 
    filter(state_label == state)
  
  state_data$animal <- factor(state_data$animal)
  
  cat("\n---", state, "---\n")
  
  # Repeated measures ANOVA
  model <- aov(percent ~ condition + Error(animal/condition), data = state_data)
  result <- summary(model)
  cat("RM-ANOVA:\n")
  print(result)
  
  # Friedman test
  mat <- state_data %>%
    select(animal, condition, percent) %>%
    pivot_wider(names_from = condition, values_from = percent) %>%
    select(-animal) %>%
    as.matrix()
  
  friedman <- friedman.test(mat)
  cat("Friedman: chi-sq =", round(friedman$statistic, 3), 
      ", p =", round(friedman$p.value, 4), "\n")
}

# 2. BOUT COUNT
cat("\n=== BOUT COUNT ===\n")

for (state_name in c("Wake", "NREM")) {
  state_data <- bout_summary_stats %>% 
    filter(state == state_name)
  
  state_data$animal <- factor(state_data$animal)
  
  cat("\n---", state_name, "---\n")
  
  # RM-ANOVA
  model <- aov(n_bouts ~ condition + Error(animal/condition), data = state_data)
  result <- summary(model)
  cat("RM-ANOVA:\n")
  print(result)
  
  # Friedman
  mat <- state_data %>%
    select(animal, condition, n_bouts) %>%
    pivot_wider(names_from = condition, values_from = n_bouts) %>%
    select(-animal) %>%
    as.matrix()
  
  friedman <- friedman.test(mat)
  cat("Friedman: chi-sq =", round(friedman$statistic, 3), 
      ", p =", round(friedman$p.value, 4), "\n")
}

# 3. BOUT DURATION
cat("\n=== BOUT DURATION ===\n")

for (state_name in c("Wake", "NREM")) {
  state_data <- bout_summary_stats %>% 
    filter(state == state_name)
  
  state_data$animal <- factor(state_data$animal)
  
  cat("\n---", state_name, "---\n")
  
  # RM-ANOVA
  model <- aov(mean_sec ~ condition + Error(animal/condition), data = state_data)
  result <- summary(model)
  cat("RM-ANOVA:\n")
  print(result)
  
  # Friedman
  mat <- state_data %>%
    select(animal, condition, mean_sec) %>%
    pivot_wider(names_from = condition, values_from = mean_sec) %>%
    select(-animal) %>%
    as.matrix()
  
  friedman <- friedman.test(mat)
  cat("Friedman: chi-sq =", round(friedman$statistic, 3), 
      ", p =", round(friedman$p.value, 4), "\n")
}

# 4. DELTA POWER
if (ANALYZE_DELTA && exists("delta_summary") && nrow(delta_summary) > 0) {
  cat("\n=== DELTA POWER ===\n")
  
  delta_summary_stats <- delta_summary
  delta_summary_stats$condition <- factor(delta_summary_stats$condition,
                                           levels = c("Baseline", "Sleep Dep", "Recovery"))
  
  for (state_name in c("NREM")) {
    state_data <- delta_summary_stats %>% 
      filter(state == state_name)
    
    state_data$animal <- factor(state_data$animal)
    
    cat("\n---", state_name, "---\n")
    
    # RM-ANOVA
    model <- aov(log_delta_power ~ condition + Error(animal/condition), data = state_data)
    result <- summary(model)
    cat("RM-ANOVA:\n")
    print(result)
    
    # Friedman
    mat <- state_data %>%
      select(animal, condition, log_delta_power) %>%
      pivot_wider(names_from = condition, values_from = log_delta_power) %>%
      select(-animal) %>%
      as.matrix()
    
    friedman <- friedman.test(mat)
    cat("Friedman: chi-sq =", round(friedman$statistic, 3), 
        ", p =", round(friedman$p.value, 4), "\n")
  }
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output saved to:", OUTPUT_PATH, "\n")
