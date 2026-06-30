# ============================================================
# POWER ANALYSIS - REPEATED MEASURES (PAIRED) DESIGN
# ============================================================
# For within-subject comparisons: BL vs SD, BL vs WO, etc.
# Answers: "How many animals do we need to detect X% change?"
#
# Includes analysis for:
#   1. Event rates
#   2. Pairwise correlations
#   3. Distance-binned correlations (close/medium/far)
# ============================================================

library(ggplot2)
library(pwr)
library(dplyr)
library(patchwork)

# ============================================================
# INPUT YOUR STATISTICS HERE
# ============================================================
# For paired designs, you need:
#   1. Mean of the measure at baseline
#   2. SD of WITHIN-ANIMAL DIFFERENCES (BL - SD), not group SD
#
# To calculate SD of differences from your data:
#   differences <- animal_BL_values - animal_SD_values
#   sd_diff <- sd(differences)
#
# If you only have group SDs, you can estimate:
#   sd_diff ~ sd_group * sqrt(2 * (1 - correlation))
# where correlation is typically 0.5-0.8 for repeated measures
# ============================================================

# Event rate statistics
event_rate_stats <- list(
  baseline_mean = 3.0,      # Mean events per ROI per epoch at BL (REPLACE)
  sd_difference = 1.5       # SD of (BL - SD) differences across animals (REPLACE)
)

# Overall correlation statistics  
correlation_stats <- list(
  baseline_mean = 0.05,     # Mean correlation at BL (REPLACE)
  sd_difference = 0.08      # SD of (BL - SD) differences across animals (REPLACE)
)

# Distance-binned correlation statistics
distance_corr_stats <- list(
  close = list(
    baseline_mean = 0.06,   # Mean correlation for close pairs at BL (REPLACE)
    sd_difference = 0.09    # SD of differences (REPLACE)
  ),
  medium = list(
    baseline_mean = 0.05,
    sd_difference = 0.08
  ),
  far = list(
    baseline_mean = 0.04,
    sd_difference = 0.07
  )
)

# Analysis parameters
effect_sizes <- c(0.20, 0.30, 0.50)  # 20%, 30%, 50% changes from baseline
alpha <- 0.05
power_target <- 0.80

# Output directory
output_dir <- "Power_Analysis"
dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# CALCULATE REQUIRED SAMPLE SIZES (PAIRED DESIGN)
# ============================================================

cat("=== POWER ANALYSIS: REPEATED MEASURES DESIGN ===\n")
cat("Design: Same animals measured at BL and SD (paired)\n")
cat("Question: How many animals needed to detect X% change?\n\n")

# Event Rate Analysis
cat("--- Event Rate ---\n")
cat("Baseline mean:", event_rate_stats$baseline_mean, "\n")
cat("SD of within-animal differences:", event_rate_stats$sd_difference, "\n\n")

event_results <- data.frame()
for (effect in effect_sizes) {
  # Expected difference = effect size * baseline mean
  expected_diff <- event_rate_stats$baseline_mean * effect
  
  # Cohen's d for paired test = mean difference / SD of differences
  cohens_d <- expected_diff / event_rate_stats$sd_difference
  
  pwr_result <- pwr.t.test(
    d = cohens_d,
    power = power_target,
    sig.level = alpha,
    type = "paired"
  )
  
  event_results <- bind_rows(event_results, data.frame(
    Effect_Pct = effect * 100,
    Expected_Diff = round(expected_diff, 2),
    Cohens_d = round(cohens_d, 3),
    N_Animals = ceiling(pwr_result$n)
  ))
}

cat("Required sample sizes (80% power, α=0.05, paired design):\n")
print(event_results)
cat("\n")

# Overall Correlation Analysis
cat("--- Overall Correlation ---\n")
cat("Baseline mean:", correlation_stats$baseline_mean, "\n")
cat("SD of within-animal differences:", correlation_stats$sd_difference, "\n\n")

corr_results <- data.frame()
for (effect in effect_sizes) {
  expected_diff <- correlation_stats$baseline_mean * effect
  cohens_d <- expected_diff / correlation_stats$sd_difference
  
  pwr_result <- pwr.t.test(
    d = cohens_d,
    power = power_target,
    sig.level = alpha,
    type = "paired"
  )
  
  corr_results <- bind_rows(corr_results, data.frame(
    Effect_Pct = effect * 100,
    Expected_Diff = round(expected_diff, 3),
    Cohens_d = round(cohens_d, 3),
    N_Animals = ceiling(pwr_result$n)
  ))
}

cat("Required sample sizes (80% power, α=0.05, paired design):\n")
print(corr_results)
cat("\n")

# Distance-Binned Correlation Analysis
cat("--- Distance-Binned Correlations ---\n")

distance_results <- data.frame()
for (dist_bin in c("close", "medium", "far")) {
  cat("\n", toupper(dist_bin), "pairs:\n", sep="")
  stats <- distance_corr_stats[[dist_bin]]
  cat("Baseline mean:", stats$baseline_mean, "\n")
  cat("SD of differences:", stats$sd_difference, "\n")
  
  for (effect in effect_sizes) {
    expected_diff <- stats$baseline_mean * effect
    cohens_d <- expected_diff / stats$sd_difference
    
    pwr_result <- pwr.t.test(
      d = cohens_d,
      power = power_target,
      sig.level = alpha,
      type = "paired"
    )
    
    distance_results <- bind_rows(distance_results, data.frame(
      Distance_Bin = dist_bin,
      Effect_Pct = effect * 100,
      Expected_Diff = round(expected_diff, 3),
      Cohens_d = round(cohens_d, 3),
      N_Animals = ceiling(pwr_result$n)
    ))
  }
}

cat("\nRequired sample sizes by distance bin (80% power, α=0.05, paired):\n")
print(distance_results)

# Save all results
write.csv(event_results, 
          file.path(output_dir, "event_rate_power_paired.csv"),
          row.names = FALSE)
write.csv(corr_results,
          file.path(output_dir, "correlation_power_paired.csv"),
          row.names = FALSE)
write.csv(distance_results,
          file.path(output_dir, "distance_correlation_power_paired.csv"),
          row.names = FALSE)

# ============================================================
# GENERATE POWER CURVES
# ============================================================

cat("\n=== Generating Figures ===\n")

generate_power_curve <- function(cohens_d, n_range = 2:30) {
  power_values <- sapply(n_range, function(n) {
    pwr.t.test(d = cohens_d, n = n, sig.level = alpha, type = "paired")$power
  })
  data.frame(n = n_range, power = power_values)
}

# ============================================================
# FIGURE 1: Event Rate Power Curves
# ============================================================

event_curves <- data.frame()
for (i in 1:nrow(event_results)) {
  curve_data <- generate_power_curve(event_results$Cohens_d[i])
  curve_data$effect <- paste0(event_results$Effect_Pct[i], "%")
  event_curves <- bind_rows(event_curves, curve_data)
}

p_event <- ggplot(event_curves, aes(x = n, y = power, color = effect)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  annotate("text", x = 25, y = 0.82, label = "80% Power", size = 3.5) +
  geom_vline(xintercept = 3, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = 4.5, y = 0.15, label = "N=3\n(this study)", 
           size = 3, color = "red") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(
    title = "Power Analysis: Event Rate (Paired Design)",
    subtitle = paste0("Baseline: ", event_rate_stats$baseline_mean, 
                      ", SD of differences: ", event_rate_stats$sd_difference),
    x = "Number of Animals",
    y = "Statistical Power",
    color = "Effect Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Event_Rate_Power_Paired.png"),
       p_event, width = 10, height = 6, dpi = 300)
cat("  ✓ Event rate power curve (paired)\n")

# ============================================================
# FIGURE 2: Correlation Power Curves
# ============================================================

corr_curves <- data.frame()
for (i in 1:nrow(corr_results)) {
  curve_data <- generate_power_curve(corr_results$Cohens_d[i])
  curve_data$effect <- paste0(corr_results$Effect_Pct[i], "%")
  corr_curves <- bind_rows(corr_curves, curve_data)
}

p_corr <- ggplot(corr_curves, aes(x = n, y = power, color = effect)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  annotate("text", x = 25, y = 0.82, label = "80% Power", size = 3.5) +
  geom_vline(xintercept = 3, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = 4.5, y = 0.15, label = "N=3\n(this study)", 
           size = 3, color = "red") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(
    title = "Power Analysis: Pairwise Correlation (Paired Design)",
    subtitle = paste0("Baseline: ", correlation_stats$baseline_mean,
                      ", SD of differences: ", correlation_stats$sd_difference),
    x = "Number of Animals",
    y = "Statistical Power",
    color = "Effect Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Correlation_Power_Paired.png"),
       p_corr, width = 10, height = 6, dpi = 300)
cat("  ✓ Correlation power curve (paired)\n")

# ============================================================
# FIGURE 3: Sample Size Comparison Bar Chart
# ============================================================

comparison_30 <- bind_rows(
  event_results %>% filter(Effect_Pct == 30) %>% mutate(Measure = "Event Rate"),
  corr_results %>% filter(Effect_Pct == 30) %>% mutate(Measure = "Correlation")
)

p_comparison <- ggplot(comparison_30, aes(x = Measure, y = N_Animals, fill = Measure)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N_Animals), vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  annotate("text", x = 1.5, y = 4, label = "This study (N=3)", color = "red", size = 4) +
  scale_fill_manual(values = c("Event Rate" = "#2E86AB", "Correlation" = "#A23B72")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Required Sample Size for 30% Effect",
    subtitle = "Paired design, 80% power, α = 0.05",
    x = NULL,
    y = "Number of Animals Required"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 13)
  )

ggsave(file.path(output_dir, "Sample_Size_Comparison_Paired.png"),
       p_comparison, width = 8, height = 6, dpi = 300)
cat("  ✓ Sample size comparison\n")

# ============================================================
# FIGURE 4: Distance Bin Comparison
# ============================================================

p_distance <- ggplot(distance_results %>% filter(Effect_Pct == 30),
                     aes(x = factor(Distance_Bin, levels = c("close", "medium", "far")), 
                         y = N_Animals, 
                         fill = Distance_Bin)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N_Animals), vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  annotate("text", x = 0.6, y = 4.5, label = "This study (N=3)", color = "red", size = 4, hjust = 0) +
  scale_fill_manual(values = c("close" = "#3A0557", "medium" = "#185FA5", "far" = "#0F6E56")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Required Sample Size by Distance Bin (30% Effect)",
    subtitle = "Paired design, 80% power, α = 0.05",
    x = "Distance Bin",
    y = "Number of Animals Required"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 13)
  )

ggsave(file.path(output_dir, "Distance_Bin_Comparison_Paired.png"),
       p_distance, width = 8, height = 6, dpi = 300)
cat("  ✓ Distance bin comparison\n")

# ============================================================
# FIGURE 5: Power at Current N (What Effects Could We Detect?)
# ============================================================

# Calculate achieved power at N=3 for various effect sizes
current_n <- 3
effect_range <- seq(0.1, 1.0, by = 0.05)

achieved_power <- data.frame()

# Event rate
for (eff in effect_range) {
  expected_diff <- event_rate_stats$baseline_mean * eff
  d <- expected_diff / event_rate_stats$sd_difference
  pwr <- pwr.t.test(d = d, n = current_n, sig.level = alpha, type = "paired")$power
  achieved_power <- bind_rows(achieved_power, data.frame(
    Measure = "Event Rate", Effect_Pct = eff * 100, Power = pwr
  ))
}

# Correlation
for (eff in effect_range) {
  expected_diff <- correlation_stats$baseline_mean * eff
  d <- expected_diff / correlation_stats$sd_difference
  pwr <- pwr.t.test(d = d, n = current_n, sig.level = alpha, type = "paired")$power
  achieved_power <- bind_rows(achieved_power, data.frame(
    Measure = "Correlation", Effect_Pct = eff * 100, Power = pwr
  ))
}

p_achieved <- ggplot(achieved_power, aes(x = Effect_Pct, y = Power, color = Measure)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  annotate("text", x = 80, y = 0.82, label = "80% Power", size = 3.5) +
  geom_vline(xintercept = 30, linetype = "dotted", color = "gray60") +
  annotate("text", x = 32, y = 0.3, label = "30%", size = 3, color = "gray40", hjust = 0) +
  scale_color_manual(values = c("Event Rate" = "#2E86AB", "Correlation" = "#A23B72")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  labs(
    title = paste0("Statistical Power at Current Sample Size (N=", current_n, ")"),
    subtitle = "What effect sizes could this study reliably detect?",
    x = "Effect Size (% Change from Baseline)",
    y = "Statistical Power",
    color = "Measure"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Achieved_Power_at_N3.png"),
       p_achieved, width = 10, height = 6, dpi = 300)
cat("  ✓ Achieved power at N=3\n")

# ============================================================
# COMBINED FIGURE
# ============================================================

p_combined <- (p_event + p_corr) / (p_comparison + p_distance) / p_achieved +
  plot_annotation(
    title = "Power Analysis Summary: Repeated Measures Design",
    subtitle = "Within-subject comparison (BL vs SD)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

ggsave(file.path(output_dir, "Power_Analysis_Combined_Paired.png"),
       p_combined, width = 16, height = 16, dpi = 300)
cat("  ✓ Combined summary figure\n")

# ============================================================
# SUMMARY TABLE
# ============================================================

summary_table <- data.frame(
  Measure = c("Event Rate", "Correlation"),
  Baseline_Mean = c(event_rate_stats$baseline_mean, correlation_stats$baseline_mean),
  SD_Difference = c(event_rate_stats$sd_difference, correlation_stats$sd_difference),
  N_for_20pct = c(
    event_results$N_Animals[event_results$Effect_Pct == 20],
    corr_results$N_Animals[corr_results$Effect_Pct == 20]
  ),
  N_for_30pct = c(
    event_results$N_Animals[event_results$Effect_Pct == 30],
    corr_results$N_Animals[corr_results$Effect_Pct == 30]
  ),
  N_for_50pct = c(
    event_results$N_Animals[event_results$Effect_Pct == 50],
    corr_results$N_Animals[corr_results$Effect_Pct == 50]
  )
)

# Add power achieved at N=3 for 30% effect
summary_table$Power_at_N3_30pct <- c(
  achieved_power$Power[achieved_power$Measure == "Event Rate" & 
                        achieved_power$Effect_Pct == 30],
  achieved_power$Power[achieved_power$Measure == "Correlation" & 
                        achieved_power$Effect_Pct == 30]
)

write.csv(summary_table,
          file.path(output_dir, "power_analysis_summary_paired.csv"),
          row.names = FALSE)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n=== SUMMARY TABLE ===\n")
print(summary_table)

cat("\n=== KEY FINDINGS ===\n")
cat("\nFor 80% power to detect a 30% change (paired design):\n")
cat("  Event Rate: N =", event_results$N_Animals[event_results$Effect_Pct == 30], "animals\n")
cat("  Correlation: N =", corr_results$N_Animals[corr_results$Effect_Pct == 30], "animals\n")
cat("\n  Distance-binned correlations:\n")
for (dbin in c("close", "medium", "far")) {
  n_req <- distance_results$N_Animals[distance_results$Distance_Bin == dbin & 
                                       distance_results$Effect_Pct == 30]
  cat("    ", tools::toTitleCase(dbin), " pairs: N =", n_req, "animals\n")
}

cat("\nPower achieved with current N=3:\n")
cat("  Event Rate (30% effect):", 
    round(achieved_power$Power[achieved_power$Measure == "Event Rate" & 
                                achieved_power$Effect_Pct == 30] * 100, 1), "%\n")
cat("  Correlation (30% effect):", 
    round(achieved_power$Power[achieved_power$Measure == "Correlation" & 
                                achieved_power$Effect_Pct == 30] * 100, 1), "%\n")

cat("\nMinimum detectable effect at 80% power with N=3:\n")
# Find effect size where power crosses 0.8
event_mde <- achieved_power %>% 
  filter(Measure == "Event Rate", Power >= 0.79) %>% 
  slice_min(Effect_Pct) %>% 
  pull(Effect_Pct)
corr_mde <- achieved_power %>% 
  filter(Measure == "Correlation", Power >= 0.79) %>% 
  slice_min(Effect_Pct) %>% 
  pull(Effect_Pct)

if (length(event_mde) > 0) {
  cat("  Event Rate:", event_mde, "% change\n")
} else {
  cat("  Event Rate: >100% change (underpowered even for large effects)\n")
}
if (length(corr_mde) > 0) {
  cat("  Correlation:", corr_mde, "% change\n")
} else {
  cat("  Correlation: >100% change (underpowered even for large effects)\n")
}

cat("\n=== INTERPRETATION ===\n")
cat("• Paired design requires fewer animals than independent groups\n")
cat("• Current N=3 is underpowered for moderate (30%) effects\n")
cat("• Would need ~", max(summary_table$N_for_30pct), " animals for adequate power\n", sep="")
cat("• Null results may reflect insufficient power, not absence of effect\n")
cat("• Individual variability observed is consistent with underpowered design\n")

cat("\nFiles saved to:", output_dir, "\n")
cat("Generated files:\n")
print(list.files(output_dir))
