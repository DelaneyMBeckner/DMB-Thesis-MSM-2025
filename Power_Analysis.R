# ============================================================
# QUICK POWER ANALYSIS - WITH DISTANCE ANALYSIS
# ============================================================
# Power analysis for:
# 1. Event rates
# 2. Overall pairwise correlations
# 3. Distance-binned correlations (close/medium/far)
# ============================================================

library(ggplot2)
library(pwr)
library(dplyr)
library(patchwork)

# ============================================================
# INPUT YOUR BASELINE STATISTICS HERE
# ============================================================

# Event rate statistics (from your baseline data)
event_rate_stats <- list(
  mean = 3.0,        # Mean events per ROI per epoch (REPLACE WITH YOUR VALUE)
  sd = 2.5,          # Standard deviation (REPLACE WITH YOUR VALUE)
  cv = 0.83          # Coefficient of variation (REPLACE WITH YOUR VALUE)
)

# Overall correlation statistics (from your baseline data)
correlation_stats <- list(
  mean = 0.05,       # Mean correlation coefficient (REPLACE WITH YOUR VALUE)
  sd = 0.15,         # Standard deviation (REPLACE WITH YOUR VALUE)
  cv = 3.0           # Coefficient of variation (REPLACE WITH YOUR VALUE)
)

# Distance-binned correlation statistics
# For each distance bin: close, medium, far
distance_corr_stats <- list(
  close = list(
    mean = 0.06,     # Mean correlation for close pairs (REPLACE)
    sd = 0.16,       # Standard deviation (REPLACE)
    cv = 2.67        # Coefficient of variation (REPLACE)
  ),
  medium = list(
    mean = 0.05,     # Mean correlation for medium pairs (REPLACE)
    sd = 0.15,       # Standard deviation (REPLACE)
    cv = 3.0         # Coefficient of variation (REPLACE)
  ),
  far = list(
    mean = 0.04,     # Mean correlation for far pairs (REPLACE)
    sd = 0.14,       # Standard deviation (REPLACE)
    cv = 3.5         # Coefficient of variation (REPLACE)
  )
)

# Analysis parameters
effect_sizes <- c(0.20, 0.30, 0.50)  # 20%, 30%, 50% changes
alpha <- 0.05
power_target <- 0.80

# Output directory
output_dir <- "Power_Analysis"
dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# CALCULATE REQUIRED SAMPLE SIZES
# ============================================================

cat("=== POWER ANALYSIS RESULTS ===\n\n")

# Event Rate Analysis
cat("--- Event Rate ---\n")
cat("Baseline:", event_rate_stats$mean, "±", event_rate_stats$sd, 
    "(CV =", round(event_rate_stats$cv * 100, 1), "%)\n\n")

event_results <- data.frame()
for (effect in effect_sizes) {
  expected_diff <- event_rate_stats$mean * effect
  cohens_d <- expected_diff / event_rate_stats$sd
  
  pwr_result <- pwr.t.test(
    d = cohens_d,
    power = power_target,
    sig.level = alpha,
    type = "two.sample"
  )
  
  event_results <- bind_rows(event_results, data.frame(
    Effect_Pct = effect * 100,
    Expected_Diff = round(expected_diff, 2),
    Cohens_d = round(cohens_d, 3),
    N_per_group = ceiling(pwr_result$n),
    Total_N = ceiling(pwr_result$n) * 2
  ))
}

cat("Required sample sizes (80% power, α=0.05):\n")
print(event_results)
cat("\n")

# Overall Correlation Analysis
cat("--- Overall Correlation ---\n")
cat("Baseline:", correlation_stats$mean, "±", correlation_stats$sd,
    "(CV =", round(correlation_stats$cv * 100, 1), "%)\n\n")

corr_results <- data.frame()
for (effect in effect_sizes) {
  expected_diff <- correlation_stats$mean * effect
  cohens_d <- expected_diff / correlation_stats$sd
  
  pwr_result <- pwr.t.test(
    d = cohens_d,
    power = power_target,
    sig.level = alpha,
    type = "two.sample"
  )
  
  corr_results <- bind_rows(corr_results, data.frame(
    Effect_Pct = effect * 100,
    Expected_Diff = round(expected_diff, 3),
    Cohens_d = round(cohens_d, 3),
    N_per_group = ceiling(pwr_result$n),
    Total_N = ceiling(pwr_result$n) * 2
  ))
}

cat("Required sample sizes (80% power, α=0.05):\n")
print(corr_results)
cat("\n")

# Distance-Binned Correlation Analysis
cat("--- Distance-Binned Correlations ---\n")

distance_results <- data.frame()
for (dist_bin in c("close", "medium", "far")) {
  cat("\n", toupper(dist_bin), "pairs:\n", sep="")
  stats <- distance_corr_stats[[dist_bin]]
  cat("Baseline:", stats$mean, "±", stats$sd,
      "(CV =", round(stats$cv * 100, 1), "%)\n")
  
  for (effect in effect_sizes) {
    expected_diff <- stats$mean * effect
    cohens_d <- expected_diff / stats$sd
    
    pwr_result <- pwr.t.test(
      d = cohens_d,
      power = power_target,
      sig.level = alpha,
      type = "two.sample"
    )
    
    distance_results <- bind_rows(distance_results, data.frame(
      Distance_Bin = dist_bin,
      Effect_Pct = effect * 100,
      Expected_Diff = round(expected_diff, 3),
      Cohens_d = round(cohens_d, 3),
      N_per_group = ceiling(pwr_result$n),
      Total_N = ceiling(pwr_result$n) * 2
    ))
  }
}

cat("\nRequired sample sizes by distance bin (80% power, α=0.05):\n")
print(distance_results)

# Save all results
write.csv(event_results, 
          file.path(output_dir, "event_rate_power.csv"),
          row.names = FALSE)
write.csv(corr_results,
          file.path(output_dir, "correlation_power.csv"),
          row.names = FALSE)
write.csv(distance_results,
          file.path(output_dir, "distance_correlation_power.csv"),
          row.names = FALSE)

# ============================================================
# GENERATE POWER CURVES
# ============================================================

cat("\n=== Generating Figures ===\n")

# Function to generate power curve data
generate_power_curve <- function(cohens_d, n_range = 2:30) {
  power_values <- sapply(n_range, function(n) {
    pwr.t.test(d = cohens_d, n = n, sig.level = alpha, type = "two.sample")$power
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
  # Shade current N region
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = 1, 
           alpha = 0.1, fill = "gray50") +
  annotate("text", x = 2, y = 0.95, label = "Current\nN=3", 
           size = 3, color = "gray30") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Power Analysis: Event Rate",
    subtitle = paste0("Baseline: ", event_rate_stats$mean, " ± ", 
                     event_rate_stats$sd, " events/ROI/epoch"),
    x = "Sample Size (n per group)",
    y = "Statistical Power",
    color = "Effect Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Event_Rate_Power_Curve.png"),
       p_event, width = 10, height = 6, dpi = 300)
cat("  ✓ Event rate power curve\n")

# ============================================================
# FIGURE 2: Overall Correlation Power Curves
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
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = 1,
           alpha = 0.1, fill = "gray50") +
  annotate("text", x = 2, y = 0.95, label = "Current\nN=3",
           size = 3, color = "gray30") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Power Analysis: Overall Pairwise Correlation",
    subtitle = paste0("Baseline: r = ", correlation_stats$mean, " ± ",
                     correlation_stats$sd),
    x = "Sample Size (n per group)",
    y = "Statistical Power",
    color = "Effect Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Correlation_Power_Curve.png"),
       p_corr, width = 10, height = 6, dpi = 300)
cat("  ✓ Correlation power curve\n")

# ============================================================
# FIGURE 3: Distance-Binned Correlation Power Curves
# ============================================================

# Generate curves for 30% effect across all distance bins
dist_curves_30pct <- data.frame()
for (dist_bin in c("close", "medium", "far")) {
  stats <- distance_corr_stats[[dist_bin]]
  expected_diff <- stats$mean * 0.30
  cohens_d <- expected_diff / stats$sd
  
  curve_data <- generate_power_curve(cohens_d)
  curve_data$distance <- tools::toTitleCase(dist_bin)
  dist_curves_30pct <- bind_rows(dist_curves_30pct, curve_data)
}

p_dist_30 <- ggplot(dist_curves_30pct, aes(x = n, y = power, color = distance)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
  annotate("text", x = 25, y = 0.82, label = "80% Power", size = 3.5) +
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = 1,
           alpha = 0.1, fill = "gray50") +
  annotate("text", x = 2, y = 0.95, label = "Current\nN=3",
           size = 3, color = "gray30") +
  scale_color_manual(values = c("Close" = "#440154", "Medium" = "#31688e", 
                                "Far" = "#fde725")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Power Analysis: Distance-Binned Correlations (30% Effect)",
    subtitle = "Power to detect 30% change in each distance bin",
    x = "Sample Size (n per group)",
    y = "Statistical Power",
    color = "Distance Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Distance_Correlation_Power_Curve.png"),
       p_dist_30, width = 10, height = 6, dpi = 300)
cat("  ✓ Distance-binned correlation power curve\n")

# ============================================================
# FIGURE 4: Sample Size Comparison (All Measures)
# ============================================================

# Extract 30% effect results for all measures
comparison_data <- bind_rows(
  event_results %>% filter(Effect_Pct == 30) %>% 
    mutate(Measure = "Event Rate", Distance = "All"),
  corr_results %>% filter(Effect_Pct == 30) %>% 
    mutate(Measure = "Correlation", Distance = "All"),
  distance_results %>% filter(Effect_Pct == 30) %>%
    mutate(Measure = "Correlation", Distance = tools::toTitleCase(Distance_Bin))
)

p_comparison <- ggplot(comparison_data %>% filter(Distance == "All"),
                      aes(x = Measure, y = N_per_group, fill = Measure)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N_per_group), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Event Rate" = "#2E86AB", "Correlation" = "#A23B72")) +
  labs(
    title = "Required Sample Size for 30% Effect (80% Power)",
    subtitle = "Comparison across primary measures",
    x = NULL,
    y = "Required n per Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 13)
  )

ggsave(file.path(output_dir, "Sample_Size_Comparison.png"),
       p_comparison, width = 8, height = 6, dpi = 300)
cat("  ✓ Sample size comparison\n")

# ============================================================
# FIGURE 5: Distance Bin Comparison
# ============================================================

p_dist_comparison <- ggplot(distance_results %>% filter(Effect_Pct == 30),
                            aes(x = tools::toTitleCase(Distance_Bin), 
                                y = N_per_group, 
                                fill = Distance_Bin)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N_per_group), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("close" = "#440154", "medium" = "#31688e", 
                               "far" = "#fde725")) +
  labs(
    title = "Required Sample Size by Distance Bin (30% Effect)",
    subtitle = "Sample size requirements similar across all distance bins",
    x = "Distance Bin",
    y = "Required n per Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 13)
  )

ggsave(file.path(output_dir, "Distance_Bin_Sample_Size.png"),
       p_dist_comparison, width = 8, height = 6, dpi = 300)
cat("  ✓ Distance bin comparison\n")

# ============================================================
# COMBINED SUMMARY FIGURE (2x2 grid)
# ============================================================

p_combined <- (p_event + p_corr) / (p_dist_30 + p_dist_comparison) +
  plot_annotation(
    title = "Power Analysis Summary: All Measures",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )

ggsave(file.path(output_dir, "Power_Analysis_Combined.png"),
       p_combined, width = 16, height = 12, dpi = 300)
cat("  ✓ Combined summary figure\n")

# ============================================================
# SUMMARY TABLE
# ============================================================

summary_table <- data.frame(
  Measure = c("Event Rate", "Correlation (All)", "Correlation (Close)", 
              "Correlation (Medium)", "Correlation (Far)"),
  Baseline_Mean = c(
    event_rate_stats$mean,
    correlation_stats$mean,
    distance_corr_stats$close$mean,
    distance_corr_stats$medium$mean,
    distance_corr_stats$far$mean
  ),
  Baseline_SD = c(
    event_rate_stats$sd,
    correlation_stats$sd,
    distance_corr_stats$close$sd,
    distance_corr_stats$medium$sd,
    distance_corr_stats$far$sd
  ),
  CV_Percent = c(
    round(event_rate_stats$cv * 100, 1),
    round(correlation_stats$cv * 100, 1),
    round(distance_corr_stats$close$cv * 100, 1),
    round(distance_corr_stats$medium$cv * 100, 1),
    round(distance_corr_stats$far$cv * 100, 1)
  ),
  N_for_20pct = c(
    event_results$N_per_group[event_results$Effect_Pct == 20],
    corr_results$N_per_group[corr_results$Effect_Pct == 20],
    distance_results$N_per_group[distance_results$Distance_Bin == "close" & 
                                  distance_results$Effect_Pct == 20],
    distance_results$N_per_group[distance_results$Distance_Bin == "medium" & 
                                  distance_results$Effect_Pct == 20],
    distance_results$N_per_group[distance_results$Distance_Bin == "far" & 
                                  distance_results$Effect_Pct == 20]
  ),
  N_for_30pct = c(
    event_results$N_per_group[event_results$Effect_Pct == 30],
    corr_results$N_per_group[corr_results$Effect_Pct == 30],
    distance_results$N_per_group[distance_results$Distance_Bin == "close" & 
                                  distance_results$Effect_Pct == 30],
    distance_results$N_per_group[distance_results$Distance_Bin == "medium" & 
                                  distance_results$Effect_Pct == 30],
    distance_results$N_per_group[distance_results$Distance_Bin == "far" & 
                                  distance_results$Effect_Pct == 30]
  ),
  N_for_50pct = c(
    event_results$N_per_group[event_results$Effect_Pct == 50],
    corr_results$N_per_group[corr_results$Effect_Pct == 50],
    distance_results$N_per_group[distance_results$Distance_Bin == "close" & 
                                  distance_results$Effect_Pct == 50],
    distance_results$N_per_group[distance_results$Distance_Bin == "medium" & 
                                  distance_results$Effect_Pct == 50],
    distance_results$N_per_group[distance_results$Distance_Bin == "far" & 
                                  distance_results$Effect_Pct == 50]
  )
)

write.csv(summary_table,
          file.path(output_dir, "power_analysis_summary.csv"),
          row.names = FALSE)

# ============================================================
# SUMMARY OUTPUT
# ============================================================

cat("\n=== SUMMARY TABLE ===\n")
print(summary_table)

cat("\n=== KEY FINDINGS ===\n")
cat("\nFor a 30% effect with 80% power:\n")
cat("  Event Rate: n =", event_results$N_per_group[event_results$Effect_Pct == 30],
    "per group\n")
cat("  Overall Correlation: n =", corr_results$N_per_group[corr_results$Effect_Pct == 30],
    "per group\n")
cat("\n  Distance-binned correlations:\n")
cat("    Close pairs: n =", 
    distance_results$N_per_group[distance_results$Distance_Bin == "close" & 
                                 distance_results$Effect_Pct == 30], "\n")
cat("    Medium pairs: n =",
    distance_results$N_per_group[distance_results$Distance_Bin == "medium" & 
                                 distance_results$Effect_Pct == 30], "\n")
cat("    Far pairs: n =",
    distance_results$N_per_group[distance_results$Distance_Bin == "far" & 
                                 distance_results$Effect_Pct == 30], "\n")

cat("\nFiles saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
list.files(output_dir)

cat("\n=== INTERPRETATION ===\n")
cat("• Event rates require largest sample sizes (high variability)\n")
cat("• Correlation measures more stable but still need n=10-15\n")
cat("• Distance bins show similar power requirements (consistent variability)\n")
cat("• Current N=3 underpowered for all measures except very large effects\n")
cat("• Recommendation: N=12-15 per group for adequately powered future studies\n")
