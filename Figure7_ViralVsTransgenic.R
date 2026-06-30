# ============================================================
# FIGURE 7: Viral vs Transgenic Correlation Magnitude
# ============================================================
# Compares mean pairwise correlation between expression types:
#   - Transgenic (CaMKIIa-Cre × jGCaMP8m): mPFCf5, mPFCf6, mPFCm9
#   - Viral (AAV-CaMKIIa-GCaMP6f): mPFCm4
#
# Input options:
#   1. CCF matrix CSVs from CCF2 output
#   2. Trajectory data CSVs from batch processing
#   3. Raw traces (computes correlations fresh)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# CONFIGURATION
# ============================================================

# Data source: "ccf_matrices", "trajectory_csv", or "raw_traces"
DATA_SOURCE <- "ccf_matrices"

# Paths
DATA_DIR <- "E:/Data_Processing/R/Data CSVs"
RESULTS_DIR <- "E:/Data_Processing/R/Results"
OUTPUT_DIR <- "E:/Data_Processing/R/Results/Figure7"

# Animals and their expression types
ANIMAL_INFO <- data.frame(
  animal = c("mPFCf5", "mPFCf6", "mPFCm9", "mPFCm4"),
  expression_type = c("Transgenic", "Transgenic", "Transgenic", "Viral"),
  stringsAsFactors = FALSE
)

condition <- c("BL", "SD", "WO")
animal = c("mPFCf5", "mPFCf6", "mPFCm9", "mPFCm4")

# Plot settings
EXPRESSION_COLORS <- c("Transgenic" = "#3498DB", "Viral" = "#E74C3C")
CONDITION_SHAPES <- c("BL" = 16, "SD" = 17, "WO" = 15)


# ============================================================
# DATA LOADING FUNCTIONS
# ============================================================

#' Load mean correlations from CCF matrix CSVs
load_from_ccf_matrices <- function(data_dir, animals, conditions) {
  
  all_correlations <- data.frame()
  
  for (animal in animal) {
    if (animal == 'mPFCf5') {
      for (condition0 in condition) {
        
        # Try different possible file paths
        ccf_filedir <- paste0(DATA_DIR, "/", animal, "_", condition0, "_NREM_ccf_matrix.csv")
  
        
        # Load CCF matrix
        ccf_matrix <- read.csv(ccf_filedir, row.names = 1, check.names = FALSE)
        ccf_matrix <- as.matrix(ccf_matrix)
        
        # Extract upper triangle (excluding diagonal)
        upper_tri <- ccf_matrix[upper.tri(ccf_matrix, diag = FALSE)]
        
        # Compute summary statistics
        mean_corr <- mean(upper_tri, na.rm = TRUE)
        median_corr <- median(upper_tri, na.rm = TRUE)
        sd_corr <- sd(upper_tri, na.rm = TRUE)
        n_pairs <- length(upper_tri)
        
        all_correlations <- rbind(all_correlations, data.frame(
          animal = animal,
          condition = condition0,
          mean_correlation = mean_corr,
          median_correlation = median_corr,
          sd_correlation = sd_corr,
          n_pairs = n_pairs
        ))
      }
    }
    if (animal == 'mPFCf6') {
       for (condition1 in condition) {
          
          # Try different possible file paths
          ccf_filedir <- paste0(DATA_DIR, "/", animal, "_", condition1, "_NREM_ccf_matrix.csv")
          
          
          # Load CCF matrix
          ccf_matrix <- read.csv(ccf_filedir, row.names = 1, check.names = FALSE)
          ccf_matrix <- as.matrix(ccf_matrix)
          
          # Extract upper triangle (excluding diagonal)
          upper_tri <- ccf_matrix[upper.tri(ccf_matrix, diag = FALSE)]
          
          # Compute summary statistics
          mean_corr <- mean(upper_tri, na.rm = TRUE)
          median_corr <- median(upper_tri, na.rm = TRUE)
          sd_corr <- sd(upper_tri, na.rm = TRUE)
          n_pairs <- length(upper_tri)
          
          all_correlations <- rbind(all_correlations, data.frame(
            animal = animal,
            condition = condition1,
            mean_correlation = mean_corr,
            median_correlation = median_corr,
            sd_correlation = sd_corr,
            n_pairs = n_pairs
          ))
      }
    }
    if (animal == 'mPFCm9') {
      for (condition2 in condition) {
        
        # Try different possible file paths
        ccf_filedir <- paste0(DATA_DIR, "/", animal, "_", condition2, "_NREM_ccf_matrix.csv")
        
        
        # Load CCF matrix
        ccf_matrix <- read.csv(ccf_filedir, row.names = 1, check.names = FALSE)
        ccf_matrix <- as.matrix(ccf_matrix)
        
        # Extract upper triangle (excluding diagonal)
        upper_tri <- ccf_matrix[upper.tri(ccf_matrix, diag = FALSE)]
        
        # Compute summary statistics
        mean_corr <- mean(upper_tri, na.rm = TRUE)
        median_corr <- median(upper_tri, na.rm = TRUE)
        sd_corr <- sd(upper_tri, na.rm = TRUE)
        n_pairs <- length(upper_tri)
        
        all_correlations <- rbind(all_correlations, data.frame(
          animal = animal,
          condition = condition2,
          mean_correlation = mean_corr,
          median_correlation = median_corr,
          sd_correlation = sd_corr,
          n_pairs = n_pairs
        ))
      }
    }
    if (animal == 'mPFCm4') {
      for (condition3 in condition) {
        
        # Try different possible file paths
        ccf_filedir <- paste0(DATA_DIR, "/", animal, "_", condition3, "_NREM_ccf_matrix.csv")
        
        
        # Load CCF matrix
        ccf_matrix <- read.csv(ccf_filedir, row.names = 1, check.names = FALSE)
        ccf_matrix <- as.matrix(ccf_matrix)
        
        # Extract upper triangle (excluding diagonal)
        upper_tri <- ccf_matrix[upper.tri(ccf_matrix, diag = FALSE)]
        
        # Compute summary statistics
        mean_corr <- mean(upper_tri, na.rm = TRUE)
        median_corr <- median(upper_tri, na.rm = TRUE)
        sd_corr <- sd(upper_tri, na.rm = TRUE)
        n_pairs <- length(upper_tri)
        
        all_correlations <- rbind(all_correlations, data.frame(
          animal = animal,
          condition = condition3,
          mean_correlation = mean_corr,
          median_correlation = median_corr,
          sd_correlation = sd_corr,
          n_pairs = n_pairs
        ))
      }
    }
  }
  
  return(all_correlations)
}


#' Load correlations from batch processing trajectory CSV
load_from_trajectory_csv <- function(results_dir, filename = "CrossAnimal_Trajectory_All.csv") {
  
  filepath <- file.path(results_dir, "between_subject", filename)
  
  if (!file.exists(filepath)) {
    # Try with suffix
    possible_files <- list.files(file.path(results_dir, "between_subject"), 
                                  pattern = "CrossAnimal_Trajectory_All.*\\.csv$",
                                  full.names = TRUE)
    if (length(possible_files) > 0) {
      filepath <- possible_files[1]
    } else {
      stop("Trajectory CSV not found: ", filepath)
    }
  }
  
  cat("Loading trajectory data from:", filepath, "\n")
  traj_data <- read.csv(filepath)
  
  # Summarize by animal and condition
  all_correlations <- traj_data %>%
    group_by(animal, condition) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )
  
  return(as.data.frame(all_correlations))
}


#' Compute correlations from raw traces
load_from_raw_traces <- function(data_dir, animals, conditions) {
  
  all_correlations <- data.frame()
  
  for (animal in animals) {
    for (condition in conditions) {
      
      traces_file <- file.path(data_dir, paste0(animal, "_", condition, "_Traces.csv"))
      
      if (!file.exists(traces_file)) {
        cat("Warning: Traces file not found:", traces_file, "\n")
        next
      }
      
      cat("Loading:", traces_file, "\n")
      
      # Load traces
      traces <- read.csv(traces_file, check.names = FALSE)
      
      # Clean up - remove Time column and header row
      colnames(traces)[1] <- "Time"
      traces <- traces[-1, ]  # Remove header row
      traces <- traces[, -1]  # Remove Time column
      
      # Convert to numeric matrix
      traces_matrix <- as.matrix(sapply(traces, as.numeric))
      
      # Compute correlation matrix
      corr_matrix <- cor(traces_matrix, use = "pairwise.complete.obs")
      
      # Extract upper triangle
      upper_tri <- corr_matrix[upper.tri(corr_matrix, diag = FALSE)]
      
      # Summary statistics
      mean_corr <- mean(upper_tri, na.rm = TRUE)
      median_corr <- median(upper_tri, na.rm = TRUE)
      sd_corr <- sd(upper_tri, na.rm = TRUE)
      n_pairs <- sum(!is.na(upper_tri))
      
      all_correlations <- rbind(all_correlations, data.frame(
        animal = animal,
        condition = condition,
        mean_correlation = mean_corr,
        median_correlation = median_corr,
        sd_correlation = sd_corr,
        n_pairs = n_pairs
      ))
    }
  }
  
  return(all_correlations)
}


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

#' Create bar plot with individual animals colored by expression type
create_barplot_by_animal <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal") %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO"),
                         labels = c("Baseline", "Sleep Dep", "Recovery")),
      expression_type = factor(expression_type, levels = c("Transgenic", "Viral"))
    )
  
  # Compute group means for reference line
  group_means <- plot_data %>%
    group_by(expression_type) %>%
    summarise(group_mean = mean(mean_correlation, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = animal, y = mean_correlation, fill = expression_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    facet_wrap(~ condition, nrow = 1) +
    scale_fill_manual(values = EXPRESSION_COLORS, name = "Expression Type") +
    labs(
      title = "Mean Pairwise Correlation by Animal and Condition",
      x = NULL,
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "top"
    )
  
  return(p)
}


#' Create grouped bar plot: expression type on x-axis, faceted by condition
create_barplot_grouped <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal") %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO"),
                         labels = c("Baseline", "Sleep Dep", "Recovery")),
      expression_type = factor(expression_type, levels = c("Transgenic", "Viral"))
    )
  
  # Summary by expression type and condition
  summary_data <- plot_data %>%
    group_by(expression_type, condition) %>%
    summarise(
      mean_corr = mean(mean_correlation, na.rm = TRUE),
      se_corr = sd(mean_correlation, na.rm = TRUE) / sqrt(n()),
      n_animals = n(),
      .groups = "drop"
    )
  
  p <- ggplot(summary_data, aes(x = expression_type, y = mean_corr, fill = expression_type)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean_corr - se_corr, ymax = mean_corr + se_corr),
                  width = 0.2) +
    geom_point(data = plot_data, aes(x = expression_type, y = mean_correlation),
               position = position_jitter(width = 0.1), size = 2, alpha = 0.7) +
    facet_wrap(~ condition, nrow = 1) +
    scale_fill_manual(values = EXPRESSION_COLORS, guide = "none") +
    labs(
      title = "Mean Pairwise Correlation: Viral vs Transgenic",
      subtitle = "Individual animals shown as points; bars show group mean ± SEM",
      x = "Expression Type",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      strip.text = element_text(size = 11, face = "bold")
    )
  
  return(p)
}


#' Create violin/box plot comparing distributions
create_violin_comparison <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal") %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO"),
                         labels = c("Baseline", "Sleep Dep", "Recovery")),
      expression_type = factor(expression_type, levels = c("Transgenic", "Viral"))
    )
  
  p <- ggplot(plot_data, aes(x = expression_type, y = mean_correlation, fill = expression_type)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
    geom_point(aes(shape = condition), size = 3, position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = EXPRESSION_COLORS, guide = "none") +
    scale_shape_manual(values = CONDITION_SHAPES, name = "Condition") +
    labs(
      title = "Correlation Magnitude by Expression Type",
      x = "Expression Type",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}


#' Create dot plot with lines connecting conditions within animals
create_dotplot_connected <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal") %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO")),
      expression_type = factor(expression_type, levels = c("Transgenic", "Viral"))
    )
  
  p <- ggplot(plot_data, aes(x = condition, y = mean_correlation, 
                              color = expression_type, group = animal)) +
    geom_line(alpha = 0.6, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = EXPRESSION_COLORS, name = "Expression Type") +
    labs(
      title = "Correlation Trajectories Across Conditions",
      x = "Condition",
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}


#' Create publication-ready figure (thesis Figure 7)
create_figure7 <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal") %>%
    mutate(
      condition = factor(condition, levels = c("BL", "SD", "WO"),
                         labels = c("Baseline", "Sleep Dep", "Recovery")),
      expression_type = factor(expression_type, levels = c("Transgenic", "Viral"))
    )
  
  # Compute fold difference
  viral_mean <- mean(plot_data$mean_correlation[plot_data$expression_type == "Viral"], na.rm = TRUE)
  trans_mean <- mean(plot_data$mean_correlation[plot_data$expression_type == "Transgenic"], na.rm = TRUE)
  fold_diff <- round(viral_mean / trans_mean, 1)
  
  p <- ggplot(plot_data, aes(x = animal, y = mean_correlation, fill = expression_type)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.7) +
    geom_point(aes(shape = condition), size = 2.5, 
               position = position_dodge(width = 0.3)) +
    scale_fill_manual(values = EXPRESSION_COLORS, name = "Expression Type") +
    scale_shape_manual(values = c("Baseline" = 16, "Sleep Dep" = 17, "Recovery" = 15),
                       name = "Condition") +
    labs(
      title = "Viral vs Transgenic Correlation Magnitude",
      subtitle = paste0("Viral expression shows ~", fold_diff, "× mean correlation (SNR equivalent, see Figure 5)"),
      x = NULL,
      y = "Mean Pairwise Correlation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(size = 10),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    guides(fill = guide_legend(order = 1),
           shape = guide_legend(order = 2))
  
  return(p)
}


# ============================================================
# STATISTICS
# ============================================================

compute_statistics <- function(data, animal_info) {
  
  plot_data <- data %>%
    left_join(animal_info, by = "animal")
  
  cat("\n============================================\n")
  cat("STATISTICAL SUMMARY: Viral vs Transgenic\n")
  cat("============================================\n\n")
  
  # Overall summary
  summary_by_type <- plot_data %>%
    group_by(expression_type) %>%
    summarise(
      n_observations = n(),
      n_animals = n_distinct(animal),
      mean_corr = mean(mean_correlation, na.rm = TRUE),
      sd_corr = sd(mean_correlation, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Overall Summary:\n")
  print(summary_by_type)
  
  # Fold difference
  viral_mean <- summary_by_type$mean_corr[summary_by_type$expression_type == "Viral"]
  trans_mean <- summary_by_type$mean_corr[summary_by_type$expression_type == "Transgenic"]
  fold_diff <- viral_mean / trans_mean
  
  cat("\nFold difference (Viral / Transgenic):", round(fold_diff, 2), "\n")
  
  # By condition
  cat("\nBy Condition:\n")
  summary_by_condition <- plot_data %>%
    group_by(expression_type, condition) %>%
    summarise(
      mean_corr = mean(mean_correlation, na.rm = TRUE),
      sd_corr = sd(mean_correlation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = expression_type, values_from = c(mean_corr, sd_corr))
  
  print(summary_by_condition)
  
  # T-test (note: very low power with n=3 vs n=1)
  cat("\nNote: Statistical tests have very low power (n=3 transgenic, n=1 viral)\n")
  cat("The difference is presented as a methodological observation, not a hypothesis test.\n")
  
  return(list(
    summary_by_type = summary_by_type,
    summary_by_condition = summary_by_condition,
    fold_difference = fold_diff
  ))
}


# ============================================================
# MAIN EXECUTION
# ============================================================

run_figure7_analysis <- function(data_source = DATA_SOURCE,
                                  data_dir = DATA_DIR,
                                  results_dir = RESULTS_DIR,
                                  output_dir = OUTPUT_DIR,
                                  animal_info = ANIMAL_INFO) {
  
  cat("============================================\n")
  cat("FIGURE 7: Viral vs Transgenic Analysis\n")
  cat("============================================\n")
  cat("Data source:", data_source, "\n\n")
  
  # Load data based on source
  if (data_source == "ccf_matrices") {
    corr_data <- load_from_ccf_matrices(data_dir, animal_info$animal, CONDITIONS)
  } else if (data_source == "trajectory_csv") {
    corr_data <- load_from_trajectory_csv(results_dir)
  } else if (data_source == "raw_traces") {
    corr_data <- load_from_raw_traces(data_dir, animal_info$animal, CONDITIONS)
  } else {
    stop("Unknown data source: ", data_source)
  }
  
  cat("\nLoaded correlation data:\n")
  print(corr_data)
  
  # Compute statistics
  stats <- compute_statistics(corr_data, animal_info)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create plots
  cat("\nCreating plots...\n")
  
  plots <- list(
    barplot_by_animal = create_barplot_by_animal(corr_data, animal_info),
    barplot_grouped = create_barplot_grouped(corr_data, animal_info),
    violin = create_violin_comparison(corr_data, animal_info),
    dotplot_connected = create_dotplot_connected(corr_data, animal_info),
    figure7 = create_figure7(corr_data, animal_info)
  )
  
  # Save plots
  cat("Saving plots to:", output_dir, "\n")
  
  ggsave(file.path(output_dir, "Figure7_ViralVsTransgenic_ByAnimal.png"),
         plots$barplot_by_animal, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(output_dir, "Figure7_ViralVsTransgenic_Grouped.png"),
         plots$barplot_grouped, width = 7.5, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "Figure7_ViralVsTransgenic_Violin.png"),
         plots$violin, width = 6, height = 5, dpi = 300)
  
  ggsave(file.path(output_dir, "Figure7_ViralVsTransgenic_Connected.png"),
         plots$dotplot_connected, width = 6, height = 5, dpi = 300)
  
  ggsave(file.path(output_dir, "Figure7_MAIN.png"),
         plots$figure7, width = 7.5, height = 6, dpi = 300)
  
  # Save data
  write.csv(corr_data, file.path(output_dir, "Figure7_correlation_data.csv"), row.names = FALSE)
  
  cat("\n============================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("Output saved to:", output_dir, "\n")
  cat("============================================\n")
  
  return(list(
    data = corr_data,
    stats = stats,
    plots = plots
  ))
}


# ============================================================
# AUTO-RUN
# ============================================================

if (interactive()) {
  cat("\n")
  cat("========================================\n")
  cat("  FIGURE 7 SCRIPT LOADED\n")
  cat("  Viral vs Transgenic Comparison\n")
  cat("========================================\n")
  cat("\nTo run analysis, execute:\n")
  cat("  results <- run_figure7_analysis()\n\n")
  cat("Data source options:\n")
  cat("  - 'ccf_matrices': Load from CCF2 output CSVs\n")
  cat("  - 'trajectory_csv': Load from batch processing output\n")
  cat("  - 'raw_traces': Compute correlations from trace files\n\n")
  cat("Example with different source:\n")
  cat("  results <- run_figure7_analysis(data_source = 'trajectory_csv')\n\n")
}
