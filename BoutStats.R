#!/usr/bin/env Rscript
# Bout_Analysis.R
# Diagnostic script to analyze bout characteristics across recordings
# Helps decide on transition analysis approach based on bout length distributions

library(data.table)

# ============================================================
# CONFIGURATION - Modify these paths as needed
# ============================================================
# Set working directory (comment out if running interactively)
# setwd("E:/Data_Processing/R")

DEFAULT_INPUT_PATH <- "./Data CSVs"
DEFAULT_OUTPUT_PATH <- "./bout_analysis_results"

# State label mapping
STATE_LABELS <- c("1" = "Wake", "2" = "NREM", "3" = "REM")

#' Identify bouts in a states dataframe
#' 
#' @param states_df Data frame with 'state' column
#' @return Data frame with bout_id, state, start_epoch, end_epoch, length
identify_bouts <- function(states_df) {
  
  if (nrow(states_df) == 0) {
    return(data.frame(
      bout_id = integer(),
      state = integer(),
      state_label = character(),
      start_epoch = integer(),
      end_epoch = integer(),
      length_epochs = integer(),
      length_seconds = numeric()
    ))
  }
  
  # Find where state changes occur
  state_changes <- c(TRUE, states_df$state[-1] != states_df$state[-nrow(states_df)])
  
  # Assign bout IDs
  bout_ids <- cumsum(state_changes)
  
  # Aggregate bout information
  bouts <- data.frame(
    epoch_num = seq_len(nrow(states_df)),
    state = states_df$state,
    bout_id = bout_ids
  )
  
  # Summarize each bout
  bout_summary <- aggregate(
    epoch_num ~ bout_id + state,
    data = bouts,
    FUN = function(x) c(start = min(x), end = max(x), length = length(x))
  )
  
  # Unpack the aggregated columns
  bout_df <- data.frame(
    bout_id = bout_summary$bout_id,
    state = bout_summary$state,
    state_label = STATE_LABELS[as.character(bout_summary$state)],
    start_epoch = bout_summary$epoch_num[, "start"],
    end_epoch = bout_summary$epoch_num[, "end"],
    length_epochs = bout_summary$epoch_num[, "length"],
    length_seconds = bout_summary$epoch_num[, "length"] * 10  # 10-second epochs
  )
  
  return(bout_df)
}


#' Create histograms of bout lengths
#' 
#' @param all_bout_data List of bout data frames from all recordings
#' @param output_path Directory for saving plots
create_bout_histograms <- function(all_bout_data, output_path = ".") {
  
  if (length(all_bout_data) == 0) {
    cat("No bout data to plot\n")
    return(NULL)
  }
  
  # Create plots directory
  plot_dir <- file.path(output_path, "bout_histograms")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # 1. Individual recording histograms
  cat("\nCreating individual recording histograms...\n")
  for (recording_id in names(all_bout_data)) {
    bout_df <- all_bout_data[[recording_id]]
    
    # Cap bout lengths at 12 for plotting
    bout_lengths_capped <- pmin(bout_df$length_epochs, 12)
    
    png(file.path(plot_dir, paste0(recording_id, "_bout_lengths.png")),
        width = 1200, height = 800, res = 120)
    
    par(mfrow = c(2, 2))
    
    # Overall histogram
    hist(bout_lengths_capped, 
         main = paste(recording_id, "- All States"),
         xlab = "Bout Length (epochs)",
         ylab = "Frequency",
         col = "skyblue",
         breaks = seq(0.5, 12.5, by = 1),
         xaxt = 'n')
    axis(1, at = 1:12, labels = c(1:11, "12+"))
    
    # Per-state histograms
    for (state_num in c(1, 2, 3)) {
      state_bouts <- bout_df[bout_df$state == state_num, ]
      if (nrow(state_bouts) > 0) {
        state_lengths_capped <- pmin(state_bouts$length_epochs, 12)
        hist(state_lengths_capped,
             main = paste(recording_id, "-", STATE_LABELS[as.character(state_num)]),
             xlab = "Bout Length (epochs)",
             ylab = "Frequency",
             col = c("coral", "lightblue", "lightgreen")[state_num],
             breaks = seq(0.5, 12.5, by = 1),
             xaxt = 'n')
        axis(1, at = 1:12, labels = c(1:11, "12+"))
      }
    }
    
    dev.off()
  }
  
  # 2. Combined histograms by day type
  cat("Creating combined histograms by day...\n")
  
  # Extract day type from recording IDs (BL, SD, WO)
  for (day_type in c("BL", "SD", "WO")) {
    # Get all recordings for this day type
    day_recordings <- names(all_bout_data)[grepl(paste0("_", day_type, "$"), names(all_bout_data))]
    
    if (length(day_recordings) == 0) next
    
    # Combine bout data from all animals on this day
    combined_bouts <- do.call(rbind, lapply(day_recordings, function(rec_id) {
      all_bout_data[[rec_id]]
    }))
    
    # Cap bout lengths at 12 for plotting
    combined_lengths_capped <- pmin(combined_bouts$length_epochs, 12)
    
    png(file.path(plot_dir, paste0("combined_", day_type, "_bout_lengths.png")),
        width = 1200, height = 800, res = 120)
    
    par(mfrow = c(2, 2))
    
    # Overall histogram for this day
    hist(combined_lengths_capped,
         main = paste(day_type, "- All States (N =", length(day_recordings), "recordings)"),
         xlab = "Bout Length (epochs)",
         ylab = "Frequency",
         col = "skyblue",
         breaks = seq(0.5, 12.5, by = 1),
         xaxt = 'n')
    axis(1, at = 1:12, labels = c(1:11, "12+"))
    abline(v = 5, col = "red", lwd = 2, lty = 2)  # Mark 5-epoch threshold
    legend("topright", legend = "5-epoch threshold", col = "red", lty = 2, lwd = 2)
    
    # Per-state histograms
    for (state_num in c(1, 2, 3)) {
      state_bouts <- combined_bouts[combined_bouts$state == state_num, ]
      if (nrow(state_bouts) > 0) {
        state_lengths_capped <- pmin(state_bouts$length_epochs, 12)
        hist(state_lengths_capped,
             main = paste(day_type, "-", STATE_LABELS[as.character(state_num)], 
                          "(n =", nrow(state_bouts), "bouts)"),
             xlab = "Bout Length (epochs)",
             ylab = "Frequency",
             col = c("coral", "lightblue", "lightgreen")[state_num],
             breaks = seq(0.5, 12.5, by = 1),
             xaxt = 'n')
        axis(1, at = 1:12, labels = c(1:11, "12+"))
        abline(v = 5, col = "red", lwd = 2, lty = 2)
        
        # Add text with percentage usable
        pct_usable <- 100 * sum(state_bouts$length_epochs >= 5) / nrow(state_bouts)
        hist_data <- hist(state_lengths_capped, plot = FALSE, 
                          breaks = seq(0.5, 12.5, by = 1))
        text(x = 9, 
             y = max(hist_data$counts) * 0.9,
             labels = sprintf("%.1f%% ≥ 5 epochs", pct_usable),
             col = "darkred", font = 2)
      }
    }
    
    dev.off()
  }
  
  cat("✓ Histograms saved to:", plot_dir, "\n")
  return(plot_dir)
}


#' Analyze a single recording's bouts
#' 
#' @param states_file Path to states CSV file
#' @param recording_id Identifier for this recording
#' @return List with bout_data and summary_stats
analyze_recording_bouts <- function(states_file, recording_id = NULL) {
  
  # Load states
  states_df <- read.csv(states_file, stringsAsFactors = FALSE)
  
  # Standardize column names
  if ("state" %in% colnames(states_df) && !"State" %in% colnames(states_df)) {
    colnames(states_df)[colnames(states_df) == "state"] <- "State"
  }
  
  if (!"State" %in% colnames(states_df)) {
    stop("No 'State' or 'state' column found in states file")
  }
  
  # Ensure State is numeric
  states_df$state <- as.numeric(states_df$State)
  
  # Identify bouts
  bout_df <- identify_bouts(states_df)
  
  if (nrow(bout_df) == 0) {
    cat("WARNING: No bouts identified in", states_file, "\n")
    return(NULL)
  }
  
  # Calculate summary statistics by state
  summary_stats <- data.frame()
  
  for (state_num in c(1, 2, 3)) {
    state_bouts <- bout_df[bout_df$state == state_num, ]
    
    if (nrow(state_bouts) > 0) {
      state_summary <- data.frame(
        recording_id = recording_id,
        state = state_num,
        state_label = STATE_LABELS[as.character(state_num)],
        total_bouts = nrow(state_bouts),
        mean_length_epochs = mean(state_bouts$length_epochs),
        median_length_epochs = median(state_bouts$length_epochs),
        min_length_epochs = min(state_bouts$length_epochs),
        max_length_epochs = max(state_bouts$length_epochs),
        sd_length_epochs = sd(state_bouts$length_epochs),
        # Count bouts by length category
        bouts_1_3_epochs = sum(state_bouts$length_epochs >= 1 & state_bouts$length_epochs <= 3),
        bouts_4_5_epochs = sum(state_bouts$length_epochs >= 4 & state_bouts$length_epochs <= 5),
        bouts_6_10_epochs = sum(state_bouts$length_epochs >= 6 & state_bouts$length_epochs <= 10),
        bouts_11_plus_epochs = sum(state_bouts$length_epochs >= 11),
        # Percentages
        pct_short_bouts = 100 * sum(state_bouts$length_epochs <= 3) / nrow(state_bouts),
        pct_usable_5epoch = 100 * sum(state_bouts$length_epochs >= 5) / nrow(state_bouts),
        pct_usable_10epoch = 100 * sum(state_bouts$length_epochs >= 10) / nrow(state_bouts)
      )
      
      summary_stats <- rbind(summary_stats, state_summary)
    }
  }
  
  return(list(
    bout_data = bout_df,
    summary_stats = summary_stats
  ))
}


#' Batch process all recordings
#' 
#' @param input_path Directory containing states files
#' @param base_names Character vector of recording identifiers
#' @param output_path Directory for output files
#' @return Combined summary data frame
batch_analyze_bouts <- function(input_path = DEFAULT_INPUT_PATH, 
                                base_names = NULL, 
                                output_path = DEFAULT_OUTPUT_PATH) {
  
  # Default base names if not provided
  if (is.null(base_names)) {
    base_names <- c(
      "mPFCm4_BL", "mPFCm4_SD", "mPFCm4_WO",
      "mPFCm9_BL", "mPFCm9_SD", "mPFCm9_WO",
      "mPFCf5_BL", "mPFCf5_SD", "mPFCf5_WO",
      "mPFCf6_BL", "mPFCf6_SD", "mPFCf6_WO"
    )
  }
  
  cat("============================================================\n")
  cat("BOUT ANALYSIS - BATCH PROCESSING\n")
  cat("============================================================\n")
  cat("Input path:", normalizePath(input_path, mustWork = FALSE), "\n")
  cat("Output path:", normalizePath(output_path, mustWork = FALSE), "\n")
  cat("Processing", length(base_names), "recordings\n\n")
  
  # Create output directory if needed
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Storage for results
  all_summaries <- data.frame()
  all_bout_data <- list()
  processed_count <- 0
  failed_count <- 0
  
  # Process each recording
  for (base in base_names) {
    cat("Processing:", base, "... ")
    
    # Construct file path
    states_file <- file.path(input_path, paste0(base, "_states_df.csv"))
    
    # Check if file exists
    if (!file.exists(states_file)) {
      cat("SKIP (file not found)\n")
      failed_count <- failed_count + 1
      next
    }
    
    # Analyze bouts
    tryCatch({
      result <- analyze_recording_bouts(states_file, recording_id = base)
      
      if (!is.null(result)) {
        all_summaries <- rbind(all_summaries, result$summary_stats)
        all_bout_data[[base]] <- result$bout_data
        processed_count <- processed_count + 1
        
        # Print quick summary
        cat("OK (")
        for (i in seq_len(nrow(result$summary_stats))) {
          row <- result$summary_stats[i, ]
          cat(row$state_label, ":", row$total_bouts, "bouts", sep = "")
          if (i < nrow(result$summary_stats)) cat(", ")
        }
        cat(")\n")
      } else {
        cat("FAILED (no bouts)\n")
        failed_count <- failed_count + 1
      }
      
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      failed_count <- failed_count + 1
    })
  }
  
  # Save combined results
  cat("\n============================================================\n")
  cat("SAVING RESULTS\n")
  cat("============================================================\n")
  
  if (nrow(all_summaries) > 0) {
    summary_file <- file.path(output_path, "bout_analysis_summary.csv")
    write.csv(all_summaries, summary_file, row.names = FALSE)
    cat("✓ Summary saved to:", summary_file, "\n")
    
    # Save individual bout data
    for (recording_id in names(all_bout_data)) {
      bout_file <- file.path(output_path, paste0(recording_id, "_bouts.csv"))
      write.csv(all_bout_data[[recording_id]], bout_file, row.names = FALSE)
    }
    cat("✓ Individual bout data saved for", length(all_bout_data), "recordings\n")
    
    # Create histograms
    create_bout_histograms(all_bout_data, output_path)
    
    # Print overall statistics
    cat("\n============================================================\n")
    cat("OVERALL STATISTICS\n")
    cat("============================================================\n")
    
    for (state_num in c(1, 2, 3)) {
      state_data <- all_summaries[all_summaries$state == state_num, ]
      if (nrow(state_data) > 0) {
        cat("\n", STATE_LABELS[as.character(state_num)], ":\n", sep = "")
        cat("  Total bouts across all recordings:", sum(state_data$total_bouts), "\n")
        cat("  Mean bout length:", round(mean(state_data$mean_length_epochs), 1), "epochs\n")
        cat("  Range:", min(state_data$min_length_epochs), "-", max(state_data$max_length_epochs), "epochs\n")
        cat("  Short bouts (1-3 epochs):", round(mean(state_data$pct_short_bouts), 1), "% on average\n")
        cat("  Usable for 5-epoch window:", round(mean(state_data$pct_usable_5epoch), 1), "% on average\n")
        cat("  Usable for 10-epoch window:", round(mean(state_data$pct_usable_10epoch), 1), "% on average\n")
      }
    }
  }
  
  cat("\n============================================================\n")
  cat("PROCESSING COMPLETE\n")
  cat("============================================================\n")
  cat("Successfully processed:", processed_count, "recordings\n")
  cat("Failed:", failed_count, "recordings\n")
  
  return(all_summaries)
}


# Main execution
if (!interactive()) {
  # Run batch analysis with default parameters
  results <- batch_analyze_bouts(
    input_path = ".",
    output_path = "./bout_analysis_results"
  )
}