# Summarize_Transition_Counts.R
# Runs Extract_Transitions on all animals/conditions, then counts windows

library(dplyr)
library(tidyr)

# ==================== CONFIGURATION ====================
# Run extraction before counting?
RUN_EXTRACTION <- TRUE  # Set to FALSE to only count existing files

# Paths
INPUT_DIR <- "E:/Data_Processing/R/Data CSVs"
RESULTS_DIR <- "E:/Data_Processing/R/Results"
EXTRACTION_SCRIPT <- "E:/Data_Processing/R/Extract_Transitions_v3.5.R"  # Path to v3.5 script

# Analysis parameters
FILE_SUFFIX <- ""  # Set to "_filtered" if using filtered analysis
WINDOW_SIZES <- c(3, 9)
EPOCH_DURATION <- 10

# Animals and conditions to process
ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9")
CONDITIONS <- c("BL", "SD", "WO")
TRANSITION_TYPES <- c("Wake2NREM", "NREM2Wake")
# =======================================================


#' Run extraction for all animals and conditions
run_batch_extraction <- function(animals = ANIMALS,
                                conditions = CONDITIONS,
                                input_dir = INPUT_DIR,
                                results_dir = RESULTS_DIR,
                                extraction_script = EXTRACTION_SCRIPT,
                                window_sizes = WINDOW_SIZES,
                                file_suffix = FILE_SUFFIX,
                                epoch_duration = EPOCH_DURATION) {
  
  cat("\n==================================================\n")
  cat("BATCH EXTRACTION: Extract_Transitions_v3_5\n")
  cat("==================================================\n")
  cat("Extraction script:", extraction_script, "\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", results_dir, "\n")
  if (file_suffix != "") cat("File suffix:", file_suffix, "\n")
  cat("==================================================\n\n")
  
  # Source the extraction script
  source(extraction_script)
  
  extraction_log <- data.frame()
  
  for (animal in animals) {
    for (condition in conditions) {
      recording_id <- paste0(animal, "_", condition)
      
      cat("\n========================================\n")
      cat("Processing:", recording_id, "\n")
      cat("========================================\n")
      
      # Build filenames
      states_file <- paste0(recording_id, "_states_df.csv")
      traces_file <- paste0(recording_id, "_Traces.csv")
      events_file <- paste0(recording_id, "_Events.csv")
      
      # Check if input files exist
      states_path <- file.path(input_dir, states_file)
      traces_path <- file.path(input_dir, traces_file)
      events_path <- file.path(input_dir, events_file)
      
      if (!file.exists(states_path)) {
        cat("  SKIPPED: States file not found:", states_file, "\n")
        extraction_log <- rbind(extraction_log, data.frame(
          Animal = animal, Condition = condition, Status = "Missing_States", stringsAsFactors = FALSE
        ))
        next
      }
      if (!file.exists(traces_path)) {
        cat("  SKIPPED: Traces file not found:", traces_file, "\n")
        extraction_log <- rbind(extraction_log, data.frame(
          Animal = animal, Condition = condition, Status = "Missing_Traces", stringsAsFactors = FALSE
        ))
        next
      }
      if (!file.exists(events_path)) {
        cat("  SKIPPED: Events file not found:", events_file, "\n")
        extraction_log <- rbind(extraction_log, data.frame(
          Animal = animal, Condition = condition, Status = "Missing_Events", stringsAsFactors = FALSE
        ))
        next
      }
      
      # Run extraction
      tryCatch({
        results <- extract_transitions_main(
          states_file = states_file,
          traces_file = traces_file,
          events_file = events_file,
          recording_id = recording_id,
          input_dir = input_dir,
          window_sizes = window_sizes,
          save_output = TRUE,
          output_dir = results_dir,
          file_suffix = file_suffix
        )
        
        extraction_log <- rbind(extraction_log, data.frame(
          Animal = animal, 
          Condition = condition, 
          Status = "Success",
          N_Transitions = nrow(results$transitions_summary),
          N_Windows = nrow(results$metadata),
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        cat("  ERROR:", e$message, "\n")
        extraction_log <<- rbind(extraction_log, data.frame(
          Animal = animal, Condition = condition, Status = paste0("Error: ", e$message), 
          stringsAsFactors = FALSE
        ))
      })
    }
  }
  
  cat("\n==================================================\n")
  cat("BATCH EXTRACTION COMPLETE\n")
  cat("==================================================\n\n")
  
  print(extraction_log)
  return(extraction_log)
}


#' Count transitions from metadata files
count_transitions <- function(results_dir = RESULTS_DIR, 
                             animals = ANIMALS,
                             conditions = CONDITIONS,
                             transition_types = TRANSITION_TYPES,
                             window_sizes = WINDOW_SIZES,
                             file_suffix = FILE_SUFFIX) {
  
  cat("==================================================\n")
  cat("Transition Window Summary\n")
  cat("==================================================\n")
  cat("Results directory:", results_dir, "\n")
  if (file_suffix != "") cat("File suffix:", file_suffix, "\n")
  cat("==================================================\n\n")
  
  summary_data <- data.frame()
  
  for (animal in animals) {
    for (condition in conditions) {
      for (trans_type in transition_types) {
        for (window_size in window_sizes) {
          
          # Build filename
          recording_id <- paste0(animal, "_", condition)
          filename <- paste0(recording_id, "_", trans_type, "_", window_size, "ep", file_suffix, "_metadata.csv")
          filepath <- file.path(results_dir, filename)
          
          # Check if file exists
          if (file.exists(filepath)) {
            # Read and count rows
            metadata <- read.csv(filepath, stringsAsFactors = FALSE)
            n_windows <- nrow(metadata)
            
            summary_data <- rbind(summary_data, data.frame(
              Animal = animal,
              Condition = condition,
              Transition_Type = trans_type,
              Window_Size = window_size,
              N_Windows = n_windows,
              stringsAsFactors = FALSE
            ))
            
            cat(sprintf("  %s: %d windows\n", filename, n_windows))
          } else {
            cat(sprintf("  %s: FILE NOT FOUND\n", filename))
          }
        }
      }
    }
  }
  
  cat("\n==================================================\n")
  cat("Summary complete!\n")
  cat("==================================================\n\n")
  
  return(summary_data)
}


#' Print formatted summary tables
print_summary <- function(summary_data) {
  
  # Overall summary by animal and condition
  cat("\n========== BY ANIMAL AND CONDITION ==========\n")
  animal_condition_summary <- summary_data %>%
    group_by(Animal, Condition) %>%
    summarise(
      Total_Windows = sum(N_Windows),
      Wake2NREM = sum(N_Windows[Transition_Type == "Wake2NREM"]),
      NREM2Wake = sum(N_Windows[Transition_Type == "NREM2Wake"]),
      .groups = "drop"
    ) %>%
    arrange(Animal, Condition)
  
  print(animal_condition_summary, n = Inf)
  
  # By window size
  cat("\n========== BY WINDOW SIZE ==========\n")
  window_summary <- summary_data %>%
    group_by(Animal, Condition, Window_Size) %>%
    summarise(
      Total_Windows = sum(N_Windows),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Window_Size,
      values_from = Total_Windows,
      names_prefix = "Window_"
    ) %>%
    arrange(Animal, Condition)
  
  print(window_summary, n = Inf)
  
  # By transition type
  cat("\n========== BY TRANSITION TYPE ==========\n")
  type_summary <- summary_data %>%
    group_by(Animal, Condition, Transition_Type) %>%
    summarise(
      Total_Windows = sum(N_Windows),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Transition_Type,
      values_from = Total_Windows
    ) %>%
    arrange(Animal, Condition)
  
  print(type_summary, n = Inf)
  
  # Grand totals
  cat("\n========== GRAND TOTALS ==========\n")
  grand_totals <- summary_data %>%
    group_by(Animal) %>%
    summarise(
      Total_All_Conditions = sum(N_Windows),
      .groups = "drop"
    )
  
  print(grand_totals, n = Inf)
  
  cat("\nOverall total windows:", sum(summary_data$N_Windows), "\n")
}


# =========================== MAIN ===========================
# Auto-run when sourced or called with Rscript
# To prevent auto-run, set: RUN_EXTRACTION <- FALSE before sourcing

# Step 1: Run batch extraction if requested
if (RUN_EXTRACTION) {
  extraction_log <- run_batch_extraction()
  
  # Save extraction log
  log_file <- file.path(RESULTS_DIR, "Extraction_Log.csv")
  write.csv(extraction_log, log_file, row.names = FALSE)
  cat("\nExtraction log saved to:", log_file, "\n")
}

# Step 2: Count transitions from metadata files
summary_data <- count_transitions()
print_summary(summary_data)

# Step 3: Save summary
output_file <- file.path(RESULTS_DIR, "Transition_Window_Summary.csv")
write.csv(summary_data, output_file, row.names = FALSE)
cat("\nSummary saved to:", output_file, "\n")
