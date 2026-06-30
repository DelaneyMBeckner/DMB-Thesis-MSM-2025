# Extract_Transitions_v3.R
# Identifies state transitions (Wake->NREM, NREM->Wake) and extracts 
# preceding epoch windows for both traces and events data
# v3: Adds file_suffix parameter for filtered analysis labeling

library(dplyr)
library(tidyr)

# ==================== CONFIGURATION ====================
SAVE_OUTPUT <- TRUE  # Set to TRUE to save CSV files
INPUT_DIR <- "E:/Data_Processing/R/Data CSVs"  # Where your raw data files are
OUTPUT_DIR <- "E:/Data_Processing/R/Results"   # Where results will be saved
EPOCH_DURATION <- 10  # seconds per epoch

# File labeling (v3)
FILE_SUFFIX <- ""  # Set to "_filtered" when using activity-filtered ROIs

# Input files (just the filenames, will be combined with INPUT_DIR)
RECORDING_ID <- "mPFCf5_BL"
STATES_FILE <- "mPFCf5_BL_states_df.csv"
TRACES_FILE <- "mPFCf5_BL_Traces.csv"
EVENTS_FILE <- "mPFCf5_BL_Events.csv"
# =======================================================


#' Load and prepare state data
load_states <- function(states_file, input_dir = INPUT_DIR) {
  states_path <- file.path(input_dir, states_file)
  cat("Loading states from:", states_path, "\n")
  states_df <- read.csv(states_path, stringsAsFactors = FALSE)
  
  # Ensure we have the right columns
  if (!"state_label" %in% colnames(states_df)) {
    if ("state" %in% colnames(states_df)) {
      state_mapping <- c("1" = "Wake", "2" = "NREM", "3" = "REM")
      states_df$state_label <- state_mapping[as.character(states_df$state)]
    } else {
      stop("Could not find 'state' or 'state_label' column in states file")
    }
  }
  
  # Ensure numeric state column
  if (!"state" %in% colnames(states_df)) {
    state_reverse_mapping <- c("Wake" = 1, "NREM" = 2, "REM" = 3)
    states_df$state <- state_reverse_mapping[states_df$state_label]
  }
  
  # Add epoch index if not present
  if (!"epoch" %in% colnames(states_df)) {
    states_df$epoch <- 0:(nrow(states_df) - 1)
  }
  
  cat("Loaded", nrow(states_df), "epochs\n")
  cat("State distribution:\n")
  print(table(states_df$state_label))
  
  return(states_df)
}


#' Identify all transitions in the recording
identify_transitions <- function(states_df) {
  cat("\nIdentifying transitions...\n")
  
  transitions <- data.frame()
  
  for (i in 2:nrow(states_df)) {
    prev_state <- states_df$state[i-1]
    curr_state <- states_df$state[i]
    
    # Check for Wakeâ†’NREM (1â†’2)
    if (prev_state == 1 && curr_state == 2) {
      transitions <- rbind(transitions, data.frame(
        transition_epoch = i - 1,  # Last epoch before transition
        transition_type = "Wake2NREM",
        from_state = "Wake",
        to_state = "NREM",
        transition_time = states_df$Time.from.Start[i]
      ))
    }
    
    # Check for NREMâ†’Wake (2â†’1)
    if (prev_state == 2 && curr_state == 1) {
      transitions <- rbind(transitions, data.frame(
        transition_epoch = i - 1,
        transition_type = "NREM2Wake",
        from_state = "NREM",
        to_state = "Wake",
        transition_time = states_df$Time.from.Start[i]
      ))
    }
  }
  
  cat("Found", nrow(transitions), "transitions:\n")
  cat("  Wakeâ†’NREM:", sum(transitions$transition_type == "Wake2NREM"), "\n")
  cat("  NREMâ†’Wake:", sum(transitions$transition_type == "NREM2Wake"), "\n")
  
  return(transitions)
}


#' Extract epoch windows before each transition
extract_transition_windows <- function(states_df, transitions, window_sizes = c(3, 9)) {
  cat("\nExtracting transition windows...\n")
  
  windows_metadata <- data.frame()
  transition_id <- 1
  
  for (idx in 1:nrow(transitions)) {
    trans_epoch <- transitions$transition_epoch[idx]
    trans_type <- transitions$transition_type[idx]
    trans_time <- transitions$transition_time[idx]
    
    for (window_size in window_sizes) {
      # Window goes from -window_size to +1 (skipping 0, which doesn't exist)
      # trans_epoch is the last epoch of the OLD state
      # trans_epoch + 1 is the first epoch of the NEW state (epoch +1)
      
      # Calculate start epoch (window_size epochs before transition)
      start_epoch <- trans_epoch - window_size + 1
      
      # Skip if too close to recording start
      if (start_epoch < 0) {
        cat("  Skipping transition at epoch", trans_epoch, 
            "- insufficient preceding epochs for", window_size, "epoch window\n")
        next
      }
      
      # End epoch is now trans_epoch + 1 (the first post-transition epoch)
      end_epoch <- trans_epoch + 1
      
      # Skip if too close to recording end
      if (end_epoch >= nrow(states_df)) {
        cat("  Skipping transition at epoch", trans_epoch, 
            "- insufficient following epochs\n")
        next
      }
      
      # Get time range
      window_start_time <- states_df$Time.from.Start[start_epoch + 1]  # +1 for R indexing
      window_end_time <- states_df$Time.from.Start[end_epoch + 1] + EPOCH_DURATION
      
      # Total epochs = window_size (pre) + 1 (post)
      total_epochs <- window_size + 1
      
      # Store metadata
      windows_metadata <- rbind(windows_metadata, data.frame(
        transition_id = transition_id,
        recording_id = RECORDING_ID,
        transition_type = trans_type,
        transition_epoch = trans_epoch,
        transition_time = trans_time,
        window_size = window_size,
        start_epoch = start_epoch,
        end_epoch = end_epoch,
        window_start_time = window_start_time,
        window_end_time = window_end_time,
        num_epochs = total_epochs,
        stringsAsFactors = FALSE
      ))
      
      transition_id <- transition_id + 1
    }
  }
  
  cat("Created", nrow(windows_metadata), "transition windows\n")
  
  return(windows_metadata)
}


#' Load and extract traces for transition windows
extract_traces <- function(traces_file, windows_metadata, input_dir = INPUT_DIR) {
  traces_path <- file.path(input_dir, traces_file)
  cat("\nLoading traces from:", traces_path, "\n")
  traces_df <- read.csv(traces_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Fix the first column name if it's empty or has special characters
  if (colnames(traces_df)[1] == "" || colnames(traces_df)[1] == "Time(s)/Cell Status") {
    colnames(traces_df)[1] <- "Time"
  }
  
  # Remove header row if first row contains text instead of numbers
  first_val <- suppressWarnings(as.numeric(traces_df[1, 1]))
  if (is.na(first_val)) {
    cat("  Removing header row from traces\n")
    traces_df <- traces_df[-1, ]
  }
  
  # Convert all to numeric
  traces_df[] <- lapply(traces_df, function(x) as.numeric(as.character(x)))
  
  # Verify Time column exists and has data
  if (!"Time" %in% colnames(traces_df)) {
    stop("Time column not found after processing. Columns: ", paste(colnames(traces_df), collapse = ", "))
  }
  if (all(is.na(traces_df$Time))) {
    stop("Time column contains only NA values after conversion")
  }
  
  cat("Loaded traces:", nrow(traces_df), "timepoints,", ncol(traces_df) - 1, "ROIs\n")
  
  # Extract traces for each window
  traces_list <- list()
  
  for (i in 1:nrow(windows_metadata)) {
    trans_id <- windows_metadata$transition_id[i]
    start_time <- windows_metadata$window_start_time[i]
    end_time <- windows_metadata$window_end_time[i]
    
    # Extract traces within time window
    window_traces <- traces_df[traces_df$Time >= start_time & traces_df$Time < end_time, ]
    
    if (nrow(window_traces) > 0) {
      traces_list[[as.character(trans_id)]] <- window_traces
      cat("  Transition", trans_id, ":", nrow(window_traces), "timepoints\n")
    } else {
      cat("  WARNING: No traces found for transition", trans_id, "\n")
    }
  }
  
  return(traces_list)
}


#' Load and extract events for transition windows
extract_events <- function(events_file, windows_metadata, input_dir = INPUT_DIR) {
  events_path <- file.path(input_dir, events_file)
  cat("\nLoading events from:", events_path, "\n")
  events_df <- read.csv(events_path, stringsAsFactors = FALSE)
  
  # Handle column name variations (R converts special chars)
  time_col <- NULL
  if ("Time..s." %in% colnames(events_df)) {
    time_col <- "Time..s."
  } else if ("Time (s)" %in% colnames(events_df)) {
    time_col <- "Time (s)"
  } else if ("Time" %in% colnames(events_df)) {
    time_col <- "Time"
  } else {
    stop("Could not find time column in events file. Available columns: ", paste(colnames(events_df), collapse = ", "))
  }
  
  cell_col <- NULL
  if ("Cell.Name" %in% colnames(events_df)) {
    cell_col <- "Cell.Name"
  } else if (" Cell Name" %in% colnames(events_df)) {
    cell_col <- " Cell Name"
  } else if ("Cell Name" %in% colnames(events_df)) {
    cell_col <- "Cell Name"
  } else {
    stop("Could not find Cell Name column in events file. Available columns: ", paste(colnames(events_df), collapse = ", "))
  }
  
  value_col <- NULL
  if ("Value" %in% colnames(events_df)) {
    value_col <- "Value"
  } else {
    stop("Could not find Value column in events file. Available columns: ", paste(colnames(events_df), collapse = ", "))
  }
  
  cat("Loaded", nrow(events_df), "events\n")
  
  # Create aggregated events dataframe
  all_events <- data.frame()
  
  for (i in 1:nrow(windows_metadata)) {
    trans_id <- windows_metadata$transition_id[i]
    trans_type <- windows_metadata$transition_type[i]
    window_size <- windows_metadata$window_size[i]
    start_time <- windows_metadata$window_start_time[i]
    end_time <- windows_metadata$window_end_time[i]
    start_epoch <- windows_metadata$start_epoch[i]
    
    # Extract events within time window
    window_events <- events_df[events_df[[time_col]] >= start_time & 
                                events_df[[time_col]] < end_time, ]
    
    if (nrow(window_events) > 0) {
      # Calculate which epoch within the window each event belongs to
      window_events$epoch_abs <- floor((window_events[[time_col]] - start_time) / EPOCH_DURATION) + start_epoch
      
      # Transition epoch is the last epoch of the old state (end_epoch - 1)
      # Pre-transition epochs get negative numbers: -window_size to -1
      # Post-transition epoch (end_epoch) gets +1
      # No epoch 0 exists
      trans_epoch <- windows_metadata$end_epoch[i] - 1  # Last epoch of old state
      
      # For pre-transition epochs (epoch_abs <= trans_epoch): 
      #   epoch_in_window = epoch_abs - trans_epoch - 1 (gives -window_size to -1)
      # For post-transition epochs (epoch_abs > trans_epoch):
      #   epoch_in_window = epoch_abs - trans_epoch (gives +1)
      window_events$epoch_in_window <- ifelse(
        window_events$epoch_abs <= trans_epoch,
        window_events$epoch_abs - trans_epoch - 1,  # Pre-transition: -window_size to -1
        window_events$epoch_abs - trans_epoch       # Post-transition: +1
      )
      
      # Add metadata columns
      window_events$transition_id <- trans_id
      window_events$transition_type <- trans_type
      window_events$window_size <- window_size
      
      # Standardize column names
      window_events$time_s <- window_events[[time_col]]
      window_events$Cell_Name <- trimws(window_events[[cell_col]])
      window_events$Value <- window_events[[value_col]]
      
      # Select only the columns we want
      window_events <- window_events[, c("transition_id", "transition_type", "window_size", 
                                          "epoch_in_window", "epoch_abs", "time_s", 
                                          "Cell_Name", "Value")]
      
      all_events <- rbind(all_events, window_events)
      
      cat("  Transition", trans_id, ":", nrow(window_events), "events\n")
    } else {
      cat("  Transition", trans_id, ": 0 events\n")
    }
  }
  
  cat("Total events extracted:", nrow(all_events), "\n")
  
  return(all_events)
}


#' Save outputs to CSV files (metadata and events only, traces never saved)
save_outputs <- function(windows_metadata, events_aggregated, recording_id, output_dir, 
                        file_suffix = "") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving outputs to:", output_dir, "\n")
  
  # Group metadata by transition type and window size
  condition_groups <- unique(windows_metadata[, c("transition_type", "window_size")])
  
  for (i in 1:nrow(condition_groups)) {
    trans_type <- condition_groups$transition_type[i]
    window_size <- condition_groups$window_size[i]
    
    # Filter metadata for this condition
    condition_meta <- windows_metadata[windows_metadata$transition_type == trans_type & 
                                        windows_metadata$window_size == window_size, ]
    
    # Save metadata (flat - directly in output_dir)
    meta_file <- file.path(output_dir, paste0(recording_id, "_", trans_type, "_", window_size, "ep", file_suffix, "_metadata.csv"))
    write.csv(condition_meta, meta_file, row.names = FALSE)
    cat("  Saved:", meta_file, "\n")
    
    # Save aggregated events for this condition (flat)
    condition_events <- events_aggregated[events_aggregated$transition_type == trans_type & 
                                           events_aggregated$window_size == window_size, ]
    if (nrow(condition_events) > 0) {
      events_file <- file.path(output_dir, paste0(recording_id, "_", trans_type, "_", window_size, "ep", file_suffix, "_events.csv"))
      write.csv(condition_events, events_file, row.names = FALSE)
      cat("  Saved:", events_file, "\n")
    }
  }
  
  cat("\nAll outputs saved to:", output_dir, "\n")
  cat("Note: Traces not saved - use metadata to extract from full traces file\n")
}


#' Main extraction function
extract_transitions_main <- function(states_file = STATES_FILE,
                                     traces_file = TRACES_FILE,
                                     events_file = EVENTS_FILE,
                                     recording_id = RECORDING_ID,
                                     input_dir = INPUT_DIR,
                                     window_sizes = c(3, 9),
                                     save_output = SAVE_OUTPUT,
                                     output_dir = OUTPUT_DIR,
                                     file_suffix = FILE_SUFFIX) {
  
  cat("==================================================\n")
  cat("Extract Transitions Analysis (v3)\n")
  cat("==================================================\n")
  cat("Recording ID:", recording_id, "\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Window sizes:", paste(window_sizes, collapse = ", "), "epochs\n")
  cat("Save output:", save_output, "\n")
  if (file_suffix != "") cat("File suffix:", file_suffix, "\n")
  cat("==================================================\n\n")
  
  # Step 1: Load states and identify transitions
  states_df <- load_states(states_file, input_dir)
  transitions <- identify_transitions(states_df)
  
  if (nrow(transitions) == 0) {
    stop("No transitions found in the recording!")
  }
  
  # Step 2: Extract windows metadata
  windows_metadata <- extract_transition_windows(states_df, transitions, window_sizes)
  
  if (nrow(windows_metadata) == 0) {
    stop("No valid transition windows could be extracted!")
  }
  
  # Step 3: Extract traces for each window
  traces_list <- extract_traces(traces_file, windows_metadata, input_dir)
  
  # Step 4: Extract and aggregate events
  events_aggregated <- extract_events(events_file, windows_metadata, input_dir)
  
  # Step 5: Save if requested
  if (save_output) {
    save_outputs(windows_metadata, events_aggregated, recording_id, output_dir, file_suffix)
  } else {
    cat("\nSAVE_OUTPUT = FALSE: Results returned but not saved to disk\n")
  }
  
  # Return all results
  results <- list(
    metadata = windows_metadata,
    traces = traces_list,
    events = events_aggregated,
    transitions_summary = transitions
  )
  
  cat("\n==================================================\n")
  cat("Extraction complete!\n")
  cat("  Transitions found:", nrow(transitions), "\n")
  cat("  Windows extracted:", nrow(windows_metadata), "\n")
  cat("  Events extracted:", nrow(events_aggregated), "\n")
  cat("==================================================\n")
  
  return(results)
}


# =========================== MAIN ===========================
# Run the extraction
if (sys.nframe() == 0) {
  results <- extract_transitions_main()
}
