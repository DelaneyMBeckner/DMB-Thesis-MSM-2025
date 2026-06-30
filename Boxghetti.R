# Fixed Neural Activity Analysis Script - Corrected Sex Extraction
# ================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)

# CORRECTED: Function to extract animal and day info from filename
parse_filename <- function(filename) {
  # Remove file extension
  base_name <- tools::file_path_sans_ext(basename(filename))
  
  # Split by underscore and extract components
  parts <- str_split(base_name, "_")[[1]]
  
  if (length(parts) >= 2) {
    # First part should be animal ID (e.g., mPFCf5, mPFCm4)
    animal <- parts[1]
    
    # Second part should be day code (BL, SD, WO)
    day_code <- parts[2]
    
    # FIXED: Proper sex extraction - look for 'f' or 'm' followed by digits
    # This handles cases like mPFCf5, mPFCf6, mPFCm4, mPFCm9
    if (str_detect(animal, "f\\d+$")) {
      # Ends with 'f' followed by digits = female
      sex <- "f"
      animal_number <- str_extract(animal, "\\d+$")
    } else if (str_detect(animal, "m\\d+$")) {
      # Ends with 'm' followed by digits = male  
      sex <- "m"
      animal_number <- str_extract(animal, "\\d+$")
    } else {
      # Fallback - try to extract any digits
      sex <- "Unknown"
      animal_number <- str_extract(animal, "\\d+")
    }
    
    # Extract experiment group (everything before the sex+number)
    experiment_group <- str_remove(animal, "[fm]\\d+$")
    
    # Convert day codes to ordered factors for proper plotting sequence
    day_order <- c("BL" = 1, "SD" = 2, "WO" = 3)
    day_numeric <- day_order[day_code]
    
    if (is.na(day_numeric)) {
      warning(paste("Unknown day code:", day_code, "in file:", filename))
      day_numeric <- 1
      day_code <- "Unknown"
    }
    
    return(list(
      animal = animal, 
      day_code = day_code,
      day_numeric = day_numeric,
      sex = sex,
      animal_number = animal_number,
      experiment_group = experiment_group
    ))
  } else {
    warning(paste("Could not parse filename:", filename))
    return(list(
      animal = "Unknown", 
      day_code = "Unknown",
      day_numeric = 1,
      sex = "Unknown",
      animal_number = "Unknown",
      experiment_group = "Unknown"
    ))
  }
}

# Function to load and process multiple CSV files with animal/day info
load_csv_files_with_metadata <- function(file_pattern = "*.csv") {
  files <- list.files(pattern = file_pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No CSV files found matching the pattern!")
  }
  
  cat("Found", length(files), "CSV files\n")
  
  # Read and combine all files with metadata
  all_data <- map_dfr(files, function(file) {
    cat("Processing:", basename(file), "\n")
    
    df <- read_csv(file, show_col_types = FALSE)
    names(df) <- trimws(names(df))
    
    # Extract animal and day info from filename
    metadata <- parse_filename(file)
    
    df$Animal <- metadata$animal
    df$Day_Code <- metadata$day_code
    df$Day <- metadata$day_numeric
    df$Sex <- metadata$sex
    df$Animal_Number <- metadata$animal_number
    df$Experiment_Group <- metadata$experiment_group
    df$File <- tools::file_path_sans_ext(basename(file))
    
    return(df)
  })
  
  # Show what was parsed
  parsed_info <- all_data %>% 
    select(File, Animal, Day_Code, Day, Sex, Animal_Number, Experiment_Group) %>% 
    distinct() %>%
    arrange(Animal, Day)
  
  cat("\nParsed file information:\n")
  print(parsed_info)
  
  # Check for potential issues
  cat("\nData validation:\n")
  cat("Unique animals found:", paste(sort(unique(all_data$Animal)), collapse = ", "), "\n")
  
  # Check sex distribution
  sex_counts <- all_data %>% 
    select(Animal, Sex) %>% 
    distinct() %>% 
    count(Sex)
  cat("Sex distribution:\n")
  print(sex_counts)
  
  return(all_data)
}

# Function to prepare data for plotting
prepare_data <- function(df) {
  cell_cols <- grep("^C\\d+$", names(df), value = TRUE)
  
  # Get unique animals and sort them for consistent ordering
  unique_animals <- sort(unique(df$Animal))
  
  df_long <- df %>%
    pivot_longer(
      cols = all_of(cell_cols),
      names_to = "Cell",
      values_to = "Event_Rate"
    ) %>%
    rename(State = `Cell Name`) %>%
    mutate(
      Cell_Num = as.numeric(str_remove(Cell, "C")),
      Animal = factor(Animal, levels = unique_animals),
      Day_Code = factor(Day_Code, levels = c("BL", "SD", "WO")),
      Day = factor(Day, levels = c(1, 2, 3), labels = c("BL", "SD", "WO")),
      Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male"))
    ) %>%
    # Remove any rows with missing Event_Rate
    filter(!is.na(Event_Rate))
  
  return(df_long)
}

# Function to generate grayscale colors for cells
generate_cell_colors <- function(n_cells) {
  # Create grayscale palette with good contrast
  if (n_cells == 1) {
    return("#404040")
  }
  # Generate grayscale values from dark to light gray
  gray_values <- seq(0.2, 0.8, length.out = n_cells)
  colors <- paste0("gray", round(gray_values * 100))
  return(colors)
}

# Function to create individual animal plots
create_individual_plots <- function(df_long) {
  unique_animals <- unique(df_long$Animal)
  
  for (animal in unique_animals) {
    animal_data <- df_long %>% filter(Animal == animal)
    
    # Get unique cells and create color palette
    unique_cells <- sort(unique(animal_data$Cell))
    n_cells <- length(unique_cells)
    cell_colors <- generate_cell_colors(n_cells)
    names(cell_colors) <- unique_cells
    
    # Create position mapping for states
    state_positions <- c("Wake" = 1, "NREM" = 2, "REM" = 3)
    animal_data <- animal_data %>%
      mutate(State_Position = state_positions[State])
    
    # State colors for box plots
    state_colors <- c("Wake" = "#B8860B", "NREM" = "#4169E1", "REM" = "#DC143C")
    animal_data <- animal_data %>%
      mutate(State_Name = case_when(
        State_Position == 1 ~ "Wake",
        State_Position == 2 ~ "NREM", 
        State_Position == 3 ~ "REM"
      ))
    
    # Create the plot for this animal
    p <- ggplot(animal_data, aes(x = State_Position, y = Event_Rate)) +
      
      # Box plots with thicker whiskers
      geom_boxplot(aes(group = State_Position, fill = State_Name),
                   alpha = 0.3, width = 0.6, 
                   outlier.alpha = 0.6, outlier.size = 1.2,
                   linewidth = 0.5, show.legend = FALSE) +
      
      # Line plots connecting states for each cell
      geom_line(aes(color = Cell, group = Cell),
                stat = "summary", fun = mean,
                linewidth = 1.0, alpha = 0.8, show.legend = FALSE) +
      
      # Points for mean values
      geom_point(aes(color = Cell, group = Cell),
                 stat = "summary", fun = mean,
                 size = 2, alpha = 0.9, show.legend = FALSE) +
      
      # Facet by Day only (since we're doing one animal at a time)
      facet_wrap(~ Day, nrow = 1, 
                 labeller = labeller(
                   Day = function(x) {
                     case_when(
                       x == "BL" ~ "Baseline",
                       x == "SD" ~ "Sleep Deprived", 
                       x == "WO" ~ "Recovery",
                       TRUE ~ as.character(x)
                     )
                   }
                 )) +
      
      # Apply colors
      scale_color_manual(values = cell_colors, name = "Cell") +
      scale_fill_manual(values = state_colors, name = "State") +
      
      # X-axis with state labels
      scale_x_continuous(
        breaks = c(1, 2, 3),
        labels = c("Wake", "NREM", "REM"),
        limits = c(0.5, 3.5)
      ) +
      
      # Y-axis
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      
      # Customize appearance
      labs(
        title = paste("Neural Activity:", animal),
        subtitle = "mPFC Recordings | Lines: Cell Activity Across States | Boxes: Event Rate Distribution",
        x = "Sleep State",
        y = "Normalized Event Rate"
      ) +
      
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none",  # Remove all legends
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3)
      )
    
    # Save individual plot
    filename <- paste0("neural_activity_", animal, ".png")
    ggsave(filename, plot = p, 
           width = 12, height = 6, 
           dpi = 300, bg = "white")
    
    cat("Individual plot saved:", filename, "\n")
  }
}

# Function to create multi-panel plot
create_multipanel_plot <- function(df_long, output_file = "neural_activity_multipanel.png") {
  
  # Calculate y-axis limits per animal (same for all days of each animal)
  y_limits <- df_long %>%
    group_by(Animal) %>%
    summarise(
      max_rate = max(Event_Rate, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      # Round up to nearest nice number
      y_max = ceiling(max_rate * 10) / 10
    )
  
  # Join y-limits back to main data
  df_plot <- df_long %>%
    left_join(y_limits, by = "Animal")
  
  # Get unique cells and create color palette
  unique_cells <- sort(unique(df_plot$Cell))
  n_cells <- length(unique_cells)
  cell_colors <- generate_cell_colors(n_cells)
  names(cell_colors) <- unique_cells
  
  # Create position mapping for states (x-axis positions)
  state_positions <- c("Wake" = 1, "NREM" = 2, "REM" = 3)
  
  # Add position information to data
  df_plot <- df_plot %>%
    mutate(State_Position = state_positions[State])
  
  # Create state colors for box plots
  state_colors <- c("Wake" = "#B8860B", "NREM" = "#4169E1", "REM" = "#DC143C")
  
  # Add state names to data for coloring
  df_plot <- df_plot %>%
    mutate(State_Name = case_when(
      State_Position == 1 ~ "Wake",
      State_Position == 2 ~ "NREM", 
      State_Position == 3 ~ "REM"
    ))
  
  # Create the faceted plot
  p <- ggplot(df_plot, aes(x = State_Position, y = Event_Rate)) +
    
    # Box plots for distribution by state (across all cells) - thicker whiskers
    geom_boxplot(aes(group = State_Position, fill = State_Name),
                 alpha = 0.3, width = 0.6, 
                 outlier.alpha = 0.6, outlier.size = 1.2,
                 linewidth = 0.5, show.legend = FALSE) +
    
    # Line plots connecting states for each cell
    geom_line(aes(color = Cell, group = Cell),
              stat = "summary", fun = mean,
              linewidth = 0.8, alpha = 0.8, show.legend = FALSE) +
    
    # Points for mean values
    geom_point(aes(color = Cell, group = Cell),
               stat = "summary", fun = mean,
               size = 1.5, alpha = 0.9, show.legend = FALSE) +
    
    # Facet by Animal (rows) and Day (columns) - UPDATED LABELING
    facet_grid(Animal ~ Day, 
               scales = "free_y",
               labeller = labeller(
                 Animal = function(x) {
                   # NEW: Create "Male 5", "Female 6" format
                   # Handle vectors by using vectorized operations
                   result <- case_when(
                     str_detect(x, "f\\d+$") ~ paste0("Female ", str_extract(x, "\\d+$")),
                     str_detect(x, "m\\d+$") ~ paste0("Male ", str_extract(x, "\\d+$")),
                     TRUE ~ as.character(x)
                   )
                   return(result)
                 },
                 Day = function(x) {
                   # Use descriptive day labels for column headers
                   case_when(
                     x == "BL" ~ "Baseline",
                     x == "SD" ~ "Sleep Deprived", 
                     x == "WO" ~ "Recovery",
                     TRUE ~ as.character(x)
                   )
                 }
               ),
               # Move row labels to left side
               switch = "y") +
    
    # Apply colors
    scale_color_manual(values = cell_colors, name = "Cell") +
    scale_fill_manual(values = state_colors, name = "State") +
    
    # X-axis with state labels
    scale_x_continuous(
      breaks = c(1, 2, 3),
      labels = c("Wake", "NREM", "REM"),
      limits = c(0.5, 3.5)
    ) +
    
    # Y-axis
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    
    # Customize appearance
    labs(
      title = "Neural Activity: Sleep Deprivation Experiment",
      subtitle = "mPFC Recordings | Lines: Cell Activity Across States | Boxes: Event Rate Distribution",
      x = "Sleep State",
      y = "Normalized Event Rate"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      strip.text.y.left = element_text(angle = 0),  # Keep row labels horizontal
      legend.position = "none",  # Remove all legends
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3)
    )
  
  # Determine figure size based on number of animals and days
  n_animals <- length(unique(df_long$Animal))
  n_days <- length(unique(df_long$Day))
  
  fig_width <- min(20, max(12, n_days * 3))
  fig_height <- min(16, max(8, n_animals * 2.5))
  
  # Save the plot
  ggsave(output_file, plot = p, 
         width = fig_width, height = fig_height, 
         dpi = 300, bg = "white")
  
  cat("Plot saved as:", output_file, "\n")
  print(p)
  return(p)
}

# Test the corrected sex extraction function
test_sex_extraction <- function() {
  test_files <- c("mPFCf5_BL_normalized_event_counts.csv",
                  "mPFCf6_SD_normalized_event_counts.csv", 
                  "mPFCm4_WO_normalized_event_counts.csv",
                  "mPFCm9_BL_normalized_event_counts.csv")
  
  cat("Testing corrected sex extraction:\n")
  cat("=================================\n")
  for (file in test_files) {
    result <- parse_filename(file)
    cat(sprintf("File: %s\n", file))
    cat(sprintf("  -> Animal: %s, Sex: %s, Number: %s, Group: %s\n\n", 
                result$animal, result$sex, result$animal_number, result$experiment_group))
  }
}

# Enhanced diagnostic function
check_data_issues <- function(df_long) {
  cat("=== ENHANCED DATA DIAGNOSTIC ===\n")
  cat("Total rows:", nrow(df_long), "\n")
  cat("Animals found:", paste(unique(df_long$Animal), collapse = ", "), "\n")
  cat("Days found:", paste(unique(df_long$Day_Code), collapse = ", "), "\n")
  cat("States found:", paste(unique(df_long$State), collapse = ", "), "\n")
  cat("Cells found:", length(unique(df_long$Cell)), "cells total\n")
  
  # Check sex distribution with more detail
  sex_summary <- df_long %>%
    select(Animal, Sex, Animal_Number) %>%
    distinct() %>%
    arrange(Animal) %>%
    mutate(Sex_Label = case_when(
      Sex == "Female" ~ "♀ Female",
      Sex == "Male" ~ "♂ Male",
      TRUE ~ as.character(Sex)
    ))
  
  cat("\nAnimal details:\n")
  print(sex_summary)
  
  # Check data completeness by animal
  completeness <- df_long %>%
    group_by(Animal, Day_Code) %>%
    summarise(
      n_states = n_distinct(State),
      n_cells = n_distinct(Cell),
      n_observations = n(),
      .groups = 'drop'
    ) %>%
    arrange(Animal, Day_Code)
  
  cat("\nData completeness by animal and day:\n")
  print(completeness)
  
  # Check for missing REM data specifically
  missing_rem <- df_long %>%
    filter(is.na(Event_Rate)) %>%
    group_by(Animal, Day_Code, State) %>%
    summarise(n_missing = n(), .groups = 'drop')
  
  if (nrow(missing_rem) > 0) {
    cat("\nMissing data found:\n")
    print(missing_rem)
  } else {
    cat("\nNo missing data found.\n")
  }
}

# Main function with better error handling
main <- function(file_pattern = "*.csv") {
  
  cat("=== Neural Activity Analysis - UPDATED VERSION ===\n\n")
  
  # Test the filename parsing first
  cat("Testing filename parsing...\n")
  test_sex_extraction()
  
  cat("Loading CSV files with pattern:", file_pattern, "\n")
  df <- load_csv_files_with_metadata(file_pattern)
  
  # Process data
  cat("\nProcessing data...\n")
  df_long <- prepare_data(df)
  
  # Run enhanced diagnostics
  check_data_issues(df_long)
  
  # Create individual plots for each animal
  cat("\nCreating individual animal plots...\n")
  create_individual_plots(df_long)
  
  # Create main multi-panel plot
  cat("\nCreating multi-panel visualization...\n")
  main_plot <- create_multipanel_plot(df_long)
  
  cat("\n=== Analysis Complete! ===\n")
  cat("Individual plots saved for each animal\n")
  cat("Main plot saved as: neural_activity_multipanel.png\n")
  cat("Current directory:", getwd(), "\n")
  
  return(list(data = df_long, main_plot = main_plot))
}

# Run the updated version
cat("=== UPDATED SCRIPT LOADED ===\n")
cat("Run: main('*normalized_event_counts.csv')\n")
cat("Or test parsing: test_sex_extraction()\n")