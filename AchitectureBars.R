# SleepArchitecture_StackedBar.R

library(ggplot2)
library(dplyr)

# ============================================================
# CONFIGURATION
# ============================================================

INPUT_PATH <- "./Data CSVs"
OUTPUT_PATH <- "./sleep_architecture_figures"

ANIMALS <- c("mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9")
CONDITIONS <- c("BL", "SD", "WO")

STATE_COLORS <- c("Wake" = "#FFFF00", "NREM" = "#0000FF", "REM" = "#FF0000")

DAY_FIRST = FALSE #if TRUE, group by day and list animals for each

# ============================================================
# LOAD AND CALCULATE PERCENTAGES
# ============================================================

all_percentages <- data.frame()

for (animal in ANIMALS) {
  for (condition in CONDITIONS) {
    
    filename <- paste0(animal, "_", condition, "_states_df.csv")
    filepath <- file.path(INPUT_PATH, filename)
    
    if (file.exists(filepath)) {
      states_df <- read.csv(filepath)
      
      total_epochs <- nrow(states_df)
      
      pct <- states_df %>%
        group_by(state_label) %>%
        summarize(n_epochs = n(), .groups = 'drop') %>%
        mutate(
          percent = 100 * n_epochs / total_epochs,
          animal = animal,
          condition = condition
        )
      
      all_percentages <- rbind(all_percentages, pct)
      cat("Loaded:", filename, "-", total_epochs, "epochs\n")
      
    } else {
      cat("Missing:", filename, "\n")
    }
  }
}

# ============================================================
# PLOT
# ============================================================

all_percentages$condition <- factor(all_percentages$condition, 
                                    levels = c("BL", "SD", "WO"),
                                    labels = c("Baseline", "Sleep Dep", "Recovery"))
all_percentages$state_label <- factor(all_percentages$state_label,
                                      levels = c("REM", "NREM", "Wake"))

if (DAY_FIRST==TRUE) {
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

if (DAY_FIRST==FALSE) {
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

# Save
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)
ggsave(file.path(OUTPUT_PATH, "sleep_architecture_stacked.png"), p, 
       width = 10, height = 5, dpi = 300)

print(p)

# ============================================================
# BOUT IDENTIFICATION
# ============================================================

identify_bouts <- function(states_df) {
  states_df <- states_df[order(states_df$`Time.from.Start`), ]
  
  # Find where state changes
  state_vec <- states_df$state
  changes <- c(TRUE, state_vec[-1] != state_vec[-length(state_vec)])
  bout_id <- cumsum(changes)
  
  # Summarize each bout
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
# CALCULATE BOUT STATS FOR ALL RECORDINGS
# ============================================================

bout_summary <- data.frame()

for (animal in ANIMALS) {
  for (condition in CONDITIONS) {
    
    filename <- paste0(animal, "_", condition, "_states_df.csv")
    filepath <- file.path(INPUT_PATH, filename)
    
    if (file.exists(filepath)) {
      states_df <- read.csv(filepath)
      bouts <- identify_bouts(states_df)
      print(head(bouts))
      
      for (state in c("Wake", "NREM")) {
        state_bouts <- bouts[bouts$state_label == state, ]
        cat(animal, condition, state, "- n_bouts:", nrow(state_bouts), "\n")
        
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
    }
  }
}

# Format as mean ± SEM
bout_table <- bout_summary %>%
  select(animal, condition, state, mean_sec, sem) %>%
  pivot_wider(names_from = state, values_from = c(mean_sec, sem))



print(bout_table)

# Save
write.csv(bout_table, file.path(OUTPUT_PATH, "bout_length_summary.csv"), row.names = FALSE)