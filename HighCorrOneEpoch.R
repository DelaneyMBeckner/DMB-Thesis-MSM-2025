# Read with headers (gets C00, C01, etc. as column names)
traces <- read.csv("E:/Data_Processing/R/Data CSVs/mPFCf5_BL_Traces.csv", stringsAsFactors = FALSE)

# Remove the status row (first data row with "undecided"s)
traces <- traces[-1, ]

# Convert to numeric (status row made them character)
traces <- traces %>% mutate(across(everything(), as.numeric))

# First column is time (named "X" from the blank header)
colnames(traces)[1] <- "Time"

# Define epoch (10 sec epochs at 10 Hz = 100 samples)
target_epoch <- 21
epoch_start <- (target_epoch - 1) * 100
epoch_end <- target_epoch * 100

# Subset to that epoch
epoch_traces <- traces %>%
  filter(row_number() > epoch_start & row_number() <= epoch_end) %>%
  select(-Time)

# Compute correlation matrix
cor_matrix <- cor(epoch_traces, use = "pairwise.complete.obs")

# Convert to long format, remove diagonal and duplicates
pair_corrs <- as.data.frame(as.table(cor_matrix)) %>%
  rename(ROI_1 = Var1, ROI_2 = Var2, correlation = Freq) %>%
  filter(as.character(ROI_1) < as.character(ROI_2)) %>%
  arrange(desc(correlation))

top_pos <- head(pair_corrs, 10)
top_neg <- tail(pair_corrs, 10)