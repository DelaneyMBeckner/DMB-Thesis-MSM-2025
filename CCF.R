setwd("C:\\Users\\TKDDM\\Downloads\\OneDrive_2025-05-02\\Data CSVs")

library(data.table)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)

animal = 'mPFCf5'
path1 = paste('.\\',animal,'\\', sep = '')
path2 = '_Traces.csv'

BLpath = paste(path1,'BL',path2, sep = '')
BLraw = read.csv(BLpath)
BL = BLraw[-c(1),-c(1)]
BL <- sapply(BL, as.numeric)

SDpath = paste(path1,'SD',path2, sep = '')
SDraw = read.csv(SDpath)
SD = SDraw[-c(1),-c(1)]
SD <- sapply(SD, as.numeric)

WOpath = paste(path1,'WO',path2, sep = '')
WOraw = read.csv(WOpath)
WO = WOraw[-c(1),-c(1)]
WO <- sapply(WO, as.numeric)

BL_heatmap = compute_pairwise_ccf_parallel(BL)
SD_heatmap = compute_pairwise_ccf_parallel(SD)
WO_heatmap = compute_pairwise_ccf_parallel(WO)

Title = paste(animal, 'Baseline Cross-correlation', sep = ' ')
BLmap = create_correlation_heatmap(BL_heatmap, title = Title)

Title = paste(animal, 'Sleep-Deprived Cross-correlation', sep = ' ')
SDmap = create_correlation_heatmap(SD_heatmap, title = Title)

Title = paste(animal, 'Recovery Cross-correlation', sep = ' ')
WOmap = create_correlation_heatmap(WO_heatmap, title = Title)

####Claude 3.7 Sonnet ####

compute_pairwise_ccf <- function(data_matrix, lag = 0, max_lag = 10) {
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix)
  
  if (is.null(series_names)) {
    series_names <- paste0("Series_", 1:n_series)
  }
  
  # Initialize the correlation matrix
  ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  rownames(ccf_matrix) <- series_names
  colnames(ccf_matrix) <- series_names
  
  # Loop through each pair of series
  for (i in 1:n_series) {
    # Auto-correlation at lag 0 is 1
    ccf_matrix[i, i] <- 1
    
    # Only process upper triangle to avoid redundant computations
    if (i < n_series) {
      for (j in (i+1):n_series) {
        # Create ts objects for just the pair we're analyzing
        series_i <- as.numeric(data_matrix[, i])
        series_j <- as.numeric(data_matrix[, j])
        
        # Calculate CCF and extract the correlation at specified lag
        ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
        lag_index <- which(ccf_result$lag == lag)
        
        if (length(lag_index) > 0) {
          correlation <- ccf_result$acf[lag_index]
          ccf_matrix[i, j] <- correlation
          ccf_matrix[j, i] <- correlation  # Matrix is symmetric for lag 0
        }
        
        # Clean up to free memory
        rm(series_i, series_j, ccf_result)
        gc()
      }
    }
    
    # Print progress
    if (i %% 5 == 0) {
      cat("Processed", i, "of", n_series, "series\n")
    }
  }
  
  return(ccf_matrix)
}

# Approach 2: Parallel processing version for multi-core systems
compute_pairwise_ccf_parallel <- function(data_matrix, lag = 0, max_lag = 10, n_cores = 4) {
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix)
  
  if (is.null(series_names)) {
    series_names <- paste0("Series_", 1:n_series)
  }
  
  # Initialize the correlation matrix
  ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  rownames(ccf_matrix) <- series_names
  colnames(ccf_matrix) <- series_names
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Create pairs for upper triangle
  pairs <- expand.grid(i = 1:n_series, j = 1:n_series)
  pairs <- pairs[pairs$i < pairs$j, ]
  
  # Compute correlations in parallel
  results <- foreach(idx = 1:nrow(pairs), .combine = rbind) %dopar% {
    i <- pairs$i[idx]
    j <- pairs$j[idx]
    
    series_i <- as.numeric(data_matrix[, i])
    series_j <- as.numeric(data_matrix[, j])
    
    ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
    lag_index <- which(ccf_result$lag == lag)
    
    if (length(lag_index) > 0) {
      correlation <- ccf_result$acf[lag_index]
      return(c(i, j, correlation))
    } else {
      return(c(i, j, NA))
    }
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Fill in the correlation matrix
  for (k in 1:nrow(results)) {
    i <- results[k, 1]
    j <- results[k, 2]
    correlation <- results[k, 3]
    
    ccf_matrix[i, j] <- correlation
    ccf_matrix[j, i] <- correlation  # Matrix is symmetric for lag 0
  }
  
  # Set diagonal to 1
  diag(ccf_matrix) <- 1
  
  return(ccf_matrix)
}

# Approach 3: Data.table efficient batch processing
compute_ccf_with_datatable <- function(data_matrix, lag = 0, max_lag = 10, batch_size = 10) {
  n_series <- ncol(data_matrix)
  series_names <- colnames(data_matrix)
  
  if (is.null(series_names)) {
    series_names <- paste0("Series_", 1:n_series)
  }
  
  # Initialize the correlation matrix
  ccf_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  rownames(ccf_matrix) <- series_names
  colnames(ccf_matrix) <- series_names
  
  # Convert to data.table for more efficient operations
  dt <- as.data.table(data_matrix)
  
  # Process in batches to control memory usage
  num_batches <- ceiling(n_series / batch_size)
  
  for (batch_i in 1:num_batches) {
    start_i <- (batch_i - 1) * batch_size + 1
    end_i <- min(batch_i * batch_size, n_series)
    
    for (batch_j in batch_i:num_batches) {
      start_j <- (batch_j - 1) * batch_size + 1
      end_j <- min(batch_j * batch_size, n_series)
      
      for (i in start_i:end_i) {
        # Auto-correlation at lag 0 is 1
        if (batch_i == batch_j) {
          ccf_matrix[i, i] <- 1
        }
        
        start_j_actual <- if (batch_i == batch_j) i + 1 else start_j
        
        if (start_j_actual <= end_j) {
          for (j in start_j_actual:end_j) {
            # Get series using data.table for better memory handling
            series_i <- dt[[i]]
            series_j <- dt[[j]]
            
            # Calculate CCF and extract correlation at lag
            ccf_result <- ccf(series_i, series_j, lag.max = max_lag, plot = FALSE)
            lag_index <- which(ccf_result$lag == lag)
            
            if (length(lag_index) > 0) {
              correlation <- ccf_result$acf[lag_index]
              ccf_matrix[i, j] <- correlation
              ccf_matrix[j, i] <- correlation  # Matrix is symmetric for lag 0
            }
          }
        }
      }
      
      # Force garbage collection after each batch
      gc()
      cat("Processed batch", batch_i, "/", batch_j, "of", num_batches, "\n")
    }
  }
  
  return(ccf_matrix)
}

create_correlation_heatmap <- function(ccf_matrix, title = "Cross-Correlation Heatmap", 
                                       triangle = "lower", show_diagonal = TRUE) {
  # Convert correlation matrix to long format for ggplot
  # Explicitly use reshape2::melt to avoid the data.table melt conflict
  melted_ccf <- reshape2::melt(ccf_matrix)
  names(melted_ccf) <- c("Series1", "Series2", "Correlation")
  
  # Filter data based on triangle parameter
  if (triangle == "upper") {
    # Keep only lower triangle (including or excluding diagonal based on show_diagonal)
    if (show_diagonal) {
      melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) >= as.numeric(melted_ccf$Series2), ]
    } else {
      melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) > as.numeric(melted_ccf$Series2), ]
    }
  } else if (triangle == "lower") {
    # Keep only upper triangle (including or excluding diagonal based on show_diagonal)
    if (show_diagonal) {
      melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) <= as.numeric(melted_ccf$Series2), ]
    } else {
      melted_ccf <- melted_ccf[as.numeric(melted_ccf$Series1) < as.numeric(melted_ccf$Series2), ]
    }
  }
  # If triangle == "full", keep all data
  
  # Create a heatmap using ggplot2
  heatmap <- ggplot(melted_ccf, aes(x = Series2, y = Series1, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "", y = "")
  
  # Only add text labels if there aren't too many variables
  if (nrow(ccf_matrix) <= 25) {
    heatmap <- heatmap + 
      geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3)
  }
  
  return(heatmap)
}