
# Debug function to inspect raw file contents
debug_file_contents <- function(file_path, n_rows = 5) {
  cat("=== DEBUGGING FILE CONTENTS ===\n")
  cat("File:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    cat("ERROR: File does not exist!\n")
    return(NULL)
  }
  
  # Read first few lines as text
  cat("\nFirst", n_rows, "lines as raw text:\n")
  lines <- readLines(file_path, n = n_rows)
  for (i in seq_along(lines)) {
    cat("Line", i, ":", lines[i], "\n")
  }
  
  # Try reading as CSV
  cat("\nReading as CSV:\n")
  tryCatch({
    df <- read.csv(file_path, nrows = n_rows, check.names = FALSE, stringsAsFactors = FALSE)
    cat("Dimensions:", nrow(df), "x", ncol(df), "\n")
    cat("Column names:", paste(colnames(df), collapse = ", "), "\n")
    cat("First row values:", paste(df[1,], collapse = ", "), "\n")
    if (nrow(df) > 1) {
      cat("Second row values:", paste(df[2,], collapse = ", "), "\n")
    }
  }, error = function(e) {
    cat("Error reading as CSV:", conditionMessage(e), "\n")
  })
  
  cat("=== END DEBUG ===\n\n")
}