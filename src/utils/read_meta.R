read_meta <- function(experiments) {
  for(exp in experiments) {
    cat("Processing experiment:", exp, "\n")
    
    # Construct file names for the current experiment
    reads_file    <- paste0(exp, "_reads.csv")
    motu_file    <- paste0(exp, "_motu.csv")
    samples_file  <- paste0(exp, "_samples.csv")
    stations_file <- paste0(exp, "_stations.csv")
    
    # Read in the data tables
    reads_data    <- read.csv(reads_file, header = TRUE, stringsAsFactors = FALSE)
    motu_data    <- read.csv(motu_file, header = TRUE, stringsAsFactors = FALSE)
    samples_data  <- read.csv(samples_file, header = TRUE, stringsAsFactors = FALSE)
    stations_data <- read.csv(stations_file, header = TRUE, stringsAsFactors = FALSE)
    
    # ---------------------------
    # Check 1:
    # Every column name in reads (except "id") should appear in samples_data$sample_id.
    # If any are missing, report and drop those columns.
    # ---------------------------
    read_cols <- colnames(reads_data)
    read_sample_cols <- setdiff(read_cols, "id")
    missing_samples <- setdiff(read_sample_cols, samples_data$sample_id)
    if(length(missing_samples) > 0) {
      cat("  WARNING: sample_id(s) not found in samples file; these columns will be dropped: ",
          paste(missing_samples, collapse = ", "), "\n")
      # Keep only the "id" column and those that appear in samples_data$sample_id
      valid_cols <- intersect(read_sample_cols, samples_data$sample_id)
      reads_data <- reads_data[, c("id", valid_cols)]
    } else {
      cat("  All read column names (except 'id') have matching sample_id in the samples file.\n")
    }
    
    # ---------------------------
    # Check 2:
    # Every "id" in reads_data must appear in motu_data$id.
    # ---------------------------
    missing_ids <- setdiff(reads_data$id, motu_data$id)
    if(length(missing_ids) > 0) {
      cat("  ERROR: The following id(s) in the reads table are missing in the motu table: ",
          paste(missing_ids, collapse = ", "), "\n")
    } else {
      cat("  All id values in reads match those in motu.\n")
    }
    
    # ---------------------------
    # Check 3:
    # Every "sample_name" in samples_data should appear in stations_data$sample_name.
    # ---------------------------
    missing_sample_names <- setdiff(samples_data$sample_name, stations_data$sample_name)
    if(length(missing_sample_names) > 0) {
      cat("  ERROR: The following sample_name value(s) in samples are missing in stations: ",
          paste(missing_sample_names, collapse = ", "), "\n")
    } else {
      cat("  All sample_name values in samples have matching entries in stations.\n")
    }
    
    # ---------------------------
    # Create the meta object for this experiment
    # ---------------------------
    meta_object <- list(
      reads    = reads_data,
      motu    = motu_data,
      samples  = samples_data,
      stations = stations_data
    )
    
    # Build a name for the meta object by appending "_meta" to the experiment name
    meta_name <- paste0(exp, "_meta")
    
    # Assign the meta object to the global environment
    assign(meta_name, meta_object, envir = .GlobalEnv)
    cat("  Created meta object:", meta_name, "\n")
    cat("Finished processing", exp, "\n\n")
  }
  
  cat("All experiments processed. Individual meta objects are now in the global environment.\n")
}
