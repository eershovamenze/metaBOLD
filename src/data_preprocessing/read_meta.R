read_meta <- function(experiments, data_dir = "data/raw") {
  # ------------------------------------------------------------
  # data_dir : folder that contains the CSV files.  Two layouts
  #            are accepted for each experiment ID <exp>:
  #   1) data_dir/<exp>_reads.csv          (flat files)
  #   2) data_dir/<exp>/<exp>_reads.csv    (each experiment in its own sub-folder)
  # ------------------------------------------------------------
  
  for (exp in experiments) {
    cat("Processing experiment:", exp, "\n")
    
    # helper that tries both layouts and stops if file not found
    fp <- function(suffix) {
      flat  <- file.path(data_dir,            paste0(exp, "_", suffix, ".csv"))
      nested <- file.path(data_dir, exp,      paste0(exp, "_", suffix, ".csv"))
      if      (file.exists(flat))   return(flat)
      else if (file.exists(nested)) return(nested)
      else stop("Cannot locate ", flat, " or ", nested)
    }
    
    # Construct full paths for the four tables
    reads_file    <- fp("reads")
    motu_file     <- fp("motu")
    samples_file  <- fp("samples")
    stations_file <- fp("stations")
    
    # ---------- read the data ----------
    reads_data    <- read.csv(reads_file,    header = TRUE, stringsAsFactors = FALSE)
    motu_data     <- read.csv(motu_file,     header = TRUE, stringsAsFactors = FALSE)
    samples_data  <- read.csv(samples_file,  header = TRUE, stringsAsFactors = FALSE)
    stations_data <- read.csv(stations_file, header = TRUE, stringsAsFactors = FALSE)
    
    # ---------- Check 1 ----------
    read_cols        <- colnames(reads_data)
    read_sample_cols <- setdiff(read_cols, "id")
    missing_samples  <- setdiff(read_sample_cols, samples_data$sample_id)
    if (length(missing_samples) > 0) {
      cat("  WARNING: sample_id(s) not found in samples file; these columns will be dropped: ",
          paste(missing_samples, collapse = ", "), "\n")
      valid_cols <- intersect(read_sample_cols, samples_data$sample_id)
      reads_data <- reads_data[, c("id", valid_cols)]
    } else {
      cat("  All read column names (except 'id') have matching sample_id in the samples file.\n")
    }
    
    # ---------- Check 2 ----------
    missing_ids <- setdiff(reads_data$id, motu_data$id)
    if (length(missing_ids) > 0) {
      cat("  ERROR: The following id(s) in the reads table are missing in the motu table: ",
          paste(missing_ids, collapse = ", "), "\n")
    } else {
      cat("  All id values in reads match those in motu.\n")
    }
    
    # ---------- Check 3 ----------
    missing_sample_names <- setdiff(samples_data$sample_name, stations_data$sample_name)
    if (length(missing_sample_names) > 0) {
      cat("  ERROR: The following sample_name value(s) in samples are missing in stations: ",
          paste(missing_sample_names, collapse = ", "), "\n")
    } else {
      cat("  All sample_name values in samples have matching entries in stations.\n")
    }
    
    # ---------- build meta object ----------
    meta_object <- list(
      reads    = reads_data,
      motu     = motu_data,
      samples  = samples_data,
      stations = stations_data
    )
    
    meta_name <- paste0(exp, "_meta")
    assign(meta_name, meta_object, envir = .GlobalEnv)
    cat("  Created meta object:", meta_name, "\n")
    cat("Finished processing", exp, "\n\n")
  }
  
  cat("All experiments processed. Individual meta objects are now in the global environment.\n")
}