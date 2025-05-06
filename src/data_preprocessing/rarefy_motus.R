#' Pool PCR replicates, rarefy, and export MOTU sequence tables
#'
#' @param meta_list   Named list of meta objects (each from read_meta()).
#'                    Example: meta_list <- read_meta(c("EXP1","EXP2"))
#' @param depth       Rarefaction depth (default 30 000).
#' @param out_dir     Folder where CSVs are written.  Will be created if missing.
#' @param parallel    Integer; number of CPU cores to use (1 = no parallel).
#' @return            Invisibly, a list with one data.frame per experiment and the rarefied matrix.
#' @export
#'
rarefy_motus <- function(meta_list,
                         depth     = 10000,
                         out_dir   = "data/processed",
                         parallel  = 1) {
  
  stopifnot(is.list(meta_list), depth > 0)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # helper that does the job for ONE experiment ---------------------------
  workhorse <- function(meta, exp_name) {
    message("• Processing ", exp_name)
    
    reads   <- meta$reads
    samples <- meta$samples
    motu    <- meta$motu
    
    # ----- pool PCR replicates -----
    pooled <- reads |>
      tidyr::pivot_longer(-id, names_to = "sample_id", values_to = "value") |>
      dplyr::left_join(samples, by = "sample_id") |>
      dplyr::group_by(id, sample_name) |>
      dplyr::summarise(value = sum(value), .groups = "drop") |>
      tidyr::pivot_wider(names_from = sample_name, values_from = value, values_fill = 0)
    
    mat <- as.matrix(pooled[,-1])
    rownames(mat) <- pooled$id
    
    # ----- rarefy -----
    if (any(colSums(mat) < depth)) {
      warning(exp_name, ": some samples have fewer than ", depth,
              " reads – they will be dropped.")
      keep <- which(colSums(mat) >= depth)
      mat  <- mat[, keep, drop = FALSE]
    }
    rare <- GUniFrac::Rarefy(t(mat), depth = depth)$otu.tab.rff
    rare <- tibble::rownames_to_column(as.data.frame(t(rare)), var = "id")
    
    # ----- save rarefied MOTU table -----
    rare_file <- file.path(out_dir, paste0(exp_name, "_rarefied_motu.csv"))
    readr::write_csv(rare, rare_file)
    message("  wrote rarefied MOTU table to ", rare_file)
    
    # ----- export MOTU sequences present after rarefaction -----
    motu_out <- dplyr::filter(motu, id %in% rare$id)
    motu_file <- file.path(out_dir, paste0(exp_name, "_motu_list.csv"))
    readr::write_csv(motu_out, motu_file)
    message("  wrote MOTU list to ", motu_file)
    
    # return the rarefied table for downstream work
    list(rarefied_motu = rare, motu_list = motu_out)
  }
  
  # run (optionally in parallel) ------------------------------------------
  exp_names <- names(meta_list)
  FUN <- function(i) workhorse(meta_list[[i]], exp_names[i])
  
  if (parallel > 1) {
    res <- pbmcapply::pbmclapply(seq_along(meta_list), FUN,
                                 mc.cores = parallel)
    names(res) <- exp_names
  } else {
    res <- setNames(lapply(seq_along(meta_list), FUN), exp_names)
  }
  invisible(res)
}


