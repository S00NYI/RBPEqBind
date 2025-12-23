#' Export Results
#'
#' Exports simulation results to JSON or CSV.
#'
#' @param results data.table of results.
#' @param output_file Output file path.
#' @param format "json" or "csv".
#' @param include_sequence Logical, whether to include sequence in output (default TRUE).
#' @importFrom jsonlite stream_out toJSON
#' @importFrom data.table fwrite
#' @export
exportResults <- function(results, output_file, format = c("json", "csv"), include_sequence = TRUE) {
  format <- match.arg(format)
  
  if (format == "json") {
    # JSONL format recommended
    # We want one line per transcript or sequence?
    # results is a long data.table (per position).
    # Group by sequence/transcript?
    # Usually results have 'seq_name' or 'transcript'.
    
    group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
    
    if (!is.null(group_col)) {
      # Grouped export
      # split by group is slow for large data.
      # better to iterate?
      # Or just dump the data.table as JSON array of objects?
      # User spec:
      # {"name":"tx001","length":200,"sequence":"...", "RBP1":[...], "RBP2":[...]}
      # This requires restructuring!
      
      # We need to pivot/flatten the per-position data into arrays.
      # This is expensive for huge data.
      
      # Helper to aggregates
      # Get RBP columns
      # non-id columns
      exclude <- c("pos", "nt", group_col, "rna_conc", grep("Conc_", names(results), value=TRUE)) # Adjust as needed
      rbp_cols <- setdiff(names(results), exclude)
      
      # Efficient approach:
      # Iterate unique transcripts
      # But iterating 100k transcripts in R is slow.
      # data.table grouping?
      
      # results[, .(RBP1 = list(RBP1), RBP2 = list(RBP2)), by = group_col]
      # This creates list columns.
      # Then toJSON.
      
      aggregated <- results[, lapply(.SD, function(x) list(x)), by = group_col, .SDcols = rbp_cols]
      
      # Add sequence if available?
      # We assume sequence is not in results data.table per pos (inefficient).
      # If 'nt' column exists, we can reconstruct sequence
      if (include_sequence && "nt" %in% names(results)) {
        seqs <- results[, .(sequence = paste(nt, collapse="")), by = group_col]
        aggregated <- merge(aggregated, seqs, by = group_col)
      }
      
      # Add length
      lens <- results[, .(length = .N), by = group_col]
      aggregated <- merge(aggregated, lens, by = group_col)
      
      # Rename group col to 'name'?
      setnames(aggregated, group_col, "name")
      
      # Write JSONL
      jsonlite::stream_out(aggregated, file(output_file), verbose = FALSE)
      
    } else {
      # No group, just dump
      jsonlite::write_json(results, output_file)
    }
    
  } else if (format == "csv") {
    data.table::fwrite(results, output_file)
  }
}

#' Export BED format
#'
#' Simple BED export for ranges exceeding threshold.
#'
#' @param results data.table.
#' @param output_file File path.
#' @param rbp RBP name to export.
#' @param threshold Threshold for score/occupancy.
#' @export
exportBed <- function(results, output_file, rbp, threshold = 0) {
  # Need pos, transcript/seq_name
  group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
  if (is.null(group_col)) stop("Results must have transcript or seq_name column")
  
  if (!rbp %in% names(results)) stop("RBP not found in results")
  
  # Filter
  bed_dt <- results[get(rbp) > threshold]
  if (nrow(bed_dt) == 0) {
    warning("No regions above threshold")
    file.create(output_file)
    return()
  }
  
  # Coordinate conversion: BED is 0-indexed, half-open.
  # pos is 1-indexed.
  # start = pos - 1
  # end = pos
  
  # Naive export: one line per position?
  # Better to merge adjacent
  # Use GenomicRanges logic if available, or manual IRanges-like reduce using data.table
  
  # Sort
  setkeyv(bed_dt, c(group_col, "pos"))
  
  # Identify jumps
  bed_dt[, g := cumsum(c(0, diff(pos) > 1)), by = group_col]
  
  # Aggregate
  regions <- bed_dt[, .(chrom = first(get(group_col)),
                        start = first(pos) - 1,
                        end = last(pos),
                        name = paste0(rbp, "_peak"),
                        score = mean(get(rbp))), by = .(get(group_col), g)]
  
  # Select BED columns: chrom, start, end, name, score
  out <- regions[, .(chrom, start, end, name, score)]
  
  data.table::fwrite(out, output_file, sep="\t", col.names=FALSE)
}


#' Import Results from JSON
#'
#' Reads previously exported simulation results from a JSON file.
#'
#' @param input_file Path to JSON file (exported by exportResults).
#' @return A data.table with the imported results.
#' @importFrom jsonlite stream_in fromJSON
#' @importFrom data.table as.data.table rbindlist
#' @export
importResults <- function(input_file) {
  if (!file.exists(input_file)) stop("File not found: ", input_file)
  
  # Try JSONL format first (stream_in)
  tryCatch({
    conn <- file(input_file, "r")
    on.exit(close(conn))
    
    # Read as JSONL (one object per line)
    raw <- jsonlite::stream_in(conn, verbose = FALSE)
    
    # If it has list columns (arrays), we need to expand them
    if (!is.null(raw) && nrow(raw) > 0) {
      # Check if any columns are lists (arrays from export)
      list_cols <- vapply(raw, is.list, logical(1))
      
      if (any(list_cols)) {
        # Need to expand list columns to long format
        # Get the RBP/metric columns that are lists
        rbp_cols <- names(list_cols)[list_cols]
        id_cols <- names(list_cols)[!list_cols]
        
        # Expand each row
        results_list <- lapply(seq_len(nrow(raw)), function(i) {
          row <- raw[i, ]
          n_pos <- length(row[[rbp_cols[1]]][[1]])
          
          expanded <- data.table::data.table(pos = seq_len(n_pos))
          
          # Add ID columns
          for (col in id_cols) {
            expanded[[col]] <- row[[col]]
          }
          
          # Add RBP columns
          for (col in rbp_cols) {
            expanded[[col]] <- unlist(row[[col]])
          }
          
          expanded
        })
        
        return(data.table::rbindlist(results_list))
      } else {
        return(data.table::as.data.table(raw))
      }
    }
    
    return(data.table::as.data.table(raw))
    
  }, error = function(e) {
    # Fall back to regular JSON
    raw <- jsonlite::fromJSON(input_file)
    return(data.table::as.data.table(raw))
  })
}


#' Create SummarizedExperiment from Results
#'
#' Converts simulation results to a SummarizedExperiment object with
#' occupancy, density, and fold-change assays.
#'
#' @param results data.table from simulateBinding or simulateBindingF.
#' @param rbp_models Optional named list of RBP models (for colData metadata).
#' @return A SummarizedExperiment object with:
#'   - assays: occupancy, density, density_fc, occupancy_fc (each n_pos × n_rbp)
#'   - rowData: pos, nt, transcript
#'   - colData: RBP name, Kd_min, Kd_max
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
makeSE <- function(results, rbp_models = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required. Install with: BiocManager::install('SummarizedExperiment')")
  }
  
  # Identify RBP columns (those without suffix that have corresponding _density)
  all_cols <- names(results)
  density_cols <- grep("_density$", all_cols, value = TRUE)
  rbp_names <- gsub("_density$", "", density_cols)
  
  if (length(rbp_names) == 0) {
    stop("No RBP columns found. Ensure results have occupancy and density columns.")
  }
  
  # Build rowData
  row_cols <- intersect(c("pos", "nt", "transcript", "seq_name"), all_cols)
  if (length(row_cols) == 0) {
    stop("No position columns (pos, nt, transcript) found in results.")
  }
  
  row_data <- S4Vectors::DataFrame(results[, ..row_cols])
  
  n_pos <- nrow(results)
  n_rbp <- length(rbp_names)
  
  # Build assay matrices
  # Each assay is n_pos rows × n_rbp cols
  make_matrix <- function(suffix) {
    cols <- paste0(rbp_names, suffix)
    missing <- setdiff(cols, all_cols)
    if (length(missing) > 0) {
      warning("Missing columns for assay: ", paste(missing, collapse = ", "))
      return(NULL)
    }
    mat <- as.matrix(results[, ..cols])
    colnames(mat) <- rbp_names
    return(mat)
  }
  
  assays <- list(
    occupancy = make_matrix(""),
    density = make_matrix("_density"),
    density_fc = make_matrix("_density_fc"),
    occupancy_fc = make_matrix("_occupancy_fc")
  )
  
  # Remove NULL assays
  assays <- assays[!vapply(assays, is.null, logical(1))]
  
  if (length(assays) == 0) {
    stop("No valid assays could be created.")
  }
  
  # Build colData (RBP metadata)
  col_data <- S4Vectors::DataFrame(
    RBP = rbp_names,
    row.names = rbp_names
  )
  
  # Add Kd ranges if rbp_models provided
  if (!is.null(rbp_models)) {
    kd_min <- vapply(rbp_names, function(r) {
      if (r %in% names(rbp_models) && "Kd" %in% names(rbp_models[[r]])) {
        min(rbp_models[[r]]$Kd, na.rm = TRUE)
      } else {
        NA_real_
      }
    }, numeric(1))
    kd_max <- vapply(rbp_names, function(r) {
      if (r %in% names(rbp_models) && "Kd" %in% names(rbp_models[[r]])) {
        max(rbp_models[[r]]$Kd, na.rm = TRUE)
      } else {
        NA_real_
      }
    }, numeric(1))
    col_data$Kd_min <- kd_min
    col_data$Kd_max <- kd_max
  }
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = row_data,
    colData = col_data
  )
  
  return(se)
}
