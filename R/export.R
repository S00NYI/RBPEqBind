#' Export Results
#'
#' Exports simulation results to JSON or CSV.
#'
#' @param results data.table of results.
#' @param output_file Output file path.
#' @param format "json" or "csv".
#' @param include_sequence Logical, include sequence in output (default TRUE).
#' @return Invisibly returns NULL. Called for side effect of writing file.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' 
#' # Export JSON
#' tmp_json <- tempfile(fileext = ".json")
#' exportResults(results, tmp_json, format = "json")
#' unlink(tmp_json)
#' 
#' # Export CSV
#' tmp_csv <- tempfile(fileext = ".csv")
#' exportResults(results, tmp_csv, format = "csv")
#' unlink(tmp_csv)
#' @importFrom jsonlite stream_out toJSON
#' @importFrom data.table fwrite
#' @export
exportResults <- function(results, output_file, format = c("json", "csv"), include_sequence = TRUE) {
  format <- match.arg(format)
  
  if (format == "json") {
    group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
    
    if (!is.null(group_col)) {
      exclude <- c("pos", "nt", group_col, "rna_conc", grep("Conc_", names(results), value = TRUE))
      rbp_cols <- setdiff(names(results), exclude)
      aggregated <- results[, lapply(.SD, function(x) list(x)), by = group_col, .SDcols = rbp_cols]
      
      if (include_sequence && "nt" %in% names(results)) {
        seqs <- results[, .(sequence = paste(nt, collapse = "")), by = group_col]
        aggregated <- merge(aggregated, seqs, by = group_col)
      }
      
      lens <- results[, .(length = .N), by = group_col]
      aggregated <- merge(aggregated, lens, by = group_col)
      setnames(aggregated, group_col, "name")
      jsonlite::stream_out(aggregated, file(output_file), verbose = FALSE)
    } else {
      jsonlite::write_json(results, output_file)
    }
  } else if (format == "csv") {
    data.table::fwrite(results, output_file)
  }
}

#' Export BED format
#'
#' Simple BED export for ranges exceeding threshold.
#' If the transcript headers contain genomic coordinates (e.g. "chr1:1000-2000(+)"),
#' it translates transcript-relative peak coordinates into genomic coordinates.
#'
#' @param results data.table.
#' @param output_file File path.
#' @param rbp RBP name to export.
#' @param threshold Threshold for score/occupancy.
#' @return Invisibly returns NULL. Called for side effect of writing BED file.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' results$transcript <- "test_transcript"
#' tmp_bed <- tempfile(fileext = ".bed")
#' exportBed(results, tmp_bed, rbp = "HH", threshold = 0.1)
#' unlink(tmp_bed)
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
  
  # Translate transcript coordinates to genomic coordinates if headers match UCSC style
  parsed_headers <- lapply(regions$chrom, parseGenomicHeader)
  is_genomic <- !vapply(parsed_headers, is.null, logical(1))
  
  if (any(is_genomic)) {
    chroms <- character(nrow(regions))
    starts <- integer(nrow(regions))
    ends <- integer(nrow(regions))
    
    for (i in seq_len(nrow(regions))) {
      p <- parsed_headers[[i]]
      if (!is.null(p)) {
        chroms[i] <- p$chrom
        if (p$strand == "-") {
          starts[i] <- p$end - regions$end[i]
          ends[i] <- p$end - regions$start[i]
        } else {
          starts[i] <- p$start + regions$start[i] - 1
          ends[i] <- p$start + regions$end[i] - 1
        }
      } else {
        chroms[i] <- regions$chrom[i]
        starts[i] <- regions$start[i]
        ends[i] <- regions$end[i]
      }
    }
    regions[, `:=`(chrom = chroms, start = starts, end = ends)]
  }
  
  # Select BED columns: chrom, start, end, name, score
  out <- regions[, .(chrom, start, end, name, score)]
  
  data.table::fwrite(out, output_file, sep="\t", col.names=FALSE)
}


#' Import Results
#'
#' Reads previously exported simulation results from a JSON or CSV file.
#'
#' @param input_file Path to JSON or CSV file (exported by exportResults).
#' @return A data.table with the imported results.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' 
#' # JSON Example
#' tmp_json <- tempfile(fileext = ".json")
#' exportResults(results, tmp_json, format = "json")
#' results_json <- importResults(tmp_json)
#' unlink(tmp_json)
#' 
#' # CSV Example
#' tmp_csv <- tempfile(fileext = ".csv")
#' exportResults(results, tmp_csv, format = "csv")
#' results_csv <- importResults(tmp_csv)
#' unlink(tmp_csv)
#' @importFrom jsonlite stream_in fromJSON
#' @importFrom data.table as.data.table rbindlist fread
#' @export
importResults <- function(input_file) {
  if (!file.exists(input_file)) stop("File not found: ", input_file)
  
  # Support CSV format directly
  ext <- tolower(tools::file_ext(input_file))
  if (ext == "csv") {
    return(data.table::fread(input_file))
  }
  
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
#' Converts simulation results to a SummarizedExperiment object.
#'
#' @param results data.table from simulateBinding or simulateBindingFasta.
#' @param rbp_models Optional named list of RBP models (for colData metadata).
#' @return A SummarizedExperiment with assays and row/col metadata.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' 
#' # Create density columns for the mock results to enable makeSE
#' results[, HH_density := HH / sum(HH)]
#' results[, HL_density := HL / sum(HL)]
#' results[, HH_density_fc := HH_density * 8]
#' results[, HL_density_fc := HL_density * 8]
#' results[, HH_occupancy_fc := HH / mean(HH)]
#' results[, HL_occupancy_fc := HL / mean(HL)]
#' 
#' se <- makeSE(results, rbp_models)
#' se
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
makeSE <- function(results, rbp_models = NULL) {
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
