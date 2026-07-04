#' Call Peaks from Binding Simulation
#'
#' Identifies peak regions based on binding occupancy.
#' If the transcript headers contain genomic coordinates (e.g. "chr1:1000-2000(+)"),
#' it appends translated genomic coordinate columns (chrom, genomic_start, genomic_end, strand).
#'
#' @param results data.table of simulation results.
#' @param rbp_primary Name of the primary RBP.
#' @param rbp_reference Name of the reference RBP (for competitive advantage) (default NULL).
#' @param method Method: "competitive" (default), "local", or "global".
#' @param threshold Threshold (ratio for competitive, occupancy for others) (default 2.0).
#' @param min_width Minimum width of peak to report (default 3).
#' @param window_size Window size for local background (default 50).
#' @return A data.table of peaks with start, end, scores, and optionally genomic coordinates.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGUACGUUUUUACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' results$transcript <- "test_transcript"
#' peaks <- callPeaks(results, rbp_primary = "HH", method = "global", threshold = 1.1)
#' head(peaks)
#' @importFrom data.table frollmean copy rleid
#' @export
callPeaks <- function(results, rbp_primary, rbp_reference = NULL, 
                       method = c("competitive", "local", "global"),
                       threshold = 2.0, min_width = 3, window_size = 50) {
  
  method <- match.arg(method)
  
  if (!rbp_primary %in% names(results)) stop("Primary RBP not found")
  
  group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
  if (is.null(group_col)) stop("No transcript column found")
  
  # Copy data to avoid modification
  dt <- copy(results)
  
  # Calculate Score
  if (method == "competitive") {
    if (is.null(rbp_reference)) stop("Reference RBP required for competitive method")
    if (!rbp_reference %in% names(results)) stop("Reference RBP not found")
    
    # Avoid div by zero
    dt[, score := get(rbp_primary) / (get(rbp_reference) + 1e-6)]
    
  } else if (method == "local") {
    # Calculate rolling mean background
    # Need to group by transcript
    dt[, bg := data.table::frollmean(get(rbp_primary), n = window_size, align = "center", fill = NA), by = group_col]
    # Handle edges
    dt[is.na(bg), bg := mean(get(rbp_primary), na.rm=TRUE), by = group_col]
    
    dt[, score := get(rbp_primary) / (bg + 1e-6)]
    
  } else if (method == "global") {
    # Global mean per transcript
    dt[, bg := mean(get(rbp_primary), na.rm=TRUE), by = group_col]
    dt[, score := get(rbp_primary) / (bg + 1e-6)]
  }
  
  # Identify peaks
  dt[, is_peak := score > threshold]
  
  if (!any(dt$is_peak, na.rm = TRUE)) {
    out <- data.table::data.table(
      placeholder = character(),
      start = integer(),
      end = integer(),
      length = integer(),
      max_score = numeric(),
      mean_score = numeric()
    )
    setnames(out, "placeholder", group_col)
    
    # Check if we should add empty genomic columns
    if (nrow(results) > 0) {
      sample_name <- results[[group_col]][1]
      parsed <- parseGenomicHeader(sample_name)
      if (!is.null(parsed)) {
        out$chrom <- character()
        out$genomic_start <- integer()
        out$genomic_end <- integer()
        out$strand <- character()
      }
    }
    return(out)
  }
  
  # Run-length encoding to merge
  dt[, cluster := rleid(is_peak), by = group_col]
  
  peaks <- dt[(is_peak), .(
    start = min(pos),
    end = max(pos),
    length = .N,
    max_score = max(score),
    mean_score = mean(score)
  ), by = c(group_col, "cluster")]
  
  peaks[, cluster := NULL]
  
  # Filter by min_width
  peaks <- peaks[length >= min_width]
  
  # Try to parse genomic coordinates from transcript name
  parsed_headers <- lapply(peaks[[group_col]], parseGenomicHeader)
  is_genomic <- !vapply(parsed_headers, is.null, logical(1))
  
  if (any(is_genomic)) {
    chroms <- character(nrow(peaks))
    g_starts <- integer(nrow(peaks))
    g_ends <- integer(nrow(peaks))
    strands <- character(nrow(peaks))
    
    for (i in seq_len(nrow(peaks))) {
      p <- parsed_headers[[i]]
      if (!is.null(p)) {
        chroms[i] <- p$chrom
        strands[i] <- p$strand
        if (p$strand == "-") {
          g_starts[i] <- p$end - peaks$end[i] + 1
          g_ends[i] <- p$end - peaks$start[i] + 1
        } else {
          g_starts[i] <- p$start + peaks$start[i] - 1
          g_ends[i] <- p$start + peaks$end[i] - 1
        }
      } else {
        chroms[i] <- NA_character_
        g_starts[i] <- NA_integer_
        g_ends[i] <- NA_integer_
        strands[i] <- NA_character_
      }
    }
    peaks[, `:=`(chrom = chroms, genomic_start = g_starts, genomic_end = g_ends, strand = strands)]
  }
  
  return(peaks)
}
