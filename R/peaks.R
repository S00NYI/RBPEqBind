#' Call Peaks from Binding Simulation
#'
#' Identifies peak regions based on binding occupancy.
#'
#' @param results data.table of simulation results.
#' @param rbp_primary Name of the primary RBP.
#' @param rbp_reference Name of the reference RBP (for competitive advantage).
#' @param method Method: "competitive", "local", or "global".
#' @param threshold Threshold (ratio for competitive, occupancy for others).
#' @param min_width Minimum width of peak to report.
#' @param window_size Window size for local background (default 50).
#' @return A data.table of peaks with start, end, scores.
#' @examples
#' \donttest{
#' peaks <- callPeaks(results, rbp_primary = "HH", method = "global", threshold = 1.5)
#' }
#' @importFrom data.table frollmean
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
    # Handle edges? fill=NA leaves NA.
    # Replace NA with global mean or 0?
    dt[is.na(bg), bg := mean(get(rbp_primary), na.rm=TRUE), by = group_col]
    
    dt[, score := get(rbp_primary) / (bg + 1e-6)]
    
  } else if (method == "global") {
    # Global mean per transcript
    dt[, bg := mean(get(rbp_primary), na.rm=TRUE), by = group_col]
    dt[, score := get(rbp_primary) / (bg + 1e-6)]
  }
  
  # Identify peaks
  dt[, is_peak := score > threshold]
  
  # Run-length encoding to merge
  # data.table way: grouping by consecutive TRUEs
  dt[, cluster := rleid(is_peak), by = group_col]
  
  peaks <- dt[is_peak == TRUE, .(
    start = min(pos),
    end = max(pos),
    length = .N,
    max_score = max(score),
    mean_score = mean(score)
  ), by = c(group_col, "cluster")]
  
  # Filter width
  peaks <- peaks[length >= min_width]
  
  # Format output
  # BED-like: start is 0-based? results usually 1-based pos.
  # Let's keep 1-based for now or follow R convention.
  # The user API spec usually implies 1-based unless exporting to BED.
  
  peaks[, cluster := NULL]
  
  return(peaks)
}
