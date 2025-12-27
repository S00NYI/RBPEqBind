#' Load RBP Model from CSV
#'
#' Reads a CSV file containing RBP binding affinities/scores.
#' The file should have a 'Motif' column and one or more RBP score columns.
#'
#' @param file Path to the CSV file.
#' @param k Optional integer. If provided, filters for k-mers of this length.
#' @param rbp Optional character vector of RBP names to load. If NULL, loads all.
#' @return A named list of data.tables, each with 'motif' and 'score' columns.
#' @examples
#' # Load sample model from package
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' raw_models <- loadModel(model_file)
#' names(raw_models)
#'
#' # Load specific RBPs
#' raw_models <- loadModel(model_file, rbp = c("HH", "HL"))
#' @importFrom data.table fread
#' @export
loadModel <- function(file, k = NULL, rbp = NULL) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  dt <- data.table::fread(file)
  
  # Find motif column (case-insensitive)
  motif_col <- grep("^motif$", names(dt), ignore.case = TRUE, value = TRUE)
  if (length(motif_col) == 0) {
    stop("Input file must contain a 'Motif' column.")
  }
  
  # Get RBP columns (everything except motif)
  all_rbps <- setdiff(names(dt), motif_col)
  
  if (length(all_rbps) == 0) {
    stop("No RBP score columns found in file.")
  }
  
  # Filter to requested RBPs
  if (!is.null(rbp)) {
    missing <- setdiff(rbp, all_rbps)
    if (length(missing) > 0) {
      stop("RBPs not found in file: ", paste(missing, collapse = ", "))
    }
    all_rbps <- rbp
  }
  
  # Filter by k-mer length
  if (!is.null(k)) {
    dt <- dt[nchar(get(motif_col)) == k]
    if (nrow(dt) == 0) {
      warning("No k-mers of length ", k, " found in model file.")
      return(list())
    }
  }
  
  # Create named list of models
  models <- lapply(all_rbps, function(rbp_name) {
    data.table::data.table(
      motif = dt[[motif_col]],
      score = dt[[rbp_name]]
    )
  })
  names(models) <- all_rbps
  
  return(models)
}


#' Set Model Affinity (Score to Kd Conversion)
#'
#' Converts a model's normalized scores to Kd values using affinity scaling.
#' Score is assumed to be proportional to affinity (Ka), so Kd = 1/Ka.
#'
#' @param models Named list of models (from loadModel), each with 'motif' and 'score'.
#' @param max_affinity Named numeric vector of max Ka per RBP, or single default (100).
#' @param min_affinity Named numeric vector of min Ka per RBP, or single default (0.00001).
#' @return A named list of data.tables, each with 'motif' and 'Kd' columns.
#' @examples
#' # Load and process model
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' raw_models <- loadModel(model_file)
#' rbp_models <- setModel(raw_models, max_affinity = 100)
#' 
#' # Check Kd values
#' head(rbp_models[[1]])
#' @export
setModel <- function(models, max_affinity = 100, min_affinity = 0.00001) {
  
  get_value <- function(vec, name, default) {
    if (length(vec) == 1 && is.null(names(vec))) {
      return(vec)
    }
    if (!is.null(names(vec)) && name %in% names(vec)) {
      return(vec[name])
    }
    return(default)
  }
  
  result <- list()
  for (name in names(models)) {
    model <- models[[name]]
    
    if (!all(c("motif", "score") %in% names(model))) {
      stop("Model ", name, " must have 'motif' and 'score' columns.")
    }
    
    dt <- data.table::copy(model)
    current_max <- get_value(max_affinity, name, 100)
    current_min <- get_value(min_affinity, name, 0.00001)
    
    # Scale score to Ka range
    if (diff(range(dt$score, na.rm = TRUE)) == 0) {
      dt[, Ka := current_max]
    } else {
      dt[, Ka := featureScale(score, MAX = current_max, MIN = current_min)]
    }
    
    # Convert to Kd
    dt[, Kd := 1 / Ka]
    
    result[[name]] <- dt[, .(motif, Kd)]
  }
  
  return(result)
}
