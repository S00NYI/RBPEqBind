#' Read RNA Accessibility From an RNAplfold File
#'
#' Reads the unpaired probabilities that \code{RNAplfold -u} writes
#' (a \code{*_lunp} file) into the matrix that the \code{accessibility}
#' argument of \code{\link{simulateBinding}} takes.
#'
#' Row \eqn{i} of the result belongs to nucleotide \eqn{i}. Column \eqn{l}
#' holds the probability that the \eqn{l} nucleotides from \eqn{i - l + 1} to
#' \eqn{i} are all unpaired. Entries where the stretch does not fit
#' (\eqn{i < l}) are \code{NA}. RBPEqBind does not fold RNA. Make the file
#' with ViennaRNA, for example
#' \code{RNAplfold -W 150 -L 100 -u 7 < transcript.fa}, with \code{-u} at
#' least as large as \code{k}. Choose the folding settings for your
#' experiment.
#'
#' @param file Path to an RNAplfold \code{*_lunp} file.
#' @return A numeric matrix with one row per nucleotide and one column per
#'   stretch length (column names \code{"1"}, \code{"2"}, ...).
#' @seealso \code{\link{simulateBinding}} for how the values are used.
#' @examples
#' # A small file in the RNAplfold format (5 nucleotides, -u 2)
#' f <- tempfile(fileext = "_lunp")
#' writeLines(c("#unpaired probabilities", " #i$\tl=1\t2\t",
#'              "1\t0.9\tNA", "2\t0.8\t0.75", "3\t0.2\t0.15",
#'              "4\t0.3\t0.1", "5\t0.95\t0.3"), f)
#' readAccessibility(f)
#' @export
readAccessibility <- function(file) {
  if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
    stop("file must be the path of an existing RNAplfold *_lunp file.")
  }
  x <- utils::read.table(file, comment.char = "#", na.strings = "NA",
                         header = FALSE, colClasses = "numeric",
                         fill = TRUE)
  if (ncol(x) < 2L || nrow(x) == 0L) {
    stop("No unpaired probabilities found in ", file, ".")
  }
  if (!all(x[[1]] == seq_len(nrow(x)))) {
    stop("The first column of ", file, " must number the nucleotides 1, 2, ...")
  }
  m <- as.matrix(x[, -1, drop = FALSE])
  # RNAplfold ends each line with a tab, which can add an empty last column
  if (all(is.na(m[, ncol(m)])) && ncol(m) > 1L) m <- m[, -ncol(m), drop = FALSE]
  dimnames(m) <- list(NULL, as.character(seq_len(ncol(m))))
  m
}


# Check an Accessibility Matrix Against a Sequence
#
# @param accessibility NULL, or a numeric matrix (or data.frame) with one row
#   per nucleotide; column l is the probability that the l nucleotides ending
#   at the row are unpaired.
# @param L Sequence length.
# @param widths Stretch lengths the model needs (k).
# @param label Text for error messages.
# @return NULL, or a numeric matrix.
checkAccessibility <- function(accessibility, L, widths, label = "accessibility") {
  if (is.null(accessibility)) return(NULL)
  if (is.data.frame(accessibility)) accessibility <- as.matrix(accessibility)
  if (!is.matrix(accessibility) || !is.numeric(accessibility)) {
    stop(label, " must be a numeric matrix with one row per nucleotide ",
         "(see readAccessibility()).")
  }
  if (nrow(accessibility) != L) {
    stop(label, " has ", nrow(accessibility), " rows, but the sequence has ",
         L, " nucleotides.")
  }
  if (ncol(accessibility) < max(widths)) {
    stop(label, " has ", ncol(accessibility), " columns, but stretches of ",
         max(widths), " nucleotides are needed (k). ",
         "Use RNAplfold -u ", max(widths), " or more.")
  }
  for (w in unique(widths)) {
    v <- accessibility[seq.int(w, L), w]
    if (anyNA(v)) {
      stop(label, ": column ", w, " has NA values in rows ", w, " to ", L, ".")
    }
    if (any(v < 0 | v > 1)) {
      stop(label, " values must be between 0 and 1.")
    }
  }
  accessibility
}

# Unpaired Probability of Each Placement
#
# @param acc Checked accessibility matrix.
# @param w Stretch length (k).
# @return Numeric vector, one value per start s = 1, ..., L - w + 1: the
#   probability that nucleotides s to s + w - 1 are all unpaired.
placementAccessibility <- function(acc, w) {
  L <- nrow(acc)
  if (L < w) return(numeric(0))
  acc[seq.int(w, L), w]
}

# Accessibility for Each Sequence of a FASTA Run
#
# @param accessibility NULL, or a list of matrices named with the sequence
#   names.
# @param seq_names Sequence names.
# @return NULL, or a list in the order of seq_names.
resolveAccessibilityList <- function(accessibility, seq_names) {
  if (is.null(accessibility)) return(NULL)
  if (!is.list(accessibility) || is.data.frame(accessibility) ||
      is.null(names(accessibility))) {
    stop("For FASTA input, accessibility must be a list of matrices named ",
         "with the sequence names.")
  }
  if (anyDuplicated(seq_names)) {
    stop("Sequence names are not unique, so accessibility cannot be matched by name.")
  }
  missing <- setdiff(seq_names, names(accessibility))
  if (length(missing) > 0L) {
    stop("accessibility needs one matrix for each sequence. Missing: ",
         paste(utils::head(missing), collapse = ", "))
  }
  accessibility[seq_names]
}
