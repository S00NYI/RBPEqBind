#' Feature Scaling
#'
#' Scales a numeric vector to a given range [MIN, MAX].
#'
#' @param X Numeric vector to scale.
#' @param MAX Maximum value of the new range.
#' @param MIN Minimum value of the new range.
#' @return Scaled numeric vector.
#' @examples
#' featureScale(1:10, MAX = 100, MIN = 1)
#' @export
featureScale <- function(X, MAX, MIN) {
  if (MAX <= MIN) {
    stop("MAX must be greater than MIN")
  }
  X_min <- min(X, na.rm = TRUE)
  X_max <- max(X, na.rm = TRUE)
  
  if (X_max == X_min) return(rep(MIN, length(X)))
  
  return(MIN + ((X - X_min) * (MAX - MIN)) / (X_max - X_min))
}

#' Generate Random RNA Sequence
#'
#' Generates a random RNA sequence with specified nucleotide frequencies.
#'
#' @param L Length of the sequence.
#' @param A Frequency of Adenine.
#' @param G Frequency of Guanine.
#' @param C Frequency of Cytosine.
#' @param U Frequency of Uracil.
#' @return A character string representing the RNA sequence.
#' @examples
#' set.seed(42)
#' generateRNA(100)
#' @export
generateRNA <- function(L, A = 0.25, G = 0.25, C = 0.25, U = 0.25) {
  if (abs(sum(c(A, G, C, U)) - 1) > 1e-9) {
    stop("Frequencies must sum to 1")
  }
  
  nts <- c("A", "G", "C", "U")
  probs <- c(A, G, C, U)
  
  seq_vec <- sample(nts, size = L, replace = TRUE, prob = probs)
  return(paste(seq_vec, collapse = ""))
}


# ==============================================================================
# K-mer Functions
# ==============================================================================

#' Extract K-mers from Sequence
#'
#' Generates a vector of all k-mers from a given sequence using a sliding window.
#'
#' @param sequence Character string of the RNA sequence.
#' @param k Integer size of the k-mer (default 5).
#' @return Character vector of k-mers.
#' @examples
#' extractKmers("ACGUACGU", k = 3)
#' @export
extractKmers <- function(sequence, k = 5) {
  if (nchar(sequence) < k) {
    return(character(0))
  }
  n <- nchar(sequence)
  starts <- 1:(n - k + 1)
  ends <- k:n
  substring(sequence, starts, ends)
}

#' Count Unique K-mers
#'
#' Returns a data.table with unique k-mers and their counts.
#'
#' @param kmers Character vector of k-mers.
#' @return data.table with columns 'motif' and 'count'.
#' @examples
#' kmers <- extractKmers("ACGUACGU", k = 3)
#' countKmers(kmers)
#' @importFrom data.table data.table
#' @export
countKmers <- function(kmers) {
  dt <- data.table::data.table(motif = kmers)
  dt[, .(count = .N), by = motif]
}
