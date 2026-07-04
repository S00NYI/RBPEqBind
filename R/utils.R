# Feature Scaling
#
# Scales a numeric vector to a given range [MIN, MAX].
#
# @param X Numeric vector to scale.
# @param MAX Maximum value of the new range.
# @param MIN Minimum value of the new range.
# @return Scaled numeric vector.
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

# Extract K-mers from Sequence
#
# Generates a vector of all k-mers from a given sequence using a sliding window.
#
# @param sequence Character string of the RNA sequence.
# @param k Integer size of the k-mer (default 5).
# @return Character vector of k-mers.
extractKmers <- function(sequence, k = 5) {
  if (nchar(sequence) < k) {
    return(character(0))
  }
  n <- nchar(sequence)
  starts <- seq_len(n - k + 1)
  ends <- k:n
  substring(sequence, starts, ends)
}

# Count Unique K-mers
#
# Returns a data.table with unique k-mers and their counts.
#
# @param kmers Character vector of k-mers.
# @return data.table with columns 'motif' and 'count'.
# @importFrom data.table data.table
countKmers <- function(kmers) {
  dt <- data.table::data.table(motif = kmers)
  dt[, .(count = .N), by = motif]
}

#' Parse Genomic Coordinates from FASTA Header
#'
#' Extracts genomic coordinates from a UCSC-style header string,
#' e.g. "chr1:1000-2000" or "chr1:1000-2000(-)".
#'
#' @param h Character string.
#' @return A list with chrom, start, end, and strand, or NULL.
parseGenomicHeader <- function(h) {
  if (is.null(h) || length(h) != 1 || is.na(h)) return(NULL)
  m <- regexec("^([^:]+):([0-9]+)-([0-9]+)(?:\\(([+-])\\))?$", h)
  reg_match <- regmatches(h, m)[[1]]
  if (length(reg_match) >= 5) {
    chrom <- reg_match[2]
    start_coord <- as.integer(reg_match[3])
    end_coord <- as.integer(reg_match[4])
    strand <- if (reg_match[5] != "") reg_match[5] else "+"
    return(list(chrom = chrom, start = start_coord, end = end_coord, strand = strand))
  }
  return(NULL)
}
