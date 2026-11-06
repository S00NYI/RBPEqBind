# Lattice Binding Model
#
# One RNA molecule is a row of L nucleotides. A bound RBP j covers w_j
# consecutive nucleotides (its footprint), and two bound proteins may not
# overlap. A protein whose footprint starts at nucleotide s has the weight
# q_j(s) = f_j * a_j(s), where f_j is the free protein concentration and
# a_j(s) = 1 / Kd of the k-mer that sets the affinity of that placement.
#
# The partition function Z sums the weights of all allowed arrangements. A
# forward pass gives Z_f(i), the sum over arrangements inside nucleotides
# 1..i, and a backward pass gives Z_b(i), the sum inside nucleotides i..L.
# The probability that a footprint of RBP j starts at s is
#   Z_f(s - 1) * q_j(s) * Z_b(s + w_j) / Z.
# Both passes are kept on the log scale, so long sequences do not overflow.

# Start Probabilities on the Lattice
#
# @param a_list List with one numeric vector per RBP: 1 / Kd for each
#   placement start s = 1, ..., L - w_j + 1 (0 means no binding).
# @param footprint Integer vector of footprints w_j, one per RBP.
# @param free Numeric vector of free protein concentrations f_j.
# @param L Sequence length.
# @return A list with `starts` (list of start-probability vectors, one per
#   RBP), `n_bound` (expected number of bound proteins of each RBP per
#   molecule) and `logZ`.
latticeStarts <- function(a_list, footprint, free, L) {
  n <- length(a_list)
  w <- as.integer(footprint)
  qa <- lapply(seq_len(n), function(j) free[j] * a_list[[j]])
  has_sites <- vapply(qa, function(x) length(x) > 0L && any(x > 0),
                      logical(1))

  # lf[i + 1] = log Z_f(i); lf[1] = log Z_f(0) = 0
  lf <- numeric(L + 1L)
  for (i in seq_len(L)) {
    acc <- 1
    base <- lf[i]
    for (j in which(has_sites)) {
      s <- i - w[j] + 1L
      if (s >= 1L) {
        q <- qa[[j]][s]
        if (q > 0) acc <- acc + exp(lf[s] - base) * q
      }
    }
    lf[i + 1L] <- base + log(acc)
  }

  # lb[i] = log Z_b(i); lb[L + 1] = log Z_b(L + 1) = 0
  lb <- numeric(L + 1L)
  for (i in rev(seq_len(L))) {
    acc <- 1
    base <- lb[i + 1L]
    for (j in which(has_sites)) {
      e <- i + w[j]
      if (e <= L + 1L) {
        q <- qa[[j]][i]
        if (q > 0) acc <- acc + exp(lb[e] - base) * q
      }
    }
    lb[i] <- base + log(acc)
  }

  logZ <- lf[L + 1L]
  starts <- lapply(seq_len(n), function(j) {
    m <- length(qa[[j]])
    if (m == 0L) return(numeric(0))
    s <- seq_len(m)
    p <- numeric(m)
    pos <- qa[[j]] > 0
    p[pos] <- exp(lf[s[pos]] + log(qa[[j]][pos]) + lb[s[pos] + w[j]] - logZ)
    p
  })
  list(starts = starts,
       n_bound = vapply(starts, sum, numeric(1)),
       logZ = logZ)
}

# Shared Protein Pool on the Lattice
#
# Solves P_j = f_j + rna_conc * <N_j>(f) for all RBPs, where <N_j> is the
# expected number of bound proteins of RBP j per molecule. Uses the same
# one-protein-at-a-time scheme as solveSharedPool().
#
# @inheritParams latticeStarts
# @param rna_conc RNA concentration.
# @param protein_totals Total protein concentrations P_j.
# @param tol Relative tolerance on the free protein levels.
# @param max_iter Maximum number of sweeps over all proteins.
# @return The result of latticeStarts() at the solution, plus `free`,
#   `iterations` and `converged`.
solveLatticePool <- function(a_list, footprint, L, rna_conc, protein_totals,
                             tol = 1e-10, max_iter = 1000L) {
  free <- as.numeric(protein_totals)
  active <- which(free > 0 & vapply(a_list, function(a) any(a > 0),
                                    logical(1)))
  converged <- length(active) == 0L || rna_conc == 0
  iter <- 0L

  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    max_change <- 0
    for (j in active) {
      p_j <- protein_totals[j]
      h <- function(u) {
        f <- free
        f[j] <- exp(u)
        f[j] + rna_conc * latticeStarts(a_list, footprint, f, L)$n_bound[j] - p_j
      }
      hi <- log(p_j)
      lo <- max(hi - log(1e18), log(.Machine$double.xmin))
      new_f <- if (h(hi) <= 0) {
        p_j
      } else if (h(lo) >= 0) {
        exp(lo)
      } else {
        exp(stats::uniroot(h, c(lo, hi), tol = tol * 1e-2)$root)
      }
      max_change <- max(max_change, abs(log(new_f) - log(free[j])))
      free[j] <- new_f
    }
    if (max_change < tol) converged <- TRUE
  }
  if (!converged) {
    warning("Lattice protein pool did not converge in ", max_iter,
            " iterations.")
  }
  out <- latticeStarts(a_list, footprint, free, L)
  c(out, list(free = free, iterations = iter, converged = converged))
}

# Footprint for Each RBP
#
# @param footprint One integer for all RBPs, or a named integer vector with
#   one value per RBP.
# @param rbp_names Names of the simulated RBPs.
# @param k K-mer size. Each footprint must be at least k.
# @return Integer vector of footprints in the order of rbp_names.
resolveFootprint <- function(footprint, rbp_names, k) {
  if (!is.numeric(footprint) || anyNA(footprint) || any(!is.finite(footprint))) {
    stop("footprint must be NULL or finite numbers.")
  }
  if (any(footprint != round(footprint))) {
    stop("footprint must be whole numbers.")
  }
  if (length(footprint) == 1L && is.null(names(footprint))) {
    w <- rep(footprint, length(rbp_names))
  } else {
    if (is.null(names(footprint)) || !all(rbp_names %in% names(footprint))) {
      stop("A footprint vector with more than one value must be named, ",
           "with one value for each RBP in protein_concs.")
    }
    w <- footprint[rbp_names]
  }
  if (any(w < k)) {
    stop("Each footprint must be at least k (", k, ").")
  }
  stats::setNames(as.integer(w), rbp_names)
}

# Coverage From Start Probabilities
#
# The probability that nucleotide x is covered by RBP j is the sum of the
# start probabilities of all placements of j that include x.
#
# @param starts Start-probability vector for one RBP.
# @param w Footprint of that RBP.
# @param L Sequence length.
# @return Numeric vector of length L.
latticeCoverage <- function(starts, w, L) {
  m <- length(starts)
  if (m == 0L) return(numeric(L))
  cs <- c(0, cumsum(starts))
  x <- seq_len(L)
  hi <- pmin(x, m)
  lo <- pmax(x - w + 1L, 1L)
  cov <- cs[hi + 1L] - cs[lo]
  cov[lo > hi] <- 0
  pmax(cov, 0)
}
