# Lattice Binding Model
#
# One RNA molecule is a row of L nucleotides. A bound RBP j covers w_j
# consecutive nucleotides (its footprint), and two bound proteins may not
# overlap. A protein whose footprint starts at nucleotide s has the weight
# q_j(s) = f_j * a_j(s), where f_j is the free protein concentration and
# a_j(s) = 1 / Kd of the k-mer that sets the affinity of that placement.
# When two proteins of the same RBP j touch (the second starts right after
# the first ends), the arrangement gets the extra factor omega_j.
#
# The partition function Z sums the weights of all allowed arrangements.
# The forward pass gives Z_f(i), the sum over arrangements inside nucleotides
# 1..i, and F_j(i), the part of Z_f(i) in which a protein of RBP j ends at i:
#   F_j(i) = q_j(s) * (Z_f(s - 1) + (omega_j - 1) * F_j(s - 1)), s = i - w_j + 1
#   Z_f(i) = Z_f(i - 1) + sum_j F_j(i)
# The backward pass gives G(i), the sum inside nucleotides i..L, and G_j(i),
# the same sum when a protein of RBP j ends at i - 1:
#   G(i)   = G(i + 1) + sum_m q_m(i) * G_m(i + w_m)
#   G_j(i) = G(i) + (omega_j - 1) * q_j(i) * G_j(i + w_j)
# The probability that a footprint of RBP j starts at s is
#   (Z_f(s - 1) + (omega_j - 1) * F_j(s - 1)) * q_j(s) * G_j(s + w_j) / Z.
# With omega_j = 1, F_j and G_j drop out. Both passes keep log Z_f and log G,
# and store F_j and G_j as ratios to them, so long sequences do not overflow.

# Start Probabilities on the Lattice
#
# @param a_list List with one numeric vector per RBP: 1 / Kd for each
#   placement start s = 1, ..., L - w_j + 1 (0 means no binding).
# @param footprint Integer vector of footprints w_j, one per RBP.
# @param free Numeric vector of free protein concentrations f_j.
# @param L Sequence length.
# @param omega Numeric vector of neighbor factors omega_j (>= 0), one per RBP.
# @return A list with `starts` (list of start-probability vectors, one per
#   RBP), `n_bound` (expected number of bound proteins of each RBP per
#   molecule) and `logZ`.
latticeStarts <- function(a_list, footprint, free, L,
                          omega = rep(1, length(a_list))) {
  n <- length(a_list)
  w <- as.integer(footprint)
  om1 <- as.numeric(omega) - 1
  coop <- any(om1 != 0)
  qa <- lapply(seq_len(n), function(j) free[j] * a_list[[j]])
  act <- which(vapply(qa, function(x) length(x) > 0L && any(x > 0),
                      logical(1)))

  # lf[i + 1] = log Z_f(i); rf[i + 1, j] = F_j(i) / Z_f(i)
  lf <- numeric(L + 1L)
  rf <- matrix(0, L + 1L, n)
  fj <- numeric(n)
  for (i in seq_len(L)) {
    acc <- 1
    base <- lf[i]
    fj[] <- 0
    for (j in act) {
      s <- i - w[j] + 1L
      if (s >= 1L) {
        q <- qa[[j]][s]
        if (q > 0) {
          x <- exp(lf[s] - base) * q
          if (coop) x <- x * (1 + om1[j] * rf[s, j])
          fj[j] <- x
          acc <- acc + x
        }
      }
    }
    lf[i + 1L] <- base + log(acc)
    if (coop) rf[i + 1L, ] <- fj / acc
  }

  # lb[i] = log G(i); ub[i, j] = G_j(i) / G(i); lb[L + 1] = 0, ub[L + 1, ] = 1
  lb <- numeric(L + 1L)
  ub <- matrix(1, L + 1L, n)
  for (i in rev(seq_len(L))) {
    acc <- 1
    base <- lb[i + 1L]
    for (j in act) {
      e <- i + w[j]
      if (e <= L + 1L) {
        q <- qa[[j]][i]
        if (q > 0) acc <- acc + exp(lb[e] - base) * q * ub[e, j]
      }
    }
    lb[i] <- base + log(acc)
    if (coop) {
      for (j in act) {
        e <- i + w[j]
        if (e <= L + 1L && om1[j] != 0) {
          q <- qa[[j]][i]
          if (q > 0) ub[i, j] <- 1 + om1[j] * q * exp(lb[e] - lb[i]) * ub[e, j]
        }
      }
    }
  }

  logZ <- lf[L + 1L]
  starts <- lapply(seq_len(n), function(j) {
    m <- length(qa[[j]])
    if (m == 0L) return(numeric(0))
    s <- seq_len(m)
    p <- numeric(m)
    pos <- qa[[j]] > 0
    sp <- s[pos]
    lp <- lf[sp] + log(qa[[j]][pos]) + lb[sp + w[j]] - logZ
    if (coop) {
      lp <- lp + log(1 + om1[j] * rf[sp, j]) + log(ub[sp + w[j], j])
    }
    p[pos] <- exp(lp)
    p
  })
  list(starts = starts,
       n_bound = vapply(starts, sum, numeric(1)),
       logZ = logZ)
}

# Shared Protein Pool on the Lattice
#
# Solves P_j = f_j + sum_t c_t * <N_j>_t(f) for all RBPs, where molecule t has
# concentration c_t and <N_j>_t is its expected number of bound proteins of
# RBP j. Each molecule is its own lattice, and all molecules draw on one
# supply of each protein. Uses the same one-protein-at-a-time scheme as
# solveSharedPool().
#
# @param molecules List with one entry per molecule. Each entry is a list
#   with `a_list` (see latticeStarts()), `L` (length) and `rna_conc`.
# @inheritParams latticeStarts
# @param protein_totals Total protein concentrations P_j.
# @param tol Relative tolerance on the free protein levels.
# @param max_iter Maximum number of sweeps over all proteins.
# @return A list with `molecules` (the result of latticeStarts() for each
#   molecule at the solution), `free`, `iterations` and `converged`.
solveLatticePoolMulti <- function(molecules, footprint, protein_totals,
                                  omega = rep(1, length(protein_totals)),
                                  tol = 1e-10, max_iter = 1000L) {
  free <- as.numeric(protein_totals)
  binds <- vapply(seq_along(free), function(j) {
    any(vapply(molecules, function(m) {
      m$rna_conc > 0 && length(m$a_list[[j]]) > 0L && any(m$a_list[[j]] > 0)
    }, logical(1)))
  }, logical(1))
  active <- which(free > 0 & binds)
  converged <- length(active) == 0L
  iter <- 0L

  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    max_change <- 0
    for (j in active) {
      p_j <- protein_totals[j]
      h <- function(u) {
        f <- free
        f[j] <- exp(u)
        bound <- 0
        for (m in molecules) {
          if (m$rna_conc > 0) {
            bound <- bound + m$rna_conc *
              latticeStarts(m$a_list, footprint, f, m$L, omega)$n_bound[j]
          }
        }
        f[j] + bound - p_j
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
  outs <- lapply(molecules, function(m) {
    latticeStarts(m$a_list, footprint, free, m$L, omega)
  })
  list(molecules = outs, free = free, iterations = iter,
       converged = converged)
}

# Shared Protein Pool on the Lattice, One Molecule
#
# solveLatticePoolMulti() for one molecule at concentration rna_conc.
#
# @inheritParams latticeStarts
# @param rna_conc RNA concentration.
# @inheritParams solveLatticePoolMulti
# @return The result of latticeStarts() at the solution, plus `free`,
#   `iterations` and `converged`.
solveLatticePool <- function(a_list, footprint, L, rna_conc, protein_totals,
                             omega = rep(1, length(a_list)),
                             tol = 1e-10, max_iter = 1000L) {
  sol <- solveLatticePoolMulti(
    list(list(a_list = a_list, L = L, rna_conc = rna_conc)),
    footprint, protein_totals, omega = omega, tol = tol, max_iter = max_iter)
  c(sol$molecules[[1]], sol[c("free", "iterations", "converged")])
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

# Neighbor Factor for Each RBP
#
# @param omega One number for all RBPs, or a named vector with one value per
#   RBP. Each value must be >= 0.
# @param rbp_names Names of the simulated RBPs.
# @return Numeric vector of omega values in the order of rbp_names.
resolveOmega <- function(omega, rbp_names) {
  if (!is.numeric(omega) || anyNA(omega) || any(!is.finite(omega)) ||
      any(omega < 0)) {
    stop("omega must be finite numbers >= 0.")
  }
  if (length(omega) == 1L && is.null(names(omega))) {
    om <- rep(omega, length(rbp_names))
  } else {
    if (is.null(names(omega)) || !all(rbp_names %in% names(omega))) {
      stop("An omega vector with more than one value must be named, ",
           "with one value for each RBP in protein_concs.")
    }
    om <- omega[rbp_names]
  }
  stats::setNames(as.numeric(om), rbp_names)
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
