# Shared Protein Pool Solver
#
# Solves the protein mass balance when all binding sites draw on one pool of
# each protein. Sites bind independently of each other (no footprint).
#
# For free protein levels f_j, the fraction of site i bound by protein j is
#   theta_ij = (f_j / Kd_ij) / (1 + sum_m f_m / Kd_im).
# The free levels satisfy, for each protein j,
#   P_j = f_j + sum_i S_i * theta_ij,
# where S_i is the concentration of site i and P_j the total protein.
#
# The solver updates one protein at a time (Gauss-Seidel). Each update is a
# root search on log(f_j), which keeps the result independent of the
# concentration unit. It repeats until no free level changes by more than
# `tol` (relative).
#
# @param site_conc Numeric vector of site concentrations S_i (>= 0).
# @param kd_mat Numeric matrix of Kd values, one row per site and one column
#   per protein. Inf means no binding.
# @param protein_totals Numeric vector of total protein concentrations P_j (>= 0).
# @param tol Relative tolerance on the free protein levels.
# @param max_iter Maximum number of sweeps over all proteins.
# @return A list with `free` (free protein for each protein), `theta` (matrix of
#   fractions bound, sites by proteins), `bound` (bound protein for each
#   protein), `iterations` and `converged`.
solveSharedPool <- function(site_conc, kd_mat, protein_totals,
                            tol = 1e-10, max_iter = 1000L) {
  kd_mat <- as.matrix(kd_mat)
  n_prot <- ncol(kd_mat)
  if (length(protein_totals) != n_prot) {
    stop("protein_totals must have one value for each column of kd_mat.")
  }
  if (length(site_conc) != nrow(kd_mat)) {
    stop("site_conc must have one value for each row of kd_mat.")
  }
  if (anyNA(site_conc) || anyNA(kd_mat) || anyNA(protein_totals)) {
    stop("site_conc, kd_mat and protein_totals must not contain NA.")
  }
  if (any(site_conc < 0)) stop("site_conc must be non-negative.")
  if (any(protein_totals < 0)) stop("protein_totals must be non-negative.")
  if (any(kd_mat <= 0)) stop("Kd values must be positive.")

  inv_kd <- 1 / kd_mat
  free <- as.numeric(protein_totals)
  q <- sweep(inv_kd, 2, free, `*`)
  active <- which(free > 0)
  converged <- length(active) == 0L || all(site_conc == 0)
  iter <- 0L

  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    max_change <- 0
    for (j in active) {
      others <- rowSums(q[, -j, drop = FALSE])
      a_j <- inv_kd[, j]
      p_j <- protein_totals[j]
      # h(u) = free + bound - total, increasing in u = log(free)
      h <- function(u) {
        f_j <- exp(u)
        x <- a_j * f_j
        f_j + sum(site_conc * x / (1 + others + x)) - p_j
      }
      hi <- log(p_j)
      lo <- max(hi - log(1e18), log(.Machine$double.xmin))
      new_f <- if (h(hi) <= 0) {
        # No binding to speak of: all protein stays free
        p_j
      } else if (h(lo) >= 0) {
        exp(lo)
      } else {
        exp(stats::uniroot(h, c(lo, hi), tol = tol * 1e-2)$root)
      }
      max_change <- max(max_change, abs(log(new_f) - log(free[j])))
      free[j] <- new_f
      q[, j] <- a_j * new_f
    }
    if (max_change < tol) converged <- TRUE
  }
  if (!converged) {
    warning("Shared protein pool did not converge in ", max_iter,
            " iterations.")
  }

  theta <- q / (1 + rowSums(q))
  bound <- colSums(site_conc * theta)
  names(free) <- colnames(kd_mat)
  names(bound) <- colnames(kd_mat)
  list(free = free, theta = theta, bound = bound,
       iterations = iter, converged = converged)
}
