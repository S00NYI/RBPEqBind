# Equilibrium Binding Solver
#
# Solves for the free RNA concentration in a competitive binding system.
#
# @param RNA_total Total RNA concentration (numeric).
# @param protein_totals Vector of total protein concentrations for N RBPs (numeric vector).
# @param Kds Vector of dissociation constants (Kd) for N RBPs (numeric vector).
# @return The equilibrium free RNA concentration.
solveEquilibrium <- function(RNA_total, protein_totals, Kds) {
  # Input validation
  if (length(protein_totals) != length(Kds)) {
    stop("Length of protein_totals and Kds must match.")
  }
  if (anyNA(protein_totals) || anyNA(Kds)) {
    stop("protein_totals and Kds must not contain NA")
  }
  if (RNA_total < 0) stop("RNA_total must be non-negative")
  if (any(protein_totals < 0)) stop("protein_totals must be non-negative")
  if (any(Kds <= 0)) stop("Kds must be positive")
  
  # If RNA is 0, free RNA is 0
  if (RNA_total == 0) return(0)
  
  # Define the objective function: f(FreeRNA) = 0
  # X = FreeRNA
  # Equation: X - RNA_total + sum( (P_tot * X) / (Kd + X) ) = 0
  obj_func <- function(X) {
    bound_sum <- sum(calcBoundProtein(X, protein_totals, Kds))
    return(X - RNA_total + bound_sum)
  }
  
  # f(0) = -RNA_total < 0 always. If f(RNA_total) <= 0 there is no interior
  # root (no protein, or all Kd infinite): all RNA is free.
  if (obj_func(RNA_total) <= 0) return(RNA_total)
  
  # Solve in log space for scale invariance across physiological concentration ranges
  hi <- log(RNA_total)
  lo <- max(hi - log(1e18), log(.Machine$double.xmin))
  if (obj_func(exp(lo)) > 0) return(exp(lo))
  
  result <- stats::uniroot(function(u) obj_func(exp(u)),
                           interval = c(lo, hi),
                           tol = .Machine$double.eps^0.5)
  return(exp(result$root))
}

# Calculate Bound Concentrations
#
# Calculates the equilibrium bound concentrations for each RBP.
#
# @param free_RNA Equilibrium free RNA concentration.
# @param protein_totals Vector of total protein concentrations.
# @param Kds Vector of dissociation constants.
# @return Vector of bound concentrations for each RBP.
calcBoundProtein <- function(free_RNA, protein_totals, Kds) {
  return((protein_totals * free_RNA) / (Kds + free_RNA))
}
