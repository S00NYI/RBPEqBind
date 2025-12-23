#' Equilibrium Binding Solver
#'
#' Solves for the free RNA concentration in a competitive binding system.
#'
#' @param RNA_total Total RNA concentration (numeric).
#' @param protein_totals Vector of total protein concentrations for N RBPs (numeric vector).
#' @param Kds Vector of dissociation constants (Kd) for N RBPs (numeric vector).
#' @return The equilibrium free RNA concentration.
#' @examples
#' # Two RBPs competing for RNA
#' free_rna <- solveEquilibrium(
#'   RNA_total = 10,
#'   protein_totals = c(100, 100),
#'   Kds = c(10, 50)
#' )
#' free_rna
#' @export
solveEquilibrium <- function(RNA_total, protein_totals, Kds) {
  # Input validation
  if (length(protein_totals) != length(Kds)) {
    stop("Length of protein_totals and Kds must match.")
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
    bound_sum <- sum((protein_totals * X) / (Kds + X))
    return(X - RNA_total + bound_sum)
  }
  
  # Root finding
  # Lower bound: 0
  # Upper bound: RNA_total (Free RNA cannot exceed Total RNA)
  # Check boundaries
  f_lower <- obj_func(0) # -RNA_total < 0
  f_upper <- obj_func(RNA_total) # sum(...) > 0, so f_upper > 0
  
  if (f_lower == 0) return(0)
  if (f_upper == 0) return(RNA_total)
  
  # Use uniroot
  result <- stats::uniroot(obj_func, interval = c(0, RNA_total), tol = .Machine$double.eps^0.25)
  return(result$root)
}

#' Calculate Bound Concentrations
#'
#' Calculates the equilibrium bound concentrations for each RBP.
#'
#' @param free_RNA Equilibrium free RNA concentration.
#' @param protein_totals Vector of total protein concentrations.
#' @param Kds Vector of dissociation constants.
#' @return Vector of bound concentrations for each RBP.
#' @examples
#' free_rna <- solveEquilibrium(10, c(100, 100), c(10, 50))
#' calcBoundProtein(free_rna, c(100, 100), c(10, 50))
#' @export
calcBoundProtein <- function(free_RNA, protein_totals, Kds) {
  return((protein_totals * free_RNA) / (Kds + free_RNA))
}
