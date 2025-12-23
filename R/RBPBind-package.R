#' RBPBind: Simulate Competitive RBP Binding to RNA
#'
#' RBPBind simulates competitive binding of multiple RNA-binding proteins (RBPs)
#' to RNA sequences using an equilibrium binding model. The package provides tools
#' for loading RBP affinity models, running binding simulations, and visualizing
#' competition patterns.
#'
#' @section Core Functions:
#' \describe{
#'   \item{Model Loading}{\code{\link{loadModel}}, \code{\link{setModel}}, \code{\link{viewModel}}}
#'   \item{Simulation}{\code{\link{simulateBinding}}, \code{\link{simulateGrid}}, 
#'         \code{\link{simulateBindingF}}, \code{\link{simulateGridF}}}
#'   \item{Visualization}{\code{\link{plotBinding}}, \code{\link{plotHeatmap}}, 
#'         \code{\link{plotGrid}}, \code{\link{plotBubble}}}
#'   \item{Export}{\code{\link{exportResults}}, \code{\link{makeSE}}}
#' }
#'
#' @section Getting Started:
#' See the package vignette for a comprehensive introduction:
#' \code{vignette("RBPBind-intro", package = "RBPBind")}
#'
#' @references
#' Yi S, Singh SS, Ye X, Krishna R, Jankowsky E, Luna JM. (2025).
#' Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.
#' bioRxiv. \url{https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2}
#'
#' @section Acknowledgments:
#' This package was developed in part with Google Antigravity.
#' For the original Python implementation of RNA-RBP interaction simulations,
#' see \url{https://github.com/S00NYI/Specificity_BITs}.
#'
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom data.table copy setkeyv setnames rleid first last
#' @importFrom pracma fsolve
#' 
#' @docType package
#' @name RBPBind
#' @keywords internal
"_PACKAGE"

