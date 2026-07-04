<p align="center">
  <img src="man/figures/logo.png" alt="RBPEqBind" height="200" />
</p>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/S00NYI/RBPEqBind)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**RBPEqBind** simulates competitive binding of multiple RNA-binding proteins (RBPs) to RNA sequences using equilibrium binding kinetics.

## Installation

### From GitHub (Development)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("S00NYI/RBPEqBind")
```

## Quick Start

```r
library(RBPEqBind)

# Load and process sample RBP models included in the package
model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
raw_models <- loadModel(model_file)
rbp_models <- setModel(raw_models, max_affinity = 100)

# Run simulation
results <- simulateBinding(
  sequence = "ACGUACGUACGUUUUUACGUACGU",
  rbp_models = rbp_models,
  protein_concs = c(HH = 100, HL = 100),
  rna_conc = 10
)

# Visualize
results$transcript <- "QuickStart"
plotBinding(results, rbp = c("HH", "HL"))
```

## Features

- **Competitive binding simulation** for multiple RBPs
- **Concentration grid sweeps** for parameter exploration
- **FASTA file support** (`simulateBindingFasta()`, `simulateGridFasta()`) for multi-sequence analysis
- **Visualization** with binding profiles and heatmaps
- **Genomic coordinate parsing** from FASTA headers for BED export and peak-calling

## Documentation

See the [package vignette](vignettes/RBPEqBind-intro.Rmd) for detailed examples covering:

1. Single sequence simulation
2. Concentration grid sweeps
3. FASTA file processing
4. Visualization options

For detailed usage, see Figure_Scripts folder in: https://github.com/S00NYI/BITS_Specificity

## Citation

If you use RBPEqBind in your research, please cite:

> Yi S, Singh SS, Ye X, Krishna R, Kothwela V, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v4

## Work in Progress

Planned developments for `RBPEqBind` include:
- **CLIP-Seq Concordance**: Directly comparing simulated occupancy with experimental CLIP-seq peaks or crosslink signal tracks.
- **RNA Structure Integration**: Adding accessibility masking using SHAPE/DMS experimental reactivity profiles or pre-computed structures.

## AI Usage Disclosure

This package's Bioconductor review remediation, API refactoring, and test suite expansion were developed in collaboration with Google Antigravity.
For the original R and python implementation of the simulation model, see the [BITS_Specificity repository](https://github.com/S00NYI/BITS_Specificity).
