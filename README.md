<p align="center">
  <img src="man/figures/logo.png" alt="RBPEqBind" height="200" />
</p>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/S00NYI/RBPEqBind)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.6.0-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**RBPEqBind** simulates competitive binding of multiple RNA-binding proteins (RBPs) to RNA sequences with an equilibrium binding model. All binding sites take protein from one pool of each RBP.

## Installation

### From Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RBPEqBind")
```

### From GitHub

<!-- Change "S00NYI/RBPEqBind@dev-v2" to "S00NYI/RBPEqBind" when dev-v2 is merged into main. -->
```r
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")

pak::pak("S00NYI/RBPEqBind@dev-v2")
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
- **FASTA file support** (`simulateBindingFasta()`, `simulateGridFasta()`) for multi-sequence analysis, with DNA or RNA input (`seq_type`)
- **Footprint** (`footprint`): a bound RBP covers a set number of nucleotides, and bound proteins may not overlap
- **One RBP pool for several RNAs** (`shared_protein = TRUE`)
- **RNA accessibility** (`accessibility`, `readAccessibility()`) from RNAplfold unpaired probabilities
- **Visualization** with binding profiles and heatmaps
- **Genomic coordinate parsing** from FASTA headers for BED export and peak-calling

## Documentation

See the [package vignette](vignettes/RBPEqBind-intro.Rmd) for detailed examples covering:

1. One RNA
2. One RNA with a concentration grid
3. Several RNAs (with the option of one RBP pool for all RNAs)
4. Several RNAs with a concentration grid

For detailed usage, see Figure_Scripts folder in: https://github.com/S00NYI/BITS_Specificity (uses an earlier version, RBPEqBind 0.99.0)

## Citation

If you use RBPEqBind in your research, please cite:

> Yi S, Singh SS, Ye X, Krishna R, Kothwela V, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v4

\* The simulation results in this manuscript were made with RBPEqBind 0.99.0. Later versions give numerically different results, but with the default settings the conclusions of those simulations do not change.

## AI Usage Disclosure

This package's Bioconductor review remediation, API refactoring, and test suite expansion were developed in collaboration with Claude Code.
For the original R and python implementation of the simulation model, see the [BITS_Specificity repository](https://github.com/S00NYI/BITS_Specificity).
