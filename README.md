# RBPBind <img src="man/figures/logo.png" align="right" height="150" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/S00NYI/RBPBind)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**RBPBind** simulates competitive binding of multiple RNA-binding proteins (RBPs) to RNA sequences using equilibrium binding kinetics.

## Installation

### From GitHub (Development)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("S00NYI/RBPBind")
```

### From Bioconductor (After Acceptance)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RBPBind")
```

## Quick Start

```r
library(RBPBind)

# Load and process RBP models
raw_models <- loadModel("model.csv")
rbp_models <- setModel(raw_models, max_affinity = 100)

# Run simulation
results <- simulateBinding(
  sequence = "ACGUACGUACGU...",
  rbp_models = rbp_models,
  protein_concs = c(RBP1 = 100, RBP2 = 100),
  rna_conc = 10
)

# Visualize
plotBinding(results, rbp = c("RBP1", "RBP2"))
```

## Features

- **Competitive binding simulation** for multiple RBPs
- **Concentration grid sweeps** for parameter exploration
- **FASTA file support** for transcriptome-wide analysis
- **Rich visualization** with binding profiles, heatmaps, and bubble plots
- **SummarizedExperiment output** for Bioconductor integration

## Documentation

See the [package vignette](vignettes/RBPBind-intro.Rmd) for detailed examples covering:

1. Single sequence simulation
2. Concentration grid sweeps
3. FASTA file processing
4. Visualization options

## Citation

If you use RBPBind in your research, please cite:

> Yi S, Singh SS, Ye X, Krishna R, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2

## Acknowledgments

This package was developed in part with **Google Antigravity**.

For the original Python implementation of RNA-RBP interaction simulations, see:
https://github.com/S00NYI/Specificity_BITs
