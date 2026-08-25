# RBPEqBind News

## Version 0.99.4

### Enhancements & Refactoring
- **Idiomatic Parallel Architecture**: Refactored `simulateGrid()`, `simulateBindingFasta()`, and `simulateGridFasta()` to adhere strictly to the `future` framework by delegating concurrency to the user-configured `future::plan()`, running sequentially by default in compliance with Bioconductor guidelines.
- **Clean Imports**: Cleaned up `jsonlite` roxygen imports in `R/export.R`, removing unused `stream_out` and `toJSON`.
- **Deterministic Visual Testing**: Implemented deterministic test fixtures and visual snapshot testing with `vdiffr` across all plotting routines.

## Version 0.99.3

### Reviewer Remediation & Numerical Hardening
- Harmonized JSON export across aggregated and unaggregated formats.
- Integrated `rtracklayer::export.bed` for BED file export.
- Replaced linear solver with scale-invariant log-space root-finder in `solveEquilibrium()`.
- Fixed whole-transcript edge background calculation in `callPeaks(method = "local")`.
- Documented raw dataset provenance in `inst/scripts/get_extdata.R`.

## Version 0.99.0

### Initial Submission Features
- **Competitive binding simulation** for N RBPs on RNA using equilibrium model
- **Concentration grid sweeps** for parameter exploration
- **FASTA file support** for transcriptome-wide binding analysis
- **Peak calling** with `callPeaks()`
- **Visualization** suite: `viewModel()`, `plotBinding()`, `plotHeatmap()`, `plotGrid()`
- **Export tools**: `exportResults()`, `exportBed()`, `makeSE()`
