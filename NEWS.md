# RBPEqBind News

## Development version (branch dev-v2, not released)

### Shared protein pool (changes numerical output)
- `simulateBinding()` now solves one mass balance for each RBP across all
  binding sites of the sequence. All sites draw on one shared supply of each
  protein, so total bound protein can no longer exceed total protein.
- The copy-number scaling of the RNA concentration (`rna_conc * count` per
  k-mer) is removed. Every copy of a k-mer now has the same binding
  probability. The copy number enters only as the amount of protein that the
  copies take from the supply.
- When protein is in large excess over RNA, results agree with 0.99.4. When
  protein is limiting, occupancies are lower than in 0.99.4.
- `simulateBinding()` now checks that `protein_concs` and `rna_conc` are
  non-negative and that `sequence` is a single string.
- The internal solvers `solveEquilibrium()` and `calcBoundProtein()` are
  replaced by the internal `solveSharedPool()`.

### Footprint lattice (new argument)
- `simulateBinding()`, `simulateGrid()`, `simulateBindingFasta()` and
  `simulateGridFasta()` gain `footprint = NULL`. With `NULL`, windows bind
  independently, as before. With a number (or a named vector, one value per
  RBP), a bound RBP covers that many nucleotides and bound proteins may not
  overlap. The model sums over all allowed arrangements on a molecule with a
  forward and a backward pass, and keeps the shared protein balance.
- With a footprint, the per-position value is the probability that the
  nucleotide is covered by the RBP. Without one, it stays the mean over the
  windows that cover the nucleotide.
- The k-mer in the middle of the footprint sets the affinity of a placement.
- The error message in `simulateGrid()` now names `simulateGridFasta()`.

### Neighbor factor (new argument)
- The same four functions gain `omega = 1`. With a footprint, two proteins of
  the same RBP whose footprints touch multiply the weight of their arrangement
  by `omega`. Values above 1 make neighbors help each other, values below 1
  make them hinder each other, and 1 has no effect. Values other than 1 need a
  footprint.
- Tests compare the lattice with brute-force enumeration and with the
  cooperative McGhee-von Hippel formula.

### Documentation
- The help page of `simulateBinding()` now describes every output column.
  It states that `_density_fc` and `_occupancy_fc` mean a fold change over a
  uniform track with one active RBP, and a ratio to the other active RBPs
  with two or more. The help page and the vignette recommend the occupancy
  columns for quantitative analysis, and two runs (with and without a
  competitor) to measure a competitor's effect. The output itself does not
  change.

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
