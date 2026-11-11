# RBPEqBind News

## Version 0.99.5

### One protein pool per RBP (changes numerical output when protein is limiting)
- `simulateBinding()` now solves one mass balance for each RBP across all
  binding windows of the sequence: total protein = free protein + bound
  protein. All windows of a sequence take protein from this one pool, so total bound
  protein can no longer exceed total protein. In 0.99.4, each k-mer was a
  site at `rna_conc` times its copy number, and the protein pool was not
  shared between different k-mers.
- Every copy of a k-mer now has the same binding probability. The copy number
  enters only as the amount of protein that the copies take from the pool.
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
- Tests compare the lattice with brute-force enumeration and with the
  McGhee-von Hippel formula.
- The error message in `simulateGrid()` now names `simulateGridFasta()`.

### One protein pool for several RNAs (new option)
- `simulateBindingFasta()` and `simulateGridFasta()` gain
  `shared_protein = FALSE`. Use it to model one protein pool that binds
  several different RNAs, for example protein titrated into a mix of
  transcripts. With `TRUE`, all sequences are in one tube and take
  protein from one pool of each RBP. The model solves one mass balance over all
  sequences, and each sequence keeps its own windows or lattice, so no
  placement crosses from one sequence to the next. With `FALSE`, each
  sequence is simulated alone, as before, and the output does not change.
- `simulateBindingFasta()` now accepts one `rna_conc` value per sequence,
  named with the sequence names or in file order. In `simulateGridFasta()`, each
  value of `rna_conc_grid` applies to every sequence.
- Tests compare the many-molecule lattice balance with enumeration, and check
  that two identical sequences at concentration c give the same result as one
  sequence alone at 2c.

### FASTA input (new argument)
- `simulateBindingFasta()` and `simulateGridFasta()` gain
  `seq_type = c("DNA", "RNA")`. With `"DNA"` (default, as before), T is
  changed to U after reading. With `"RNA"`, the file can hold U.
- A letter outside the chosen alphabet (for example U in a file read as DNA)
  now stops the function with an error. Before, Biostrings removed such
  letters with only a warning, which shortened the sequence.

### RNA accessibility (new argument and function)
- `simulateBinding()`, `simulateGrid()`, `simulateBindingFasta()` and
  `simulateGridFasta()` gain `accessibility = NULL`. It takes a matrix with
  one row per nucleotide, where column l is the probability that the l
  nucleotides ending there are all unpaired (the layout of an RNAplfold
  `-u` file). The FASTA functions take a list of such matrices, named with
  the sequence names. A protein binds only an unpaired k-mer window. The
  binding weight of each window (with a footprint: of each placement) is
  multiplied by a, the probability that the k-mer is unpaired. This assumes
  that the window switches between paired and unpaired and that a bound
  protein holds it unpaired. Kd does not change, and a applies to every
  RBP. With `NULL`, the output does not change.
- New `readAccessibility()` reads an RNAplfold `*_lunp` file. RBPEqBind does
  not fold RNA and does not depend on ViennaRNA.
- Tests check that accessibility of 1 gives the same result as `NULL`,
  compare single sites with the three-state formula, check the mass balance,
  and check the FASTA and shared-protein paths.

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
