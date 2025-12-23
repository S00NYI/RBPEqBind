# RBPBind News

## Version 0.99.0

### New Features
- Initial Bioconductor submission
- **Competitive binding simulation** for N RBPs using equilibrium model
- **Concentration grid sweeps** for parameter exploration
- **FASTA file support** for transcriptome-wide analysis
- **SummarizedExperiment output** for Bioconductor integration

### Visualization
- `plotBinding()` - Position-wise binding profiles
- `plotHeatmap()` - Multi-RBP heatmaps
- `plotGrid()` - Concentration grid visualizations
- `plotBubble()` - Split-dot bubble plots for competitive analysis
- `viewModel()` - Affinity distribution histograms

### Export Functions
- JSON/CSV export with `exportResults()`
- BED format export with `exportBed()`
- SummarizedExperiment creation with `makeSE()`

### Peak Calling
- `callPeaks()` with competitive, local, and global methods
