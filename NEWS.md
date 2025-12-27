# RBPEqBind News

## Version 0.99.0

### New Features
- Initial Bioconductor submission
- **Competitive binding simulation** for N RBPs on RNA using equilibrium model
- **Concentration grid sweeps** for parameter exploration
- **FASTA file support** for transcriptome-wide binding analysis

### Visualization
- `viewModel()` - Affinity distribution histograms
- `plotBinding()` - Position-wise binding profiles
- `plotHeatmap()` - Multi-RBP heatmaps
- `plotGrid()` - Concentration grid visualizations
- `plotBubble()` - Split-dot bubble plots for competitive analysis

### Export Functions
- JSON/CSV export with `exportResults()`
- BED format export with `exportBed()
- SummarizedExperiment creation with `makeSE()`

### Future Implementations
- `callPeaks()`
