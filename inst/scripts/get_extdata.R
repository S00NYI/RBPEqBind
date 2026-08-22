# ==============================================================================
# Origin and Description of Raw Data in inst/extdata/
# ==============================================================================
#
# This script documents the source and generation of example datasets included 
# in the RBPEqBind package, per Bioconductor guidelines:
# https://contributions.bioconductor.org/data.html#raw-data-and-the-instextdata-directory
#
# 1. model_RBP.csv
# ----------------
# Description:
#   A synthetic 5-mer RBP affinity table containing relative binding scores 
#   for 4 model RBPs:
#     - HH: High specificity, High affinity (strongly prefers poly-U/poly-G)
#     - HL: High specificity, Low affinity
#     - LH: Low specificity, High affinity
#     - LL: Low specificity, Low affinity (non-specific baseline)
#   Total rows: 1024 (all possible 4^5 RNA 5-mers).
#
# Generation:
#   Generated synthetically using the model generation routines described in:
#   Yi S, et al. (2025). Inherent Specificity and Mutational Sensitivity as 
#   Quantitative Metrics for RBP Binding. bioRxiv 10.1101/2025.03.28.646018.
#
# 2. test_transcripts.fa
# ----------------------
# Description:
#   A FASTA file containing 4 synthetic transcript sequences (TX1 to TX4)
#   with embedded high-affinity motifs and UCSC-style coordinate headers:
#     - >chr1:1000-2000(+) TX1 (synthetic 1000 nt transcript)
#     - >chr2:5000-6000(-) TX2 (synthetic 1000 nt transcript, negative strand)
#     - >chr3:100-300(+) TX3 (synthetic 200 nt transcript)
#     - >chr4:10-50(+) TX4 (synthetic 40 nt transcript)
#
# Purpose:
#   Used in automated testing and package examples to verify multi-sequence 
#   processing, FASTA parsing, BED coordinate translation, and export functions.
# ==============================================================================
