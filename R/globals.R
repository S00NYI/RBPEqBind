# Global variables for data.table NSE
# This prevents R CMD check NOTEs about "no visible binding for global variable"

utils::globalVariables(c(
  # data.table special symbols
  ".", ".N", ".SD", ":=",
  
  # Column names used in data.table expressions
  "..cols", "..row_cols",
  "bg", "chrom", "cluster", "end", "first", "g", "group",
  "height", "is_peak", "Ka", "Kd", "label", "last", "motif",
  "name", "nt", "pos", "RBP", "rleid", "rna_conc", "score",
  "start", "value", "Value", "width", "x", "y",
  
  # data.table functions used without namespace
  "copy", "setkeyv", "setnames"
))
