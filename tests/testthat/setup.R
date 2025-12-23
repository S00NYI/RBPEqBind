# Setup file for testthat
# Loads the package before running tests

library(data.table)
library(ggplot2)

# Load all package functions
devtools::load_all()
