# Test Simulation Functions
# Tests for simulateBinding, simulateGrid, and FASTA variants

test_that("simulateBinding returns correct structure", {
  # Create mock models
  models <- list(
    RBP1 = data.table::data.table(
      motif = c("AA", "UU", "CC", "GG"),
      Kd = c(10, 100, 50, 50)
    )
  )
  
  seq <- "AAUUCCGG"
  result <- simulateBinding(seq, models, c(RBP1 = 100), rna_conc = 10, k = 2)
  
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), nchar(seq))
  expect_true("pos" %in% names(result))
  expect_true("nt" %in% names(result))
  expect_true("RBP1" %in% names(result))
  expect_true("RBP1_density" %in% names(result))
  expect_true("RBP1_density_fc" %in% names(result))
  expect_true("RBP1_occupancy_fc" %in% names(result))
})

test_that("simulateBinding occupancy is between 0 and 1", {
  models <- list(
    RBP1 = data.table::data.table(
      motif = c("AA", "UU"),
      Kd = c(10, 100)
    )
  )
  
  result <- simulateBinding("AAAUUUAAA", models, c(RBP1 = 100), rna_conc = 10, k = 2)
  
  expect_true(all(result$RBP1 >= 0))
  expect_true(all(result$RBP1 <= 1))
})

test_that("simulateBinding handles multiple RBPs", {
  models <- list(
    HH = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100)),
    LL = data.table::data.table(motif = c("AA", "UU"), Kd = c(100, 10))
  )
  
  result <- simulateBinding("AAAUUU", models, c(HH = 100, LL = 100), rna_conc = 10, k = 2)
  
  expect_true("HH" %in% names(result))
  expect_true("LL" %in% names(result))
  expect_true("HH_density_fc" %in% names(result))
  expect_true("LL_occupancy_fc" %in% names(result))
})

test_that("simulateBinding high affinity has higher occupancy", {
  models <- list(
    HighAff = data.table::data.table(motif = "AA", Kd = 1),
    LowAff = data.table::data.table(motif = "AA", Kd = 1000)
  )
  
  result <- simulateBinding("AAAA", models, c(HighAff = 100, LowAff = 100), rna_conc = 10, k = 2)
  
  expect_true(mean(result$HighAff) > mean(result$LowAff))
})

test_that("simulateBinding errors on invalid sequence type", {
  models <- list(RBP1 = data.table::data.table(motif = "AA", Kd = 10))
  
  expect_error(simulateBinding(c("AA", "UU"), models, c(RBP1 = 100), rna_conc = 10, k = 2))
})

test_that("simulateGrid returns correct structure", {
  models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100))
  )
  
  seq <- "AAAUUU"
  result <- simulateGrid(
    seq, 
    models, 
    protein_conc_grid = list(RBP1 = c(10, 100)),
    rna_conc_grid = c(10),
    k = 2,
    parallel = FALSE
  )
  
  expect_s3_class(result, "data.table")
  expect_true("Conc_RBP1" %in% names(result))
  expect_true("rna_conc" %in% names(result))
  
  # Should have 2 concentration combos * 6 positions = 12 rows
  expect_equal(nrow(result), 12)
})

test_that("simulateBindingF processes FASTA files", {
  skip_if_not(file.exists("../../data/test_transcripts.fa"))
  
  models <- list(
    RBP1 = data.table::data.table(motif = c("AAAAA", "UUUUU"), Kd = c(10, 100))
  )
  
  result <- simulateBindingF(
    "../../data/test_transcripts.fa",
    models,
    c(RBP1 = 100),
    rna_conc = 10,
    k = 5,
    parallel = FALSE
  )
  
  expect_s3_class(result, "data.table")
  expect_true("transcript" %in% names(result))
})

test_that("density_fc is density * L for single RBP", {
  models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100))
  )
  
  seq <- "AAUUAA"  # Length 6
  result <- simulateBinding(seq, models, c(RBP1 = 100), rna_conc = 10, k = 2)
  
  L <- nrow(result)
  # density_fc should be density * L for single RBP
  expected_fc <- result$RBP1_density * L
  
  expect_equal(result$RBP1_density_fc, expected_fc, tolerance = 1e-10)
})

test_that("density_fc uses sum(others) for multiple active RBPs", {
  models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100)),
    RBP2 = data.table::data.table(motif = c("AA", "UU"), Kd = c(100, 10))
  )
  
  # Both RBPs are active (conc > 0)
  result <- simulateBinding("AAUUAA", models, c(RBP1 = 100, RBP2 = 50), rna_conc = 10, k = 2)
  
  # For RBP1, fc = RBP1_density / RBP2_density
  expected_fc <- result$RBP1_density / result$RBP2_density
  expected_fc[is.infinite(expected_fc)] <- NA_real_
  
  expect_equal(result$RBP1_density_fc, expected_fc, tolerance = 1e-10)
})
