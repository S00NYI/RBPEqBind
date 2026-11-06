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

test_that("simulateBinding errors when k does not match model motif length", {
  models <- list(
    RBP1 = data.table::data.table(
      motif = c("AAAAA", "UUUUU"),
      Kd = c(10, 100)
    )
  )
  expect_error(simulateBinding("AAAUUU", models, c(RBP1 = 100), k = 4), "uses 5-mers but k = 4")
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

test_that("every copy of a motif gets the same occupancy as a single-copy motif of equal Kd", {
  # AAAAA appears many times, CCCCC once, both with the same Kd. With one
  # shared protein supply, both see the same free protein, so their windows
  # must have the same occupancy, however many copies exist.
  # All other 5-mers share one weak Kd, so a nucleotide's occupancy depends
  # only on how its covering windows bind. The protein supply (50) is below
  # the AAAAA sites (20 copies x rna_conc 5 = 100), so the pool matters here.
  b <- c("A", "C", "G", "U")
  all5 <- do.call(paste0, expand.grid(b, b, b, b, b, stringsAsFactors = FALSE))
  kd <- rep(1e6, length(all5))
  kd[all5 %in% c("AAAAA", "CCCCC")] <- 10
  models <- list(RBP1 = data.table::data.table(motif = all5, Kd = kd))
  seq <- paste0("G", strrep("AAAAAG", 20), "CCCCCGGGGG")
  res <- simulateBinding(seq, models, c(RBP1 = 50), rna_conc = 5, k = 5)
  # Compare the middle nucleotide of an interior AAAAA block and of CCCCC.
  # Each is covered by five windows, of which one is the motif itself.
  a_mid <- 1 + 6 * 10 + 3
  c_mid <- regexpr("CCCCC", seq)[1] + 2
  expect_equal(substr(seq, a_mid - 2, a_mid + 2), "AAAAA")
  expect_equal(res$RBP1[a_mid], res$RBP1[c_mid], tolerance = 1e-10)
  expect_gt(res$RBP1[c_mid], 0.01)
})

test_that("simulateBinding occupancy falls as RNA rises at fixed protein", {
  models <- list(RBP1 = data.table::data.table(motif = c("AA", "AU", "UA", "UU"),
                                               Kd = c(10, 50, 50, 100)))
  occ <- vapply(c(1, 10, 100, 1000), function(r) {
    mean(simulateBinding("AAUUAAUU", models, c(RBP1 = 100), rna_conc = r, k = 2)$RBP1)
  }, numeric(1))
  expect_true(all(diff(occ) < 0))
})

test_that("simulateBinding at protein excess matches the independent-site formula", {
  models <- list(RBP1 = data.table::data.table(motif = c("AA", "AU", "UA", "UU"),
                                               Kd = c(10, 50, 50, 100)))
  seq <- "AAUU"
  res <- simulateBinding(seq, models, c(RBP1 = 100), rna_conc = 1e-9, k = 2)
  theta <- function(kd) (100 / kd) / (1 + 100 / kd)
  w <- theta(c(10, 50, 100))            # windows AA, AU, UU
  expected <- c(w[1], mean(w[1:2]), mean(w[2:3]), w[3])
  expect_equal(res$RBP1, expected, tolerance = 1e-8)
})

test_that("simulateBinding rejects negative concentrations", {
  models <- list(RBP1 = data.table::data.table(motif = "AA", Kd = 10))
  expect_error(simulateBinding("AAA", models, c(RBP1 = -1), k = 2), "non-negative")
  expect_error(simulateBinding("AAA", models, c(RBP1 = 1), rna_conc = -1, k = 2), "non-negative")
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
    k = 2
  )
  
  expect_s3_class(result, "data.table")
  expect_true("Conc_RBP1" %in% names(result))
  expect_true("rna_conc" %in% names(result))
  
  # Should have 2 concentration combos * 6 positions = 12 rows
  expect_equal(nrow(result), 12)
})

test_that("simulateBindingFasta processes FASTA files", {
  fasta_file <- system.file("extdata", "test_transcripts.fa", package = "RBPEqBind")
  skip_if(fasta_file == "")
  
  models <- list(
    RBP1 = data.table::data.table(motif = c("AAAAA", "UUUUU"), Kd = c(10, 100))
  )
  
  # Test single file
  result <- simulateBindingFasta(
    fasta_file,
    models,
    c(RBP1 = 100),
    rna_conc = 10,
    k = 5
  )
  
  expect_s3_class(result, "data.table")
  expect_true("transcript" %in% names(result))
  
  # Test multiple files
  result_multi <- simulateBindingFasta(
    c(fasta_file, fasta_file),
    models,
    c(RBP1 = 100),
    rna_conc = 10,
    k = 5
  )
  expect_s3_class(result_multi, "data.table")
  expect_equal(nrow(result_multi), 2 * nrow(result))
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
