# Test Export Functions
# Tests for exportResults, import_results, and makeSE

test_that("makeSE creates SummarizedExperiment with correct structure", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Create mock simulation result
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = runif(10, 0, 0.5),
    RBP2 = runif(10, 0, 0.5)
  )
  mock_result[, RBP1_density := RBP1 / sum(RBP1)]
  mock_result[, RBP2_density := RBP2 / sum(RBP2)]
  mock_result[, RBP1_density_fc := RBP1_density * 10]
  mock_result[, RBP2_density_fc := RBP2_density * 10]
  mock_result[, RBP1_occupancy_fc := RBP1 / mean(RBP1)]
  mock_result[, RBP2_occupancy_fc := RBP2 / mean(RBP2)]
  
  se <- makeSE(mock_result)
  
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(nrow(se), 10)
  expect_equal(ncol(se), 2)
  
  # Check assays
  assay_names <- names(SummarizedExperiment::assays(se))
  expect_true("occupancy" %in% assay_names)
  expect_true("density" %in% assay_names)
  expect_true("density_fc" %in% assay_names)
  expect_true("occupancy_fc" %in% assay_names)
})

test_that("makeSE rowData contains position info", {
  skip_if_not_installed("SummarizedExperiment")
  
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = runif(5)
  )
  mock_result[, RBP1_density := RBP1 / sum(RBP1)]
  mock_result[, RBP1_density_fc := RBP1_density * 5]
  mock_result[, RBP1_occupancy_fc := RBP1 / mean(RBP1)]
  
  se <- makeSE(mock_result)
  
  rd <- SummarizedExperiment::rowData(se)
  expect_true("pos" %in% names(rd))
  expect_true("nt" %in% names(rd))
  expect_equal(rd$pos, 1:5)
})

test_that("makeSE colData contains RBP info", {
  skip_if_not_installed("SummarizedExperiment")
  
  mock_result <- data.table::data.table(
    pos = 1:3,
    nt = c("A", "C", "G"),
    Alpha = runif(3),
    Beta = runif(3)
  )
  mock_result[, Alpha_density := Alpha / sum(Alpha)]
  mock_result[, Beta_density := Beta / sum(Beta)]
  mock_result[, Alpha_density_fc := Alpha_density * 3]
  mock_result[, Beta_density_fc := Beta_density * 3]
  mock_result[, Alpha_occupancy_fc := Alpha / mean(Alpha)]
  mock_result[, Beta_occupancy_fc := Beta / mean(Beta)]
  
  se <- makeSE(mock_result)
  
  cd <- SummarizedExperiment::colData(se)
  expect_true("RBP" %in% names(cd))
  expect_equal(cd$RBP, c("Alpha", "Beta"))
})

test_that("makeSE includes Kd ranges when models provided", {
  skip_if_not_installed("SummarizedExperiment")
  
  mock_result <- data.table::data.table(
    pos = 1:3,
    nt = c("A", "C", "G"),
    RBP1 = runif(3)
  )
  mock_result[, RBP1_density := RBP1 / sum(RBP1)]
  mock_result[, RBP1_density_fc := RBP1_density * 3]
  mock_result[, RBP1_occupancy_fc := RBP1 / mean(RBP1)]
  
  mock_models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(0.1, 100))
  )
  
  se <- makeSE(mock_result, rbp_models = mock_models)
  
  cd <- SummarizedExperiment::colData(se)
  expect_true("Kd_min" %in% names(cd))
  expect_true("Kd_max" %in% names(cd))
  expect_equal(as.numeric(cd$Kd_min[1]), 0.1)
  expect_equal(as.numeric(cd$Kd_max[1]), 100)
})

test_that("exportResults creates JSON file", {
  mock_result <- data.table::data.table(
    pos = 1:3,
    nt = c("A", "C", "G"),
    RBP1 = c(0.1, 0.2, 0.3),
    transcript = "TX1"
  )
  
  tmp_file <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_file))
  
  exportResults(mock_result, tmp_file, format = "json")
  
  expect_true(file.exists(tmp_file))
  expect_gt(file.size(tmp_file), 0)
})

test_that("exportResults creates CSV file", {
  mock_result <- data.table::data.table(
    pos = 1:3,
    nt = c("A", "C", "G"),
    RBP1 = c(0.1, 0.2, 0.3)
  )
  
  tmp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_file))
  
  exportResults(mock_result, tmp_file, format = "csv")
  
  expect_true(file.exists(tmp_file))
  
  # Read back and verify
  read_back <- data.table::fread(tmp_file)
  expect_equal(nrow(read_back), 3)
})
