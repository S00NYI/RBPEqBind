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

test_that("importResults imports JSON and CSV files correctly", {
  mock_result <- data.table::data.table(
    pos = 1:3,
    nt = c("A", "C", "G"),
    RBP1 = c(0.1, 0.2, 0.3),
    transcript = "TX1"
  )
  
  # Test JSON Import
  tmp_json <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_json), add = TRUE)
  exportResults(mock_result, tmp_json, format = "json")
  imported_json <- importResults(tmp_json)
  expect_s3_class(imported_json, "data.table")
  expect_equal(imported_json$pos, mock_result$pos)
  expect_equal(imported_json$RBP1, mock_result$RBP1)
  
  # Test CSV Import
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)
  exportResults(mock_result, tmp_csv, format = "csv")
  imported_csv <- importResults(tmp_csv)
  expect_s3_class(imported_csv, "data.table")
  expect_equal(imported_csv$pos, mock_result$pos)
  expect_equal(imported_csv$RBP1, mock_result$RBP1)
})

test_that("exportBed handles genomic coordinates correctly", {
  # Mock genomic sequence header
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
    transcript = "chr1:1000-2000(+)"
  )
  
  tmp_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp_bed), add = TRUE)
  
  # Export with threshold 0.25 (should keep pos 3, 4, 5 -> start=2, end=5 relative)
  # For genomic chr1:1000-2000(+):
  # BED start = 1000 + 3 - 2 = 1001
  # BED end = 1000 + 5 - 1 = 1004
  exportBed(mock_result, tmp_bed, rbp = "RBP1", threshold = 0.25)
  
  expect_true(file.exists(tmp_bed))
  bed_content <- data.table::fread(tmp_bed, header = FALSE)
  expect_equal(nrow(bed_content), 1)
  expect_equal(bed_content$V1[1], "chr1")
  expect_equal(bed_content$V2[1], 1001)
  expect_equal(bed_content$V3[1], 1004)
  
  # Test with negative strand: chr1:1000-2000(-)
  mock_result_neg <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
    transcript = "chr1:1000-2000(-)"
  )
  tmp_bed_neg <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp_bed_neg), add = TRUE)
  
  # Export with threshold 0.25 (keeps pos 3, 4, 5)
  # For negative strand:
  # BED start = header_end - pos_end = 2000 - 5 = 1995
  # BED end = header_end - pos_start + 1 = 2000 - 3 + 1 = 1998
  exportBed(mock_result_neg, tmp_bed_neg, rbp = "RBP1", threshold = 0.25)
  
  expect_true(file.exists(tmp_bed_neg))
  bed_content_neg <- data.table::fread(tmp_bed_neg, header = FALSE)
  expect_equal(nrow(bed_content_neg), 1)
  expect_equal(bed_content_neg$V1[1], "chr1")
  expect_equal(bed_content_neg$V2[1], 1995)
  expect_equal(bed_content_neg$V3[1], 1998)
})
