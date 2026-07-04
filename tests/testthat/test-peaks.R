# Test Peak Calling Functions
# Tests for callPeaks

test_that("callPeaks returns correct structure and values", {
  # Mock simulation results
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = c(0.1, 0.1, 0.9, 0.9, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1),
    RBP2 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    transcript = "TX1"
  )
  
  # 1. Global method
  peaks_global <- callPeaks(mock_result, rbp_primary = "RBP1", method = "global", threshold = 1.5, min_width = 2)
  expect_s3_class(peaks_global, "data.table")
  expect_equal(nrow(peaks_global), 1)
  expect_equal(peaks_global$start[1], 3)
  expect_equal(peaks_global$end[1], 5)
  expect_equal(peaks_global$length[1], 3)
  
  # 2. Competitive method
  peaks_comp <- callPeaks(mock_result, rbp_primary = "RBP1", rbp_reference = "RBP2", method = "competitive", threshold = 5.0, min_width = 2)
  expect_s3_class(peaks_comp, "data.table")
  expect_equal(nrow(peaks_comp), 1)
  expect_equal(peaks_comp$start[1], 3)
  
  # 3. Local method
  peaks_local <- callPeaks(mock_result, rbp_primary = "RBP1", method = "local", threshold = 1.5, min_width = 2, window_size = 4)
  expect_s3_class(peaks_local, "data.table")
})

test_that("callPeaks handles no peaks gracefully", {
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
    transcript = "TX1"
  )
  
  # Should run with no warnings and return empty data.table
  expect_warning(
    peaks <- callPeaks(mock_result, rbp_primary = "RBP1", method = "global", threshold = 5.0),
    regexp = NA  # Assert no warning is thrown
  )
  
  expect_s3_class(peaks, "data.table")
  expect_equal(nrow(peaks), 0)
  expect_true("start" %in% names(peaks))
  expect_true("transcript" %in% names(peaks))
})

test_that("callPeaks parses genomic coordinates from headers", {
  # Positive strand
  mock_result_pos <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.1, 0.9, 0.9, 0.1),
    transcript = "chr1:1000-2000(+)"
  )
  
  peaks_pos <- callPeaks(mock_result_pos, rbp_primary = "RBP1", method = "global", threshold = 1.5, min_width = 2)
  
  expect_true("chrom" %in% names(peaks_pos))
  expect_true("genomic_start" %in% names(peaks_pos))
  expect_true("genomic_end" %in% names(peaks_pos))
  expect_true("strand" %in% names(peaks_pos))
  
  expect_equal(peaks_pos$chrom[1], "chr1")
  expect_equal(peaks_pos$strand[1], "+")
  # pos 3 to 4. Genomic start = 1000 + 3 - 1 = 1002, end = 1000 + 4 - 1 = 1003
  expect_equal(peaks_pos$genomic_start[1], 1002)
  expect_equal(peaks_pos$genomic_end[1], 1003)
  
  # Negative strand
  mock_result_neg <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.1, 0.9, 0.9, 0.1),
    transcript = "chr1:1000-2000(-)"
  )
  
  peaks_neg <- callPeaks(mock_result_neg, rbp_primary = "RBP1", method = "global", threshold = 1.5, min_width = 2)
  
  expect_equal(peaks_neg$chrom[1], "chr1")
  expect_equal(peaks_neg$strand[1], "-")
  # pos 3 to 4. Negative strand:
  # Genomic start = 2000 - 4 + 1 = 1997
  # Genomic end = 2000 - 3 + 1 = 1998
  expect_equal(peaks_neg$genomic_start[1], 1997)
  expect_equal(peaks_neg$genomic_end[1], 1998)
})
