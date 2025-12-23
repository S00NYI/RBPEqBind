# Test Utility Functions
# Tests for featureScale, generateRNA, extractKmers, countKmers

test_that("featureScale scales to correct range", {
  input <- c(0, 0.5, 1)
  result <- featureScale(input, MAX = 100, MIN = 0)
  
  expect_equal(min(result), 0)
  expect_equal(max(result), 100)
  expect_equal(result[2], 50)
})

test_that("featureScale handles constant input", {
  input <- c(5, 5, 5)
  result <- featureScale(input, MAX = 100, MIN = 0)
  
  # All same value - should return MAX (or handle gracefully)
  expect_true(all(is.finite(result)))
})

test_that("generateRNA returns correct length", {
  result <- generateRNA(100)
  
  expect_equal(nchar(result), 100)
})
  
test_that("generateRNA only contains valid nucleotides", {
  result <- generateRNA(1000)
  valid_nts <- c("A", "C", "G", "U")
  
  chars <- strsplit(result, "")[[1]]
  expect_true(all(chars %in% valid_nts))
})

test_that("generateRNA respects frequency parameters", {
  # Generate with only A
  result <- generateRNA(100, A = 1, G = 0, C = 0, U = 0)
  
  chars <- strsplit(result, "")[[1]]
  expect_true(all(chars == "A"))
})

test_that("extractKmers returns correct count", {
  seq <- "ACGUA"  # Length 5
  k <- 3
  # K-mers: ACG, CGU, GUA = 3 k-mers (length - k + 1)
  
  result <- extractKmers(seq, k = k)
  
  expect_length(result, 3)
  expect_equal(result[1], "ACG")
  expect_equal(result[2], "CGU")
  expect_equal(result[3], "GUA")
})

test_that("extractKmers handles edge cases", {
  # Sequence same length as k
  result <- extractKmers("ACGUA", k = 5)
  expect_length(result, 1)
  expect_equal(result, "ACGUA")
  
  # k larger than sequence
  result <- extractKmers("ACG", k = 5)
  expect_length(result, 0)
})

test_that("countKmers returns data.table with correct columns", {
  kmers <- c("AA", "AA", "UU", "AA", "CC")
  
  result <- countKmers(kmers)
  
  expect_s3_class(result, "data.table")
  expect_true("motif" %in% names(result))
  expect_true("count" %in% names(result))
})

test_that("countKmers counts correctly", {
  kmers <- c("AA", "AA", "AA", "UU", "CC")
  
  result <- countKmers(kmers)
  
  aa_count <- result[motif == "AA", count]
  expect_equal(aa_count, 3)
  
  uu_count <- result[motif == "UU", count]
  expect_equal(uu_count, 1)
})
