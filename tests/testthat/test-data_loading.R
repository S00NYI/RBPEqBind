# Test Model Loading Functions
# Tests for loadModel and setModel

test_that("loadModel reads multi-RBP CSV correctly", {
  # Test with actual data file
  skip_if_not(file.exists("../../data/model_RBP.csv"))
  
  result <- loadModel("../../data/model_RBP.csv")
  
  expect_type(result, "list")
  expect_true(length(result) >= 1)
  
  # Check structure of each RBP model
  for (rbp in names(result)) {
    expect_true("motif" %in% names(result[[rbp]]))
    expect_true("score" %in% names(result[[rbp]]))
  }
})

test_that("loadModel filters by k-mer length", {
  skip_if_not(file.exists("../../data/model_RBP.csv"))
  
  result <- loadModel("../../data/model_RBP.csv", k = 5)
  
  for (rbp in names(result)) {
    motif_lengths <- nchar(result[[rbp]]$motif)
    expect_true(all(motif_lengths == 5))
  }
})

test_that("loadModel filters by RBP names", {
  skip_if_not(file.exists("../../data/model_RBP.csv"))
  
  result <- loadModel("../../data/model_RBP.csv", rbp = c("HH", "LL"))
  
  expect_equal(length(result), 2)
  expect_true("HH" %in% names(result))
  expect_true("LL" %in% names(result))
})

test_that("loadModel errors on missing file", {
  expect_error(loadModel("nonexistent_file.csv"))
})

test_that("setModel converts score to Kd correctly", {
  # Create mock raw model
  mock_raw <- list(
    RBP1 = data.table::data.table(
      motif = c("AAAAA", "UUUUU"),
      score = c(1.0, 0.0)  # Max and min scores
    )
  )
  
  result <- setModel(mock_raw, max_affinity = 100, min_affinity = 0.01)
  
  expect_type(result, "list")
  expect_true("Kd" %in% names(result$RBP1))
  expect_true("motif" %in% names(result$RBP1))
  
  # High score should have low Kd (high affinity)
  expect_true(result$RBP1$Kd[1] < result$RBP1$Kd[2])
})

test_that("setModel accepts per-RBP affinity values", {
  mock_raw <- list(
    HH = data.table::data.table(motif = "AAAAA", score = 0.5),
    LL = data.table::data.table(motif = "UUUUU", score = 0.5)
  )
  
  result <- setModel(mock_raw, 
                      max_affinity = c(HH = 100, LL = 10),
                      min_affinity = c(HH = 0.01, LL = 0.001))
  
  expect_true("Kd" %in% names(result$HH))
  expect_true("Kd" %in% names(result$LL))
})
