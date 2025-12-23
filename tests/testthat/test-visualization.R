# Test Visualization Functions
# Tests for plotBinding, plotHeatmap, plotGrid, plotBubble, viewModel

test_that("plotBinding returns ggplot object", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = runif(10, 0, 0.5),
    RBP1_density = runif(10, 0.05, 0.15),
    RBP1_density_fc = runif(10, 0.5, 1.5),
    RBP1_occupancy_fc = runif(10, 0.5, 1.5),
    transcript = "TX1"
  )
  
  p <- plotBinding(mock_result, rbp = "RBP1", transcript = "TX1")
  
  expect_s3_class(p, "ggplot")
})

test_that("plotBinding accepts all metric types", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = runif(5),
    RBP1_density = runif(5),
    RBP1_density_fc = runif(5),
    RBP1_occupancy_fc = runif(5),
    transcript = "TX1"
  )
  
  for (metric in c("occupancy", "density", "density_fc", "occupancy_fc")) {
    p <- plotBinding(mock_result, rbp = "RBP1", transcript = "TX1", metric = metric)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plotHeatmap returns ggplot object", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = runif(10),
    RBP2 = runif(10),
    RBP1_density = runif(10),
    RBP2_density = runif(10),
    transcript = "TX1"
  )
  
  p <- plotHeatmap(mock_result, transcript = "TX1")
  
  expect_s3_class(p, "ggplot")
})

test_that("plotHeatmap respects zlim parameter", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = runif(5),
    RBP1_density = runif(5),
    transcript = "TX1"
  )
  
  p <- plotHeatmap(mock_result, transcript = "TX1", metric = "occupancy", zlim = c(0, 1))
  
  expect_s3_class(p, "ggplot")
})

test_that("viewModel returns ggplot for raw models", {
  skip_if_not_installed("ggplot2")
  
  mock_raw <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU", "CC"), score = c(0.5, 0.8, 0.3))
  )
  
  p <- viewModel(mock_raw, rbp = "RBP1")
  
  expect_s3_class(p, "ggplot")
})

test_that("viewModel returns ggplot for processed models", {
  skip_if_not_installed("ggplot2")
  
  mock_processed <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU", "CC"), Kd = c(10, 100, 50))
  )
  
  p <- viewModel(mock_processed, rbp = "RBP1", metric = "Kd")
  expect_s3_class(p, "ggplot")
  
  p <- viewModel(mock_processed, rbp = "RBP1", metric = "Ka")
  expect_s3_class(p, "ggplot")
})

test_that("viewModel handles multiple RBPs", {
  skip_if_not_installed("ggplot2")
  
  mock_models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100)),
    RBP2 = data.table::data.table(motif = c("AA", "UU"), Kd = c(50, 50))
  )
  
  p <- viewModel(mock_models, rbp = c("RBP1", "RBP2"))
  
  expect_s3_class(p, "ggplot")
})
