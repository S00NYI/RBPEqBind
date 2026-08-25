# Test Visualization Functions
# Tests for plotBinding, plotHeatmap, plotGrid, viewModel

test_that("plotBinding returns ggplot object and matches snapshot", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = c(0.1, 0.2, 0.4, 0.3, 0.1, 0.05, 0.05, 0.1, 0.2, 0.3),
    RBP1_density = c(0.05, 0.1, 0.2, 0.15, 0.05, 0.025, 0.025, 0.05, 0.1, 0.15),
    RBP1_density_fc = c(0.5, 1.0, 2.0, 1.5, 0.5, 0.25, 0.25, 0.5, 1.0, 1.5),
    RBP1_occupancy_fc = c(0.5, 1.0, 2.0, 1.5, 0.5, 0.25, 0.25, 0.5, 1.0, 1.5),
    transcript = "TX1"
  )
  
  p <- plotBinding(mock_result, rbp = "RBP1", transcript = "TX1")
  expect_s3_class(p, "ggplot")
  
  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("plotBinding-basic", p)
  }
})

test_that("plotBinding accepts all metric types", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.2, 0.4, 0.3, 0.1),
    RBP1_density = c(0.1, 0.2, 0.4, 0.3, 0.1) / 1.1,
    RBP1_density_fc = c(0.5, 1.0, 2.0, 1.5, 0.5),
    RBP1_occupancy_fc = c(0.5, 1.0, 2.0, 1.5, 0.5),
    transcript = "TX1"
  )
  
  for (metric in c("occupancy", "density", "density_fc", "occupancy_fc")) {
    p <- plotBinding(mock_result, rbp = "RBP1", transcript = "TX1", metric = metric)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plotHeatmap returns ggplot object and matches snapshot", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:10,
    nt = c("A", "C", "G", "U", "A", "C", "G", "U", "A", "C"),
    RBP1 = c(0.1, 0.2, 0.4, 0.3, 0.1, 0.05, 0.05, 0.1, 0.2, 0.3),
    RBP2 = c(0.3, 0.2, 0.1, 0.05, 0.05, 0.1, 0.2, 0.4, 0.3, 0.1),
    RBP1_density = c(0.05, 0.1, 0.2, 0.15, 0.05, 0.025, 0.025, 0.05, 0.1, 0.15),
    RBP2_density = c(0.15, 0.1, 0.05, 0.025, 0.025, 0.05, 0.1, 0.2, 0.15, 0.05),
    transcript = "TX1"
  )
  
  p <- plotHeatmap(mock_result, transcript = "TX1")
  expect_s3_class(p, "ggplot")
  
  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("plotHeatmap-basic", p)
  }
})

test_that("plotHeatmap respects zlim parameter", {
  skip_if_not_installed("ggplot2")
  
  mock_result <- data.table::data.table(
    pos = 1:5,
    nt = c("A", "C", "G", "U", "A"),
    RBP1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
    RBP1_density = c(0.1, 0.2, 0.3, 0.4, 0.5) / 1.5,
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

test_that("viewModel returns ggplot and matches snapshot for processed models", {
  skip_if_not_installed("ggplot2")
  
  mock_processed <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU", "CC"), Kd = c(10, 100, 50))
  )
  
  p <- viewModel(mock_processed, rbp = "RBP1", metric = "Kd")
  expect_s3_class(p, "ggplot")
  
  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("viewModel-Kd", p)
  }
  
  p2 <- viewModel(mock_processed, rbp = "RBP1", metric = "Ka")
  expect_s3_class(p2, "ggplot")
})

test_that("viewModel handles multiple RBPs with Okabe-Ito palette", {
  skip_if_not_installed("ggplot2")
  
  mock_models <- list(
    RBP1 = data.table::data.table(motif = c("AA", "UU"), Kd = c(10, 100)),
    RBP2 = data.table::data.table(motif = c("AA", "UU"), Kd = c(50, 50))
  )
  
  p <- viewModel(mock_models, rbp = c("RBP1", "RBP2"))
  expect_s3_class(p, "ggplot")
  
  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("viewModel-multiRBP", p)
  }
})

test_that("plotGrid returns ggplot object and matches snapshot", {
  skip_if_not_installed("ggplot2")
  
  mock_grid <- data.table::data.table(
    pos = rep(1:5, 4),
    nt = rep(c("A", "C", "G", "U", "A"), 4),
    Conc_RBP1 = rep(c(10, 100), each = 10),
    Conc_RBP2 = rep(c(10, 100), 2, each = 5),
    rna_conc = 10,
    RBP1 = c(0.1, 0.2, 0.3, 0.2, 0.1,  0.2, 0.3, 0.4, 0.3, 0.2,
             0.3, 0.4, 0.5, 0.4, 0.3,  0.4, 0.5, 0.6, 0.5, 0.4),
    RBP2 = c(0.2, 0.3, 0.2, 0.1, 0.1,  0.3, 0.4, 0.3, 0.2, 0.1,
             0.1, 0.2, 0.3, 0.2, 0.1,  0.2, 0.3, 0.4, 0.3, 0.2)
  )
  mock_grid[, RBP1_density := RBP1 / sum(RBP1), by = .(Conc_RBP1, Conc_RBP2)]
  mock_grid[, RBP1_density_fc := RBP1_density * 5]
  
  p <- plotGrid(mock_grid, rbp1 = "RBP1", rbp2 = "RBP2", roi_range = c(1, 5), metric = "density")
  expect_s3_class(p, "ggplot")
  
  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("plotGrid-density", p)
  }
})
