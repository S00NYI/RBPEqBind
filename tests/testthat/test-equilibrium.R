# Test Equilibrium Functions
# Tests for solveEquilibrium and calcBoundProtein

test_that("solveEquilibrium returns positive value", {
  result <- solveEquilibrium(
    RNA_total = 10,
    protein_totals = c(P1 = 100),
    Kds = c(P1 = 10)
  )
  
  expect_true(result > 0)
  expect_true(result <= 10)  # Free RNA can't exceed total
})

test_that("solveEquilibrium returns lower free RNA with high affinity", {
  # High affinity (low Kd) should bind more, leaving less free RNA
  result_high_aff <- solveEquilibrium(
    RNA_total = 10,
    protein_totals = c(P1 = 100),
    Kds = c(P1 = 1)  # Low Kd = high affinity
  )
  
  result_low_aff <- solveEquilibrium(
    RNA_total = 10,
    protein_totals = c(P1 = 100),
    Kds = c(P1 = 1000)  # High Kd = low affinity
  )
  
  expect_true(result_high_aff < result_low_aff)
})

test_that("solveEquilibrium handles multiple proteins", {
  result <- solveEquilibrium(
    RNA_total = 10,
    protein_totals = c(P1 = 100, P2 = 100),
    Kds = c(P1 = 10, P2 = 100)
  )
  
  expect_true(result > 0)
  expect_true(is.numeric(result))
})

test_that("calcBoundProtein returns named vector", {
  result <- calcBoundProtein(
    free_RNA = 5,
    protein_totals = c(P1 = 100, P2 = 50),
    Kds = c(P1 = 10, P2 = 100)
  )
  
  expect_type(result, "double")
  expect_equal(names(result), c("P1", "P2"))
})

test_that("calcBoundProtein values are sensible", {
  free_RNA <- 5
  P_total <- 100
  Kd <- 10
  
  result <- calcBoundProtein(
    free_RNA = free_RNA,
    protein_totals = c(P1 = P_total),
    Kds = c(P1 = Kd)
  )
  
  # Bound protein should be less than total
  expect_true(result["P1"] >= 0)
  expect_true(result["P1"] <= P_total)
  
  # Manual calculation: Bound = P_total * free_RNA / (Kd + free_RNA)
  expected <- P_total * free_RNA / (Kd + free_RNA)
  expect_equal(as.numeric(result["P1"]), expected, tolerance = 1e-6)
})

test_that("calcBoundProtein high affinity binds more", {
  result <- calcBoundProtein(
    free_RNA = 10,
    protein_totals = c(HighAff = 100, LowAff = 100),
    Kds = c(HighAff = 1, LowAff = 1000)
  )
  
  expect_true(result["HighAff"] > result["LowAff"])
})
