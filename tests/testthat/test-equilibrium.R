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

test_that("equilibrium conserves RNA mass across scales", {
  set.seed(1)
  for (i in seq_len(100)) {
    n <- sample(seq_len(5), 1)
    Rt <- 10^runif(1, -6, 4)
    P <- 10^runif(n, -2, 4)
    K <- 10^runif(n, -4, 3)
    X <- solveEquilibrium(Rt, P, K)
    B <- calcBoundProtein(X, P, K)
    expect_lt(abs(X + sum(B) - Rt) / Rt, 1e-7)
  }
})

test_that("occupancies never exceed the simplex", {
  set.seed(2)
  for (i in seq_len(100)) {
    n <- sample(seq_len(5), 1)
    Rt <- 10^runif(1, -6, 4)
    P <- 10^runif(n, -2, 4)
    K <- 10^runif(n, -4, 3)
    X <- solveEquilibrium(Rt, P, K)
    expect_lte(sum(calcBoundProtein(X, P, K)) / Rt, 1 + 1e-9)
  }
})

test_that("results are invariant to a change of concentration units", {
  ref <- calcBoundProtein(solveEquilibrium(10, 100, 10), 100, 10) / 10
  for (s in c(1e-3, 1e-6, 1e-9, 1e-12)) {
    got <- calcBoundProtein(solveEquilibrium(10 * s, 100 * s, 10 * s), 100 * s, 10 * s) / (10 * s)
    expect_equal(as.numeric(got), as.numeric(ref), tolerance = 1e-6)
  }
})

test_that("single-RBP case matches the stable closed form", {
  stable <- function(Rt, P, K) {
    b <- K + P - Rt
    2 * K * Rt / (b + sqrt(b^2 + 4 * K * Rt))
  }
  for (p in list(c(10, 100, 10), c(10, 1e4, 0.01),
                 c(1e-3, 1e-3, 1e-3), c(1, 1e6, 1e-4))) {
    expect_equal(solveEquilibrium(p[1], p[2], p[3]),
                 stable(p[1], p[2], p[3]), tolerance = 1e-6)
  }
})

test_that("solveEquilibrium edge cases produce expected roots", {
  expect_equal(solveEquilibrium(10, c(100, 100), c(Inf, Inf)), 10)
  expect_equal(solveEquilibrium(10, c(0, 0), c(1, 1)), 10)
  expect_equal(solveEquilibrium(0, c(100, 100), c(1, 1)), 0)
  expect_true(solveEquilibrium(1e-12, c(100), c(0.01)) < 1e-12)
})
