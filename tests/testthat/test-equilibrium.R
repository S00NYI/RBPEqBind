# Test Shared Protein Pool Solver
# Tests for solveSharedPool

random_case <- function(n_sites, n_prot) {
  kd <- matrix(10^runif(n_sites * n_prot, -3, 3), n_sites, n_prot)
  colnames(kd) <- paste0("P", seq_len(n_prot))
  list(S = 10^runif(n_sites, -3, 3), kd = kd, P = 10^runif(n_prot, -2, 4))
}

test_that("free plus bound equals total for every protein", {
  set.seed(1)
  for (i in seq_len(100)) {
    cs <- random_case(sample(seq_len(50), 1), sample(seq_len(4), 1))
    sol <- solveSharedPool(cs$S, cs$kd, cs$P)
    expect_true(sol$converged)
    expect_lt(max(abs(sol$free + sol$bound - cs$P) / cs$P), 1e-8)
  }
})

test_that("bound protein never exceeds total protein", {
  set.seed(2)
  for (i in seq_len(100)) {
    cs <- random_case(sample(seq_len(50), 1), sample(seq_len(4), 1))
    sol <- solveSharedPool(cs$S, cs$kd, cs$P)
    expect_true(all(sol$bound <= cs$P * (1 + 1e-9)))
    expect_true(all(sol$free >= 0))
  }
})

test_that("fractions bound stay in [0, 1] and sum to at most 1 per site", {
  set.seed(3)
  cs <- random_case(40, 3)
  sol <- solveSharedPool(cs$S, cs$kd, cs$P)
  expect_true(all(sol$theta >= 0 & sol$theta <= 1))
  expect_true(all(rowSums(sol$theta) <= 1 + 1e-12))
})

test_that("one site and one protein match the exact quadratic solution", {
  # Complex C solves C^2 - (P + S + K) C + P S = 0 (smaller root)
  exact <- function(S, P, K) {
    b <- P + S + K
    2 * P * S / (b + sqrt(b^2 - 4 * P * S))
  }
  for (p in list(c(10, 100, 10), c(10, 1e4, 0.01), c(1e-3, 1e-3, 1e-3),
                 c(1, 1e6, 1e-4), c(1000, 100, 50))) {
    sol <- solveSharedPool(p[1], matrix(p[3], dimnames = list(NULL, "P")), p[2])
    expect_equal(as.numeric(sol$theta[1, 1]), exact(p[1], p[2], p[3]) / p[1],
                 tolerance = 1e-8)
  }
})

test_that("protein in large excess gives the independent-site closed form", {
  kd <- matrix(c(5, 50, 500, 20, 2, 200), 3, 2,
               dimnames = list(NULL, c("A", "B")))
  P <- c(A = 100, B = 300)
  sol <- solveSharedPool(rep(1e-9, 3), kd, P)
  q <- sweep(1 / kd, 2, P, `*`)
  expect_equal(sol$theta, q / (1 + rowSums(q)), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("results do not depend on the concentration unit", {
  set.seed(4)
  cs <- random_case(20, 2)
  ref <- solveSharedPool(cs$S, cs$kd, cs$P)$theta
  for (s in c(1e-3, 1e-6, 1e-9, 1e3)) {
    got <- solveSharedPool(cs$S * s, cs$kd * s, cs$P * s)$theta
    expect_equal(got, ref, tolerance = 1e-7)
  }
})

test_that("more sites leave less free protein for each site", {
  kd <- matrix(10, dimnames = list(NULL, "P"))
  occ <- vapply(c(1, 10, 100, 1000), function(S) {
    solveSharedPool(S, kd, 100)$theta[1, 1]
  }, numeric(1))
  expect_true(all(diff(occ) < 0))
})

test_that("a competitor lowers the occupancy of a protein at shared sites", {
  kd <- matrix(c(10, 10), 1, 2, dimnames = list(NULL, c("A", "B")))
  alone <- solveSharedPool(5, kd, c(100, 0))$theta[1, "A"]
  with_b <- solveSharedPool(5, kd, c(100, 100))$theta[1, "A"]
  expect_lt(with_b, alone)
})

test_that("edge cases give the expected values", {
  kd <- matrix(c(10, 20), 2, 1, dimnames = list(NULL, "P"))
  none <- solveSharedPool(c(1, 1), kd, 0)
  expect_equal(none$free, c(P = 0))
  expect_true(all(none$theta == 0))

  no_sites <- solveSharedPool(c(0, 0), kd, 50)
  expect_equal(no_sites$free, c(P = 50))

  no_binding <- solveSharedPool(c(1, 1), matrix(Inf, 2, 1, dimnames = list(NULL, "P")), 50)
  expect_true(all(no_binding$theta == 0))
  expect_equal(no_binding$free, c(P = 50), tolerance = 1e-10)
})

test_that("invalid input stops with an error", {
  kd <- matrix(10, dimnames = list(NULL, "P"))
  expect_error(solveSharedPool(-1, kd, 10), "non-negative")
  expect_error(solveSharedPool(1, kd, -10), "non-negative")
  expect_error(solveSharedPool(1, matrix(0, dimnames = list(NULL, "P")), 10), "positive")
  expect_error(solveSharedPool(c(1, 2), kd, 10), "one value for each row")
})
