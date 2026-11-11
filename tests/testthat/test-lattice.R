# Test Lattice Binding Model
# Tests for latticeStarts, solveLatticePool, latticeCoverage, resolveFootprint
# and simulateBinding(footprint = ...)

# Reference: enumerate every allowed arrangement of non-overlapping proteins.
brute_force <- function(q_list, w, L) {
  n <- length(q_list)
  starts <- lapply(q_list, function(q) numeric(length(q)))
  Z <- 0
  visit <- function(i, weight, placed) {
    if (i > L) {
      Z <<- Z + weight
      for (pl in placed) {
        starts[[pl[1]]][pl[2]] <<- starts[[pl[1]]][pl[2]] + weight
      }
      return(invisible(NULL))
    }
    visit(i + 1L, weight, placed)
    for (j in seq_len(n)) {
      if (i <= length(q_list[[j]])) {
        visit(i + w[j], weight * q_list[[j]][i], c(placed, list(c(j, i))))
      }
    }
  }
  visit(1L, 1, list())
  list(Z = Z, starts = lapply(starts, function(s) s / Z))
}

test_that("one RBP matches brute-force enumeration", {
  set.seed(10)
  for (L in c(5, 9, 12)) {
    w <- 3L
    q <- 10^runif(L - w + 1, -1, 1)
    lat <- latticeStarts(list(q), w, 1, L)
    ref <- brute_force(list(q), w, L)
    expect_equal(lat$logZ, log(ref$Z), tolerance = 1e-12)
    expect_equal(lat$starts[[1]], ref$starts[[1]], tolerance = 1e-12)
  }
})

test_that("two RBPs with different footprints match brute-force enumeration", {
  set.seed(11)
  L <- 13L
  w <- c(3L, 5L)
  q <- list(10^runif(L - w[1] + 1, -1, 1), 10^runif(L - w[2] + 1, -1, 1))
  lat <- latticeStarts(q, w, c(1, 1), L)
  ref <- brute_force(q, w, L)
  expect_equal(lat$logZ, log(ref$Z), tolerance = 1e-12)
  expect_equal(lat$starts[[1]], ref$starts[[1]], tolerance = 1e-12)
  expect_equal(lat$starts[[2]], ref$starts[[2]], tolerance = 1e-12)
})

test_that("a single placement gives q / (1 + q)", {
  for (q in c(0.01, 1, 50)) {
    lat <- latticeStarts(list(q), 7L, 1, 7L)
    expect_equal(lat$starts[[1]], q / (1 + q), tolerance = 1e-14)
  }
})

test_that("a long uniform lattice matches McGhee-von Hippel", {
  # Bound ligands per lattice unit (nu) for footprint n and weight q solve
  # nu = q * (1 - n nu)^n / (1 - (n - 1) nu)^(n - 1)
  # (McGhee and von Hippel, J Mol Biol 86:469, 1974).
  mvh <- function(q, n) {
    f <- function(nu) nu - q * (1 - n * nu)^n / (1 - (n - 1) * nu)^(n - 1)
    stats::uniroot(f, c(0, 1 / n - 1e-12), tol = 1e-14)$root
  }
  L <- 3000L
  for (n in c(1L, 5L, 11L)) {
    for (q in c(0.05, 1, 20)) {
      lat <- latticeStarts(list(rep(q, L - n + 1)), n, 1, L)
      mid <- lat$starts[[1]][1400:1600]
      expect_equal(mean(mid), mvh(q, n), tolerance = 1e-8)
    }
  }
})

test_that("coverage stays in [0, 1] and summed over RBPs stays at most 1", {
  set.seed(12)
  L <- 200L
  w <- c(4L, 9L)
  q <- list(10^runif(L - w[1] + 1, -2, 3), 10^runif(L - w[2] + 1, -2, 3))
  lat <- latticeStarts(q, w, c(1, 1), L)
  cov <- cbind(latticeCoverage(lat$starts[[1]], w[1], L),
               latticeCoverage(lat$starts[[2]], w[2], L))
  expect_true(all(cov >= 0 & cov <= 1 + 1e-12))
  expect_true(all(rowSums(cov) <= 1 + 1e-12))
})

test_that("long sequences with strong binding do not overflow", {
  L <- 20000L
  lat <- latticeStarts(list(rep(1e3, L - 6)), 7L, 1, L)
  expect_true(is.finite(lat$logZ))
  cov <- latticeCoverage(lat$starts[[1]], 7L, L)
  expect_true(all(is.finite(cov)))
  expect_true(all(cov <= 1 + 1e-12))
  expect_lte(lat$n_bound, L / 7)
})

test_that("the lattice protein pool satisfies the mass balance", {
  set.seed(13)
  L <- 120L
  w <- c(5L, 7L)
  a <- list(10^runif(L - w[1] + 1, -3, 0), 10^runif(L - w[2] + 1, -3, 0))
  P <- c(50, 20)
  for (rna in c(0.1, 5, 50)) {
    sol <- solveLatticePool(a, w, L, rna, P)
    expect_true(sol$converged)
    expect_lt(max(abs(sol$free + rna * sol$n_bound - P) / P), 1e-8)
  }
})

test_that("resolveFootprint checks its input", {
  expect_equal(resolveFootprint(7, c("A", "B"), 5), c(A = 7L, B = 7L))
  expect_equal(resolveFootprint(c(B = 9, A = 6), c("A", "B"), 5), c(A = 6L, B = 9L))
  expect_error(resolveFootprint(4, "A", 5), "at least k")
  expect_error(resolveFootprint(c(6, 7), c("A", "B"), 5), "named")
  expect_error(resolveFootprint(6.5, "A", 5), "whole numbers")
})

test_that("simulateBinding with a footprint equal to the sequence length gives q / (1 + q)", {
  models <- list(RBP1 = data.table::data.table(motif = c("AAAAA", "UUUUU"), Kd = c(10, 100)))
  res <- simulateBinding("AAAAA", models, c(RBP1 = 30), rna_conc = 1e-9, k = 5,
                         footprint = 5)
  expect_equal(res$RBP1, rep(3 / 4, 5), tolerance = 1e-8)
})

test_that("a footprint limits how many proteins fit on a U-tract", {
  b <- c("A", "C", "G", "U")
  all5 <- do.call(paste0, expand.grid(b, b, b, b, b, stringsAsFactors = FALSE))
  models <- list(RBP1 = data.table::data.table(motif = all5,
                                               Kd = ifelse(all5 == "UUUUU", 1e-3, 1e6)))
  seq <- paste0("GGG", strrep("U", 20), "GGG")
  indep <- simulateBinding(seq, models, c(RBP1 = 1e4), rna_conc = 1e-9, k = 5)
  latt <- simulateBinding(seq, models, c(RBP1 = 1e4), rna_conc = 1e-9, k = 5,
                          footprint = 10)
  # Without a footprint every U-window is near 1. With a 10-nt footprint at
  # most 2 proteins fit on 20 nt, so covered nucleotides sum to at most 20.
  # Positions 8-19 are covered only by UUUUU windows
  expect_gt(min(indep$RBP1[8:19]), 0.99)
  expect_lte(sum(latt$RBP1), 20 + 1e-8)
  expect_true(all(latt$RBP1 <= 1 + 1e-12))
})

test_that("footprints can differ between RBPs and total coverage stays at most 1", {
  models <- list(
    A = data.table::data.table(motif = c("AA", "AU", "UA", "UU"), Kd = c(1, 10, 10, 100)),
    B = data.table::data.table(motif = c("AA", "AU", "UA", "UU"), Kd = c(100, 10, 10, 1))
  )
  res <- simulateBinding("AAUUAAUUAAUUAAUU", models, c(A = 50, B = 50),
                         rna_conc = 5, k = 2, footprint = c(A = 3, B = 4))
  expect_true(all(res$A + res$B <= 1 + 1e-12))
  expect_true(all(res$A >= 0 & res$B >= 0))
})

test_that("lattice results do not depend on the concentration unit", {
  models <- list(RBP1 = data.table::data.table(motif = c("AA", "AU", "UA", "UU"),
                                               Kd = c(10, 50, 50, 100)))
  scaled <- function(s) {
    m <- list(RBP1 = data.table::copy(models$RBP1))
    m$RBP1$Kd <- m$RBP1$Kd * s
    simulateBinding("AAUUAAUUAA", m, c(RBP1 = 100 * s), rna_conc = 5 * s,
                    k = 2, footprint = 3)$RBP1
  }
  ref <- scaled(1)
  for (s in c(1e-3, 1e3)) expect_equal(scaled(s), ref, tolerance = 1e-7)
})

test_that("simulateGrid passes the footprint on", {
  models <- list(RBP1 = data.table::data.table(motif = c("AA", "AU", "UA", "UU"),
                                               Kd = c(10, 50, 50, 100)))
  g <- simulateGrid("AAUUAAUU", models, protein_conc_grid = list(RBP1 = c(10, 100)),
                    rna_conc_grid = 1, k = 2, footprint = 3)
  one <- simulateBinding("AAUUAAUU", models, c(RBP1 = 100), rna_conc = 1, k = 2,
                         footprint = 3)
  expect_equal(g[Conc_RBP1 == 100]$RBP1, one$RBP1)
})
