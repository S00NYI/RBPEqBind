# Test RNA Accessibility
# Tests for readAccessibility(), checkAccessibility() and the accessibility
# argument of simulateBinding(), simulateGrid(), simulateBindingFasta() and
# simulateGridFasta()

all3 <- do.call(paste0, expand.grid(c("A", "C", "G", "U"), c("A", "C", "G", "U"),
                                    c("A", "C", "G", "U"), stringsAsFactors = FALSE))
set.seed(7)
acc_models <- list(
  R1 = data.table::data.table(motif = all3, Kd = 10^runif(64, 0, 3)),
  R2 = data.table::data.table(motif = all3, Kd = 10^runif(64, 0, 3)))
acc_seq <- function(n) paste(sample(c("A", "C", "G", "U"), n, TRUE), collapse = "")
write_fa <- function(seqs) {
  f <- tempfile(fileext = ".fa")
  writeLines(paste0(">", names(seqs), "\n", gsub("U", "T", seqs)), f)
  f
}
# Accessibility from per-nucleotide unpaired probabilities, assuming the
# nucleotides are independent (only to build test input)
acc_from_nt <- function(p, u) {
  L <- length(p)
  m <- matrix(NA_real_, L, u)
  for (l in seq_len(u)) for (i in seq.int(l, L)) m[i, l] <- prod(p[(i - l + 1):i])
  m
}
cols <- c("R1", "R2")

test_that("readAccessibility() reads an RNAplfold _lunp file", {
  f <- tempfile(fileext = "_lunp")
  writeLines(c("#unpaired probabilities", " #i$\tl=1\t2\t3\t",
               "1\t0.9\tNA\tNA\t", "2\t0.8\t0.75\tNA\t", "3\t0.2\t0.15\t0.1\t",
               "4\t0.3\t0.1\t0.05\t"), f)
  m <- readAccessibility(f)
  expect_equal(dim(m), c(4L, 3L))
  expect_equal(colnames(m), c("1", "2", "3"))
  expect_equal(unname(m[2, 2]), 0.75)
  expect_true(is.na(m[1, 2]))
  expect_error(readAccessibility(tempfile()), "existing")
})

test_that("checkAccessibility() rejects bad input", {
  m <- acc_from_nt(rep(0.5, 10), 4)
  expect_null(RBPEqBind:::checkAccessibility(NULL, 10, 3))
  expect_error(RBPEqBind:::checkAccessibility(m, 11, 3), "rows")
  expect_error(RBPEqBind:::checkAccessibility(m, 10, 5), "columns")
  bad <- m; bad[6, 3] <- NA
  expect_error(RBPEqBind:::checkAccessibility(bad, 10, 3), "NA")
  bad <- m; bad[6, 3] <- 1.2
  expect_error(RBPEqBind:::checkAccessibility(bad, 10, 3), "between 0 and 1")
  expect_error(RBPEqBind:::checkAccessibility("x", 10, 3), "numeric matrix")
  # NA where the stretch does not fit is allowed
  expect_identical(RBPEqBind:::checkAccessibility(m, 10, 4), m)
})

test_that("accessibility of 1 everywhere gives the same result as NULL", {
  s <- acc_seq(80)
  ones <- matrix(1, nchar(s), 9)
  for (fp in list(NULL, 5L)) {
    ref <- simulateBinding(s, acc_models, c(R1 = 30, R2 = 10), rna_conc = 5,
                           k = 3, footprint = fp)
    got <- simulateBinding(s, acc_models, c(R1 = 30, R2 = 10), rna_conc = 5,
                           k = 3, footprint = fp, accessibility = ones)
    expect_equal(as.data.frame(got), as.data.frame(ref), tolerance = 1e-8)
  }
})

test_that("one window without a footprint follows the three-state formula", {
  # One k-mer, no RNA depletion (rna_conc = 0): occupancy = x / (1 + x),
  # x = a * P / Kd, with a = probability that the k-mer is unpaired.
  s <- "ACG"
  kd <- acc_models$R1$Kd[acc_models$R1$motif == s]
  for (a in c(1, 0.3, 0)) {
    m <- matrix(c(NA, NA, a, NA, NA, NA, NA, NA, a), 3, 3)
    m[, 1] <- 1; m[2:3, 2] <- 1
    r <- simulateBinding(s, acc_models["R1"], c(R1 = 50), rna_conc = 0, k = 3,
                         accessibility = m)
    x <- a * 50 / kd
    expect_equal(r$R1, rep(x / (1 + x), 3), tolerance = 1e-10)
  }
})

test_that("one placement on the lattice uses a of the middle k-mer", {
  s <- "ACGUA"   # footprint 5: one placement; the middle 3-mer (CGU) sets Kd
  kd <- acc_models$R1$Kd[acc_models$R1$motif == "CGU"]
  m <- matrix(1, 5, 5); m[upper.tri(m)] <- NA
  m[4, 3] <- 0.25   # nucleotides 2 to 4 (the middle 3-mer)
  m[5, 5] <- 0.1    # nucleotides 1 to 5 (the whole footprint): not used
  r <- simulateBinding(s, acc_models["R1"], c(R1 = 50), rna_conc = 0, k = 3,
                       footprint = 5L, accessibility = m)
  x <- 0.25 * 50 / kd
  expect_equal(r$R1, rep(x / (1 + x), 5), tolerance = 1e-10)
})

test_that("bound protein follows the mass balance", {
  # One site, one protein: occupancy theta = a f / (Kd + a f), and
  # P = f + c * theta.
  P <- 50; rna <- 20; a <- 0.3
  kd <- acc_models$R1$Kd[acc_models$R1$motif == "CGU"]
  f <- stats::uniroot(function(f) f + rna * a * f / (kd + a * f) - P, c(0, P),
                      tol = 1e-14)$root
  want <- a * f / (kd + a * f)
  m3 <- matrix(1, 3, 3); m3[upper.tri(m3)] <- NA; m3[3, 3] <- a
  r <- simulateBinding("CGU", acc_models["R1"], c(R1 = P), rna_conc = rna,
                       k = 3, accessibility = m3)
  expect_equal(r$R1, rep(want, 3), tolerance = 1e-8)
  m5 <- matrix(1, 5, 5); m5[upper.tri(m5)] <- NA; m5[4, 3] <- a
  r5 <- simulateBinding("ACGUA", acc_models["R1"], c(R1 = P), rna_conc = rna,
                        k = 3, footprint = 5L, accessibility = m5)
  expect_equal(r5$R1, rep(want, 5), tolerance = 1e-8)
})

test_that("a closed nucleotide is never covered without a footprint", {
  s <- acc_seq(40)
  p <- rep(0.9, 40); p[20] <- 0
  m <- acc_from_nt(p, 5)
  # Every window that contains nucleotide 20 has a = 0
  r0 <- simulateBinding(s, acc_models, c(R1 = 100, R2 = 100), rna_conc = 1,
                        k = 3, accessibility = m)
  expect_equal(r0$R1[20], 0)
  expect_equal(r0$R2[20], 0)
})

test_that("simulateGrid() passes accessibility on", {
  s <- acc_seq(30)
  m <- acc_from_nt(runif(30, 0.2, 1), 5)
  g <- simulateGrid(s, acc_models, list(R1 = c(10, 100), R2 = 20), rna_conc_grid = 1,
                    k = 3, footprint = 5L, accessibility = m)
  one <- simulateBinding(s, acc_models, c(R1 = 100, R2 = 20), rna_conc = 1, k = 3,
                         footprint = 5L, accessibility = m)
  expect_equal(g[Conc_R1 == 100, R1], one$R1, tolerance = 1e-10)
})

test_that("FASTA runs match simulateBinding() per sequence and by name", {
  seqs <- c(tA = acc_seq(50), tB = acc_seq(35))
  fa <- write_fa(seqs)
  acc <- list(tB = acc_from_nt(runif(35, 0.1, 1), 5),
              tA = acc_from_nt(runif(50, 0.1, 1), 5))
  for (fp in list(NULL, 5L)) {
    res <- simulateBindingFasta(fa, acc_models, c(R1 = 20, R2 = 40), rna_conc = 2,
                                k = 3, footprint = fp, accessibility = acc)
    for (t in names(seqs)) {
      one <- simulateBinding(seqs[[t]], acc_models, c(R1 = 20, R2 = 40), rna_conc = 2,
                             k = 3, footprint = fp, accessibility = acc[[t]])
      expect_equal(res[transcript == t, R1], one$R1, tolerance = 1e-10)
    }
    # One sequence with shared protein pool equals it alone
    fa1 <- write_fa(seqs["tA"])
    alone <- simulateBindingFasta(fa1, acc_models, c(R1 = 20, R2 = 40), rna_conc = 2,
                                  k = 3, footprint = fp, accessibility = acc["tA"])
    shared <- simulateBindingFasta(fa1, acc_models, c(R1 = 20, R2 = 40), rna_conc = 2,
                                   k = 3, footprint = fp, accessibility = acc["tA"],
                                   shared_protein = TRUE)
    expect_equal(as.data.frame(shared), as.data.frame(alone), tolerance = 1e-8)
  }
  expect_error(simulateBindingFasta(fa, acc_models, c(R1 = 20), k = 3,
                                    accessibility = acc["tA"]), "Missing: tB")
  expect_error(simulateBindingFasta(fa, acc_models, c(R1 = 20), k = 3,
                                    accessibility = acc[[1]]), "list of matrices")
})

test_that("shared protein pool with accessibility of 1 gives the same result as NULL", {
  seqs <- c(tA = acc_seq(50), tB = acc_seq(35))
  fa <- write_fa(seqs)
  ones <- lapply(seqs, function(s) matrix(1, nchar(s), 5))
  for (fp in list(NULL, 5L)) {
    ref <- simulateBindingFasta(fa, acc_models, c(R1 = 20, R2 = 40), rna_conc = 3,
                                k = 3, footprint = fp, shared_protein = TRUE)
    got <- simulateBindingFasta(fa, acc_models, c(R1 = 20, R2 = 40), rna_conc = 3,
                                k = 3, footprint = fp, shared_protein = TRUE,
                                accessibility = ones)
    expect_equal(as.data.frame(got), as.data.frame(ref), tolerance = 1e-8)
    g <- simulateGridFasta(fa, acc_models, list(R1 = 20, R2 = 40), rna_conc_grid = 3,
                           k = 3, footprint = fp, shared_protein = TRUE,
                           accessibility = ones)
    expect_equal(g$R1, ref$R1, tolerance = 1e-8)
  }
})
