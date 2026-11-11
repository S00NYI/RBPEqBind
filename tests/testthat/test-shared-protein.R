# Test Shared Protein Supply Over Many Sequences
# Tests for simulateBindingFasta(shared_protein = TRUE), simulateGridFasta(
# shared_protein = TRUE), solveLatticePoolMulti and resolveRnaConc

write_fasta <- function(seqs) {
  f <- tempfile(fileext = ".fa")
  writeLines(paste0(">", names(seqs), "\n", gsub("U", "T", seqs)), f)
  f
}

# Complete 3-mer models for two RBPs, so no motif is missing
all3 <- do.call(paste0, expand.grid(c("A", "C", "G", "U"), c("A", "C", "G", "U"),
                                    c("A", "C", "G", "U"), stringsAsFactors = FALSE))
set.seed(42)
models3 <- list(
  R1 = data.table::data.table(motif = all3, Kd = 10^runif(64, 0, 3)),
  R2 = data.table::data.table(motif = all3, Kd = 10^runif(64, 0, 3)))
rand_seq <- function(n) paste(sample(c("A", "C", "G", "U"), n, TRUE), collapse = "")
seqs3 <- c(tA = rand_seq(60), tB = rand_seq(35), tC = rand_seq(90))
occ_cols <- c("R1", "R2")

test_that("one sequence gives the same result with and without shared_protein", {
  fa <- write_fasta(seqs3[1])
  for (fp in list(NULL, 5L)) {
    alone <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40), rna_conc = 2,
                                  k = 3, footprint = fp)
    shared <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40), rna_conc = 2,
                                   k = 3, footprint = fp,
                                   shared_protein = TRUE)
    expect_equal(as.data.frame(shared), as.data.frame(alone), tolerance = 1e-8)
  }
})

test_that("two identical sequences at c each equal one sequence alone at 2c", {
  # Both copies see the same free protein, and together they hold the same
  # sites as one copy at twice the concentration.
  s <- seqs3[["tA"]]
  fa <- write_fasta(c(copy1 = s, copy2 = s))
  for (fp in list(NULL, 5L)) {
    shared <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40), rna_conc = 1.5,
                                   k = 3, footprint = fp,
                                   shared_protein = TRUE)
    alone2c <- simulateBinding(s, models3, c(R1 = 20, R2 = 40), rna_conc = 3,
                               k = 3, footprint = fp)
    for (cp in c("copy1", "copy2")) {
      expect_equal(shared[transcript == cp, R1], alone2c$R1, tolerance = 1e-8)
      expect_equal(shared[transcript == cp, R2], alone2c$R2, tolerance = 1e-8)
    }
  }
})

test_that("a shared protein pool lowers occupancy for one RBP at limiting protein", {
  fa <- write_fasta(seqs3)
  m1 <- models3["R1"]
  for (fp in list(NULL, 5L)) {
    alone <- simulateBindingFasta(fa, m1, c(R1 = 5), rna_conc = 1, k = 3,
                                  footprint = fp)
    shared <- simulateBindingFasta(fa, m1, c(R1 = 5), rna_conc = 1, k = 3,
                                   footprint = fp, shared_protein = TRUE)
    expect_equal(shared[, .(pos, transcript)], alone[, .(pos, transcript)])
    expect_true(all(shared$R1 <= alone$R1 + 1e-12))
    expect_lt(mean(shared$R1 / alone$R1), 0.95)
  }
})

test_that("the lattice mass balance over many molecules matches enumeration", {
  # Enumerate every arrangement of one RBP on each small molecule to get its
  # expected number of bound proteins, then solve f + sum c_t N_t(f) = P.
  enum_n_bound <- function(q, w, L) {
    tot <- 0
    Z <- 0
    visit <- function(i, weight, n) {
      if (i > L) {
        Z <<- Z + weight
        tot <<- tot + weight * n
        return(invisible(NULL))
      }
      visit(i + 1L, weight, n)
      if (i <= length(q)) visit(i + w, weight * q[i], n + 1L)
    }
    visit(1L, 1, 0L)
    tot / Z
  }
  set.seed(7)
  w <- 3L
  Ls <- c(8L, 11L, 6L)
  conc <- c(2, 0.5, 4)
  a <- lapply(Ls, function(L) 10^runif(L - w + 1, -2, 0))
  P <- 6
  mols <- lapply(seq_along(Ls), function(t) list(a_list = list(a[[t]]), L = Ls[t],
                                                 rna_conc = conc[t]))
  sol <- solveLatticePoolMulti(mols, w, P)
  expect_true(sol$converged)
  h <- function(u) {
    f <- exp(u)
    f + sum(vapply(seq_along(Ls), function(t) conc[t] * enum_n_bound(f * a[[t]], w, Ls[t]),
                   numeric(1))) - P
  }
  f_ref <- exp(uniroot(h, c(log(1e-9), log(P)), tol = 1e-14)$root)
  expect_equal(sol$free, f_ref, tolerance = 1e-8)
  for (t in seq_along(Ls)) {
    expect_equal(sol$molecules[[t]]$n_bound,
                 enum_n_bound(f_ref * a[[t]], w, Ls[t]), tolerance = 1e-8)
  }
})

test_that("solveLatticePool equals solveLatticePoolMulti with one molecule", {
  set.seed(3)
  a <- list(10^runif(20, -2, 0), 10^runif(18, -2, 0))
  one <- solveLatticePool(a, c(5L, 7L), 24L, rna_conc = 2, protein_totals = c(3, 8))
  multi <- solveLatticePoolMulti(list(list(a_list = a, L = 24L, rna_conc = 2)),
                                 c(5L, 7L), c(3, 8))
  expect_identical(one$free, multi$free)
  expect_identical(one$starts, multi$molecules[[1]]$starts)
})

test_that("rna_conc can be set for each sequence, by name or in order", {
  fa <- write_fasta(seqs3)
  by_order <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40),
                                   rna_conc = c(1, 3, 0.5), k = 3, footprint = 5L,
                                   shared_protein = TRUE)
  by_name <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40),
                                  rna_conc = c(tC = 0.5, tA = 1, tB = 3), k = 3,
                                  footprint = 5L, shared_protein = TRUE)
  expect_equal(by_name, by_order)
  # Without a shared protein pool, each sequence uses its own value
  sep <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40),
                              rna_conc = c(tC = 0.5, tA = 1, tB = 3), k = 3)
  ref <- simulateBinding(seqs3[["tB"]], models3, c(R1 = 20, R2 = 40),
                         rna_conc = 3, k = 3)
  expect_equal(sep[transcript == "tB", R1], ref$R1)
})

test_that("bad rna_conc and shared_protein values give clear errors", {
  fa <- write_fasta(seqs3)
  m <- c(R1 = 20, R2 = 40)
  expect_error(simulateBindingFasta(fa, models3, m, rna_conc = c(1, 2), k = 3),
               "one value per sequence")
  expect_error(simulateBindingFasta(fa, models3, m, rna_conc = c(tA = 1, tB = 2), k = 3),
               "Missing: tC")
  expect_error(simulateBindingFasta(fa, models3, m, rna_conc = -1, k = 3),
               "non-negative")
  expect_error(simulateBindingFasta(fa, models3, m, k = 3, shared_protein = NA),
               "TRUE or FALSE")
  expect_error(simulateGridFasta(fa, models3, list(R1 = 1), 1, k = 3,
                                 shared_protein = "yes"), "TRUE or FALSE")
})

test_that("a sequence shorter than k is dropped with a warning", {
  fa <- write_fasta(c(seqs3[1], tiny = "AC"))
  expect_warning(
    res <- simulateBindingFasta(fa, models3, c(R1 = 20, R2 = 40), rna_conc = 1,
                                k = 3, shared_protein = TRUE),
    "returned NULL for tiny")
  expect_setequal(unique(res$transcript), "tA")
})

test_that("simulateGridFasta with a shared protein pool matches simulateBindingFasta", {
  fa <- write_fasta(seqs3)
  grid <- suppressMessages(simulateGridFasta(
    fa, models3, protein_conc_grid = list(R1 = c(5, 50), R2 = 10),
    rna_conc_grid = c(0.5, 2), k = 3, footprint = 5L, shared_protein = TRUE))
  expect_equal(nrow(grid), 2 * 2 * sum(nchar(seqs3)))
  one <- simulateBindingFasta(fa, models3, c(R1 = 50, R2 = 10), rna_conc = 2,
                              k = 3, footprint = 5L, shared_protein = TRUE)
  pt <- grid[Conc_R1 == 50 & rna_conc == 2]
  expect_equal(pt$R1, one$R1)
  expect_equal(pt$transcript, one$transcript)
})

test_that("seq_type reads DNA and RNA FASTA files to the same result", {
  dna <- write_fasta(seqs3)
  rna <- tempfile(fileext = ".fa")
  writeLines(paste0(">", names(seqs3), "\n", seqs3), rna)
  a <- simulateBindingFasta(dna, models3, c(R1 = 20, R2 = 40), rna_conc = 2, k = 3)
  b <- simulateBindingFasta(rna, models3, c(R1 = 20, R2 = 40), rna_conc = 2, k = 3,
                            seq_type = "RNA")
  expect_equal(a, b)
  g1 <- simulateGridFasta(dna, models3, list(R1 = 20, R2 = 40), rna_conc_grid = 2, k = 3)
  g2 <- simulateGridFasta(rna, models3, list(R1 = 20, R2 = 40), rna_conc_grid = 2, k = 3,
                          seq_type = "RNA")
  expect_equal(g1, g2)
  # A letter outside the chosen alphabet stops with an error, not a silent drop
  expect_error(simulateBindingFasta(rna, models3, c(R1 = 20, R2 = 40), k = 3),
               "seq_type = \"RNA\"")
  expect_error(simulateBindingFasta(dna, models3, c(R1 = 20, R2 = 40), k = 3,
                                    seq_type = "RNA"), "seq_type = \"DNA\"")
  expect_error(simulateBindingFasta(dna, models3, c(R1 = 20, R2 = 40), k = 3,
                                    seq_type = "rna"))
})
