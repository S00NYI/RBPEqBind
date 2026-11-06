#' Simulate Competitive Binding
#'
#' Simulates competitive binding of multiple RBPs to a single RNA sequence.
#'
#' Each k-mer window of the sequence is a binding site. All sites draw on one
#' shared supply of each protein. For free protein concentrations \eqn{f_j},
#' the fraction of site \eqn{i} bound by RBP \eqn{j} is
#' \deqn{\theta_{ij} = \frac{f_j / K_{d,ij}}{1 + \sum_m f_m / K_{d,im}}.}
#' The free concentrations satisfy the mass balance
#' \eqn{P_j = f_j + c_{RNA} \sum_i \theta_{ij}}, summed over all windows, where
#' \eqn{P_j} is the total concentration of RBP \eqn{j}. Sites bind
#' independently of each other. The occupancy at a nucleotide is the mean of
#' \eqn{\theta} over the windows that cover it.
#'
#' With \code{footprint}, the RNA is a lattice instead. A bound RBP covers
#' \code{footprint} consecutive nucleotides, and two bound proteins may not
#' overlap. A placement has the weight \eqn{f_j / K_d}, where \eqn{K_d} belongs
#' to the k-mer in the middle of the footprint. The model sums over all
#' allowed arrangements on one molecule (forward and backward passes) and
#' solves the same shared protein balance, with bound protein equal to
#' \code{rna_conc} times the expected number of bound proteins per molecule.
#' The occupancy at a nucleotide is then the probability that the nucleotide
#' is covered by the RBP.
#'
#' With \code{omega}, two proteins of the same RBP whose footprints touch
#' (the second starts right after the first ends) give their arrangement the
#' extra factor \code{omega}. Values above 1 make neighbors help each other,
#' values below 1 make them hinder each other, and 1 means no effect.
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of data.frames/data.tables with 'motif' and 'Kd'.
#' @param protein_concs Named numeric vector of total protein concentrations.
#' @param rna_conc Concentration of the RNA molecule (default 10). Each window
#'   is one site at this concentration.
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) for independent windows, or the
#'   number of nucleotides a bound RBP covers. Either one integer for all RBPs
#'   or a named integer vector with one value per RBP. Each value must be at
#'   least \code{k}.
#' @param omega Neighbor factor (default 1, no effect). Either one number for
#'   all RBPs or a named vector with one value per RBP. Each value must be
#'   \code{>= 0}. Values other than 1 need a \code{footprint}.
#' @return A data.table with one row per nucleotide and these columns:
#'   \describe{
#'     \item{\code{pos}, \code{nt}}{Position and nucleotide.}
#'     \item{\code{<RBP>}}{Occupancy, from 0 to 1. Without a footprint: the
#'       mean binding probability of the windows that cover the nucleotide.
#'       With a footprint: the probability that the nucleotide is covered by
#'       the RBP. Use this column for quantitative analysis.}
#'     \item{\code{<RBP>_density}}{Occupancy scaled to sum to 1 over the
#'       sequence.}
#'     \item{\code{<RBP>_density_fc}, \code{<RBP>_occupancy_fc}}{The meaning
#'       depends on the number of active RBPs (those with a concentration
#'       above 0). With one active RBP: fold change over a uniform track
#'       (density times sequence length, or occupancy over its mean). With two
#'       or more: the value of this RBP divided by the sum of the other active
#'       RBPs at the same nucleotide, \code{NA} where that sum is 0. Do not
#'       compare these columns between runs with different numbers of active
#'       RBPs.}
#'   }
#'   To measure the effect of a competitor, run \code{simulateBinding()} with
#'   and without it, and compare the \code{<RBP>} occupancy columns.
#' @examples
#' # Load and process model
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' raw_models <- loadModel(model_file, rbp = c("HH", "HL"))
#' rbp_models <- setModel(raw_models, max_affinity = 100)
#'
#' # Simulate binding on a short sequence
#' seq <- "ACGUACGUACGUACGUACGU"
#' results <- simulateBinding(seq, rbp_models, c(HH = 100, HL = 100), rna_conc = 10)
#' head(results)
#' @importFrom data.table as.data.table setkey rbindlist
#' @export
simulateBinding <- function(sequence, rbp_models, protein_concs, rna_conc = 10,
                            k = 5, footprint = NULL, omega = 1) {
  

  if (!is.character(sequence) || length(sequence) != 1L) {
    stop("sequence must be a single character string.")
  }
  if (!is.numeric(rna_conc) || length(rna_conc) != 1L || is.na(rna_conc) || rna_conc < 0) {
    stop("rna_conc must be a single non-negative number.")
  }
  om_vec <- checkProteinInputs(rbp_models, protein_concs, footprint, omega)
  rbp_names <- names(protein_concs)

  sites <- prepareKmerSites(sequence, rbp_models, rbp_names, k)
  if (is.null(sites)) return(NULL)
  p0_vec <- protein_concs[rbp_names]

  if (is.null(footprint)) {
    # Independent windows. All windows draw on one supply of each protein.
    # Each unique k-mer is a site species at rna_conc times its copy number.
    pool <- solveSharedPool(rna_conc * sites$count, sites$kd_mat, p0_vec)
    final_mat <- windowTrack(pool$theta[sites$kmer_index, , drop = FALSE],
                             sites$L, k)
  } else {
    # Lattice: a bound RBP covers `w` nucleotides and bound proteins may not
    # overlap. The k-mer in the middle of the footprint sets the affinity.
    w_vec <- resolveFootprint(footprint, rbp_names, k)
    pool <- solveLatticePool(latticeAffinity(sites, w_vec, k), w_vec,
                             sites$L, rna_conc, p0_vec, omega = om_vec)
    final_mat <- latticeTrack(pool$starts, w_vec, sites$L)
  }
  buildTrackTable(sequence, final_mat, rbp_names, protein_concs)
}

# ==============================================================================
# Internal steps shared by simulateBinding() and the shared-supply runs
# ==============================================================================

# Result Table for One Sequence
#
# @param sequence RNA sequence string.
# @param final_mat Occupancy matrix, one row per nucleotide and one column per
#   RBP.
# @param rbp_names Names of the simulated RBPs.
# @param protein_concs Named protein concentrations (sets the active RBPs).
# @return The data.table that simulateBinding() returns.
buildTrackTable <- function(sequence, final_mat, rbp_names, protein_concs) {
  L <- nchar(sequence)
  seq_chars <- strsplit(sequence, "")[[1]]
  
  res <- data.table::data.table(
    pos = seq_len(L),
    nt = seq_chars
  )
  
  res <- cbind(res, data.table::as.data.table(final_mat))
  colnames(res)[3:ncol(res)] <- rbp_names
  
  # Add density column for each RBP
  for (rbp in rbp_names) {
    occ_col <- res[[rbp]]
    sum_occ <- sum(occ_col, na.rm = TRUE)
    res[[paste0(rbp, "_density")]] <- if (sum_occ > 0) occ_col / sum_occ else 0
  }
  
  # Calculate FC based on number of active RBPs
  # Active RBPs are determined by protein_concs > 0
  active_rbps <- names(protein_concs)[protein_concs > 0]
  n_active <- length(active_rbps)
  L <- nrow(res)  # sequence length
  
  for (rbp in rbp_names) {
    occ_col <- res[[rbp]]
    dens_col <- res[[paste0(rbp, "_density")]]
    
    if (n_active <= 1) {
      # Case 1: Single RBP - fold change over uniform distribution
      # fc = density / (1/L) = density * L
      # For occupancy: fc = occupancy / mean(occupancy)
      mean_occ <- mean(occ_col, na.rm = TRUE)
      res[[paste0(rbp, "_density_fc")]] <- dens_col * L
      res[[paste0(rbp, "_occupancy_fc")]] <- if (mean_occ > 0) occ_col / mean_occ else NA_real_
    } else {
      # Case 2: Multiple RBPs - ratio vs sum of other RBPs
      # Only consider other ACTIVE RBPs
      other_active <- setdiff(active_rbps, rbp)
      
      if (length(other_active) > 0) {
        # Sum of other active RBP densities per position
        other_dens_sum <- rowSums(res[, paste0(other_active, "_density"), with = FALSE], na.rm = TRUE)
        other_occ_sum <- rowSums(res[, other_active, with = FALSE], na.rm = TRUE)
        
        # fc = self / sum(others), NA if denominator is 0
        res[[paste0(rbp, "_density_fc")]] <- ifelse(other_dens_sum > 0, dens_col / other_dens_sum, NA_real_)
        res[[paste0(rbp, "_occupancy_fc")]] <- ifelse(other_occ_sum > 0, occ_col / other_occ_sum, NA_real_)
      } else {
        # Edge case: this RBP is the only active one (shouldn't happen in n_active > 1)
        res[[paste0(rbp, "_density_fc")]] <- dens_col * L
        mean_occ <- mean(occ_col, na.rm = TRUE)
        res[[paste0(rbp, "_occupancy_fc")]] <- if (mean_occ > 0) occ_col / mean_occ else NA_real_
      }
    }
  }
  
  return(res)
}


# Check the Protein Inputs of a Simulation
#
# @return The omega value of each RBP (see resolveOmega()).
checkProteinInputs <- function(rbp_models, protein_concs, footprint, omega) {
  if (!all(names(protein_concs) %in% names(rbp_models))) {
    stop("All RBPs in protein_concs must have a corresponding model in rbp_models")
  }
  if (!is.numeric(protein_concs) || anyNA(protein_concs) || any(protein_concs < 0)) {
    stop("protein_concs must be non-negative numbers.")
  }
  om_vec <- resolveOmega(omega, names(protein_concs))
  if (is.null(footprint) && any(om_vec != 1)) {
    stop("omega acts between neighboring proteins and needs a footprint.")
  }
  om_vec
}

# K-mer Sites of One Sequence
#
# @param sequence RNA sequence string.
# @param rbp_models Named list of RBP models.
# @param rbp_names Names of the simulated RBPs.
# @param k K-mer size.
# @return NULL (with a warning) when the sequence has no k-mers. Otherwise a
#   list with `kmer_index` (the unique k-mer of each window), `count` (copies
#   of each unique k-mer), `kd_mat` (Kd of each unique k-mer, one column per
#   RBP) and `L` (sequence length).
prepareKmerSites <- function(sequence, rbp_models, rbp_names, k) {
  kmers <- extractKmers(sequence, k)
  if (length(kmers) == 0) {
    warning("No k-mers found for sequence length ", nchar(sequence), " k = ", k)
    return(NULL)
  }

  binding_params <- countKmers(kmers)
  if (nrow(binding_params) == 0) {
    warning("Unique k-mers count is 0")
    return(NULL)
  }

  for (rbp in rbp_names) {
    model <- rbp_models[[rbp]]
    if (!is.data.frame(model) || !all(c("motif", "Kd") %in% colnames(model))) {
      stop("Model for ", rbp, " must be a dataframe with 'motif' and 'Kd' columns")
    }

    model_k <- unique(nchar(model$motif))
    if (length(model_k) != 1L || model_k != k) {
      stop("Model for ", rbp, " uses ", paste(model_k, collapse = "/"),
           "-mers but k = ", k, ".")
    }

    match_idx <- match(binding_params$motif, model$motif)
    if (anyNA(match_idx)) {
      missing <- head(binding_params$motif[is.na(match_idx)])
      warning("Some motifs in sequence not found in model for ", rbp, ": ", paste(missing, collapse=","), "... using Inf Kd (no binding)")
      kds <- rep(Inf, nrow(binding_params))
      kds[!is.na(match_idx)] <- model$Kd[match_idx[!is.na(match_idx)]]
    } else {
      kds <- model$Kd[match_idx]
    }

    binding_params[[paste0("Kd_", rbp)]] <- kds
  }

  kd_mat <- do.call(cbind, lapply(rbp_names, function(rbp) {
    binding_params[[paste0("Kd_", rbp)]]
  }))
  colnames(kd_mat) <- rbp_names
  list(kmer_index = match(kmers, binding_params$motif),
       count = binding_params$count, kd_mat = kd_mat, L = nchar(sequence))
}

# Nucleotide Occupancy From Window Occupancy
#
# @param window_occupancy Matrix of window occupancies, one row per window and
#   one column per RBP.
# @param L Sequence length.
# @param k K-mer size.
# @return Matrix with one row per nucleotide: the mean occupancy of the
#   windows that cover it.
windowTrack <- function(window_occupancy, L, k) {
  W <- nrow(window_occupancy)
  n_rbps <- ncol(window_occupancy)
  # Mean over the windows that cover each nucleotide (cumsum trick)
  diff_mat <- matrix(0, nrow = L + 1, ncol = n_rbps)
  diff_mat[seq_len(W), ] <- diff_mat[seq_len(W), ] + window_occupancy
  diff_mat[(k + 1):(L + 1), ] <- diff_mat[(k + 1):(L + 1), ] - window_occupancy
  pos_sums <- apply(diff_mat, 2, cumsum)[seq_len(L), , drop = FALSE]

  count_diff <- numeric(L + 1)
  count_diff[seq_len(W)] <- count_diff[seq_len(W)] + 1
  count_diff[(k + 1):(L + 1)] <- count_diff[(k + 1):(L + 1)] - 1
  pos_counts <- cumsum(count_diff)[seq_len(L)]
  pos_counts[pos_counts == 0] <- 1

  pos_sums / pos_counts
}

# Placement Affinities on the Lattice
#
# @param sites Result of prepareKmerSites().
# @param w_vec Footprint of each RBP.
# @param k K-mer size.
# @return List with one vector per RBP: 1 / Kd for each placement start. The
#   k-mer in the middle of the footprint sets the affinity.
latticeAffinity <- function(sites, w_vec, k) {
  L <- sites$L
  lapply(seq_along(w_vec), function(j) {
    kd_win <- sites$kd_mat[sites$kmer_index, j]
    m <- L - w_vec[j] + 1L
    if (m < 1L) return(numeric(0))
    offset <- (w_vec[j] - k) %/% 2L
    1 / kd_win[seq_len(m) + offset]
  })
}

# Nucleotide Occupancy on the Lattice
#
# @param starts List of start-probability vectors, one per RBP.
# @param w_vec Footprint of each RBP.
# @param L Sequence length.
# @return Matrix with one row per nucleotide: the probability that the
#   nucleotide is covered by each RBP.
latticeTrack <- function(starts, w_vec, L) {
  n_rbps <- length(w_vec)
  final_mat <- vapply(seq_len(n_rbps), function(j) {
    latticeCoverage(starts[[j]], w_vec[j], L)
  }, numeric(L))
  matrix(final_mat, nrow = L, ncol = n_rbps)
}

#' Simulate Concentration Grid
#'
#' Runs simulations across a grid of protein/RNA concentrations.
#' Concurrency is controlled by the user's active \code{\link[future]{plan}} (default is sequential).
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration.
#' @param rna_conc_grid Numeric vector of RNA concentrations to sweep.
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) or footprint sizes; see
#'   \code{\link{simulateBinding}}.
#' @param omega Neighbor factor (default 1); see \code{\link{simulateBinding}}.
#' @return A data.table with results for all grid combinations.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' grid <- simulateGrid("ACGUACGU", rbp_models,
#'   protein_conc_grid = list(HH = c(10, 100), HL = c(10, 100)),
#'   rna_conc_grid = 10)
#' head(grid)
#' @importFrom future.apply future_lapply
#' @export
simulateGrid <- function(sequence, rbp_models, protein_conc_grid, rna_conc_grid,
                         k = 5, footprint = NULL, omega = 1) {
  
  # Validate: sequence must be a single string
  if (!is.character(sequence) || length(sequence) != 1) {
    stop("sequence must be a single character string. For FASTA files, use simulateGridFasta().")
  }
  
  # Expand protein concentration grid
  prot_grid <- do.call(expand.grid, protein_conc_grid)
  
  # Create full grid: RNA conc × protein concentrations
  grid_args <- list(
    rna_conc = rna_conc_grid,
    ProteinGridRow = seq_len(nrow(prot_grid))
  )
  full_grid <- do.call(expand.grid, grid_args)
  n_jobs <- nrow(full_grid)
  
  # Job function for one grid point
  run_job <- function(i) {
    row <- full_grid[i, ]
    rna_c <- row$rna_conc
    
    # Extract protein concentrations
    prot_row_idx <- row$ProteinGridRow
    prot_concs_df <- prot_grid[prot_row_idx, , drop = FALSE]
    prot_concs <- unlist(prot_concs_df)
    names(prot_concs) <- colnames(prot_grid)
    
    # Run simulation
    res <- simulateBinding(sequence, rbp_models, prot_concs, rna_conc = rna_c,
                           k = k, footprint = footprint, omega = omega)
    
    if (is.null(res)) return(NULL)
    
    # Add metadata
    res$rna_conc <- rna_c
    for (p in names(prot_concs)) {
      res[[paste0("Conc_", p)]] <- prot_concs[[p]]
    }
    
    return(res)
  }
  
  # Execute via future.apply (respects user-configured future::plan)
  results_list <- future.apply::future_lapply(seq_len(n_jobs), run_job, 
                                              future.seed = TRUE, 
                                              future.packages = c("data.table"))
  
  # Combine results
  final_dt <- data.table::rbindlist(results_list)
  return(final_dt)
}


# ==============================================================================
# FASTA Input Wrappers
# ==============================================================================

#' Simulate Binding from FASTA File(s)
#'
#' Simulates competitive binding for all sequences in one or more FASTA files.
#' Concurrency is controlled by the user's active \code{\link[future]{plan}} (default is sequential).
#'
#' By default (\code{shared_supply = FALSE}) each sequence is simulated alone,
#' with the full protein supply, as with \code{\link{simulateBinding}}. With
#' \code{shared_supply = TRUE}, all sequences are in one tube and draw on one
#' supply of each protein. The model then solves one mass balance,
#' \eqn{P_j = f_j + \sum_t c_t \langle N_j \rangle_t}, where \eqn{c_t} is the
#' concentration of sequence \eqn{t} and \eqn{\langle N_j \rangle_t} is its
#' expected number of bound proteins of RBP \eqn{j} (summed over windows
#' without a footprint). Each sequence keeps its own windows or lattice, so no
#' placement crosses from one sequence to the next. All sequences see the
#' same free protein level \eqn{f_j}. This is the setup of an experiment that
#' titrates protein into a mix of RNAs, and it assumes that the mix is well
#' mixed and at equilibrium.
#'
#' @param fasta_file Character vector of path(s) to FASTA files.
#' @param rbp_models Named list of RBP models.
#' @param protein_concs Named numeric vector of protein concentrations.
#' @param rna_conc RNA concentration (default 10). Either one value for all
#'   sequences, or one value per sequence: a vector named with the sequence
#'   names, or an unnamed vector in file order.
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) or footprint sizes; see
#'   \code{\link{simulateBinding}}.
#' @param omega Neighbor factor (default 1); see \code{\link{simulateBinding}}.
#' @param shared_supply If \code{TRUE}, all sequences draw on one supply of
#'   each protein. If \code{FALSE} (default), each sequence is simulated alone.
#' @return A data.table with combined results, including a transcript column.
#'   The per-sequence columns are the same as for
#'   \code{\link{simulateBinding}}. The density and fold-change columns are
#'   computed within each sequence.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' fasta_file <- system.file("extdata", "test_transcripts.fa", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBindingFasta(fasta_file, rbp_models, c(HH = 100, HL = 100))
#' head(results)
#'
#' # All transcripts in one tube, sharing one supply of each protein
#' pooled <- simulateBindingFasta(fasta_file, rbp_models, c(HH = 100, HL = 100),
#'   rna_conc = 1, shared_supply = TRUE)
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table rbindlist
#' @importFrom future.apply future_lapply
#' @export
simulateBindingFasta <- function(fasta_file, rbp_models, protein_concs, rna_conc = 10,
                                 k = 5, footprint = NULL, omega = 1,
                                 shared_supply = FALSE) {
  checkSharedSupply(shared_supply)
  seqs <- readFastaSequences(fasta_file)
  rna <- resolveRnaConc(rna_conc, names(seqs))

  if (shared_supply) {
    return(simulateSharedSupply(seqs, rbp_models, protein_concs, rna, k,
                                footprint = footprint, omega = omega))
  }

  run_seq <- function(i) {
    res <- simulateBinding(seqs[[i]], rbp_models, protein_concs, rna[[i]], k,
                           footprint = footprint, omega = omega)
    if (is.null(res)) {
      warning("Simulation returned NULL for ", names(seqs)[i])
      return(NULL)
    }
    res$transcript <- names(seqs)[i]
    return(res)
  }

  # Execute via future.apply (respects user-configured future::plan)
  results_list <- future.apply::future_lapply(seq_along(seqs), run_seq,
                                              future.seed = TRUE,
                                              future.packages = c("data.table", "Biostrings"))

  data.table::rbindlist(results_list)
}


#' Simulate Grid from FASTA File(s)
#'
#' Simulates competitive binding for all sequences in one or more FASTA files across a concentration grid.
#' Concurrency is controlled by the user's active \code{\link[future]{plan}} (default is sequential).
#'
#' With \code{shared_supply = TRUE}, each grid point puts all sequences in
#' one tube, each at the RNA concentration of that grid point, and they draw
#' on one supply of each protein (see \code{\link{simulateBindingFasta}}).
#' The work is then split over grid points instead of sequences.
#'
#' @param fasta_file Character vector of path(s) to FASTA files.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration.
#' @param rna_conc_grid Numeric vector of RNA concentrations to sweep. Each
#'   value applies to every sequence.
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) or footprint sizes; see
#'   \code{\link{simulateBinding}}.
#' @param omega Neighbor factor (default 1); see \code{\link{simulateBinding}}.
#' @param shared_supply If \code{TRUE}, all sequences draw on one supply of
#'   each protein. If \code{FALSE} (default), each sequence is simulated alone.
#' @return A data.table with combined results, including a transcript column.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' fasta_file <- system.file("extdata", "test_transcripts.fa", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateGridFasta(
#'   fasta_file = fasta_file,
#'   rbp_models = rbp_models,
#'   protein_conc_grid = list(HH = c(10, 100), HL = c(10, 100)),
#'   rna_conc_grid = 10,
#'   k = 5
#' )
#' head(results)
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table rbindlist
#' @importFrom future.apply future_lapply
#' @export
simulateGridFasta <- function(fasta_file, rbp_models, protein_conc_grid, rna_conc_grid,
                              k = 5, footprint = NULL, omega = 1,
                              shared_supply = FALSE) {
  checkSharedSupply(shared_supply)
  seqs <- readFastaSequences(fasta_file)
  n_seqs <- length(seqs)

  if (shared_supply) {
    prot_grid <- do.call(expand.grid, protein_conc_grid)
    full_grid <- expand.grid(rna_conc = rna_conc_grid,
                             ProteinGridRow = seq_len(nrow(prot_grid)))
    message("Running grid simulation on ", n_seqs,
            " sequences with one shared protein supply...")
    run_point <- function(i) {
      rna_c <- full_grid$rna_conc[i]
      prot_concs <- unlist(prot_grid[full_grid$ProteinGridRow[i], , drop = FALSE])
      names(prot_concs) <- colnames(prot_grid)
      res <- simulateSharedSupply(seqs, rbp_models, prot_concs,
                                  rep(rna_c, n_seqs), k,
                                  footprint = footprint, omega = omega)
      if (is.null(res)) return(NULL)
      res$rna_conc <- rna_c
      for (p in names(prot_concs)) {
        res[[paste0("Conc_", p)]] <- prot_concs[[p]]
      }
      res
    }
    results_list <- future.apply::future_lapply(seq_len(nrow(full_grid)), run_point,
                                                future.seed = TRUE,
                                                future.packages = c("data.table"))
    return(data.table::rbindlist(results_list))
  }

  message("Running grid simulation on ", n_seqs, " sequences...")

  run_seq_grid <- function(i) {
    res <- simulateGrid(seqs[[i]], rbp_models, protein_conc_grid, rna_conc_grid,
                        k = k, footprint = footprint, omega = omega)
    if (!is.null(res)) {
      res$transcript <- names(seqs)[i]
    }
    return(res)
  }

  # Execute via future.apply (respects user-configured future::plan)
  results_list <- future.apply::future_lapply(seq_len(n_seqs), run_seq_grid,
                                              future.seed = TRUE,
                                              future.packages = c("data.table", "Biostrings"))

  data.table::rbindlist(results_list)
}


# ==============================================================================
# Shared protein supply over many sequences
# ==============================================================================

# Read Sequences From FASTA File(s)
#
# @param fasta_file Character vector of FASTA paths.
# @return Named character vector of RNA sequences (T changed to U).
readFastaSequences <- function(fasta_file) {
  for (f in fasta_file) {
    if (!file.exists(f)) stop("FASTA file not found: ", f)
  }
  seqs_dna <- do.call(c, lapply(fasta_file, Biostrings::readDNAStringSet))
  seq_names <- names(seqs_dna)
  if (is.null(seq_names) || any(seq_names == "")) {
    warning("Sequence names missing, using Tx indices")
    seq_names <- paste0("Tx", seq_along(seqs_dna))
  }
  seqs <- gsub("T", "U", as.character(seqs_dna), fixed = TRUE)
  stats::setNames(unname(seqs), seq_names)
}

# Check the shared_supply Argument
checkSharedSupply <- function(shared_supply) {
  if (!is.logical(shared_supply) || length(shared_supply) != 1L ||
      is.na(shared_supply)) {
    stop("shared_supply must be TRUE or FALSE.")
  }
}

# RNA Concentration of Each Sequence
#
# @param rna_conc One value, a vector named with the sequence names, or an
#   unnamed vector with one value per sequence in order.
# @param seq_names Sequence names.
# @return Numeric vector in the order of seq_names.
resolveRnaConc <- function(rna_conc, seq_names) {
  if (!is.numeric(rna_conc) || anyNA(rna_conc) || any(rna_conc < 0)) {
    stop("rna_conc must be non-negative numbers.")
  }
  n <- length(seq_names)
  if (length(rna_conc) == 1L && is.null(names(rna_conc))) {
    return(rep(as.numeric(rna_conc), n))
  }
  if (!is.null(names(rna_conc))) {
    if (anyDuplicated(seq_names)) {
      stop("Sequence names are not unique, so rna_conc cannot be matched by name.")
    }
    if (!all(seq_names %in% names(rna_conc))) {
      stop("A named rna_conc must have one value for each sequence. Missing: ",
           paste(head(setdiff(seq_names, names(rna_conc))), collapse = ", "))
    }
    return(as.numeric(rna_conc[seq_names]))
  }
  if (length(rna_conc) != n) {
    stop("An unnamed rna_conc must have one value or one value per sequence (",
         n, ").")
  }
  as.numeric(rna_conc)
}

# Simulate Many Sequences With One Shared Protein Supply
#
# @param seqs Named character vector of RNA sequences.
# @param rbp_models Named list of RBP models.
# @param protein_concs Named protein concentrations.
# @param rna Numeric vector of RNA concentrations, one per sequence.
# @param k K-mer size.
# @param footprint,omega See simulateBinding().
# @return A data.table as from simulateBindingFasta(), or NULL when no
#   sequence has a k-mer.
simulateSharedSupply <- function(seqs, rbp_models, protein_concs, rna, k,
                                 footprint = NULL, omega = 1) {
  om_vec <- checkProteinInputs(rbp_models, protein_concs, footprint, omega)
  rbp_names <- names(protein_concs)
  p0_vec <- protein_concs[rbp_names]

  sites <- lapply(seqs, prepareKmerSites, rbp_models = rbp_models,
                  rbp_names = rbp_names, k = k)
  keep <- which(!vapply(sites, is.null, logical(1)))
  for (i in setdiff(seq_along(seqs), keep)) {
    warning("Simulation returned NULL for ", names(seqs)[i])
  }
  if (length(keep) == 0L) return(NULL)

  if (is.null(footprint)) {
    # Every unique k-mer of every sequence is one site species, at the
    # concentration of its sequence times its copy number.
    site_conc <- unlist(lapply(keep, function(i) rna[i] * sites[[i]]$count))
    kd_mat <- do.call(rbind, lapply(keep, function(i) sites[[i]]$kd_mat))
    pool <- solveSharedPool(site_conc, kd_mat, p0_vec)
    n_unique <- vapply(keep, function(i) length(sites[[i]]$count), integer(1))
    first_row <- cumsum(c(0L, utils::head(n_unique, -1L)))
    tracks <- lapply(seq_along(keep), function(t) {
      s <- sites[[keep[t]]]
      windowTrack(pool$theta[first_row[t] + s$kmer_index, , drop = FALSE],
                  s$L, k)
    })
  } else {
    w_vec <- resolveFootprint(footprint, rbp_names, k)
    molecules <- lapply(keep, function(i) {
      list(a_list = latticeAffinity(sites[[i]], w_vec, k), L = sites[[i]]$L,
           rna_conc = rna[i])
    })
    pool <- solveLatticePoolMulti(molecules, w_vec, p0_vec, omega = om_vec)
    tracks <- lapply(seq_along(keep), function(t) {
      latticeTrack(pool$molecules[[t]]$starts, w_vec, sites[[keep[t]]]$L)
    })
  }

  data.table::rbindlist(lapply(seq_along(keep), function(t) {
    res <- buildTrackTable(seqs[[keep[t]]], tracks[[t]], rbp_names, protein_concs)
    res$transcript <- names(seqs)[keep[t]]
    res
  }))
}
