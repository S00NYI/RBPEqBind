#' Simulate Competitive Binding
#'
#' Simulates competitive binding of multiple RBPs to a single RNA sequence.
#'
#' Each k-mer window of the sequence is a binding site. All sites draw on one
#' shared protein pool of each protein. For free protein concentrations \eqn{f_j},
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
#' With \code{accessibility}, RNA structure can hide binding sites. For a
#' k-mer window that starts at \eqn{s}, let \eqn{a(s)} be the probability
#' that its \code{k} nucleotides are all unpaired. A protein binds only an
#' unpaired window, with the same \eqn{K_d}. The model assumes that each RNA
#' molecule switches between paired and unpaired at that window, and that a
#' bound protein holds the window unpaired. Then the window has three states,
#' with the weights \eqn{1 - a} (paired), \eqn{a} (unpaired, free) and
#' \eqn{a f_j / K_{d,j}} (unpaired, bound by RBP \eqn{j}), and
#' \deqn{\theta_{ij} = \frac{a f_j / K_{d,ij}}{1 + \sum_m a f_m / K_{d,im}}.}
#' So the binding weight of the window, \eqn{f_j / K_d}, is multiplied by
#' \eqn{a(s)}. \eqn{K_d} itself does not change. With a footprint, the weight
#' of each placement is multiplied by \eqn{a} of the k-mer in the middle of
#' the footprint (the k-mer that sets the affinity). The same \eqn{a(s)}
#' applies to every RBP. The model does not include the effect of a bound
#' protein on the structure of other windows. RBPEqBind does not fold RNA:
#' compute the probabilities with a folding program such as ViennaRNA
#' \code{RNAplfold -u} and read them with \code{\link{readAccessibility}}.
#' For a model in which binding and folding are computed together, see
#' RBPBind (Gaither et al. 2022, J Mol Biol).
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
#' @param accessibility \code{NULL} (default, RNA fully accessible) or a
#'   numeric matrix with one row per nucleotide of \code{sequence}. Column
#'   \eqn{l} is the probability that the \eqn{l} nucleotides ending at that
#'   row are all unpaired (the layout of an RNAplfold \code{*_lunp} file; see
#'   \code{\link{readAccessibility}}). It needs at least \code{k} columns.
#'   Values must be between 0 and 1.
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
                            k = 5, footprint = NULL,
                            accessibility = NULL) {
  

  if (!is.character(sequence) || length(sequence) != 1L) {
    stop("sequence must be a single character string.")
  }
  if (!is.numeric(rna_conc) || length(rna_conc) != 1L || is.na(rna_conc) || rna_conc < 0) {
    stop("rna_conc must be a single non-negative number.")
  }
  checkProteinInputs(rbp_models, protein_concs)
  rbp_names <- names(protein_concs)
  w_vec <- if (is.null(footprint)) NULL else resolveFootprint(footprint, rbp_names, k)
  acc <- checkAccessibility(accessibility, nchar(sequence), k)

  sites <- prepareKmerSites(sequence, rbp_models, rbp_names, k)
  if (is.null(sites)) return(NULL)
  p0_vec <- protein_concs[rbp_names]

  if (is.null(footprint)) {
    # Independent windows. All windows draw on one pool of each protein.
    ws <- windowSites(sites, rna_conc, acc, k)
    pool <- solveSharedPool(ws$site_conc, ws$kd_mat, p0_vec)
    final_mat <- windowTrack(pool$theta[ws$index, , drop = FALSE], sites$L, k)
  } else {
    # Lattice: a bound RBP covers `w` nucleotides and bound proteins may not
    # overlap. The k-mer in the middle of the footprint sets the affinity.
    pool <- solveLatticePool(latticeAffinity(sites, w_vec, k, acc), w_vec,
                             sites$L, rna_conc, p0_vec)
    final_mat <- latticeTrack(pool$starts, w_vec, sites$L)
  }
  buildTrackTable(sequence, final_mat, rbp_names, protein_concs)
}

# Window Sites for the Model Without a Footprint
#
# Without accessibility, each unique k-mer is one site species at rna_conc
# times its copy number. With accessibility, each window is its own site
# species at rna_conc, with Kd / a(s): the binding weight of the window is
# multiplied by a(s), the probability that the window is unpaired.
#
# @param sites Result of prepareKmerSites().
# @param rna_conc RNA concentration.
# @param acc NULL or a checked accessibility matrix.
# @param k K-mer size.
# @return List with `site_conc`, `kd_mat` and `index` (the site row of each
#   window).
windowSites <- function(sites, rna_conc, acc, k) {
  if (is.null(acc)) {
    return(list(site_conc = rna_conc * sites$count, kd_mat = sites$kd_mat,
                index = sites$kmer_index))
  }
  n_win <- length(sites$kmer_index)
  list(site_conc = rep(rna_conc, n_win),
       kd_mat = sites$kd_mat[sites$kmer_index, , drop = FALSE] /
         placementAccessibility(acc, k),
       index = seq_len(n_win))
}

# ==============================================================================
# Internal steps shared by simulateBinding() and the shared-protein runs
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
# @return NULL (invisibly); stops on invalid input.
checkProteinInputs <- function(rbp_models, protein_concs) {
  if (!all(names(protein_concs) %in% names(rbp_models))) {
    stop("All RBPs in protein_concs must have a corresponding model in rbp_models")
  }
  if (!is.numeric(protein_concs) || anyNA(protein_concs) || any(protein_concs < 0)) {
    stop("protein_concs must be non-negative numbers.")
  }
  invisible(NULL)
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
# @param acc NULL, or a checked accessibility matrix. The weight of each
#   placement is then multiplied by the probability that the k-mer in the
#   middle of the footprint is unpaired.
# @return List with one vector per RBP: 1 / Kd for each placement start
#   (times the unpaired probability). The k-mer in the middle of the
#   footprint sets the affinity.
latticeAffinity <- function(sites, w_vec, k, acc = NULL) {
  L <- sites$L
  a_k <- if (is.null(acc)) NULL else placementAccessibility(acc, k)
  lapply(seq_along(w_vec), function(j) {
    kd_win <- sites$kd_mat[sites$kmer_index, j]
    m <- L - w_vec[j] + 1L
    if (m < 1L) return(numeric(0))
    idx <- seq_len(m) + (w_vec[j] - k) %/% 2L
    a <- 1 / kd_win[idx]
    if (!is.null(a_k)) a <- a * a_k[idx]
    a
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
#' @param accessibility \code{NULL} (default) or an accessibility matrix for
#'   \code{sequence}; see \code{\link{simulateBinding}}.
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
                         k = 5, footprint = NULL,
                         accessibility = NULL) {
  
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
                           k = k, footprint = footprint,
                           accessibility = accessibility)
    
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
#' By default (\code{shared_protein = FALSE}) each sequence is simulated alone,
#' with the full protein amount, as with \code{\link{simulateBinding}}. With
#' \code{shared_protein = TRUE}, all sequences are in one tube and draw on one
#' pool of each protein. The model then solves one mass balance,
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
#' @param shared_protein If \code{TRUE}, all sequences draw on one pool of
#'   each protein. If \code{FALSE} (default), each sequence is simulated alone.
#' @param accessibility \code{NULL} (default) or a list with one
#'   accessibility matrix per sequence, named with the sequence names (the
#'   FASTA headers); see \code{\link{simulateBinding}}.
#' @param seq_type Alphabet of the FASTA file: \code{"DNA"} (default; T is
#'   changed to U after reading) or \code{"RNA"} (the file holds U). A letter
#'   that does not belong to the chosen alphabet (for example U in a DNA file)
#'   stops the function with an error.
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
#' # All transcripts in one tube, sharing one pool of each protein
#' pooled <- simulateBindingFasta(fasta_file, rbp_models, c(HH = 100, HL = 100),
#'   rna_conc = 1, shared_protein = TRUE)
#' @importFrom Biostrings readDNAStringSet readRNAStringSet
#' @importFrom data.table rbindlist
#' @importFrom future.apply future_lapply
#' @export
simulateBindingFasta <- function(fasta_file, rbp_models, protein_concs, rna_conc = 10,
                                 k = 5, footprint = NULL,
                                 shared_protein = FALSE, accessibility = NULL,
                                 seq_type = c("DNA", "RNA")) {
  checkSharedProtein(shared_protein)
  seqs <- readFastaSequences(fasta_file, match.arg(seq_type))
  rna <- resolveRnaConc(rna_conc, names(seqs))
  acc_list <- resolveAccessibilityList(accessibility, names(seqs))

  if (shared_protein) {
    return(simulateSharedSupply(seqs, rbp_models, protein_concs, rna, k,
                                footprint = footprint,
                                accessibility = acc_list))
  }

  run_seq <- function(i) {
    res <- simulateBinding(seqs[[i]], rbp_models, protein_concs, rna[[i]], k,
                           footprint = footprint,
                           accessibility = acc_list[[i]])
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
#' With \code{shared_protein = TRUE}, each grid point puts all sequences in
#' one tube, each at the RNA concentration of that grid point, and they draw
#' on one pool of each protein (see \code{\link{simulateBindingFasta}}).
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
#' @param shared_protein If \code{TRUE}, all sequences draw on one pool of
#'   each protein. If \code{FALSE} (default), each sequence is simulated alone.
#' @param accessibility \code{NULL} (default) or a list with one
#'   accessibility matrix per sequence, named with the sequence names (the
#'   FASTA headers); see \code{\link{simulateBinding}}.
#' @param seq_type Alphabet of the FASTA file: \code{"DNA"} (default; T is
#'   changed to U after reading) or \code{"RNA"} (the file holds U). A letter
#'   that does not belong to the chosen alphabet (for example U in a DNA file)
#'   stops the function with an error.
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
#' @importFrom Biostrings readDNAStringSet readRNAStringSet
#' @importFrom data.table rbindlist
#' @importFrom future.apply future_lapply
#' @export
simulateGridFasta <- function(fasta_file, rbp_models, protein_conc_grid, rna_conc_grid,
                              k = 5, footprint = NULL,
                              shared_protein = FALSE, accessibility = NULL,
                              seq_type = c("DNA", "RNA")) {
  checkSharedProtein(shared_protein)
  seqs <- readFastaSequences(fasta_file, match.arg(seq_type))
  n_seqs <- length(seqs)
  acc_list <- resolveAccessibilityList(accessibility, names(seqs))

  if (shared_protein) {
    prot_grid <- do.call(expand.grid, protein_conc_grid)
    full_grid <- expand.grid(rna_conc = rna_conc_grid,
                             ProteinGridRow = seq_len(nrow(prot_grid)))
    message("Running grid simulation on ", n_seqs,
            " sequences with one shared protein pool...")
    run_point <- function(i) {
      rna_c <- full_grid$rna_conc[i]
      prot_concs <- unlist(prot_grid[full_grid$ProteinGridRow[i], , drop = FALSE])
      names(prot_concs) <- colnames(prot_grid)
      res <- simulateSharedSupply(seqs, rbp_models, prot_concs,
                                  rep(rna_c, n_seqs), k,
                                  footprint = footprint,
                                  accessibility = acc_list)
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
                        k = k, footprint = footprint,
                        accessibility = acc_list[[i]])
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
# Shared protein pool over many sequences
# ==============================================================================

# Read Sequences From FASTA File(s)
#
# @param fasta_file Character vector of FASTA paths.
# @param seq_type "DNA" (T is changed to U) or "RNA".
# @return Named character vector of RNA sequences.
readFastaSequences <- function(fasta_file, seq_type = "DNA") {
  for (f in fasta_file) {
    if (!file.exists(f)) stop("FASTA file not found: ", f)
  }
  reader <- if (seq_type == "RNA") Biostrings::readRNAStringSet else
    Biostrings::readDNAStringSet
  other <- if (seq_type == "RNA") "DNA" else "RNA"
  # Biostrings drops letters outside the alphabet with only a warning. Stop
  # instead, so that no nucleotide is lost without notice.
  seq_set <- withCallingHandlers(
    do.call(c, lapply(fasta_file, reader)),
    warning = function(w) {
      if (grepl("invalid one-letter sequence codes", conditionMessage(w))) {
        stop(conditionMessage(w), ". The file has letters that are not ",
             seq_type, ". Use seq_type = \"", other, "\" if the file holds ",
             other, ".", call. = FALSE)
      }
    })
  seq_names <- names(seq_set)
  if (is.null(seq_names) || any(seq_names == "")) {
    warning("Sequence names missing, using Tx indices")
    seq_names <- paste0("Tx", seq_along(seq_set))
  }
  seqs <- gsub("T", "U", as.character(seq_set), fixed = TRUE)
  stats::setNames(unname(seqs), seq_names)
}

# Check the shared_protein Argument
checkSharedProtein <- function(shared_protein) {
  if (!is.logical(shared_protein) || length(shared_protein) != 1L ||
      is.na(shared_protein)) {
    stop("shared_protein must be TRUE or FALSE.")
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
# @param footprint See simulateBinding().
# @param accessibility NULL, or a list of accessibility matrices in the order
#   of seqs.
# @return A data.table as from simulateBindingFasta(), or NULL when no
#   sequence has a k-mer.
simulateSharedSupply <- function(seqs, rbp_models, protein_concs, rna, k,
                                 footprint = NULL,
                                 accessibility = NULL) {
  checkProteinInputs(rbp_models, protein_concs)
  rbp_names <- names(protein_concs)
  p0_vec <- protein_concs[rbp_names]
  w_vec <- if (is.null(footprint)) NULL else resolveFootprint(footprint, rbp_names, k)
  acc <- lapply(seq_along(seqs), function(i) {
    checkAccessibility(accessibility[[i]], nchar(seqs[[i]]), k,
                       label = paste0("accessibility for ", names(seqs)[i]))
  })

  sites <- lapply(seqs, prepareKmerSites, rbp_models = rbp_models,
                  rbp_names = rbp_names, k = k)
  keep <- which(!vapply(sites, is.null, logical(1)))
  for (i in setdiff(seq_along(seqs), keep)) {
    warning("Simulation returned NULL for ", names(seqs)[i])
  }
  if (length(keep) == 0L) return(NULL)

  if (is.null(footprint)) {
    # Every unique k-mer of every sequence is one site species, at the
    # concentration of its sequence times its copy number (with
    # accessibility: every window, see windowSites()).
    ws <- lapply(keep, function(i) windowSites(sites[[i]], rna[i], acc[[i]], k))
    site_conc <- unlist(lapply(ws, `[[`, "site_conc"))
    kd_mat <- do.call(rbind, lapply(ws, `[[`, "kd_mat"))
    pool <- solveSharedPool(site_conc, kd_mat, p0_vec)
    n_rows <- vapply(ws, function(x) length(x$site_conc), integer(1))
    first_row <- cumsum(c(0L, utils::head(n_rows, -1L)))
    tracks <- lapply(seq_along(keep), function(t) {
      windowTrack(pool$theta[first_row[t] + ws[[t]]$index, , drop = FALSE],
                  sites[[keep[t]]]$L, k)
    })
  } else {
    molecules <- lapply(keep, function(i) {
      list(a_list = latticeAffinity(sites[[i]], w_vec, k, acc[[i]]),
           L = sites[[i]]$L, rna_conc = rna[i])
    })
    pool <- solveLatticePoolMulti(molecules, w_vec, p0_vec)
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
