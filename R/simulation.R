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
#' @return A data.table containing per-position binding probabilities/occupancies.
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
  if (!all(names(protein_concs) %in% names(rbp_models))) {
    stop("All RBPs in protein_concs must have a corresponding model in rbp_models")
  }
  if (!is.numeric(protein_concs) || anyNA(protein_concs) || any(protein_concs < 0)) {
    stop("protein_concs must be non-negative numbers.")
  }
  if (!is.numeric(rna_conc) || length(rna_conc) != 1L || is.na(rna_conc) || rna_conc < 0) {
    stop("rna_conc must be a single non-negative number.")
  }
  rbp_names <- names(protein_concs)
  n_rbps <- length(rbp_names)
  om_vec <- resolveOmega(omega, rbp_names)
  if (is.null(footprint) && any(om_vec != 1)) {
    stop("omega acts between neighboring proteins and needs a footprint.")
  }
  

  kmers <- extractKmers(sequence, k)
  if (length(kmers) == 0) {
    warning("No k-mers found for sequence length ", nchar(sequence), " k = ", k)
    return(NULL)
  }
  
  unique_kmers <- countKmers(kmers)
  if (nrow(unique_kmers) == 0) {
      warning("Unique k-mers count is 0")
      return(NULL)
  }
  
  binding_params <- unique_kmers
  
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
  
  p0_vec <- protein_concs[rbp_names]
  kmer_indices <- match(kmers, binding_params$motif)
  L <- nchar(sequence)
  W <- length(kmers)

  if (is.null(footprint)) {
    # Independent windows. All windows draw on one supply of each protein.
    # Each unique k-mer is a site species at rna_conc times its copy number.
    kd_mat <- do.call(cbind, lapply(rbp_names, function(rbp) {
      binding_params[[paste0("Kd_", rbp)]]
    }))
    colnames(kd_mat) <- rbp_names
    site_conc <- rna_conc * binding_params$count

    pool <- solveSharedPool(site_conc, kd_mat, p0_vec)
    window_occupancy <- pool$theta[kmer_indices, , drop = FALSE]

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

    final_mat <- pos_sums / pos_counts
  } else {
    # Lattice: a bound RBP covers `w` nucleotides and bound proteins may not
    # overlap. The k-mer in the middle of the footprint sets the affinity.
    w_vec <- resolveFootprint(footprint, rbp_names, k)
    a_list <- lapply(seq_len(n_rbps), function(j) {
      kd_win <- binding_params[[paste0("Kd_", rbp_names[j])]][kmer_indices]
      m <- L - w_vec[j] + 1L
      if (m < 1L) return(numeric(0))
      offset <- (w_vec[j] - k) %/% 2L
      1 / kd_win[seq_len(m) + offset]
    })
    pool <- solveLatticePool(a_list, w_vec, L, rna_conc, p0_vec, omega = om_vec)
    final_mat <- vapply(seq_len(n_rbps), function(j) {
      latticeCoverage(pool$starts[[j]], w_vec[j], L)
    }, numeric(L))
    final_mat <- matrix(final_mat, nrow = L, ncol = n_rbps)
  }
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
#' @param fasta_file Character vector of path(s) to FASTA files.
#' @param rbp_models Named list of RBP models.
#' @param protein_concs Named numeric vector of protein concentrations.
#' @param rna_conc RNA concentration (default 10).
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) or footprint sizes; see
#'   \code{\link{simulateBinding}}.
#' @param omega Neighbor factor (default 1); see \code{\link{simulateBinding}}.
#' @return A data.table with combined results, including a transcript column.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' fasta_file <- system.file("extdata", "test_transcripts.fa", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBindingFasta(fasta_file, rbp_models, c(HH = 100, HL = 100))
#' head(results)
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table rbindlist
#' @importFrom future.apply future_lapply
#' @export
simulateBindingFasta <- function(fasta_file, rbp_models, protein_concs, rna_conc = 10,
                                 k = 5, footprint = NULL, omega = 1) {
  
  # Input validation
  for (f in fasta_file) {
    if (!file.exists(f)) stop("FASTA file not found: ", f)
  }
  
  # Read sequences from all files and combine them
  seqs_dna_list <- lapply(fasta_file, Biostrings::readDNAStringSet)
  seqs_dna <- do.call(c, seqs_dna_list)
  n_seqs <- length(seqs_dna)
  
  seq_names <- names(seqs_dna)
  if (is.null(seq_names) || any(seq_names == "")) {
    warning("Sequence names missing, using Tx indices")
    seq_names <- paste0("Tx", seq_len(n_seqs))
  }
  
  get_seq_str <- function(i) {
    seq_str <- as.character(seqs_dna[[i]])
    gsub("T", "U", seq_str, fixed = TRUE)
  }
  
  run_seq <- function(i) {
    seq_str <- get_seq_str(i)
    res <- simulateBinding(seq_str, rbp_models, protein_concs, rna_conc, k,
                           footprint = footprint, omega = omega)
    if (is.null(res)) {
      warning("Simulation returned NULL for ", seq_names[i])
      return(NULL)
    }
    res$transcript <- seq_names[i]
    return(res)
  }
  
  # Execute via future.apply (respects user-configured future::plan)
  results_list <- future.apply::future_lapply(seq_len(n_seqs), run_seq, 
                                              future.seed = TRUE, 
                                              future.packages = c("data.table", "Biostrings"))
  
  data.table::rbindlist(results_list)
}


#' Simulate Grid from FASTA File(s)
#'
#' Simulates competitive binding for all sequences in one or more FASTA files across a concentration grid.
#' Concurrency is controlled by the user's active \code{\link[future]{plan}} (default is sequential).
#'
#' @param fasta_file Character vector of path(s) to FASTA files.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration.
#' @param rna_conc_grid Numeric vector of RNA concentrations to sweep.
#' @param k K-mer size (default 5).
#' @param footprint \code{NULL} (default) or footprint sizes; see
#'   \code{\link{simulateBinding}}.
#' @param omega Neighbor factor (default 1); see \code{\link{simulateBinding}}.
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
                              k = 5, footprint = NULL, omega = 1) {
  
  # Input validation
  for (f in fasta_file) {
    if (!file.exists(f)) stop("FASTA file not found: ", f)
  }
  
  # Read sequences from all files and combine them
  seqs_dna_list <- lapply(fasta_file, Biostrings::readDNAStringSet)
  seqs_dna <- do.call(c, seqs_dna_list)
  n_seqs <- length(seqs_dna)
  
  seq_names <- names(seqs_dna)
  if (is.null(seq_names) || any(seq_names == "")) {
    warning("Sequence names missing, using Tx indices")
    seq_names <- paste0("Tx", seq_len(n_seqs))
  }
  
  message("Running grid simulation on ", n_seqs, " sequences...")
  
  run_seq_grid <- function(i) {
    seq_str <- as.character(seqs_dna[[i]])
    seq_str <- gsub("T", "U", seq_str, fixed = TRUE)
    
    res <- simulateGrid(seq_str, rbp_models, protein_conc_grid, rna_conc_grid,
                        k = k, footprint = footprint, omega = omega)
    if (!is.null(res)) {
      res$transcript <- seq_names[i]
    }
    return(res)
  }
  
  # Execute via future.apply (respects user-configured future::plan)
  results_list <- future.apply::future_lapply(seq_len(n_seqs), run_seq_grid, 
                                              future.seed = TRUE, 
                                              future.packages = c("data.table", "Biostrings"))
  
  data.table::rbindlist(results_list)
}
