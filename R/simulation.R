#' Simulate Competitive Binding
#'
#' Simulates competitive binding of multiple RBPs to a single RNA sequence.
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of data.frames/data.tables with 'motif' and 'Kd'.
#' @param protein_concs Named numeric vector of total protein concentrations.
#' @param rna_conc Total RNA concentration (default 10).
#' @param k K-mer size (default 5).
#' @return A data.table containing per-position binding probabilities/occupancies.
#' @examples
#' # Load and process model
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPBind")
#' raw_models <- loadModel(model_file, rbp = c("HH", "HL"))
#' rbp_models <- setModel(raw_models, max_affinity = 100)
#'
#' # Simulate binding on a short sequence
#' seq <- "ACGUACGUACGUACGUACGU"
#' results <- simulateBinding(seq, rbp_models, c(HH = 100, HL = 100), rna_conc = 10)
#' head(results)
#' @importFrom data.table as.data.table setkey rbindlist
#' @export
simulateBinding <- function(sequence, rbp_models, protein_concs, rna_conc = 10, k = 5) {
  

  if (!all(names(protein_concs) %in% names(rbp_models))) {
    stop("All RBPs in protein_concs must have a corresponding model in rbp_models")
  }
  rbp_names <- names(protein_concs)
  n_rbps <- length(rbp_names)
  

  kmers <- extractKmers(sequence, k)
  if (length(kmers) == 0) {
    warning(paste("No k-mers found for sequence length", nchar(sequence), "k =", k))
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
      stop(paste("Model for", rbp, "must be a dataframe with 'motif' and 'Kd' columns"))
    }

    match_idx <- match(binding_params$motif, model$motif)
    if (any(is.na(match_idx))) {
      missing <- head(binding_params$motif[is.na(match_idx)])
      warning(paste("Some motifs in sequence not found in model for", rbp, ":", paste(missing, collapse=","), "... using Inf Kd (no binding)"))
      kds <- rep(Inf, nrow(binding_params))
      kds[!is.na(match_idx)] <- model$Kd[match_idx[!is.na(match_idx)]]
    } else {
      kds <- model$Kd[match_idx]
    }
    
    binding_params[[paste0("Kd_", rbp)]] <- kds
  }
  
  binding_params$FreeRNA <- 0.0
  for (rbp in rbp_names) {
    binding_params[[paste0("Bound_", rbp)]] <- 0.0
  }
  
  p0_vec <- protein_concs[rbp_names]
  n_unique <- nrow(binding_params)
  res_Req <- numeric(n_unique)
  res_Bound <- matrix(0, nrow = n_unique, ncol = n_rbps)
  colnames(res_Bound) <- rbp_names
  
  for (i in seq_len(n_unique)) {
    eR0 <- rna_conc * binding_params$count[i]
    current_kds <- numeric(n_rbps)
    for (j in seq_len(n_rbps)) {
      current_kds[j] <- binding_params[[paste0("Kd_", rbp_names[j])]][i]
    }
    
    X_eq <- solveEquilibrium(eR0, p0_vec, current_kds)
    bound_eq <- calcBoundProtein(X_eq, p0_vec, current_kds)
    
    if (eR0 > 0) {
      occ_vec <- bound_eq / eR0
    } else {
      occ_vec <- rep(0, n_rbps)
    }
    res_Bound[i, ] <- occ_vec
  }
  
  binding_res <- cbind(binding_params, res_Bound)
  kmer_indices <- match(kmers, binding_res$motif)
  
  L <- nchar(sequence)
  W <- length(kmers)
  window_occupancy <- res_Bound[kmer_indices, , drop = FALSE]
  
  # Use cumsum trick for efficient position mapping
  diff_mat <- matrix(0, nrow = L + 1, ncol = n_rbps)
  diff_mat[1:W, ] <- diff_mat[1:W, ] + window_occupancy
  diff_mat[(k + 1):(L + 1), ] <- diff_mat[(k + 1):(L + 1), ] - window_occupancy
  pos_sums <- apply(diff_mat, 2, cumsum)[1:L, , drop = FALSE]
  
  count_diff <- numeric(L + 1)
  count_diff[1:W] <- count_diff[1:W] + 1
  count_diff[(k + 1):(L + 1)] <- count_diff[(k + 1):(L + 1)] - 1
  pos_counts <- cumsum(count_diff)[1:L]
  pos_counts[pos_counts == 0] <- 1
  
  final_mat <- pos_sums / pos_counts
  seq_chars <- strsplit(sequence, "")[[1]]
  
  res <- data.table::data.table(
    pos = 1:L,
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
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration.
#' @param rna_conc_grid Numeric vector of RNA concentrations to sweep.
#' @param k K-mer size (default 5).
#' @param parallel Logical, whether to use parallel processing.
#' @param n_cores Integer number of cores.
#' @return A data.table with results for all grid combinations.
#' @examples
#' \donttest{
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' grid <- simulateGrid("ACGUACGU", rbp_models,
#'   protein_conc_grid = list(HH = c(10, 100), HL = c(10, 100)),
#'   rna_conc_grid = 10, parallel = FALSE)
#' }
#' @importFrom parallel mclapply detectCores
#' @export
simulateGrid <- function(sequence, rbp_models, protein_conc_grid, rna_conc_grid, 
                          k = 5, parallel = TRUE, n_cores = NULL) {
  
  # Validate: sequence must be a single string
  if (!is.character(sequence) || length(sequence) != 1) {
    stop("sequence must be a single character string. For FASTA files, use simulateGridF().")
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
  
  # Parallel setup
  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    if (n_cores < 1) n_cores <- 1
  } else {
    n_cores <- 1
  }
  
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
    res <- simulateBinding(sequence, rbp_models, prot_concs, rna_conc = rna_c, k = k)
    
    if (is.null(res)) return(NULL)
    
    # Add metadata
    res$rna_conc <- rna_c
    for (p in names(prot_concs)) {
      res[[paste0("Conc_", p)]] <- prot_concs[[p]]
    }
    
    return(res)
  }
  
  # Execute
  if (n_cores > 1 && .Platform$OS.type != "windows") {
    results_list <- parallel::mclapply(seq_len(n_jobs), run_job, mc.cores = n_cores)
  } else if (n_cores > 1 && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, requireNamespace("data.table", quietly = TRUE))
    parallel::clusterExport(cl, varlist = c("sequence", "rbp_models", "full_grid", "prot_grid", 
                                             "simulateBinding", "solveEquilibrium", 
                                             "calcBoundProtein", "extractKmers", "countKmers"), 
                            envir = environment())
    results_list <- parallel::parLapply(cl, seq_len(n_jobs), run_job)
    parallel::stopCluster(cl)
  } else {
    results_list <- lapply(seq_len(n_jobs), run_job)
  }
  
  # Combine results
  final_dt <- data.table::rbindlist(results_list)
  return(final_dt)
}


# ==============================================================================
# FASTA Input Wrappers
# ==============================================================================

#' Simulate Binding from FASTA File
#'
#' Simulates binding for all sequences in a FASTA file.
#'
#' @param fasta_file Path to FASTA file.
#' @param rbp_models Named list of RBP models.
#' @param protein_concs Named numeric vector of protein concentrations.
#' @param rna_conc RNA concentration (default 10).
#' @param k K-mer size (default 5).
#' @param parallel Logical.
#' @param n_cores Integer number of cores.
#' @return A data.table with combined results, including transcript column.
#' @examples
#' \donttest{
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPBind")
#' fasta_file <- system.file("extdata", "test_transcripts.fa", package = "RBPBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBindingF(fasta_file, rbp_models, c(HH = 100, HL = 100))
#' }
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table rbindlist
#' @export
simulateBindingF <- function(fasta_file, rbp_models, protein_concs, rna_conc = 10, 
                               k = 5, parallel = TRUE, n_cores = NULL) {
  
  if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)
  
  seqs_dna <- Biostrings::readDNAStringSet(fasta_file)
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
  
  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    if (n_cores < 1) n_cores <- 1
  } else {
    n_cores <- 1
  }
  
  run_seq <- function(i) {
    seq_str <- get_seq_str(i)
    res <- simulateBinding(seq_str, rbp_models, protein_concs, rna_conc, k)
    if (is.null(res)) {
      warning(paste("Simulation returned NULL for", seq_names[i]))
      return(NULL)
    }
    res$transcript <- seq_names[i]
    return(res)
  }
  
  if (n_cores > 1 && .Platform$OS.type != "windows") {
    results_list <- parallel::mclapply(seq_len(n_seqs), run_seq, mc.cores = n_cores)
  } else if (n_cores > 1 && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      requireNamespace("data.table", quietly = TRUE)
      requireNamespace("Biostrings", quietly = TRUE)
    })
    parallel::clusterExport(cl, varlist = c("seqs_dna", "seq_names", "rbp_models", 
                                             "protein_concs", "rna_conc", "k",
                                             "simulateBinding", "get_seq_str"),
                            envir = environment())
    results_list <- parallel::parLapply(cl, seq_len(n_seqs), run_seq)
    parallel::stopCluster(cl)
  } else {
    results_list <- lapply(seq_len(n_seqs), run_seq)
  }
  
  data.table::rbindlist(results_list)
}


#' Simulate Grid from FASTA File
#'
#' Simulates binding for all sequences in a FASTA file across a concentration grid.
#'
#' @param fasta_file Path to FASTA file.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration.
#' @param rna_conc_grid Numeric vector of RNA concentrations.
#' @param k K-mer size (default 5).
#' @param parallel Logical.
#' @param n_cores Integer number of cores.
#' @return A data.table with combined results, including transcript column.
#' @examples
#' \donttest{
#' # See vignette for full FASTA grid example
#' }
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table rbindlist
#' @export
simulateGridF <- function(fasta_file, rbp_models, protein_conc_grid, rna_conc_grid,
                            k = 5, parallel = TRUE, n_cores = NULL) {
  
  if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)
  
  seqs_dna <- Biostrings::readDNAStringSet(fasta_file)
  n_seqs <- length(seqs_dna)
  
  seq_names <- names(seqs_dna)
  if (is.null(seq_names) || any(seq_names == "")) {
    warning("Sequence names missing, using Tx indices")
    seq_names <- paste0("Tx", seq_len(n_seqs))
  }
  
  message("Running grid simulation on ", n_seqs, " sequences...")
  
  # Parallelize at sequence level, not inside simulateGrid
  # to avoid nested parallelization
  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    if (n_cores < 1) n_cores <- 1
  } else {
    n_cores <- 1
  }
  
  run_seq_grid <- function(i) {
    seq_str <- as.character(seqs_dna[[i]])
    seq_str <- gsub("T", "U", seq_str, fixed = TRUE)
    
    # IMPORTANT: parallel=FALSE here to avoid nested parallelization
    res <- simulateGrid(seq_str, rbp_models, protein_conc_grid, rna_conc_grid,
                         k = k, parallel = FALSE)
    if (!is.null(res)) {
      res$transcript <- seq_names[i]
    }
    return(res)
  }
  
  if (n_cores > 1 && .Platform$OS.type != "windows") {
    results_list <- parallel::mclapply(seq_len(n_seqs), run_seq_grid, mc.cores = n_cores)
  } else {
    results_list <- lapply(seq_len(n_seqs), run_seq_grid)
  }
  
  data.table::rbindlist(results_list)
}
