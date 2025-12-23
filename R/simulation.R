#' Simulate Competitive Binding
#'
#' Simulates competitive binding of multiple RBPs to a single RNA sequence.
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of data.frames/data.tables. Each must have 'motif' and 'Kd' columns. The names of the list must match names in protein_concs.
#' @param protein_concs Named numeric vector of total protein concentrations.
#' @param rna_conc Total RNA concentration (default 10).
#' @param k K-mer size (default 5).
#' @return A data.table containing per-position binding probabilities/occupancies.
#' @importFrom data.table as.data.table setkey rbindlist
#' @export
simulateBinding <- function(sequence, rbp_models, protein_concs, rna_conc = 10, k = 5) {
  
  # 1. Validation
  if (!all(names(protein_concs) %in% names(rbp_models))) {
    stop("All RBPs in protein_concs must have a corresponding model in rbp_models")
  }
  rbp_names <- names(protein_concs)
  n_rbps <- length(rbp_names)
  
  # 2. Extract and count k-mers
  kmers <- extractKmers(sequence, k)
  if (length(kmers) == 0) {
    warning(paste("No k-mers found for sequence length", nchar(sequence), "k =", k))
    return(NULL)
  }
  
  # Get unique k-mers and their counts (effective concentration multiplier)
  # using countKmers from kmer.R (assuming package is loaded or access via ::)
  # For implementation inside package, we can call countKmers directly
  unique_kmers <- countKmers(kmers) # columns: motif, count
  if (nrow(unique_kmers) == 0) {
      warning("Unique k-mers count is 0")
      return(NULL)
  }
  
  # 3. Model Lookup
  # We need a table: motif | Kd_RBP1 | Kd_RBP2 ...
  # Initialize with unique motifs
  binding_params <- unique_kmers
  
  for (rbp in rbp_names) {
    model <- rbp_models[[rbp]]
    if (!is.data.frame(model) || !all(c("motif", "Kd") %in% colnames(model))) {
      stop(paste("Model for", rbp, "must be a dataframe with 'motif' and 'Kd' columns"))
    }
    # Join
    # We use match/merge. data.table merge is fast
    # Assuming valid Kds for all motifs. If missing, default to Inf? Or Error?
    # For now, let's assume models are complete for 4^k space or at least cover the sequence
    
    # Efficient lookup using match
    match_idx <- match(binding_params$motif, model$motif)
    if (any(is.na(match_idx))) {
      # Identify missing motifs
      missing <- head(binding_params$motif[is.na(match_idx)])
      warning(paste("Some motifs in sequence not found in model for", rbp, ":", paste(missing, collapse=","), "... using Inf Kd (no binding)"))
      kds <- rep(Inf, nrow(binding_params))
      kds[!is.na(match_idx)] <- model$Kd[match_idx[!is.na(match_idx)]]
    } else {
      kds <- model$Kd[match_idx]
    }
    
    binding_params[[paste0("Kd_", rbp)]] <- kds
  }
  
  # 4. Solve Equilibrium
  # We will iterate over rows of binding_params
  # Vectorization of solveEquilibrium is not trivial because uniroot is not vectorized.
  # But we can apply it per unique k-mer.
  
  # Prepare result columns
  binding_params$FreeRNA <- 0.0
  for (rbp in rbp_names) {
    binding_params[[paste0("Bound_", rbp)]] <- 0.0
  }
  
  # This loop can be parallelized if needed, but for typical k-mer counts (4^5=1024) it's fast sequentially
  # For larger k or many unique, parallel might be needed?
  # 4^7 = 16384. Sequential is likely fine (16k calls to uniroot takes < 1s usually).
  
  # Extract vectors for speed
  p0_vec <- protein_concs[rbp_names]
  
  # Loop
  n_unique <- nrow(binding_params)
  
  # Pre-allocate result vectors
  res_Req <- numeric(n_unique)
  res_Bound <- matrix(0, nrow=n_unique, ncol=n_rbps)
  colnames(res_Bound) <- rbp_names
  
  for (i in seq_len(n_unique)) {
    # Effective R0 for this k-mer species = rna_conc * count
    # Wait, is R0 the single molecule concentration or bulk?
    # In the python script: eR0 = R0 * appearance.
    # This implies we are solving for the equilibrium of *this specific k-mer species* in the tube,
    # treating all instances of the k-mer as a pool of identical independent sites.
    # Yes, standard approximation.
    
    eR0 <- rna_conc * binding_params$count[i]
    
    # Get Kds for this row
    # Columns are Kd_RBP1, Kd_RBP2...
    # We can access by name
    current_kds <- numeric(n_rbps)
    for (j in seq_len(n_rbps)) {
      current_kds[j] <- binding_params[[paste0("Kd_", rbp_names[j])]][i]
    }
    
    # Solve
    # Assuming solveEquilibrium is exported from local package
    X_eq <- solveEquilibrium(eR0, p0_vec, current_kds)
    
    # Calc bound
    bound_eq <- calcBoundProtein(X_eq, p0_vec, current_kds)
    
    # Store
    # Note: Python code divides by appearance because we want "per instance" occupancy?
    # Python: simulation.loc[..., 'Req'] = Req/appearance
    # Yes. The solver gives Total Free RNA (bulk).
    # If we have N sites, and X is free RNA, then X/N is meaningless?
    # Wait.
    # eR0 is Total Concentration of this site type.
    # X_eq is Free Concentration of this site type.
    # Bound_eq is Bound Concentration (total bound to this site type).
    # Probability of *a single site* being bound = Bound_eq / eR0 ?
    # Or Bound_eq / count?
    # If concentration is defined such that 1 nM = X molecules/vol.
    # If we have `count` sites, the max occupancy sums to `count * conc_per_site`.
    # Actually, usually fractional occupancy is Bound / Total.
    # So Bound_i / eR0 is the fractional occupancy.
    # This represents the probability that any single instance of the k-mer is bound by RBP i.
    
    # Important: Python returns `Req/appearance`.
    # Req from allBind is Free RNA concentration.
    # So `Req/appearance` is roughly Free amount per site? No.
    # If eR0 = R0 * appearance, then `Req` (Free) / `appearance` scales it back to R0 scale?
    # Let's check Python again:
    # eR0 = R0 * appearance
    # Req = allBind(..., eR0, ...) # Returns free RNA concentration X
    # simulation... = Req/appearance.
    #
    # If we want probability of being bound, we usually want Bound / Total.
    # Bound_total = eR0 - Req.
    # Bound_i = ...
    #
    # The Python script visualization uses `Req` (Free?) and `HHReq` (Bound HH).
    # `simulation_per_pos['HH/R0'] = (simulation_per_pos['HHReq'])/R0`
    # This implies HHReq is a concentration (comparable to R0).
    # Since `HHReq` in simulation was `HHReq/appearance`, and `HHReq` (before div) was `(HH0 * Req)/(Kd + Req)`.
    # Wait. `(P0 * X) / (Kd + X)` is the Bound Protein concentration.
    # So `HHReq` (before div) is the total concentration of HH bound to this k-mer species.
    # Dividing by `appearance` gives the concentration "per site instance" (if we imagine splitting the volume?).
    # Effectively, `HHReq_per_instance = HHReq_total / count`.
    # And since `eR0 = R0 * count`, `HHReq_total / eR0` would be fractional occupancy.
    # `HHReq_total / count` = `(frac * eR0) / count` = `frac * R0`.
    # So yes, the value stored is `Occupancy * R0`.
    # To get Probability (0-1), we divide by R0 later.
    
    # So:
    # 1. Solve for X (Free RNA total for this k-mer type).
    # 2. Calc Bound_i (Bound RBP i total for this k-mer type).
    # 3. Store `Bound_i / count` (which is Occupancy * R0).
    
    # Actually, we should just return Occupancy (0-1).
    # Occupancy_i = Bound_i / eR0.
    # But for "Per Position Aggregation", if we sum probabilities, we get expected number of bound proteins?
    # Python sums `Req` (which is `Bound/count`).
    # `simulation_per_pos.loc[...] += simulation.Req[idx]`
    # And then divides by `update_counts` (overlap depth).
    # So it computes average occupancy over the window.
    
    # So, I will store `Occupancy_i = Bound_total_i / eR0`.
    # Wait, check edge case if eR0 is 0? (Shouldn't happen if count > 0 and rna_conc > 0).
    
    if (eR0 > 0) {
      occ_vec <- bound_eq / eR0
    } else {
      occ_vec <- rep(0, n_rbps)
    }
    
    res_Bound[i, ] <- occ_vec
  }
  
  # Add results to binding_params
  binding_res <- cbind(binding_params, res_Bound)
  
  # 5. Map back to positions
  # Re-construct full kmer vector
  # Need to map 'kmers' (from step 2) to 'binding_res' rows.
  # match(kmers, binding_res$motif)
  
  kmer_indices <- match(kmers, binding_res$motif)
  
  # We need a per-position array of length L = nchar(sequence)
  L <- nchar(sequence)
  # Initialize outcome matrix: L rows, N_rbps columns
  pos_mat <- matrix(0, nrow=L, ncol=n_rbps)
  # Also need to track counts to average
  count_vec <- numeric(L)
  
  # We can loop through the k-mers (1 to W)
  W <- length(kmers) # L - k + 1
  
  # Vectorized update is hard in R without a loop or C++ (Rcpp).
  # But loop over W (length of seq) is slow in R?
  # L ~ 100-10000. 10k loop is fine. 100k might be slow.
  # Python does a loop over W.
  #
  # Optimization: Use `GIntervalTree` or simply working with Rle?
  # Or just a loop. For sim_transcriptome, we process many small sequences.
  # Optimized approach:
  # Create a matrix of occupancy for all windows: W x N_rbps
  # This corresponds to `binding_res[kmer_indices, ]`.
  
  window_occupancy <- res_Bound[kmer_indices, , drop=FALSE] # Matrix W x N
  
  # We need to add window_occupancy[i] to pos_mat[i : i+k-1]
  # This is a convolution-like operation.
  # Can be done with `filter`?
  # Actually, `count_vec` is simple:
  # count_vec is how many k-mers cover position p.
  # For p in 1..L: min(p, k, L-p+1) (roughly).
  #
  # For the sum of occupancies:
  # value at pos p = Sum_{i s.t. window i covers p} (occupancy of window i)
  # Window i covers p if i <= p < i+k.
  # => i in (p-k+1) : p.
  # So value[p] is sum of sequential W values.
  # This is a moving sum (rolling sum)!
  # We want `roll sum` of `window_occupancy`.
  # But mapped back to correct indices.
  #
  # Wait.
  # Window 1 covers 1..k.
  # Window 2 covers 2..k+1.
  # ...
  # This is exactly running sum of stream of length W, giving output of length L?
  # No, running sum of length k?
  # Yes because each position is covered by k windows (max).
  # So we pad the `window_occupancy` with zeros (k-1 zeros at start/end?) and run sum?
  #
  # Let's consider 1 RBP.
  # Occ vector `v` of length W.
  # Pos 1: v[1]
  # Pos 2: v[1] + v[2]
  # ...
  # Pos k: v[1] + ... + v[k]
  # Pos k+1: v[2] + ... + v[k+1]
  #
  # This is `stats::filter(v, rep(1, k), sides=1)`?
  # Need to be careful with alignment.
  #
  # Simple loop is safest for correctness first.
  # For 4 RBP, loop 1:W.
  
  for (i in seq_len(W)) {
    # positions i to i+k-1
    # Add scalar to block
    idx_range <- i:(i+k-1)
    
    # For each RBP
    # This inner addition might be the bottleneck.
    # pos_mat[idx_range, ] <- pos_mat[idx_range, ] + window_occupancy[i, ]
    # In R, matrix subset assignment is not super fast.
    
    # Alternative: use `rle` or `cumsum` trick.
    # Cumsum trick:
    # Add v at start, subtract at end+1. Then cumsum.
    # pos 1: +v[1]
    # pos 1+k: -v[1]
    # ...
    # cumsum will carry the value.
    
    # Efficient algorithm for 1D array updates:
    # diff_array = zeros(L + 1)
    # for i in 1:W:
    #   val = occ[i]
    #   diff_array[i] += val
    #   diff_array[i+k] -= val
    # result = cumsum(diff_array)[1:L]
    
    vals <- window_occupancy[i, ]
    pos_mat[i, ] <- pos_mat[i, ] + vals # Actually this is wrong, I need to add to DIFF matrix
  }
  
  # Implementing cumsum trick for N columns:
  # matrix extended: L+1 rows.
  diff_mat <- matrix(0, nrow=L+1, ncol=n_rbps)
  
  # We can construct the diff vector directly without loop?
  # We have value V_i starting at i and ending at i+k (exclusive).
  # So add V_i at index i.
  # Subtract V_i at index i+k.
  
  # i goes 1..W.
  # Add at 1..W.
  # Subtract at (1..W) + k => (1+k)..(W+k) = (k+1)..(L+1).
  
  # So, diff_mat[1:W, ] += window_occupancy
  # diff_mat[(k+1):(L+1), ] -= window_occupancy
  
  # Only careful about dimensions.
  # window_occupancy is W x N.
  
  diff_mat[1:W, ] <- diff_mat[1:W, ] + window_occupancy
  diff_mat[(k+1):(L+1), ] <- diff_mat[(k+1):(L+1), ] - window_occupancy
  
  # Cumsum columns
  pos_sums <- apply(diff_mat, 2, cumsum)[1:L, , drop=FALSE]
  
  # Counts
  # Same logic for counts (val=1)
  count_diff <- numeric(L+1)
  count_diff[1:W] <- count_diff[1:W] + 1
  count_diff[(k+1):(L+1)] <- count_diff[(k+1):(L+1)] - 1
  pos_counts <- cumsum(count_diff)[1:L]
  # Avoid div by zero
  pos_counts[pos_counts == 0] <- 1
  
  # Average
  final_mat <- pos_sums / pos_counts
  
  # Build result data.table
  # columns: pos, nt, RBP1, RBP2...
  
  # Get nucleotide at each pos
  # sequence split
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
#' Runs simulations across a grid of protein and RNA concentrations for a single sequence.
#'
#' @param sequence RNA sequence string.
#' @param rbp_models Named list of RBP models.
#' @param protein_conc_grid List of numeric vectors for each RBP concentration to sweep.
#' @param rna_conc_grid Numeric vector of RNA concentrations to sweep.
#' @param k K-mer size (default 5).
#' @param parallel Logical, whether to use parallel processing.
#' @param n_cores Integer number of cores. Defaults to detectCores() - 1.
#' @return A data.table with results for all grid combinations.
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
    parallel::clusterEvalQ(cl, library(data.table))
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
#' Simulates binding for all sequences in a FASTA file with fixed concentrations.
#'
#' @param fasta_file Path to FASTA file.
#' @param rbp_models Named list of RBP models.
#' @param protein_concs Named numeric vector of protein concentrations.
#' @param rna_conc RNA concentration (default 10).
#' @param k K-mer size (default 5).
#' @param parallel Logical.
#' @param n_cores Integer number of cores.
#' @return A data.table with combined results, including `transcript` column.
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
    parallel::clusterEvalQ(cl, { library(data.table); library(Biostrings) })
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
#' @return A data.table with combined results, including `transcript` column.
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
