#' Plot Binding Profile
#'
#' Plots the binding profile for one or more RBPs.
#'
#' @param results data.table of simulation results.
#' @param rbp Character vector of RBP names to plot.
#' @param transcript Optional transcript name to filter (default NULL).
#' @param metric Type of metric: "occupancy" (default), "density", "density_fc", or "occupancy_fc".
#' @param ylim Optional y-axis limits (default NULL).
#' @param xlim Optional x-axis limits (default NULL).
#' @param xaxis_type Type of x-axis labels: "nucleotide" (default), "position", or "both".
#' @param rna_conc Optional numeric. Filter by RNA concentration (default NULL).
#' @param protein_conc Optional named numeric vector. Filter by protein concentration (default NULL).
#' @param window Optional integer. Applies K-mer sliding window average (default NULL).
#' @return A ggplot object.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' results$transcript <- "test_transcript"
#' plotBinding(results, rbp = c("HH", "HL"))
#' @export
plotBinding <- function(results, rbp, transcript = NULL, metric = c("occupancy", "density", "density_fc", "occupancy_fc"), 
                         ylim = NULL, xlim = NULL, xaxis_type = c("nucleotide", "position", "both"),
                         rna_conc = NULL, protein_conc = NULL, window = NULL) {
  
  metric <- match.arg(metric)
  xaxis_type <- match.arg(xaxis_type)
  
  # Determine which columns to use based on metric
  if (metric == "occupancy") {
    col_suffix <- ""
    y_label <- "Occupancy (Probability)"
  } else if (metric == "density") {
    col_suffix <- "_density"
    y_label <- "Density"
  } else if (metric == "density_fc") {
    col_suffix <- "_density_fc"
    y_label <- "Density Fold Enrichment"
  } else {
    col_suffix <- "_occupancy_fc"
    y_label <- "Occupancy Fold Enrichment"
  }
  
  # Build column names to extract
  value_cols <- paste0(rbp, col_suffix)
  
  # Check columns exist
  missing_cols <- setdiff(value_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Columns not found in results: ", paste(missing_cols, collapse=", "), 
         ". Available: ", paste(names(results), collapse=", "))
  }
  
  data <- copy(results)
  
  # Filters
  if (!is.null(transcript)) {
    group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
    if (!is.null(group_col)) {
       data <- data[get(group_col) == transcript]
    }
  }
  
  if (!is.null(rna_conc)) {
     target_rna <- rna_conc
     if ("rna_conc" %in% names(data)) {
        data <- data[rna_conc == target_rna]
     }
  }
  
  if (!is.null(protein_conc)) {
     for (p_name in names(protein_conc)) {
        col <- paste0("Conc_", p_name)
        val <- protein_conc[[p_name]]
        if (col %in% names(data)) {
           data <- data[get(col) == val]
        }
     }
  }
  
  if (!is.null(xlim)) {
    data <- data[pos >= xlim[1] & pos <= xlim[2]]
  }
  
  # Handle window: create K-mer averaged data
  if (!is.null(window) && window > 1 && "nt" %in% names(data)) {
    seq_dt <- unique(data[, .(pos, nt)])
    setkey(seq_dt, pos)
    full_seq <- paste(seq_dt$nt, collapse = "")
    n_pos <- nrow(seq_dt)
    
    valid_starts <- seq_dt$pos[seq_len(n_pos - window + 1)]
    
    windowed_list <- list()
    for (i in seq_along(rbp)) {
      r <- rbp[i]
      col <- value_cols[i]
      for (j in seq_along(valid_starts)) {
        start_pos <- valid_starts[j]
        end_pos <- start_pos + window - 1
        subset_data <- data[pos >= start_pos & pos <= end_pos]
        avg_val <- mean(subset_data[[col]], na.rm = TRUE)
        kmer <- substr(full_seq, j, j + window - 1)
        windowed_list[[length(windowed_list) + 1]] <- data.table::data.table(
          pos = start_pos,
          RBP = r,
          Value = avg_val,
          kmer = kmer
        )
      }
    }
    melted <- rbindlist(windowed_list)
    kmer_labels <- unique(melted[, .(pos, kmer)])
    setkey(kmer_labels, pos)
    
  } else {
    # Standard melt - create column mapping for renaming
    melted <- data.table::melt(data, id.vars = "pos", measure.vars = value_cols, 
                               variable.name = "RBP", value.name = "Value")
    # Rename RBP values to remove suffix
    if (col_suffix != "") {
      melted[, RBP := gsub(col_suffix, "", RBP)]
    }
    kmer_labels <- NULL
  }
  
  title_text <- "Binding Profile"
  if (!is.null(window) && window > 1) {
    title_text <- paste0(title_text, " (", window, "-mer window)")
  }
  
  p <- ggplot2::ggplot(melted, ggplot2::aes(x = pos, y = Value, color = RBP)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = title_text, x = "Position", y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 14, face = 'bold'),
      legend.text = ggplot2::element_text(size = 14)
    )
  
  if (!is.null(ylim)) {
    p <- p + ggplot2::ylim(ylim)
  }
  
  # Handle X-axis labels
  if (!is.null(window) && window > 1 && !is.null(kmer_labels)) {
    p <- p + ggplot2::scale_x_continuous(breaks = kmer_labels$pos, labels = kmer_labels$kmer) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
  } else if ("nt" %in% names(data)) {
    labels_dt <- unique(data[, .(pos, nt)])
    setkey(labels_dt, pos)
    
    if (xaxis_type == "position") {
      r <- range(labels_dt$pos)
      start_10 <- ceiling(r[1]/10)*10
      end_10 <- floor(r[2]/10)*10
      if (start_10 <= end_10) {
        brks <- seq(from = start_10, to = end_10, by = 10)
        p <- p + ggplot2::scale_x_continuous(breaks = brks)
      }
    } else if (xaxis_type == "nucleotide") {
      p <- p + ggplot2::scale_x_continuous(breaks = labels_dt$pos, labels = labels_dt$nt)
    } else if (xaxis_type == "both") {
      brks <- labels_dt$pos
      lbls <- labels_dt$nt
      tens_idx <- which(brks %% 10 == 0)
      lbls[tens_idx] <- paste0(lbls[tens_idx], "\n", brks[tens_idx])
      p <- p + ggplot2::scale_x_continuous(breaks = brks, labels = lbls)
    }
  }
    
  return(p)
}

#' Plot Binding Heatmap
#'
#' Plots a heatmap of binding for multiple RBPs.
#'
#' @param results data.table of simulation results.
#' @param rbps Optional vector of RBP names to include (default NULL, infers all).
#' @param transcript Transcript name to filter (required).
#' @param metric Type of metric: "occupancy" (default), "density", "density_fc", or "occupancy_fc".
#' @param zlim Optional numeric vector for color scale limits (default NULL).
#' @param xlim Optional x-axis limits (default NULL).
#' @param xaxis_type Type of x-axis labels: "nucleotide" (default), "position", or "both".
#' @param window Optional integer for K-mer sliding window average (default NULL).
#' @return A ggplot object.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' results$transcript <- "test_transcript"
#' plotHeatmap(results, transcript = "test_transcript")
#' @export
plotHeatmap <- function(results, transcript, rbps = NULL, xlim = NULL, 
                         xaxis_type = c("nucleotide", "position", "both"),
                         metric = c("occupancy", "density", "density_fc", "occupancy_fc"),
                         window = NULL, zlim = NULL) {
  
  xaxis_type <- match.arg(xaxis_type)
  metric <- match.arg(metric)
  
  # Infer RBPs if not provided
  if (is.null(rbps)) {
    meta <- c("pos", "nt", "transcript", "seq_name", "rna_conc", grep("Conc_|_density|_fc", names(results), value=TRUE))
    rbps <- setdiff(names(results), meta)
    rbps <- rbps[!grepl("_density|_fc", rbps)]
  }
  
  # Determine columns based on metric
  if (metric == "occupancy") {
    col_suffix <- ""
    legend_label <- "Occupancy"
  } else if (metric == "density") {
    col_suffix <- "_density"
    legend_label <- "Density"
  } else if (metric == "density_fc") {
    col_suffix <- "_density_fc"
    legend_label <- "Density FC"
  } else {
    col_suffix <- "_occupancy_fc"
    legend_label <- "Occupancy FC"
  }
  
  value_cols <- paste0(rbps, col_suffix)
  
  group_col <- if ("transcript" %in% names(results)) "transcript" else if ("seq_name" %in% names(results)) "seq_name" else NULL
  if (is.null(group_col)) stop("No transcript column found")
  
  data <- results[get(group_col) == transcript]
  
  if (!is.null(xlim)) {
    data <- data[pos >= xlim[1] & pos <= xlim[2]]
  }
  
  # Handle window
  if (!is.null(window) && window > 1 && "nt" %in% names(data)) {
    seq_dt <- unique(data[, .(pos, nt)])
    setkey(seq_dt, pos)
    full_seq <- paste(seq_dt$nt, collapse = "")
    n_pos <- nrow(seq_dt)
    valid_starts <- seq_dt$pos[seq_len(n_pos - window + 1)]
    
    windowed_list <- list()
    for (i in seq_along(rbps)) {
      r <- rbps[i]
      col <- value_cols[i]
      for (j in seq_along(valid_starts)) {
        start_pos <- valid_starts[j]
        end_pos <- start_pos + window - 1
        subset_data <- data[pos >= start_pos & pos <= end_pos]
        avg_val <- mean(subset_data[[col]], na.rm = TRUE)
        kmer <- substr(full_seq, j, j + window - 1)
        windowed_list[[length(windowed_list) + 1]] <- data.table::data.table(
          pos = start_pos, RBP = r, Value = avg_val, kmer = kmer
        )
      }
    }
    melted <- rbindlist(windowed_list)
    kmer_labels <- unique(melted[, .(pos, kmer)])
    setkey(kmer_labels, pos)
  } else {
    melted <- data.table::melt(data, id.vars = "pos", measure.vars = value_cols, 
                               variable.name = "RBP", value.name = "Value")
    if (col_suffix != "") melted[, RBP := gsub(col_suffix, "", RBP)]
    kmer_labels <- NULL
  }
  
  title_suffix <- if (!is.null(window) && window > 1) paste0(" (", window, "-mer)") else ""
  
  p <- ggplot2::ggplot(melted, ggplot2::aes(x = pos, y = RBP, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = legend_label, limits = zlim) +
    ggplot2::labs(title = paste0("Binding Heatmap: ", transcript, title_suffix), x = "Position", y = "RBP") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 14, face = 'bold'))

  if (!is.null(kmer_labels)) {
    p <- p + ggplot2::scale_x_continuous(breaks = kmer_labels$pos, labels = kmer_labels$kmer, expand = c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
  } else if ("nt" %in% names(data)) {
    labels_dt <- unique(data[, .(pos, nt)])
    if (xaxis_type == "nucleotide") {
      p <- p + ggplot2::scale_x_continuous(breaks = labels_dt$pos, labels = labels_dt$nt, expand = c(0,0))
    } else if (xaxis_type == "both") {
      lbls <- labels_dt$nt
      tens_idx <- which(labels_dt$pos %% 10 == 0)
      lbls[tens_idx] <- paste0(lbls[tens_idx], "\n", labels_dt$pos[tens_idx])
      p <- p + ggplot2::scale_x_continuous(breaks = labels_dt$pos, labels = lbls, expand = c(0,0))
    }
  }
  
  return(p)
}

#' Plot Competition Grid
#'
#' Visualizes binding across a concentration grid.
#'
#' @param results data.table from simulateGrid.
#' @param rbp1 Name of first RBP (subdivides each cell).
#' @param rbp2 Name of second RBP (X-axis).
#' @param rbp1_concs Optional vector of RBP1 concentrations to include (default NULL).
#' @param roi_range Numeric vector (start, end) for region of interest.
#' @param metric "occupancy" (default), "density", "density_fc", or "occupancy_fc".
#' @return A ggplot object.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateGrid(
#'   sequence = "ACGUACGU",
#'   rbp_models = rbp_models,
#'   protein_conc_grid = list(HH = c(10, 100), HL = c(10, 100)),
#'   rna_conc_grid = 10,
#'   parallel = FALSE
#' )
#' plotGrid(results, rbp1 = "HH", rbp2 = "HL", roi_range = c(1, 8))
#' @export
plotGrid <- function(results, rbp1, rbp2, roi_range, metric = c("occupancy", "density", "density_fc", "occupancy_fc"),
                      rbp1_concs = NULL) {
  
  metric <- match.arg(metric)
  
  # Determine which column to use
  if (metric == "occupancy") {
    value_col <- rbp1
    metric_title <- "Occupancy"
  } else if (metric == "density") {
    value_col <- paste0(rbp1, "_density")
    metric_title <- "Density"
  } else if (metric == "density_fc") {
    value_col <- paste0(rbp1, "_density_fc")
    metric_title <- "Density FC"
  } else {
    value_col <- paste0(rbp1, "_occupancy_fc")
    metric_title <- "Occupancy FC"
  }
  
  # Filter by ROI
  dt <- results[pos >= roi_range[1] & pos <= roi_range[2]]
  
  conc_col_1 <- paste0("Conc_", rbp1)
  conc_col_2 <- paste0("Conc_", rbp2)
  
  if (!all(c(conc_col_1, conc_col_2, "rna_conc", value_col) %in% names(dt))) {
    stop("Required columns not found. Looking for: ", value_col)
  }
  
  # Filter by rbp1_concs if specified
  if (!is.null(rbp1_concs)) {
    dt <- dt[get(conc_col_1) %in% rbp1_concs]
  }
  
  # Aggregate: average the pre-calculated metric values within ROI
  agg <- dt[, .(Value = mean(get(value_col), na.rm=TRUE)), 
            by = c("rna_conc", conc_col_1, conc_col_2)]
  
  rna_vals <- sort(unique(agg$rna_conc))
  rbp1_vals <- sort(unique(agg[[conc_col_1]]))
  rbp2_vals <- sort(unique(agg[[conc_col_2]]))
  
  n_rna <- length(rna_vals)
  n_rbp1 <- length(rbp1_vals)
  n_rbp2 <- length(rbp2_vals)
  
  # Create tile coordinates
  tile_data <- data.table::data.table()
  
  for (i in seq_along(rna_vals)) {
    for (j in seq_along(rbp2_vals)) {
      for (k in seq_along(rbp1_vals)) {
        val <- agg[rna_conc == rna_vals[i] & 
                   get(conc_col_1) == rbp1_vals[k] &
                   get(conc_col_2) == rbp2_vals[j], Value]
        if (length(val) == 0) val <- NA
        
        tile_width <- 1 / n_rbp1
        x_center <- j + (k - 0.5) * tile_width - 0.5
        y_center <- i
        
        tile_data <- rbind(tile_data, data.table::data.table(
          x = x_center, y = y_center,
          width = tile_width * 0.9, height = 0.9,
          value = val, rbp1_conc = rbp1_vals[k], rbp2_conc = rbp2_vals[j], rna = rna_vals[i]
        ))
      }
    }
  }
  
  p <- ggplot2::ggplot(tile_data, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(ggplot2::aes(width = width, height = height), color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(name = metric_title, na.value = "gray90") +
    ggplot2::scale_x_continuous(breaks = seq_len(n_rbp2), labels = rbp2_vals, 
                                 name = paste("Concentration", rbp2, "(nM)"), expand = c(0.05, 0.05)) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_rna), labels = rna_vals, 
                                 name = "RNA Concentration (nM)", expand = c(0.05, 0.05)) +
    ggplot2::labs(
      title = paste("Competition Grid:", rbp1, "vs", rbp2),
      subtitle = paste("ROI:", roi_range[1], "-", roi_range[2], "| Metric:", metric_title,
                       "\nSubdivisions:", paste(rbp1_vals, collapse=", "))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14, face = 'bold'),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}

#' View Model Affinity Distributions
#'
#' Creates a histogram of Kd or Ka distributions for RBP models.
#'
#' @param rbp_models Named list of RBP model data.tables (from setModel).
#' @param rbp Character vector of RBP names to plot (default NULL, plots all).
#' @param metric Display metric: "Kd" (default) or "Ka" (affinity = 1/Kd).
#' @param bins Number of histogram bins (default 50).
#' @param colors Named vector of colors per RBP (default NULL, uses Okabe-Ito palette).
#' @param alpha Transparency of histogram bars (default 0.5).
#' @param log_scale Logical, use log10 scale for x-axis (default FALSE).
#' @return A ggplot object.
#' @examples
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file))
#' viewModel(rbp_models, metric = "Kd")
#' @export
viewModel <- function(rbp_models, rbp = NULL, metric = c("Kd", "Ka"),
                         bins = 50, colors = NULL, alpha = 0.5, log_scale = FALSE) {
  
  metric <- match.arg(metric)
  
  # Select RBPs to plot
  if (is.null(rbp)) {
    rbp <- names(rbp_models)
  }
  missing_rbp <- setdiff(rbp, names(rbp_models))
  if (length(missing_rbp) > 0) {
    stop("RBPs not found in models: ", paste(missing_rbp, collapse = ", "))
  }
  
  # Build combined data
  plot_data <- data.table::rbindlist(lapply(rbp, function(r) {
    model <- rbp_models[[r]]
    
    if ("Kd" %in% names(model)) {
      val <- if (metric == "Kd") model$Kd else 1 / model$Kd
    } else if ("score" %in% names(model)) {
      # Raw model - score is already normalized relative affinity
      val <- model$score
    } else {
      stop("Model ", r, " has neither 'Kd' nor 'score' column.")
    }
    
    data.table::data.table(RBP = r, value = val)
  }))
  
  # Set RBP factor order
  plot_data[, RBP := factor(RBP, levels = rbp)]
  
  # Default colors if not provided
  if (is.null(colors)) {
    # Use Okabe-Ito colorblind-friendly palette (Okabe & Ito, 2008)
    default_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    colors <- setNames(rep_len(default_palette, length(rbp)), rbp)
  }
  
  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = RBP)) +
    ggplot2::geom_histogram(bins = bins, alpha = alpha, position = "identity") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = paste(metric, "Distribution by RBP"),
      x = if (metric == "Kd") "Kd (nM)" else "Ka (1/nM)",
      y = "Count",
      fill = "RBP"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "right"
    )
  
  if (log_scale) {
    p <- p + ggplot2::scale_x_log10()
  }
  
  return(p)
}
