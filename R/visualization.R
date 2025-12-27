#' Plot Binding Profile
#'
#' Plots the binding profile for one or more RBPs.
#'
#' @param results data.table of simulation results.
#' @param rbp Character vector of RBP names to plot.
#' @param transcript Optional transcript name to filter.
#' @param metric Type of metric: "occupancy", "density", "density_fc", or "occupancy_fc".
#' @param ylim Optional y-axis limits.
#' @param xlim Optional x-axis limits.
#' @param xaxis_type Type of x-axis labels: "nucleotide", "position", or "both".
#' @param rna_conc Optional numeric. Filter by RNA concentration.
#' @param protein_conc Optional named numeric vector. Filter by protein concentration.
#' @param window Optional integer. Applies K-mer sliding window average.
#' @return A ggplot object.
#' @examples
#' if (FALSE) {
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file, rbp = c("HH", "HL")))
#' results <- simulateBinding("ACGUACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' plotBinding(results, rbp = c("HH", "HL"))
#' }
#' @export
plotBinding <- function(results, rbp, transcript = NULL, metric = c("occupancy", "density", "density_fc", "occupancy_fc"), 
                         ylim = NULL, xlim = NULL, xaxis_type = c("nucleotide", "position", "both"),
                         rna_conc = NULL, protein_conc = NULL, window = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work. Please install it.")
  }
  
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
    
    valid_starts <- seq_dt$pos[1:(n_pos - window + 1)]
    
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
        windowed_list[[length(windowed_list) + 1]] <- data.table(
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
#' @param rbps Optional vector of RBP names.
#' @param transcript Transcript name to filter.
#' @param metric Type of metric: "occupancy", "density", "density_fc", or "occupancy_fc".
#' @param zlim Optional numeric vector for color scale limits.
#' @param xlim Optional x-axis limits.
#' @param xaxis_type Type of x-axis labels: "nucleotide", "position", or "both".
#' @param window Optional integer for K-mer sliding window average.
#' @return A ggplot object.
#' @examples
#' if (FALSE) {
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file))
#' results <- simulateBinding("ACGUACGUACGU", rbp_models, c(HH = 100, HL = 100))
#' plotHeatmap(results, transcript = NULL)
#' }
#' @export
plotHeatmap <- function(results, transcript, rbps = NULL, xlim = NULL, 
                         xaxis_type = c("nucleotide", "position", "both"),
                         metric = c("occupancy", "density", "density_fc", "occupancy_fc"),
                         window = NULL, zlim = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work. Please install it.")
  }
  
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
    valid_starts <- seq_dt$pos[1:(n_pos - window + 1)]
    
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
        windowed_list[[length(windowed_list) + 1]] <- data.table(
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
#' @param rbp1_concs Optional vector of RBP1 concentrations to include.
#' @param roi_range Numeric vector (start, end) for region of interest.
#' @param metric "occupancy", "density", "density_fc", or "occupancy_fc".
#' @return A ggplot object.
#' @examples
#' if (FALSE) {
#' # See vignette for full grid simulation example
#' plotGrid(grid_results, rbp1 = "HH", rbp2 = "HL", roi_range = c(1, 20))
#' }
#' @export
plotGrid <- function(results, rbp1, rbp2, roi_range, metric = c("occupancy", "density", "density_fc", "occupancy_fc"),
                      rbp1_concs = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  
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
  tile_data <- data.table()
  
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
        
        tile_data <- rbind(tile_data, data.table(
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

#' Plot Competition Bubble (Split-Dot)
#'
#' Visualizes competitive binding between two RBPs using split-dot bubbles.
#'
#' @param results data.table from simulateGrid.
#' @param rbp_x Name of first RBP (one half of split dot).
#' @param rbp_y Name of second RBP (other half of split dot).
#' @param roi_range Numeric vector (start, end) for region of interest.
#' @param rna_conc RNA concentration to filter.
#' @param protein_conc Named numeric vector to filter other protein concentrations.
#' @param metric "occupancy", "density", "density_fc", or "occupancy_fc".
#' @param scale_factor Scaling factor for bubble size (default 0.4).
#' @param colors Named vector of colors for the RBPs.
#' @param legend_breaks Numeric vector of values to show in size legend.
#' @return A ggplot object.
#' @examples
#' if (FALSE) {
#' # See vignette for bubble plot example
#' plotBubble(grid_results, rbp_x = "HH", rbp_y = "HL", roi_range = c(1, 20))
#' }
#' @export
plotBubble <- function(results, rbp_x, rbp_y, roi_range, rna_conc = NULL, 
                        protein_conc = NULL,
                        metric = c("occupancy", "density", "density_fc", "occupancy_fc"),
                        scale_factor = 0.4, colors = NULL,
                        legend_breaks = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  
  metric <- match.arg(metric)
  
  # Determine value columns based on metric
  if (metric == "occupancy") {
    col_x <- rbp_x
    col_y <- rbp_y
    metric_title <- "Occupancy"
  } else if (metric == "density") {
    col_x <- paste0(rbp_x, "_density")
    col_y <- paste0(rbp_y, "_density")
    metric_title <- "Density"
  } else if (metric == "density_fc") {
    col_x <- paste0(rbp_x, "_density_fc")
    col_y <- paste0(rbp_y, "_density_fc")
    metric_title <- "Density FC"
  } else {
    col_x <- paste0(rbp_x, "_occupancy_fc")
    col_y <- paste0(rbp_y, "_occupancy_fc")
    metric_title <- "Occupancy FC"
  }
  
  # Filter by ROI
  dt <- results[pos >= roi_range[1] & pos <= roi_range[2]]
  
  # Filter by RNA concentration
  if (!is.null(rna_conc)) {
    target_rna <- rna_conc
    dt <- dt[rna_conc == target_rna]
  } else if ("rna_conc" %in% names(dt)) {
    if (length(unique(dt$rna_conc)) > 1) {
      stop("Multiple RNA concentrations found. Please specify 'rna_conc'.")
    }
  }
  
  # Filter by other protein concentrations
  if (!is.null(protein_conc)) {
    for (p_name in names(protein_conc)) {
      col <- paste0("Conc_", p_name)
      val <- protein_conc[[p_name]]
      if (col %in% names(dt)) {
        dt <- dt[get(col) == val]
      }
    }
  }
  
  # Get concentration columns for X and Y axes
  conc_col_x <- paste0("Conc_", rbp_x)
  conc_col_y <- paste0("Conc_", rbp_y)
  
  if (!all(c(conc_col_x, conc_col_y, col_x, col_y) %in% names(dt))) {
    stop("Required columns not found. Need: ", conc_col_x, " ", conc_col_y, " ", col_x, " ", col_y)
  }
  
  # Aggregate: average metric values within ROI
  agg <- dt[, .(
    Val_X = mean(get(col_x), na.rm=TRUE),
    Val_Y = mean(get(col_y), na.rm=TRUE)
  ), by = c(conc_col_x, conc_col_y)]
  
  # Get unique concentrations
  x_vals <- sort(unique(agg[[conc_col_x]]))
  y_vals <- sort(unique(agg[[conc_col_y]]))
  
  # Map to grid indices
  agg$x_idx <- match(agg[[conc_col_x]], x_vals)
  agg$y_idx <- match(agg[[conc_col_y]], y_vals)
  
  # Create arc data for split dots
  make_arc <- function(cx, cy, r, start_rad, end_rad, n_points=20) {
    if (r < 1e-6) return(NULL)
    theta <- seq(start_rad, end_rad, length.out = n_points)
    data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
  }
  
  # Build polygon data
  poly_data <- data.table()
  counter <- 1
  
  # Get max values for scaling
  max_val <- max(c(agg$Val_X, agg$Val_Y), na.rm=TRUE)
  if (max_val <= 0) max_val <- 1
  
  for (i in seq_len(nrow(agg))) {
    cx <- agg$x_idx[i]
    cy <- agg$y_idx[i]
    
    # Scale radius by value (handle NA values)
    val_x <- agg$Val_X[i]
    val_y <- agg$Val_Y[i]
    r_x <- if (!is.na(val_x) && val_x > 0) sqrt(val_x / max_val) * scale_factor else 0
    r_y <- if (!is.na(val_y) && val_y > 0) sqrt(val_y / max_val) * scale_factor else 0
    
    # Bottom-right arc (RBP X): -pi/4 to 3*pi/4
    if (r_x > 1e-6) {
      arc_x <- make_arc(cx, cy, r_x, -pi/4, 3*pi/4)
      arc_x <- rbind(data.frame(x=cx, y=cy), arc_x)
      arc_x$group <- counter
      arc_x$RBP <- rbp_x
      arc_x$value <- val_x
      poly_data <- rbind(poly_data, arc_x)
      counter <- counter + 1
    }
    
    # Top-left arc (RBP Y): 3*pi/4 to 7*pi/4
    if (r_y > 1e-6) {
      arc_y <- make_arc(cx, cy, r_y, 3*pi/4, 7*pi/4)
      arc_y <- rbind(data.frame(x=cx, y=cy), arc_y)
      arc_y$group <- counter
      arc_y$RBP <- rbp_y
      arc_y$value <- val_y
      poly_data <- rbind(poly_data, arc_y)
      counter <- counter + 1
    }
  }
  
  # Default colors
  if (is.null(colors)) {
    colors <- c("skyblue", "salmon")
    names(colors) <- c(rbp_x, rbp_y)
  }
  
  # Extract ROI sequence for subtitle
  if ("nt" %in% names(results)) {
    seq_dt <- unique(results[pos >= roi_range[1] & pos <= roi_range[2], .(pos, nt)])
    setkey(seq_dt, pos)
    roi_seq <- paste(seq_dt$nt, collapse = "")
  } else {
    roi_seq <- ""
  }
  
  # Create size legend data - position to the right of the main plot
  n_x <- length(x_vals)
  n_y <- length(y_vals)
  
  # Use custom legend breaks if provided, otherwise default to 25/50/75/100%
  if (!is.null(legend_breaks)) {
    legend_vals <- legend_breaks
  } else {
    legend_vals <- c(0.25, 0.5, 0.75, 1.0) * max_val
  }
  
  legend_radii <- sqrt(legend_vals / max_val) * scale_factor
  max_radius <- max(legend_radii)
  
  # Position legend left edge (where title starts)
  legend_left_edge <- n_x + 1.2  # Start beyond plot edge with buffer
  
  # Position legend circles within the y-axis range
  legend_y_spacing <- max(1, (n_y - 1) / (length(legend_vals) + 1))
  legend_y_positions <- 1 + legend_y_spacing * (seq_along(legend_vals))
  
  # Build circles - left edge of each circle aligns with legend_left_edge
  legend_circles <- data.table()
  for (i in seq_along(legend_vals)) {
    theta <- seq(0, 2*pi, length.out = 30)
    # Center x = left_edge + radius, so left edge of circle = left_edge
    circle_center_x <- legend_left_edge + legend_radii[i]
    circle <- data.frame(
      x = circle_center_x + legend_radii[i] * cos(theta),
      y = legend_y_positions[i] + legend_radii[i] * sin(theta),
      group = 1000 + i,
      RBP = "legend",
      value = legend_vals[i]
    )
    legend_circles <- rbind(legend_circles, circle)
  }
  
  # Legend labels - position to the right of the largest circle
  legend_labels <- data.frame(
    x = legend_left_edge + 2 * max_radius + 0.3,  # Right of largest circle
    y = legend_y_positions,
    label = round(legend_vals, 3)
  )
  
  # Position title above the legend circles, left-aligned with left edge
  legend_title_y <- max(legend_y_positions) + legend_y_spacing * 0.6
  legend_title_df <- data.frame(
    x = legend_left_edge,  # Aligned with left edge of circles
    y = legend_title_y,
    label = metric_title
  )
  
  # Build subtitle
  subtitle_parts <- c(roi_seq, paste("\nROI:", roi_range[1], "-", roi_range[2]))
  if (!is.null(rna_conc)) subtitle_parts <- c(subtitle_parts, paste("| RNA:", rna_conc, "nM"))
  if (!is.null(protein_conc)) {
    fixed_str <- paste(names(protein_conc), "=", protein_conc, collapse = ", ")
    subtitle_parts <- c(subtitle_parts, paste("| Fixed:", fixed_str))
  }
  
  # Combine main plot data with legend circles
  all_poly <- rbind(poly_data, legend_circles, fill = TRUE)
  
  # Plot: use coord_fixed to maintain aspect ratio (circles stay circular)
  p <- ggplot2::ggplot(all_poly, ggplot2::aes(x = x, y = y, fill = RBP, group = group)) +
    ggplot2::geom_polygon(alpha = 0.8, color = "black", linewidth = 0.2) +
    # Color scale: include legend color but filter it from the color legend display
    ggplot2::scale_fill_manual(
      values = c(colors, "legend" = "gray80"),
      breaks = c(rbp_x, rbp_y),
      name = "RBP"
    ) +
    # Add legend value labels
    ggplot2::geom_text(data = legend_labels, 
                       ggplot2::aes(x = x, y = y, label = label),
                       hjust = 0, size = 3, inherit.aes = FALSE) +
    # Add legend title
    ggplot2::geom_text(data = legend_title_df, 
                       ggplot2::aes(x = x, y = y, label = label),
                       fontface = "bold", size = 3.5, hjust = 0, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_x), labels = x_vals,
                                 name = paste("Concentration", rbp_x, "(nM)")) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_y), labels = y_vals,
                                 name = paste("Concentration", rbp_y, "(nM)")) +
    # coord_fixed maintains aspect ratio (circles stay circular), clip="off" allows legend outside
    ggplot2::coord_fixed(ratio = 1, xlim = c(0.5, n_x + 0.5), ylim = c(0.5, n_y + 0.5), clip = "off") +
    ggplot2::labs(
      title = paste("Competitive Binding:", rbp_x, "vs", rbp_y, "-", metric_title),
      subtitle = paste(subtitle_parts, collapse = " ")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14, face = 'bold'),
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 150, b = 10, l = 10)
    )
  
  return(p)
}


#' View Model Affinity Distributions
#'
#' Creates a histogram of Kd or Ka distributions for RBP models.
#'
#' @param rbp_models Named list of RBP model data.tables (from setModel).
#' @param rbp Character vector of RBP names to plot. If NULL, plots all.
#' @param metric Display metric: "Kd" or "Ka" (affinity = 1/Kd).
#' @param bins Number of histogram bins (default 50).
#' @param colors Named vector of colors per RBP. If NULL, uses default palette.
#' @param alpha Transparency of histogram bars (default 0.5).
#' @param log_scale Logical, use log10 scale for x-axis (default FALSE).
#' @return A ggplot object.
#' @examples
#' if (FALSE) {
#' model_file <- system.file("extdata", "model_RBP.csv", package = "RBPEqBind")
#' rbp_models <- setModel(loadModel(model_file))
#' viewModel(rbp_models, metric = "Kd")
#' }
#' @export
viewModel <- function(rbp_models, rbp = NULL, metric = c("Kd", "Ka"),
                         bins = 50, colors = NULL, alpha = 0.5, log_scale = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work. Please install it.")
  }
  
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
    
    # Handle both raw models (score) and processed models (Kd)
    if ("Kd" %in% names(model)) {
      # Processed model
      if (metric == "Ka") {
        val <- 1 / model$Kd
      } else {
        val <- model$Kd
      }
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
    default_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                         "#FFFF33", "#A65628", "#F781BF", "#999999")
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
