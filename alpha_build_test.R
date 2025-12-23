# Alpha Build Test Script for RBPBind
# This script is a standalone test for the RBPBind package using real model data.
# It demonstrates loading models from CSV, simulating binding, and visualizing results.

# ==============================================================================
# 1. Setup and Library Loading
# ==============================================================================
# Ensure required packages are installed
required_pkgs <- c("devtools", "data.table", "ggplot2", "parallel", "Biostrings")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing required package:", pkg))
    install.packages(pkg)
  }
}


# Load libraries explicitly
library(devtools)
library(data.table)
library(ggplot2)
library(parallel)
# Biostrings is loaded via namespace usually, but good to have if we use DNAStringSet
if (requireNamespace("Biostrings", quietly = TRUE)) {
  library(Biostrings)
}

# Load RBPBind package from current directory
message("Loading RBPBind package...")
devtools::load_all(".")

# ==============================================================================
# 2. Data Loading
# ==============================================================================
# Load RBP models from consolidated CSV file
model_file <- "data/model_RBP.csv"
message("Loading RBP models from ", model_file, "...")

# Step 1: Load raw models (score column per RBP)
raw_models <- loadModel(model_file)
print("Raw models loaded:")
print(names(raw_models))
message("Score ranges (raw): HH=", round(range(raw_models$HH$score)[1], 3), "-", 
        round(range(raw_models$HH$score)[2], 3))

p_model_raw <- viewModel(
  rbp_models = raw_models,
  rbp = c("HH", "HL", "LH", "LL"),
  bins = 100,
  alpha = 0.4
)
print(p_model_raw)

# Step 2: Set affinity ranges (Ka → Kd conversion)
# High affinity = low Kd, Low affinity = high Kd
max_affinities <- c(
  "HH" = 100,  # Ka_max = 100 -> Kd_min = 0.01 nM
  "HL" = 100,
  "LH" = 10,   # Ka_max = 10 -> Kd_min = 0.1 nM
  "LL" = 10
)

min_affinities <- c(
  "HH" = 0.001,  # Ka_min = 0.001 -> Kd_max = 1000 nM
  "HL" = 0.001,
  "LH" = 0.001,
  "LL" = 0.001
)

rbp_models <- setModel(raw_models, 
                         max_affinity = max_affinities, 
                         min_affinity = min_affinities)

# Check Kd ranges
message("Checking Kd ranges (after affinity setting)...")
message("HH Kd range: ", round(min(rbp_models$HH$Kd), 2), " - ", round(max(rbp_models$HH$Kd), 0))
message("LH Kd range: ", round(min(rbp_models$LH$Kd), 2), " - ", round(max(rbp_models$LH$Kd), 0))

# --- Visualize Model Affinity Distributions ---
message("\nVisualizing model affinity distributions...")
p_model <- viewModel(
  rbp_models = rbp_models,
  rbp = c("HH", "HL", "LH", "LL"),
  metric = "Ka",
  bins = 100,
  alpha = 0.4
)
print(p_model)

# Custom colors example
p_model_custom <- viewModel(
  rbp_models = rbp_models,
  rbp = c("LL"),
  metric = "Ka",
  bins = 50,
  colors = c("LL" = "#377EB8"),
  alpha = 0.6
)
print(p_model_custom)

# ==============================================================================
# 1. Simulate Binding on Random Sequence (Scenario 1)
# ==============================================================================
message("Using reference RNA sequence...")
# Sequence from original Python script for comparison
seq <- "UAGCGGUGCGAUUGGCCCGUGGACCGCGUUUUUGCACUCAUCGUUUCGCACUAAGUACAUAUAGUUGCGACAAAGCCGCUUAUGAGUUGGGGGUAUAUUC"

# Define concentrations (nM)
# Competitive binding of all 4 RBPs
prot_concs <- c("HH" = 100, "HL" = 100, "LH" = 100, "LL" = 100)
rna_conc <- 10.0 # nM

message("Running simulation (k=5)...")
k_size <- nchar(rbp_models$HH$motif[1])
message("Detected k-mer size: ", k_size)

res <- simulateBinding(
  sequence = seq,
  rbp_models = rbp_models,
  protein_concs = prot_concs,
  rna_conc = rna_conc,
  k = k_size
)

message("Simulation complete.")
print(head(res))

# --- 1a. Visualization (Single Simulation) ---
message("Creating binding profile plots (Scenario 1)...")

# 1. Occupancy (Absolute Probability)
# Plot all RBPs overlaid with Y-axis limited to 0-1, X-axis Zoomed (20-60), and Both Labels
message("Plotting Occupancy (All RBPs, Zoomed)...")
p1 <- plotBinding(res, rbp = c("HH", "HL", "LH", "LL"), metric = "occupancy_fc", 
                   transcript = "Custom RNA", ylim = c(0, 2), xlim = c(20, 60), xaxis_type = "both",
                   window = 5)
print(p1)

# 2. Density (Normalized Shape)
# This mimics the original Python script visualization
message("Plotting Density (All RBPs)...")
p2 <- plotBinding(res, rbp = c("HH", "HL", "LH", "LL"), metric = "density_fc", 
                   transcript = "Custom RNA", window = 5)
print(p2)

# 3. heatmap
message("Creating heatmap (Zoomed)...")
# Add transcript column required for heatmap
res$transcript <- "Custom RNA"

p3 <- plotHeatmap(res, transcript = "Custom RNA", xlim = c(20, 60), 
                   xaxis_type = "both", metric = "occupancy_fc")
print(p3)

# 4. heatmap with density_fc metric
message("Creating heatmap with FC over Mean...")
p4 <- plotHeatmap(res, transcript = "Custom RNA", xaxis_type = "both", metric = 'density_fc', window = 5)
print(p4)

# --- 1b. Export (Single Simulation) ---
message("Testing Export for Scenario 1...")
exportResults(res, output_file = "test_res_1.json", format = "json")
message("Exported to test_res_1.json")

# --- 1c. SummarizedExperiment Conversion ---
message("Converting to SummarizedExperiment...")
se <- makeSE(res, rbp_models = rbp_models)
print(se)
message("SE assays: ", paste(names(SummarizedExperiment::assays(se)), collapse = ", "))
message("SE colData: ", paste(colnames(SummarizedExperiment::colData(se)), collapse = ", "))

# ==============================================================================
# 2. Concentration Grid Sweep (Scenario 2)
# ==============================================================================
message("Running grid simulation (Concentration Sweep)...")
# ==============================================================================
# METRIC INTERPRETATION:
# - occupancy: Raw binding probability (0-1) at each position.
# - density: Distribution of binding along sequence (sum = 1 per RBP).
# - density_fc: Fold-change over mean (>1 = enriched, <1 = depleted).
# ==============================================================================

# Grid setup (user-specified concentrations)
prot_grid <- list(
  HH = c(0, 10, 25, 50, 100, 250, 500, 1000), 
  HL = c(0, 10, 25, 50, 100, 250, 500, 1000),
  LH = c(0, 10, 25, 50, 100, 250, 500, 1000),
  LL = c(0, 10, 25, 50, 100, 250, 500, 1000)
)
rna_grid <- c(0, 10, 25, 50, 100, 250, 500, 1000)

message("Running Grid Simulation (4 RBPs)...")
res_grid <- simulateGrid(
  sequence = seq,  # Single sequence string
  rbp_models = rbp_models,
  protein_conc_grid = prot_grid,
  rna_conc_grid = rna_grid,
  k = k_size,
  parallel = TRUE
)
message("Grid simulation complete. Rows: ", nrow(res_grid))

# --- 2a. Binding Profile (with filters) ---
message("Visualizing Binding Profile (Scenario 2)...")
p_prof <- plotBinding(
  results = res_grid,
  rbp = c("HH", "HL", "LH", "LL"),
  rna_conc = 10,
  protein_conc = c(HH=10, HL=10, LH=1000, LL=1000),
  metric = "density",
  # ylim = c(0, 0.2)
)
print(p_prof)

p_bubble <- plotBubble(
  results = res_grid,
  rbp_x = "HH",
  rbp_y = "LL",
  roi_range = c(29, 33),
  rna_conc = 10,
  protein_conc = c(HL=0, LH=0),
  metric = "density_fc",
  colors = c("HH" = "skyblue", "LL" = "salmon"),
)
print(p_bubble)

# --- 2b. Separate 2-RBP Grid for plotGrid ---
message("Running 2-RBP Grid Simulation (HH vs LL)...")
prot_grid_2rbp <- list(
  HH = c(0, 10, 50, 100, 500, 1000), 
  LL = c(0, 10, 50, 100, 500, 1000)
)
rna_grid_2rbp <- c(0, 10, 50, 100, 500, 1000)

# Need to filter models to only HH and LL
rbp_models_2rbp <- rbp_models[c("HH", "LL")]

res_grid_2rbp <- simulateGrid(
  sequence = seq,  # Single sequence string
  rbp_models = rbp_models_2rbp,
  protein_conc_grid = prot_grid_2rbp,
  rna_conc_grid = rna_grid_2rbp,
  k = k_size,
  parallel = TRUE
)
message("2-RBP Grid simulation complete. Rows: ", nrow(res_grid_2rbp))

# --- 2c. Competition Grid Plot ---
message("Visualizing Competition Grid (HH vs LL)...")
# plotGrid: Y-axis = RNA conc, X-axis = LL conc, subdivided by HH conc
p_grid <- plotGrid(
  results = res_grid_2rbp,
  rbp1 = "HH",  # Subdivisions within each cell
  rbp2 = "LL",  # X-axis
  rbp1_concs = c(10, 50, 100, 500, 1000),
  roi_range = c(29, 33),
  metric = "density"
)
print(p_grid)

p_prof <- plotBinding(
  results = res_grid_2rbp,
  rbp = c("HH", "LL"),
  rna_conc = 10,
  protein_conc = c(HH=10, LL=1000),
  metric = "density_fc",
  # ylim = c(0, 0.1)
)
print(p_prof)

# --- 2d. Export (Grid) ---
message("Testing Export for Scenario 2...")
exportResults(res_grid, output_file = "test_res_2.json", format = "json")
message("Exported to test_res_2.json")

# ==============================================================================
# 3. FASTA + Fixed Concentrations (Scenario 3 - RNACompete)
# ==============================================================================
message("\n=== SCENARIO 3: FASTA + Fixed Concentrations (RNACompete) ===")

# Load RNACompete model using unified API
rnacompete_file <- "data/rnacompete_z.csv"
message("Loading RNACompete data from ", rnacompete_file, "...")

raw_models_rc <- loadModel(rnacompete_file)
message("Loaded ", length(raw_models_rc), " RBPs: ", paste(names(raw_models_rc), collapse = ", "))

# Set affinity ranges (Ka → Kd)
# U2AF2:  Ka_max = 1/100  = 0.01,   Ka_min = 0.00001
# PTBP1:  Ka_max = 1/300  = 0.0033, Ka_min = 0.00001
# HNRNPC: Ka_max = 1/700  = 0.0014, Ka_min = 0.00001
rbp_models_rc <- setModel(
  raw_models_rc,
  max_affinity = c(U2AF2 = 1/100, PTBP1 = 1/100, HNRNPC = 1/100),
  min_affinity = 0.00001
)

message("RNACompete models created:")
message("  U2AF2 Kd range: ", round(min(rbp_models_rc$U2AF2$Kd), 0), " - ", round(max(rbp_models_rc$U2AF2$Kd), 0))
message("  PTBP1 Kd range: ", round(min(rbp_models_rc$PTBP1$Kd), 0), " - ", round(max(rbp_models_rc$PTBP1$Kd), 0))
message("  HNRNPC Kd range: ", round(min(rbp_models_rc$HNRNPC$Kd), 0), " - ", round(max(rbp_models_rc$HNRNPC$Kd), 0))

# Run FASTA simulation with fixed concentrations (from legacy script)
fasta_file <- "data/test_transcripts.fa"
message("\nRunning simulateBindingF on ", fasta_file, "...")
res_fasta <- simulateBindingF(
  fasta_file = fasta_file,
  rbp_models = rbp_models_rc,
  protein_concs = c(U2AF2 = 500, PTBP1 = 200, HNRNPC = 200),
  rna_conc = 10,
  k = 7
)
message("Simulation complete. Transcripts: ", uniqueN(res_fasta$transcript))
message("Positions per transcript: ", nrow(res_fasta) / uniqueN(res_fasta$transcript))

# Visualization: Focus on PTBP2 transcript, ROI = positions 1851-2051
ptbp2_data <- res_fasta[grepl("PTBP2", transcript)]
message("\nPTBP2 data rows: ", nrow(ptbp2_data))

message("Creating line plot for PTBP2 ROI (1851-2051)...")
p_sc3 <- plotBinding(
  results = ptbp2_data,
  rbp = c("U2AF2", "PTBP1", "HNRNPC"),
  xlim = c(1851, 2051),
  metric = "density_fc",
  window = 7
)
print(p_sc3)

message("Creating heatmap for PTBP2 ROI (1851-2051)...")
p_sc3_heatmap <- plotHeatmap(ptbp2_data, transcript = unique(ptbp2_data$transcript)[1], 
                               xaxis_type = "both", xlim = c(1851, 2051), 
                               metric = "density_fc")
print(p_sc3_heatmap)

# ==============================================================================
# 4. FASTA + Grid Sweep (Scenario 4 - RNACompete)
# ==============================================================================
message("\n=== SCENARIO 4: FASTA + Grid Sweep (RNACompete) ===")

# Grid sweep with 3 RBPs (matching legacy 6_1 script)
message("Running simulateGridF (3 RBPs)...")
res_grid_fasta <- simulateGridF(
  fasta_file = fasta_file,
  rbp_models = rbp_models_rc[c("U2AF2", "PTBP1", "HNRNPC")],
  protein_conc_grid = list(
    U2AF2 = c(0, 100, 200, 500),
    PTBP1 = c(0, 100, 200, 500),
    HNRNPC = c(0, 100, 200, 500)
  ),
  rna_conc_grid = c(10),
  k = 7,
  parallel = TRUE
)
message("Grid simulation complete. Rows: ", nrow(res_grid_fasta))

# Filter to PTBP2 transcript
ptbp2_grid <- res_grid_fasta[grepl("PTBP2", transcript)]
tx_name <- unique(ptbp2_grid$transcript)[1]

# --- Scenario 4a: U2AF2 alone (200 nM) ---
message("\nCreating heatmap: U2AF2 alone (500 nM)...")
grid_4a <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 0 & Conc_HNRNPC == 0]
p_sc4a <- plotHeatmap(grid_4a, transcript = tx_name, xlim = c(1851, 2051), 
                        xaxis_type = "both", metric = "density_fc",
                       window = 7,
                       zlim = c(0, 3.5))
print(p_sc4a)

# --- Scenario 4b: U2AF2 + PTBP1 (500, 200 nM) ---
message("Creating heatmap: U2AF2 + PTBP1 (500, 200 nM)...")
grid_4b <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 200 & Conc_HNRNPC == 0]
p_sc4b <- plotHeatmap(grid_4b, transcript = tx_name, xlim = c(1851, 2051), 
                        xaxis_type = "both", metric = "density_fc",
                       window = 7,
                       zlim = c(0, 3.5))
print(p_sc4b)

# --- Scenario 4c: U2AF2 + HNRNPC (500, 200 nM) ---
message("Creating heatmap: U2AF2 + HNRNPC (500, 200 nM)...")
grid_4c <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 0 & Conc_HNRNPC == 200]
p_sc4c <- plotHeatmap(grid_4c, transcript = tx_name, xlim = c(1851, 2051), 
                        xaxis_type = "both", metric = "density_fc",
                       window = 7,
                       zlim = c(0, 3.5))
print(p_sc4c)




# --- Scenario 4a: U2AF2 alone (200 nM) ---
message("\nCreating line plot: U2AF2 alone (500 nM)...")
grid_4a <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 0 & Conc_HNRNPC == 0]
p_line4a <- plotBinding(grid_4a, rbp = "U2AF2", xlim = c(1851, 2051), 
                         metric = "density_fc", window = 5, ylim = c(0, 4))
print(p_line4a)

# --- Scenario 4b: U2AF2 + PTBP1 (500, 200 nM) ---
message("Creating line plot: U2AF2 + PTBP1 (500, 200 nM)...")
grid_4b <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 200 & Conc_HNRNPC == 0]
p_line4b <- plotBinding(grid_4b, rbp = c("U2AF2", "PTBP1"), xlim = c(1851, 2051), 
                         metric = "density_fc", window = 5, ylim = c(0, 4))
print(p_line4b)

# --- Scenario 4c: U2AF2 + HNRNPC (200, 200 nM) ---
message("Creating line plot: U2AF2 + HNRNPC (500, 200 nM)...")
grid_4c <- ptbp2_grid[Conc_U2AF2 == 500 & Conc_PTBP1 == 0 & Conc_HNRNPC == 200]
p_line4c <- plotBinding(grid_4c, rbp = c("U2AF2", "HNRNPC"), xlim = c(1851, 2051), 
                          metric = "density_fc", window = 5, ylim = c(0, 4))
print(p_line4c)

message("\n=== Alpha Test Complete! ===")
