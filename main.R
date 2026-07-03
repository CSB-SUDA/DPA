#============================================================
# DPA workflow
#
# Step 1. Construct state-specific weighted PPI networks and Detect network modules
# Step 2. Identify module rewiring events
# Step 3. Calculate degree perturbation (DP) scores
#============================================================

curr_path <- "F:/GC_Stage/GC_Stage.lxy/Write/code"
setwd(curr_path)

##------------------------------------------------------------
## Load functions
##------------------------------------------------------------
source("R/global.R")
source("R/load_string.R")
source("R/build_state_network.R")
source("R/detect_modules.R")
source("R/module_rewiring.R")
source("R/calculate_DP.R")

##------------------------------------------------------------
## Output directory
##------------------------------------------------------------
outdir <- "result"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

##------------------------------------------------------------
## Load expression matrix
##------------------------------------------------------------
expr <- read.csv(
  "F:/GC_Stage/GC_Stage.lxy/Experiment/0-FromLXY/3-网络构建-buildNetwork/matrixSelected.csv",
  row.names = 1,
  check.names = FALSE
)

##------------------------------------------------------------
## Load phenotype information
##------------------------------------------------------------
phenotype <- read.csv(
  "F:/GC_Stage/GC_Stage.lxy/Experiment/0-FromLXY/3-网络构建-buildNetwork/pheotype.csv",
  row.names = 1
)

#phenotype$sampleID <- make.names(phenotype$sampleID)

##------------------------------------------------------------
## Load STRING protein-protein interaction network
##------------------------------------------------------------
ppi <- load_string(
  ppi_dir = "data",
  score_cutoff = 400
)

#============================================================
# Step 1. (1) Construct state-specific weighted PPI networks
#============================================================

network <- build_state_network(
  expr = expr,
  phenotype = phenotype,
  ppi = ppi
)

write.csv(
  network,
  file.path(outdir, "State-specific_weighted_PPI_network.csv"),
  row.names = FALSE,
  quote = FALSE
)

#============================================================
# Step 1. (2) Detect network modules
#============================================================

module_res <- detect_modules(network)

gene_module <- module_res$gene_module

write.csv(
  gene_module,
  file.path(outdir, "Gene_module_table.csv"),
  row.names = FALSE,
  quote = FALSE
)

module_change <- module_res$module_change
write.csv(
  module_change,
  file.path(outdir, "Module_change_&_conservation_score.csv"),
  row.names = FALSE,
  quote = FALSE
)

#============================================================
# Step 2. Identify module rewiring during disease progression
#============================================================

rewired_res <- identify_module_rewiring(
  gene_module,
  network
)

rewired_df <- rewired_res$rewired_edges_summary

write.csv(
  rewired_df,
  file.path(outdir, "Module_rewiring_table.csv"),
  row.names = FALSE,
  quote = FALSE
)

#============================================================
# Step 3. Calculate dynamic perturbation (DP) scores
#============================================================

## Gene DP
nodeDP <- calculate_node_DP(
  rewired_df,
  gene_module
)

write.csv(
  nodeDP,
  file.path(outdir, "DP_node.csv"),
  row.names = FALSE,
  quote = FALSE
)

## Module DP
moduleDP <- calculate_module_DP(
  rewired_df,
  gene_module
)

## Save module DP
write.csv(
  moduleDP,
  file.path(outdir,"DP_module.csv"),
  row.names = FALSE,
  quote = FALSE
)
