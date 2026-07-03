# DPA
Dynamic Perturbation Analysis (DPA) is an R framework for identifying dynamic network perturbations during disease progression using state-specific weighted protein-protein interaction (PPI) networks.

The workflow integrates:

- Construction of state-specific weighted PPI networks and identification of network modules
- Module rewiring analysis
- Degree of Perturbation (DP) score calculation

---
## Workflow
![Workflow](Workflow.png)

---

## Installation

Clone this repository

```r
git clone https://github.com/yourname/DPA.git
```

Install required packages

```r
install.packages(c(
    "igraph",
    "progress"
))
```

---

## Input

### Expression matrix

Rows

- genes

Columns

- samples

Example

| Gene | Sample1 | Sample2 | ... |
|------|----------|----------|-----|

---

### Phenotype table

| sampleID | phase |
|----------|-------|
| S1 | N |
| S2 | N |
| S3 | I |
| ... | ... |

---

### STRING network

Required columns

```
protein1
protein2
```

---

## Usage
Run the complete DPA workflow using the `main.R` script:

```r
source("main.R")
```

The workflow consists of three major steps:

### Load functions
```r
source("R/global.R")
source("R/load_string.R")
source("R/build_state_network.R")
source("R/detect_modules.R")
source("R/module_rewiring.R")
source("R/calculate_DP.R")
```

### Step 1. Construct State-specific Networks and Detect Modules

Build state-specific weighted PPI networks:

```r
network <- build_state_network(
    expr,
    phenotype,
    ppi
)
```

Detect network modules and evaluate module conservation:

```r
module_res <- detect_modules(network)
```

**Outputs**

- State-specific weighted PPI network
- Gene–module assignment table
- Module conservation scores

---

### Step 2. Identify module rewiring

Identify rewiring events between adjacent disease states:

```r
rewired_res <- identify_module_rewiring(
    gene_module,
    network
)
```

**Output**

- Module rewiring table


---

### Step 3. Calculate DP scores

Calculate node-level DP scores:

```r
nodeDP <- calculate_node_DP(
    rewired_df,
    gene_module
)
```

Calculate module-level DP scores:

```r
moduleDP <- calculate_module_DP(
    rewired_df,
    gene_module
)
```

**Outputs**

- Gene DP scores
- Module DP scores

---

## Output

| File | Description |
|------|-------------|
| `State-specific_weighted_PPI_network.csv` | State-specific weighted PPI network |
| `Gene_module_table.csv` | Gene-to-module assignment table |
| `Module_change_&_conservation_score.csv` | Module conservation scores across disease states |
| `Module_rewiring_table.csv` | Rewired interactions between adjacent states |
| `DP_node.csv` | Node-level Degree of Perturbation (DP) scores |
| `DP_module.csv` | Module-level Degree of Perturbation (DP) scores |

---

## Citation

If you use **DPA** in your work, please cite our manuscript (coming soon).

```
