
# TF_lib

An R helper library for discovering transcription-factor (TF) regulators of any hub-gene set and for comparing their expression in two bulk RNA-seq contrasts.

---

## Overview  

Hub genes identified from network or pathway analyses often highlight core biological processes. By uncovering the TFs that govern these genes, researchers can expose upstream regulatory programs and potential therapeutic entry points.  
`tf_enrichment_lib.R` wraps this end-to-end workflow—from CHEA3 queries to publication-ready figures—into a single, lightweight script.

---
## Scripts in this repo

TF_Analysis.R

---

## Clone the repository  

Choose one of the three common methods below, then move into the project folder:

```bash
# HTTPS
git clone https://github.com/Abtint/TF-Analysis-for-Hub-Genes.git

# SSH (requires a configured SSH key)
git clone git@github.com:Abtint/TF-Analysis-for-Hub-Genes.git

# GitHub CLI
gh repo clone Abtint/TF-Analysis-for-Hub-Genes

cd TF-Analysis-for-Hub-Genes
````

---

## Installation

```r
install.packages(c("httr", "jsonlite", "data.table", "ggplot2"))
# Then, inside your R session:
source("tf_enrichment_lib.R")
```

---

## Quick start

```r
# 1. Supply a character vector of hub genes
hub_genes <- c("GENE1", "GENE2", "GENE3")   # replace with your genes

# 2. Load two differential-expression (DE) tables as data.table objects
deg_a <- data.table::fread("conditionA_vs_control.csv")   # first contrast
deg_b <- data.table::fread("conditionB_vs_control.csv")   # second contrast

# 3. Run the complete workflow
results <- tf_enrichment_workflow(
             hub_genes   = hub_genes,
             deg1_dt     = deg_a,
             deg2_dt     = deg_b,
             deg1_name   = "Condition A",
             deg2_name   = "Condition B",
             logfc1_col  = "log2FC",   # column with log-fold-changes in deg_a
             logfc2_col  = "log2FC"    # column with log-fold-changes in deg_b
           )

# 4. Explore results
results$top_tfs        # top-ranked regulators from CHEA3
results$shared_compare # side-by-side statistics for shared TFs
results$plot           # ggplot object for quick visualisation
```

---

## Function glossary

| Function                    | Purpose                                                               |
| --------------------------- | --------------------------------------------------------------------- |
| `run_chea3()`               | Submit hub-gene list to the CHEA3 API and return enrichment table     |
| `select_top_tfs()`          | Keep the strongest regulators by percentage rank                      |
| `intersect_degs_with_tfs()` | Merge TF list with a DE table and filter by adjusted P value          |
| `compare_two_conditions()`  | Evaluate shared TFs across two contrasts and flag stronger regulation |
| `plot_shared_tfs()`         | Horizontal bar plot of log-fold-changes for shared TFs                |
| `tf_enrichment_workflow()`  | One-call wrapper that chains everything together                      |

---

## Input requirements

* **Hub genes** – character vector of official gene symbols
* **Differential-expression tables** – `data.table` (or `data.frame`) with at least:

  * gene-symbol column
  * log-fold-change column
  * adjusted P-value column

---

## Output

* Data tables for raw enrichment, top-ranked TFs, significant TFs per condition, and shared comparison
* Bar plot as a `ggplot` object

---

## License

MIT

---

## Citation

Tondar, A. 2025. https://github.com/Abtint/TF-Analysis-for-Hub-Genes. 

```
```


