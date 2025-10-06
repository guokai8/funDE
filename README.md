# funDE <img src="man/figures/logo.png" align="right" height="139" />

> **Fun**ctional **D**ifferential **E**xpression Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/guokai8/funDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/guokai8/funDE/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/guokai8/funDE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/guokai8/funDE)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

## Overview

**funDE** makes differential expression analysis more fun(ctional)! Analyze your RNA-seq data at multiple biological levels:

🧬 **Genes** → Individual gene expression  
👨‍👩‍👧‍👦 **Gene Families** → Coordinated regulation of related genes  
🛤️ **Pathways** → Biological process activity  
🔬 **Enrichment** → Gene set over-representation  

## Why funDE?

Traditional gene-level analysis can miss important biology:

| Problem | funDE Solution |
|---------|----------------|
| Subtle coordinated changes | ✅ Pathway-level detection |
| Gene family compensation | ✅ Family aggregation |
| Lost in multiple testing | ✅ Reduced testing burden |
| Hard to interpret | ✅ Functional interpretation |

## Installation

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("guokai8/funDE")

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "limma", "GSVA", "AUCell", "fgsea"))

# Optional: for pathway databases
install.packages("msigdbr")
```

## Quick Start

```r
library(funDE)

# Load your data
data("example_counts")      # Gene × sample count matrix
data("example_metadata")    # Sample information

# 1. Gene-level analysis
gene_results <- analyze_genes(
  counts = example_counts,
  metadata = example_metadata,
  design = ~ condition,
  contrast = c("condition", "treated", "control")
)

# View results
head(gene_results$results)
summary(gene_results)

# 2. Pathway analysis
pathway_results <- analyze_pathways(
  counts = example_counts,
  metadata = example_metadata,
  design = ~ condition,
  contrast = c("condition", "treated", "control"),
  pathways = "Hallmark",              # 50 well-defined pathways
  scoring_method = "GSVA"
)

# 3. Enrichment analysis
sig_genes <- gene_results$results %>%
  filter(padj < 0.05) %>%
  pull(gene)

enrichment <- test_enrichment(
  genes = sig_genes,
  universe = rownames(gene_results$results),
  gene_sets = get_pathways("KEGG", species = "human")
)

# 4. Visualize
plot_enrichment(enrichment, top_n = 20)

plot_pathway_heatmap(
  pathway_results,
  metadata = example_metadata,
  annotation_cols = "condition"
)

# 5. Compare levels
plot_multilevel(
  results_list = list(
    genes = gene_results,
    pathways = pathway_results
  ),
  plot_type = "concordance"
)
```

## Features

### 🔧 Core Analysis
- **Multi-level DE**: Genes, families, pathways, regions
- **Flexible backends**: DESeq2, edgeR, limma
- **Pathway scoring**: GSVA, AUCell, fgsea

### 📊 Rich Visualizations
- Multi-level comparison plots
- Pathway activity heatmaps
- Interactive enrichment plots
- Coordinated regulation views

### 🗃️ Comprehensive Databases
- MSigDB collections (Hallmark, KEGG, Reactome, GO)
- Gene family annotations (HGNC)
- Custom pathway support

### 🧪 Single-cell Support
- Pseudobulk aggregation
- Cell type-specific analysis
- Multi-sample integration

## Documentation

- [Introduction](vignettes/introduction.Rmd)
- [Gene Family Analysis](vignettes/gene_family_analysis.Rmd)
- [Pathway Analysis](vignettes/pathway_analysis.Rmd)
- [Single-cell Workflow](vignettes/scrna_workflow.Rmd)

## Citation

If you use funDE in your research, please cite:

```
# Citation will be added upon publication
```

## License

MIT License. See [LICENSE](LICENSE) for details.
