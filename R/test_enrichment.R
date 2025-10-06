#' Test Gene Set Enrichment
#'
#' @description
#' Performs enrichment analysis on a list of genes using hypergeometric test
#' or Fisher's exact test. This replaces clusterProfiler functionality with
#' a simpler, more transparent implementation.
#'
#' @param genes Character vector of gene symbols or IDs (e.g., DE genes)
#' @param universe Character vector of all genes in the analysis (background)
#' @param gene_sets Named list of gene sets (e.g., from get_pathways())
#' @param method Character, statistical test: "hypergeometric" (default) or "fisher"
#' @param min_set_size Minimum genes in a set (default: 5)
#' @param max_set_size Maximum genes in a set (default: 500)
#' @param p_adjust_method Method for p-value adjustment: "BH" (default), "bonferroni", "fdr"
#'
#' @return Data.frame with columns:
#'   \item{gene_set}{Name of gene set}
#'   \item{set_size}{Total genes in set}
#'   \item{overlap}{Number of genes overlapping with input}
#'   \item{genes}{Comma-separated list of overlapping genes}
#'   \item{pvalue}{Raw p-value}
#'   \item{padj}{Adjusted p-value}
#'   \item{odds_ratio}{Enrichment odds ratio}
#'   \item{expected}{Expected overlap by chance}
#'   \item{enrichment}{Observed to Expected ratio}
#'
#' @details
#' **Statistical Methods:**
#'
#' *Hypergeometric Test:*
#' - Tests over-representation of gene set in input genes
#' - More conservative, appropriate for most cases
#' - Formula: P(X >= k) where k is observed overlap
#'
#' *Fisher's Exact Test:*
#' - Two-sided test
#' - Tests both enrichment and depletion
#' - Provides odds ratio estimate
#'
#' **Biological Interpretation:**
#' - padj < 0.05: Significant enrichment
#' - enrichment > 1: Over-represented (more than expected by chance)
#' - enrichment < 1: Under-represented
#' - Large overlap with small set_size: Highly specific enrichment
#'
#' @examples
#' # Get significant genes
#' gene_results <- analyze_genes(...)
#' sig_genes <- gene_results$results %>%
#'   filter(padj < 0.05) %>%
#'   pull(gene)
#'
#' # Get pathways
#' pathways <- get_pathways("Hallmark", species = "human")
#'
#' # Test enrichment
#' enrichment <- test_enrichment(
#'   genes = sig_genes,
#'   universe = rownames(gene_results$results),
#'   gene_sets = pathways,
#'   method = "hypergeometric"
#' )
#'
#' # View top enriched pathways
#' head(enrichment)
#'
#' # Visualize
#' plot_enrichment(enrichment, top_n = 20)
#'
#' @seealso \code{\link{get_pathways}}, \code{\link{plot_enrichment}}
#' @export
test_enrichment <- function(genes,
                           universe,
                           gene_sets,
                           method = "hypergeometric",
                           min_set_size = 5,
                           max_set_size = 500,
                           p_adjust_method = "BH") {
  
  # Validate inputs
  if (!is.character(genes)) stop("genes must be a character vector")
  if (!is.character(universe)) stop("universe must be a character vector")
  if (!is.list(gene_sets)) stop("gene_sets must be a named list")
  
  # Remove duplicates
  genes <- unique(genes)
  universe <- unique(universe)
  
  # Filter gene sets by size
  gene_sets <- gene_sets[sapply(gene_sets, length) >= min_set_size]
  gene_sets <- gene_sets[sapply(gene_sets, length) <= max_set_size]
  
  if (length(gene_sets) == 0) {
    stop("No gene sets remain after filtering by size")
  }
  
  # Initialize results
  results <- data.frame(
    gene_set = character(),
    set_size = integer(),
    overlap = integer(),
    genes = character(),
    pvalue = numeric(),
    odds_ratio = numeric(),
    expected = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Test each gene set
  for (set_name in names(gene_sets)) {
    set_genes <- intersect(gene_sets[[set_name]], universe)
    overlap_genes <- intersect(genes, set_genes)
    
    k <- length(overlap_genes)  # observed overlap
    m <- length(set_genes)      # genes in set (in universe)
    n <- length(universe) - m   # genes not in set
    q <- length(genes)          # genes drawn
    
    # Calculate p-value based on method
    if (method == "hypergeometric") {
      # One-sided test (enrichment)
      pval <- phyper(k - 1, m, n, q, lower.tail = FALSE)
      
      # Odds ratio approximation
      odds_ratio <- (k * (n - q + k)) / ((q - k) * (m - k))
      
    } else if (method == "fisher") {
      # Create contingency table
      cont_table <- matrix(
        c(k, q - k, m - k, n - q + k),
        nrow = 2,
        dimnames = list(
          c("In_Set", "Not_In_Set"),
          c("In_Genes", "Not_In_Genes")
        )
      )
      
      fisher_result <- fisher.test(cont_table, alternative = "greater")
      pval <- fisher_result$p.value
      odds_ratio <- fisher_result$estimate
      
    } else {
      stop("method must be 'hypergeometric' or 'fisher'")
    }
    
    # Expected overlap by chance
    expected <- (m * q) / length(universe)
    
    # Store results
    results <- rbind(results, data.frame(
      gene_set = set_name,
      set_size = m,
      overlap = k,
      genes = paste(overlap_genes, collapse = ","),
      pvalue = pval,
      odds_ratio = odds_ratio,
      expected = expected,
      stringsAsFactors = FALSE
    ))
  }
  
  # Adjust p-values
  results$padj <- p.adjust(results$pvalue, method = p_adjust_method)
  
  # Calculate enrichment ratio
  results$enrichment <- results$overlap / results$expected
  
  # Sort by p-value
  results <- results[order(results$pvalue), ]
  rownames(results) <- NULL
  
  # Add metadata
  attr(results, "method") <- method
  attr(results, "n_genes") <- length(genes)
  attr(results, "n_universe") <- length(universe)
  attr(results, "n_gene_sets") <- nrow(results)
  
  class(results) <- c("functionalDE_enrichment", "data.frame")
  
  return(results)
}


#' Plot Enrichment Results
#'
#' @description
#' Visualizes gene set enrichment analysis results
#'
#' @param enrichment_results Output from test_enrichment()
#' @param top_n Number of top gene sets to display (default: 20)
#' @param plot_type Type of plot: "bar" (default), "dot", "lollipop"
#' @param color_by What to color by: "padj", "enrichment", "overlap"
#' @param order_by What to order by: "padj" (default), "enrichment", "overlap"
#' @param pval_threshold Show only gene sets with padj below this (default: 0.05)
#' @param wrap_labels Wrap long gene set names (default: TRUE)
#' @param label_width Maximum label width (default: 40)
#'
#' @return A ggplot2 object
#'
#' @examples
#' enrichment <- test_enrichment(...)
#' 
#' # Bar plot
#' plot_enrichment(enrichment, top_n = 15, plot_type = "bar")
#'
#' # Dot plot colored by enrichment
#' plot_enrichment(enrichment, plot_type = "dot", color_by = "enrichment")
#'
#' @export
plot_enrichment <- function(enrichment_results,
                           top_n = 20,
                           plot_type = "bar",
                           color_by = "padj",
                           order_by = "padj",
                           pval_threshold = 0.05,
                           wrap_labels = TRUE,
                           label_width = 40) {
  
  # Filter and select top results
  plot_data <- enrichment_results %>%
    filter(padj < pval_threshold) %>%
    arrange(!!sym(order_by)) %>%
    head(top_n) %>%
    mutate(
      gene_set = factor(gene_set, levels = rev(gene_set)),
      neg_log10_padj = -log10(padj)
    )
  
  if (nrow(plot_data) == 0) {
    stop("No gene sets meet the p-value threshold")
  }
  
  # Wrap labels if requested
  if (wrap_labels) {
    plot_data$gene_set <- stringr::str_wrap(plot_data$gene_set, width = label_width)
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = gene_set))
  
  # Add geoms based on plot type
  if (plot_type == "bar") {
    p <- p +
      geom_col(aes(y = neg_log10_padj, fill = !!sym(color_by))) +
      coord_flip() +
      labs(y = "-log10(adjusted p-value)", x = "Gene Set")
    
  } else if (plot_type == "dot") {
    p <- p +
      geom_point(aes(y = neg_log10_padj, 
                     size = overlap, 
                     color = !!sym(color_by))) +
      coord_flip() +
      labs(y = "-log10(adjusted p-value)", 
           x = "Gene Set",
           size = "Overlap",
           color = color_by)
    
  } else if (plot_type == "lollipop") {
    p <- p +
      geom_segment(aes(xend = gene_set, y = 0, yend = neg_log10_padj),
                   color = "grey60") +
      geom_point(aes(y = neg_log10_padj, color = !!sym(color_by)),
                size = 3) +
      coord_flip() +
      labs(y = "-log10(adjusted p-value)", 
           x = "Gene Set",
           color = color_by)
  }
  
  # Theme
  p <- p +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}