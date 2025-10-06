#' Calculate Pathway Activity Scores
#'
#' @description
#' Calculates pathway activity scores from gene expression data using various
#' scoring methods. This is used internally by analyze_pathways() but can also
#' be used standalone for custom analyses.
#'
#' @param expression_matrix A matrix of expression values (genes x samples)
#' @param gene_sets Named list of gene sets (pathways)
#' @param method Scoring method: "GSVA" (default), "ssGSEA", "AUCell", "mean", "median", "maxmean"
#' @param min_set_size Minimum genes required in a gene set (default: 1)
#' @param max_set_size Maximum genes allowed in a gene set (default: Inf)
#' @param parallel Logical, use parallel processing (default: FALSE)
#' @param ... Additional arguments passed to specific scoring methods
#'
#' @return Matrix of pathway scores (pathways x samples)
#'
#' @details
#' **Scoring Methods:**
#'
#' *GSVA (Gene Set Variation Analysis):*
#' - Non-parametric, rank-based method
#' - Estimates variation of pathway activity over sample population
#' - Good for detecting subtle coordinated changes
#' - Values roughly centered around 0
#'
#' *ssGSEA (single-sample Gene Set Enrichment Analysis):*
#' - Sample-level enrichment scores
#' - Based on Kolmogorov-Smirnov statistic
#' - Values between -1 and 1
#' - Good for comparing pathway activity across samples
#'
#' *AUCell (Area Under the Curve):*
#' - Area under the recovery curve
#' - Robust to dropouts (good for single-cell)
#' - Values between 0 and 1
#' - Less sensitive to outliers
#'
#' *Mean Expression:*
#' - Simple average expression of pathway genes
#' - Fast and interpretable
#' - Good for initial exploration
#' - Values in same scale as input expression
#'
#' *Median Expression:*
#' - Median expression of pathway genes
#' - More robust to outliers than mean
#' - Good for noisy data
#'
#' *MaxMean:*
#' - Mean of top 50% expressed genes in pathway
#' - Reduces noise from lowly expressed genes
#' - Good compromise between mean and maximum
#'
#' @examples
#' # Get pathway definitions
#' pathways <- get_pathways("Hallmark", species = "human")
#'
#' # Calculate GSVA scores
#' gsva_scores <- score_pathway_activity(
#'   expression_matrix = log2(counts + 1),
#'   gene_sets = pathways,
#'   method = "GSVA"
#' )
#'
#' # Calculate AUCell scores (good for single-cell)
#' aucell_scores <- score_pathway_activity(
#'   expression_matrix = log2(counts + 1),
#'   gene_sets = pathways,
#'   method = "AUCell"
#' )
#'
#' # Simple mean expression scores
#' mean_scores <- score_pathway_activity(
#'   expression_matrix = log2(counts + 1),
#'   gene_sets = pathways,
#'   method = "mean"
#' )
#'
#' @seealso \code{\link{analyze_pathways}}, \code{\link{get_pathways}}
#' @export
score_pathway_activity <- function(expression_matrix,
                                  gene_sets,
                                  method = "GSVA",
                                  min_set_size = 1,
                                  max_set_size = Inf,
                                  parallel = FALSE,
                                  ...) {
  
  # Validate inputs
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data.frame")
  }
  
  if (is.data.frame(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }
  
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list")
  }
  
  # Filter gene sets by size
  set_sizes <- sapply(gene_sets, length)
  gene_sets <- gene_sets[set_sizes >= min_set_size & set_sizes <= max_set_size]
  
  if (length(gene_sets) == 0) {
    stop("No gene sets remain after size filtering")
  }
  
  # Filter gene sets to only include genes in expression matrix
  available_genes <- rownames(expression_matrix)
  gene_sets <- lapply(gene_sets, function(genes) {
    intersect(genes, available_genes)
  })
  
  # Remove empty gene sets
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  
  if (length(gene_sets) == 0) {
    stop("No gene sets have genes in expression matrix")
  }
  
  message("Scoring ", length(gene_sets), " pathways using ", method, " method")
  
  # Calculate scores based on method
  if (method == "GSVA") {
    scores <- calculate_gsva_scores(expression_matrix, gene_sets, parallel, ...)
    
  } else if (method == "ssGSEA") {
    scores <- calculate_ssgsea_scores(expression_matrix, gene_sets, parallel, ...)
    
  } else if (method == "AUCell") {
    scores <- calculate_aucell_scores(expression_matrix, gene_sets, ...)
    
  } else if (method == "mean") {
    scores <- calculate_mean_scores(expression_matrix, gene_sets)
    
  } else if (method == "median") {
    scores <- calculate_median_scores(expression_matrix, gene_sets)
    
  } else if (method == "maxmean") {
    scores <- calculate_maxmean_scores(expression_matrix, gene_sets)
    
  } else {
    stop("Unsupported scoring method: ", method, 
         "\nAvailable methods: GSVA, ssGSEA, AUCell, mean, median, maxmean")
  }
  
  # Add attributes
  attr(scores, "method") <- method
  attr(scores, "n_pathways") <- nrow(scores)
  attr(scores, "n_samples") <- ncol(scores)
  attr(scores, "pathway_sizes") <- sapply(gene_sets, length)
  
  return(scores)
}

#' Calculate GSVA scores
#' @keywords internal
calculate_gsva_scores <- function(expression_matrix, gene_sets, parallel, ...) {
  check_package("GSVA", "GSVA scoring")
  
  if (parallel) {
    BiocParallel::register(BiocParallel::MulticoreParam())
  }
  
  # Handle different GSVA versions
  scores <- tryCatch({
    # Try new GSVA interface (v1.50+)
    param <- GSVA::gsvaParam(expression_matrix, gene_sets)
    GSVA::gsva(param, verbose = FALSE, ...)
  }, error = function(e) {
    # Fallback to old GSVA interface  
    GSVA::gsva(
      expr = expression_matrix,
      gset.idx.list = gene_sets,
      method = "gsva",
      parallel.sz = if (parallel) 0 else 1,
      verbose = FALSE,
      ...
    )
  })
  
  return(scores)
}

#' Calculate ssGSEA scores
#' @keywords internal
calculate_ssgsea_scores <- function(expression_matrix, gene_sets, parallel, ...) {
  check_package("GSVA", "ssGSEA scoring")
  
  if (parallel) {
    BiocParallel::register(BiocParallel::MulticoreParam())
  }
  
  # Handle different GSVA versions
  scores <- tryCatch({
    # Try new GSVA interface (v1.50+)
    param <- GSVA::ssgseaParam(expression_matrix, gene_sets)
    GSVA::gsva(param, verbose = FALSE, ...)
  }, error = function(e) {
    # Fallback to old GSVA interface
    GSVA::gsva(
      expr = expression_matrix,
      gset.idx.list = gene_sets,
      method = "ssgsea",
      parallel.sz = if (parallel) 0 else 1,
      verbose = FALSE,
      ...
    )
  })
  
  return(scores)
}

#' Calculate AUCell scores
#' @keywords internal
calculate_aucell_scores <- function(expression_matrix, gene_sets, ...) {
  check_package("AUCell", "AUCell scoring")
  
  # Build gene rankings
  rankings <- AUCell::AUCell_buildRankings(expression_matrix, ...)
  
  # Calculate AUC scores
  auc_scores <- AUCell::AUCell_calcAUC(gene_sets, rankings, ...)
  
  # Extract scores matrix
  scores <- AUCell::getAUC(auc_scores)
  
  return(scores)
}

#' Calculate mean expression scores
#' @keywords internal
calculate_mean_scores <- function(expression_matrix, gene_sets) {
  n_pathways <- length(gene_sets)
  n_samples <- ncol(expression_matrix)
  
  scores <- matrix(0, nrow = n_pathways, ncol = n_samples)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expression_matrix)
  
  for (i in seq_along(gene_sets)) {
    pathway_genes <- gene_sets[[i]]
    
    if (length(pathway_genes) == 1) {
      scores[i, ] <- expression_matrix[pathway_genes, ]
    } else {
      scores[i, ] <- colMeans(expression_matrix[pathway_genes, , drop = FALSE])
    }
  }
  
  return(scores)
}

#' Calculate median expression scores
#' @keywords internal
calculate_median_scores <- function(expression_matrix, gene_sets) {
  n_pathways <- length(gene_sets)
  n_samples <- ncol(expression_matrix)
  
  scores <- matrix(0, nrow = n_pathways, ncol = n_samples)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expression_matrix)
  
  for (i in seq_along(gene_sets)) {
    pathway_genes <- gene_sets[[i]]
    
    if (length(pathway_genes) == 1) {
      scores[i, ] <- expression_matrix[pathway_genes, ]
    } else {
      scores[i, ] <- apply(expression_matrix[pathway_genes, , drop = FALSE], 2, median)
    }
  }
  
  return(scores)
}

#' Calculate maxmean scores (mean of top 50% genes)
#' @param expression_matrix Expression matrix (genes x samples)
#' @param gene_sets Named list of gene sets
#' @return Matrix of pathway scores
#' @keywords internal
#' @noRd
calculate_maxmean_scores <- function(expression_matrix, gene_sets) {
  n_pathways <- length(gene_sets)
  n_samples <- ncol(expression_matrix)
  
  scores <- matrix(0, nrow = n_pathways, ncol = n_samples)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expression_matrix)
  
  for (i in seq_along(gene_sets)) {
    pathway_genes <- gene_sets[[i]]
    
    if (length(pathway_genes) == 1) {
      scores[i, ] <- expression_matrix[pathway_genes, ]
    } else {
      pathway_expr <- expression_matrix[pathway_genes, , drop = FALSE]
      n_top <- ceiling(length(pathway_genes) / 2)
      
      # For each sample, take mean of top 50% expressed genes
      for (j in seq_len(n_samples)) {
        sample_expr <- pathway_expr[, j]
        top_genes <- order(sample_expr, decreasing = TRUE)[1:n_top]
        scores[i, j] <- mean(sample_expr[top_genes])
      }
    }
  }
  
  return(scores)
}

#' Compare pathway scoring methods
#'
#' @description
#' Compares different pathway scoring methods on the same data to help
#' users choose the most appropriate method for their analysis.
#'
#' @param expression_matrix Expression matrix (genes x samples)
#' @param gene_sets Named list of gene sets
#' @param methods Character vector of methods to compare (default: c("GSVA", "ssGSEA", "mean"))
#' @param sample_pathways Number of pathways to use for comparison (default: 10)
#' @param plot Logical, create comparison plot (default: TRUE)
#'
#' @return List containing scores from each method and correlation matrix
#'
#' @examples
#' pathways <- get_pathways("Hallmark", species = "human")
#' comparison <- compare_scoring_methods(
#'   expression_matrix = log2(counts + 1),
#'   gene_sets = pathways,
#'   methods = c("GSVA", "ssGSEA", "mean", "AUCell")
#' )
#'
#' @export
compare_scoring_methods <- function(expression_matrix,
                                   gene_sets,
                                   methods = c("GSVA", "ssGSEA", "mean"),
                                   sample_pathways = 10,
                                   plot = TRUE) {
  
  # Sample pathways if too many
  if (length(gene_sets) > sample_pathways) {
    set.seed(42)  # For reproducibility
    sampled_pathways <- sample(gene_sets, sample_pathways)
  } else {
    sampled_pathways <- gene_sets
  }
  
  # Calculate scores for each method
  all_scores <- list()
  
  for (method in methods) {
    message("Calculating scores using ", method)
    
    tryCatch({
      scores <- score_pathway_activity(
        expression_matrix = expression_matrix,
        gene_sets = sampled_pathways,
        method = method
      )
      all_scores[[method]] <- scores
    }, error = function(e) {
      warning("Failed to calculate scores for method ", method, ": ", e$message)
    })
  }
  
  # Calculate correlations between methods
  if (length(all_scores) >= 2) {
    method_names <- names(all_scores)
    cor_matrix <- matrix(1, nrow = length(method_names), ncol = length(method_names))
    rownames(cor_matrix) <- colnames(cor_matrix) <- method_names
    
    for (i in 1:(length(method_names) - 1)) {
      for (j in (i + 1):length(method_names)) {
        method1 <- method_names[i]
        method2 <- method_names[j]
        
        # Calculate correlation for overlapping pathways
        common_pathways <- intersect(rownames(all_scores[[method1]]), 
                                    rownames(all_scores[[method2]]))
        
        if (length(common_pathways) > 0) {
          scores1 <- as.vector(all_scores[[method1]][common_pathways, ])
          scores2 <- as.vector(all_scores[[method2]][common_pathways, ])
          
          cor_val <- cor(scores1, scores2, use = "complete.obs")
          cor_matrix[i, j] <- cor_matrix[j, i] <- cor_val
        }
      }
    }
  } else {
    cor_matrix <- NULL
  }
  
  result <- list(
    scores = all_scores,
    correlations = cor_matrix,
    pathways_used = names(sampled_pathways)
  )
  
  if (plot && !is.null(cor_matrix)) {
    # Simple correlation heatmap
    print("Method Correlations:")
    print(round(cor_matrix, 3))
  }
  
  return(result)
}