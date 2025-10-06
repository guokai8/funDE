#' Filter low count genes and features
#' @param counts Count matrix
#' @param min_count Minimum count threshold
#' @param min_samples Minimum number of samples
#' @return Filtered count matrix
#' @keywords internal
filter_low_counts <- function(counts, min_count = 10, min_samples = 2) {
  keep <- rowSums(counts >= min_count) >= min_samples
  
  if (sum(keep) == 0) {
    stop("No features pass filtering criteria")
  }
  
  message("Filtered from ", nrow(counts), " to ", sum(keep), " features")
  
  return(counts[keep, , drop = FALSE])
}

#' Calculate size factors using DESeq2 method
#' @param counts Count matrix
#' @return Size factors
#' @keywords internal
calculate_size_factors <- function(counts) {
  # Geometric mean of each gene across samples
  geo_means <- exp(rowMeans(log(counts + 1)))
  
  # Exclude genes with zero geometric mean
  geo_means[geo_means == 0] <- NA
  
  # Calculate size factors for each sample
  size_factors <- apply(counts, 2, function(cnts) {
    median((cnts / geo_means)[is.finite(geo_means) & geo_means > 0], na.rm = TRUE)
  })
  
  # Replace any zero or NA size factors with 1
  size_factors[is.na(size_factors) | size_factors <= 0] <- 1
  
  return(size_factors)
}

#' Normalize counts using size factors
#' @param counts Count matrix
#' @param size_factors Size factors
#' @return Normalized counts
#' @keywords internal
normalize_counts <- function(counts, size_factors = NULL) {
  if (is.null(size_factors)) {
    size_factors <- calculate_size_factors(counts)
  }
  
  normalized <- sweep(counts, 2, size_factors, "/")
  
  return(normalized)
}

#' Perform variance stabilizing transformation
#' @param counts Count matrix
#' @param method Method: "log2", "rlog", "vst"
#' @return Transformed matrix
#' @keywords internal
transform_counts <- function(counts, method = "log2") {
  if (method == "log2") {
    # Simple log2 transformation with pseudocount
    transformed <- log2(counts + 1)
    
  } else if (method == "rlog") {
    check_package("DESeq2", "rlog transformation")
    
    # Create DESeq2 object for rlog
    coldata <- data.frame(condition = factor(rep("A", ncol(counts))))
    rownames(coldata) <- colnames(counts)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = coldata,
      design = ~ 1
    )
    
    transformed <- DESeq2::rlog(dds, blind = TRUE)
    transformed <- SummarizedExperiment::assay(transformed)
    
  } else if (method == "vst") {
    check_package("DESeq2", "VST transformation")
    
    # Create DESeq2 object for VST
    coldata <- data.frame(condition = factor(rep("A", ncol(counts))))
    rownames(coldata) <- colnames(counts)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = coldata,
      design = ~ 1
    )
    
    transformed <- DESeq2::vst(dds, blind = TRUE)
    transformed <- SummarizedExperiment::assay(transformed)
    
  } else {
    stop("Unsupported transformation method: ", method)
  }
  
  return(transformed)
}

#' Calculate correlation between two vectors
#' @param x First vector
#' @param y Second vector
#' @param method Correlation method
#' @return Correlation coefficient
#' @keywords internal
calculate_correlation <- function(x, y, method = "pearson") {
  if (length(x) != length(y)) {
    stop("Vectors must have same length")
  }
  
  # Remove NAs
  complete_cases <- complete.cases(x, y)
  if (sum(complete_cases) < 3) {
    return(NA)
  }
  
  cor(x[complete_cases], y[complete_cases], method = method)
}

#' Calculate robust statistics
#' @param x Numeric vector
#' @return List of robust statistics
#' @keywords internal
robust_stats <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(list(
      median = NA,
      mad = NA,
      q25 = NA,
      q75 = NA,
      iqr = NA
    ))
  }
  
  list(
    median = median(x),
    mad = mad(x),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75),
    iqr = IQR(x)
  )
}

#' Aggregate genes to functional units
#' @param counts Count matrix (genes x samples)
#' @param grouping Named vector mapping genes to groups
#' @param method Aggregation method: "mean", "median", "sum"
#' @return Aggregated count matrix
#' @keywords internal
aggregate_genes <- function(counts, grouping, method = "mean") {
  # Match genes
  common_genes <- intersect(rownames(counts), names(grouping))
  
  if (length(common_genes) == 0) {
    stop("No genes found in both counts and grouping")
  }
  
  counts_subset <- counts[common_genes, , drop = FALSE]
  grouping_subset <- grouping[common_genes]
  
  # Aggregate by group
  groups <- unique(grouping_subset)
  aggregated <- matrix(0, nrow = length(groups), ncol = ncol(counts))
  rownames(aggregated) <- groups
  colnames(aggregated) <- colnames(counts)
  
  for (group in groups) {
    group_genes <- names(grouping_subset)[grouping_subset == group]
    
    if (length(group_genes) == 1) {
      aggregated[group, ] <- counts_subset[group_genes, ]
    } else {
      group_counts <- counts_subset[group_genes, , drop = FALSE]
      
      if (method == "mean") {
        aggregated[group, ] <- round(colMeans(group_counts))
      } else if (method == "median") {
        aggregated[group, ] <- round(apply(group_counts, 2, median))
      } else if (method == "sum") {
        aggregated[group, ] <- colSums(group_counts)
      } else {
        stop("Unsupported aggregation method: ", method)
      }
    }
  }
  
  # Ensure integer matrix for DESeq2 compatibility
  mode(aggregated) <- "integer"
  
  return(aggregated)
}

#' Calculate pathway activity scores
#' @param expression_matrix Expression matrix (genes x samples)
#' @param gene_sets Named list of gene sets
#' @param method Scoring method: "mean", "gsva", "ssgsea"
#' @return Pathway scores matrix
#' @keywords internal
calculate_pathway_scores <- function(expression_matrix, gene_sets, method = "mean") {
  # Ensure gene sets contain genes in expression matrix
  available_genes <- rownames(expression_matrix)
  gene_sets <- lapply(gene_sets, function(genes) {
    intersect(genes, available_genes)
  })
  
  # Remove empty gene sets
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  
  if (length(gene_sets) == 0) {
    stop("No gene sets have genes in expression matrix")
  }
  
  if (method == "mean") {
    # Simple mean expression
    scores <- matrix(0, nrow = length(gene_sets), ncol = ncol(expression_matrix))
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
    
  } else if (method == "gsva") {
    check_package("GSVA", "GSVA scoring")
    scores <- GSVA::gsva(expression_matrix, gene_sets, method = "gsva")
    
  } else if (method == "ssgsea") {
    check_package("GSVA", "ssGSEA scoring")
    scores <- GSVA::gsva(expression_matrix, gene_sets, method = "ssgsea")
    
  } else {
    stop("Unsupported scoring method: ", method)
  }
  
  return(scores)
}