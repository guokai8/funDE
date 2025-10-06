#' Gene-level Differential Expression Analysis
#'
#' @description
#' Performs differential expression analysis at the gene level using DESeq2,
#' edgeR, or limma. This is the foundation for all higher-level analyses in funDE.
#'
#' @param counts A matrix of raw counts (genes x samples)
#' @param metadata A data.frame containing sample information with rownames matching colnames of counts
#' @param design A formula specifying the design matrix (e.g., ~ condition)
#' @param contrast Character vector of length 3: c("variable", "level1", "level2") where level1 is compared to level2
#' @param method Method for DE analysis: "DESeq2" (default), "edgeR", or "limma"
#' @param filter_low_counts Logical, whether to filter lowly expressed genes (default: TRUE)
#' @param min_count Minimum count threshold for filtering (default: 10)
#' @param min_samples Minimum number of samples that must have min_count (default: 2)
#' @param alpha Significance threshold for adjusted p-values (default: 0.05)
#' @param lfc_threshold Log2 fold change threshold for significance (default: 0)
#' @param transform_method Method for count transformation: "log2", "rlog", "vst" (default: "log2")
#' @param parallel Logical, use parallel processing when available (default: FALSE)
#' @param ... Additional arguments passed to the DE method
#'
#' @return A list of class "functionalDE_result" containing:
#'   \item{results}{Data.frame with DE statistics for each gene}
#'   \item{counts_raw}{Original count matrix}
#'   \item{counts_filtered}{Filtered count matrix}
#'   \item{counts_normalized}{Normalized count matrix}
#'   \item{counts_transformed}{Transformed count matrix for visualization}
#'   \item{metadata}{Sample metadata used}
#'   \item{design}{Design formula used}
#'   \item{contrast}{Contrast specification used}
#'   \item{method}{Method used for analysis}
#'   \item{params}{List of parameters used}
#'   \item{summary}{Summary statistics}
#'
#' @details
#' **Methods:**
#'
#' *DESeq2:*
#' - Recommended for most RNA-seq datasets
#' - Handles count data with negative binomial model
#' - Built-in normalization and shrinkage
#' - Best for datasets with replicates
#'
#' *edgeR:*
#' - Good for datasets with fewer replicates
#' - Quasi-likelihood framework
#' - Robust dispersion estimation
#' - Fast and memory efficient
#'
#' *limma:*
#' - Uses voom transformation for count data
#' - Linear modeling framework
#' - Good for complex experimental designs
#' - Handles missing values well
#'
#' **Results Interpretation:**
#' - log2FoldChange: log2(level1/level2), positive = upregulated in level1
#' - padj: Benjamini-Hochberg adjusted p-values
#' - baseMean: Average normalized expression across all samples
#' - stat: Test statistic from the chosen method
#'
#' @examples
#' # Basic analysis with DESeq2
#' results <- analyze_genes(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control")
#' )
#'
#' # View results
#' head(results$results)
#' summary(results)
#'
#' # Get significant genes
#' sig_genes <- results$results %>%
#'   filter(padj < 0.05, abs(log2FoldChange) > 1)
#'
#' # With edgeR and different parameters
#' results_edger <- analyze_genes(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   method = "edgeR",
#'   min_count = 5,
#'   lfc_threshold = 0.5
#' )
#'
#' @seealso \code{\link{analyze_families}}, \code{\link{analyze_pathways}}
#' @export
analyze_genes <- function(counts,
                         metadata,
                         design,
                         contrast = NULL,
                         method = "DESeq2",
                         filter_low_counts = TRUE,
                         min_count = 10,
                         min_samples = 2,
                         alpha = 0.05,
                         lfc_threshold = 0,
                         transform_method = "log2",
                         parallel = FALSE,
                         ...) {
  
  # Validate inputs
  counts <- validate_counts(counts)
  metadata <- validate_metadata(metadata, counts)
  design <- validate_design(design, metadata)
  contrast <- validate_contrast(contrast, metadata)
  
  # Store original data
  counts_raw <- counts
  
  # Filter low count genes
  if (filter_low_counts) {
    counts <- filter_low_counts(counts, min_count, min_samples)
  }
  
  # Ensure metadata order matches counts
  metadata <- metadata[colnames(counts), , drop = FALSE]
  
  message("Analyzing ", nrow(counts), " genes across ", ncol(counts), " samples using ", method)
  
  # Perform DE analysis based on method
  if (method == "DESeq2") {
    results <- run_deseq2(counts, metadata, design, contrast, alpha, lfc_threshold, parallel, ...)
    
  } else if (method == "edgeR") {
    results <- run_edger(counts, metadata, design, contrast, alpha, lfc_threshold, ...)
    
  } else if (method == "limma") {
    results <- run_limma(counts, metadata, design, contrast, alpha, lfc_threshold, ...)
    
  } else {
    stop("Unsupported method: ", method, ". Choose from: DESeq2, edgeR, limma")
  }
  
  # Add gene names if missing
  if (!"gene" %in% colnames(results)) {
    results$gene <- rownames(results)
  }
  
  # Reorder columns
  standard_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  available_cols <- intersect(standard_cols, colnames(results))
  other_cols <- setdiff(colnames(results), standard_cols)
  results <- results[, c(available_cols, other_cols)]
  
  # Calculate normalized counts
  counts_normalized <- normalize_counts(counts)
  
  # Transform counts for visualization
  counts_transformed <- transform_counts(counts_normalized, method = transform_method)
  
  # Create summary
  summary_stats <- list(
    n_genes_total = nrow(counts_raw),
    n_genes_analyzed = nrow(counts),
    n_samples = ncol(counts),
    n_significant = sum(results$padj < alpha, na.rm = TRUE),
    n_upregulated = sum(results$padj < alpha & results$log2FoldChange > lfc_threshold, na.rm = TRUE),
    n_downregulated = sum(results$padj < alpha & results$log2FoldChange < -lfc_threshold, na.rm = TRUE),
    median_baseMean = median(results$baseMean, na.rm = TRUE),
    median_lfc = median(abs(results$log2FoldChange[results$padj < alpha]), na.rm = TRUE)
  )
  
  # Create result object
  result <- list(
    results = results,
    counts_raw = counts_raw,
    counts_filtered = counts,
    counts_normalized = counts_normalized,
    counts_transformed = counts_transformed,
    metadata = metadata,
    design = design,
    contrast = contrast,
    method = method,
    params = list(
      filter_low_counts = filter_low_counts,
      min_count = min_count,
      min_samples = min_samples,
      alpha = alpha,
      lfc_threshold = lfc_threshold,
      transform_method = transform_method,
      level = "gene"
    ),
    summary = summary_stats
  )
  
  class(result) <- c("functionalDE_result", "list")
  
  return(result)
}

#' Run DESeq2 analysis
#' @keywords internal
run_deseq2 <- function(counts, metadata, design, contrast, alpha, lfc_threshold, parallel, ...) {
  check_package("DESeq2", "DESeq2 analysis")
  
  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design
  )
  
  # Set parallel if requested
  if (parallel) {
    BiocParallel::register(BiocParallel::MulticoreParam())
  }
  
  # Run DESeq2 pipeline
  dds <- DESeq2::DESeq(dds, parallel = parallel, ...)
  
  # Extract results
  if (is.null(contrast)) {
    res <- DESeq2::results(dds, alpha = alpha, lfcThreshold = lfc_threshold)
  } else {
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha, lfcThreshold = lfc_threshold)
  }
  
  # Convert to data.frame
  results_df <- as.data.frame(res)
  results_df$gene <- rownames(results_df)
  
  return(results_df)
}

#' Run edgeR analysis
#' @keywords internal
run_edger <- function(counts, metadata, design, contrast, alpha, lfc_threshold, ...) {
  check_package("edgeR", "edgeR analysis")
  
  # Create DGEList
  dge <- edgeR::DGEList(counts = counts, samples = metadata)
  
  # Calculate normalization factors
  dge <- edgeR::calcNormFactors(dge)
  
  # Create design matrix
  design_matrix <- model.matrix(design, data = metadata)
  
  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design_matrix)
  
  # Fit GLM
  fit <- edgeR::glmQLFit(dge, design_matrix, ...)
  
  # Perform test
  if (is.null(contrast)) {
    # Use last coefficient
    test <- edgeR::glmQLFTest(fit, coef = ncol(design_matrix))
  } else {
    # Make contrast
    contrast_vector <- rep(0, ncol(design_matrix))
    var_name <- contrast[1]
    level1 <- contrast[2]
    level2 <- contrast[3]
    
    # Find contrast coefficients
    coef_names <- colnames(design_matrix)
    level1_coef <- grep(paste0(var_name, level1), coef_names)
    level2_coef <- grep(paste0(var_name, level2), coef_names)
    
    if (length(level1_coef) > 0) contrast_vector[level1_coef] <- 1
    if (length(level2_coef) > 0) contrast_vector[level2_coef] <- -1
    
    test <- edgeR::glmQLFTest(fit, contrast = contrast_vector)
  }
  
  # Extract results
  results_df <- edgeR::topTags(test, n = Inf, adjust.method = "BH")$table
  
  # Standardize column names
  colnames(results_df)[colnames(results_df) == "logFC"] <- "log2FoldChange"
  colnames(results_df)[colnames(results_df) == "logCPM"] <- "baseMean"
  colnames(results_df)[colnames(results_df) == "LR"] <- "stat"
  colnames(results_df)[colnames(results_df) == "PValue"] <- "pvalue"
  colnames(results_df)[colnames(results_df) == "FDR"] <- "padj"
  
  # Add missing columns
  if (!"lfcSE" %in% colnames(results_df)) {
    results_df$lfcSE <- NA
  }
  
  results_df$gene <- rownames(results_df)
  
  return(results_df)
}

#' Run limma analysis
#' @keywords internal
run_limma <- function(counts, metadata, design, contrast, alpha, lfc_threshold, ...) {
  check_package("limma", "limma analysis")
  
  # Create design matrix
  design_matrix <- model.matrix(design, data = metadata)
  
  # Transform counts using voom
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design_matrix, ...)
  
  # Fit linear model
  fit <- limma::lmFit(v, design_matrix)
  
  # Apply empirical Bayes
  fit <- limma::eBayes(fit)
  
  # Extract results
  if (is.null(contrast)) {
    # Use last coefficient
    results_df <- limma::topTable(fit, coef = ncol(design_matrix), number = Inf, adjust.method = "BH")
  } else {
    # Make contrast
    contrast_matrix <- limma::makeContrasts(
      contrasts = paste0(contrast[1], contrast[2], "-", contrast[1], contrast[3]),
      levels = design_matrix
    )
    
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    
    results_df <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
  }
  
  # Standardize column names
  colnames(results_df)[colnames(results_df) == "logFC"] <- "log2FoldChange"
  colnames(results_df)[colnames(results_df) == "AveExpr"] <- "baseMean"
  colnames(results_df)[colnames(results_df) == "t"] <- "stat"
  colnames(results_df)[colnames(results_df) == "P.Value"] <- "pvalue"
  colnames(results_df)[colnames(results_df) == "adj.P.Val"] <- "padj"
  colnames(results_df)[colnames(results_df) == "B"] <- "lfcSE"  # Not exactly SE, but similar
  
  results_df$gene <- rownames(results_df)
  
  return(results_df)
}

#' Print method for functionalDE_result
#' @param x A functionalDE_result object
#' @param ... Additional arguments (ignored)
#' @export
print.functionalDE_result <- function(x, ...) {
  cat("funDE Results:", x$params$level, "level\n")
  
  # Handle method printing for different levels
  if (is.list(x$method)) {
    # For pathway/family level with scoring and DE methods
    cat("Scoring method:", x$method$scoring, "\n")
    cat("DE method:", x$method$de, "\n")
  } else {
    # For gene level with single method
    cat("Method:", x$method, "\n")
  }
  
  # Use appropriate label based on level
  feature_label <- switch(x$params$level,
    "gene" = "Genes",
    "pathway" = "Pathways", 
    "family" = "Gene families",
    "transcript" = "Transcripts",
    "Features"
  )
  
  cat(feature_label, "analyzed:", nrow(x$results), "\n")
  cat("Significant", tolower(feature_label), "(padj <", x$params$alpha, "):", x$summary$n_significant, "\n")
  cat("  Upregulated:", x$summary$n_upregulated, "\n")
  cat("  Downregulated:", x$summary$n_downregulated, "\n")
  
  if (!is.null(x$contrast)) {
    cat("Contrast:", paste(x$contrast, collapse = " "), "\n")
  }
}

#' Summary method for functionalDE_result
#' @param object A functionalDE_result object
#' @param ... Additional arguments (ignored)
#' @export
summary.functionalDE_result <- function(object, ...) {
  print(object)
  cat("\nTop 10 significant genes:\n")
  
  sig_genes <- object$results %>%
    filter(padj < object$params$alpha) %>%
    arrange(padj) %>%
    head(10)
  
  if (nrow(sig_genes) > 0) {
    print(sig_genes[, c("gene", "log2FoldChange", "padj")])
  } else {
    cat("No significant genes found\n")
  }
}