#' Pathway-level Differential Expression Analysis
#'
#' @description
#' Performs differential expression analysis at the pathway level by first
#' calculating pathway activity scores, then testing for differential activity
#' between conditions. This approach can detect coordinated changes in gene
#' expression that might be missed at the individual gene level.
#'
#' @param counts A matrix of raw counts (genes x samples)
#' @param metadata A data.frame containing sample information
#' @param design A formula specifying the design matrix
#' @param contrast Character vector of length 3: c("variable", "level1", "level2")
#' @param pathways Character database name (see get_pathways()) or named list of gene sets
#' @param species Character, "human" or "mouse" (used when pathways is a database name)
#' @param scoring_method Method for pathway activity scoring: "GSVA" (default), "ssGSEA", "AUCell", "mean"
#' @param de_method Method for differential testing: "limma" (default), "t.test", "wilcox"
#' @param min_pathway_size Minimum genes required in pathway (default: 5)
#' @param max_pathway_size Maximum genes allowed in pathway (default: 500)
#' @param transform_method Count transformation before scoring: "log2", "vst", "rlog"
#' @param alpha Significance threshold (default: 0.05)
#' @param ... Additional arguments passed to scoring method
#'
#' @return A list of class "functionalDE_result" containing:
#'   \item{results}{Data.frame with pathway-level DE statistics including:
#'     - pathway: Pathway name
#'     - n_genes: Number of genes in pathway found in data
#'     - genes: Comma-separated list of genes in pathway (present in data)
#'     - log2FoldChange: Pathway-level fold change
#'     - pvalue: Raw p-value
#'     - padj: Adjusted p-value}
#'   \item{pathway_scores}{Matrix of pathway activity scores}
#'   \item{pathways_used}{List of pathways used in analysis}
#'   \item{counts_transformed}{Transformed count matrix used}
#'   \item{metadata}{Sample metadata}
#'   \item{method}{Methods used (scoring and DE)}
#'   \item{params}{Parameters used}
#'   \item{summary}{Summary statistics}
#'
#' @details
#' **Scoring Methods:**
#'
#' *GSVA (Gene Set Variation Analysis):*
#' - Non-parametric, rank-based method
#' - Good for detecting pathway-level changes
#' - Accounts for correlation structure
#' - Recommended for most analyses
#'
#' *ssGSEA (single-sample GSEA):*
#' - Single-sample pathway enrichment scores
#' - Normalized enrichment scores per sample
#' - Good for comparing across samples
#'
#' *AUCell:*
#' - Area Under the Curve approach
#' - Robust to outliers
#' - Good for single-cell data
#'
#' *Mean Expression:*
#' - Simple average expression of pathway genes
#' - Fast and interpretable
#' - Good for initial exploration
#'
#' **Biological Interpretation:**
#' - Positive log2FoldChange: Higher pathway activity in level1 vs level2
#' - Pathway activity represents coordinated expression changes
#' - Significant pathways suggest functional changes between conditions
#'
#' @examples
#' # Basic pathway analysis with Hallmark pathways
#' pathway_results <- analyze_pathways(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   pathways = "Hallmark",
#'   scoring_method = "GSVA"
#' )
#'
#' # View results
#' head(pathway_results$results)
#' summary(pathway_results)
#'
#' # With custom pathways and different scoring
#' custom_pathways <- get_pathways("KEGG", species = "mouse")
#' pathway_results <- analyze_pathways(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   pathways = custom_pathways,
#'   scoring_method = "AUCell",
#'   min_pathway_size = 10
#' )
#'
#' # Visualize pathway activity
#' plot_pathway_heatmap(pathway_results, top_n = 20)
#'
#' @seealso \code{\link{get_pathways}}, \code{\link{score_pathway_activity}}, \code{\link{plot_pathway_heatmap}}
#' @export
analyze_pathways <- function(counts,
                            metadata,
                            design,
                            contrast = NULL,
                            pathways = "Hallmark",
                            species = "human",
                            scoring_method = "GSVA",
                            de_method = "limma",
                            min_pathway_size = 5,
                            max_pathway_size = 500,
                            transform_method = "log2",
                            alpha = 0.05,
                            ...) {
  
  # Validate inputs
  counts <- validate_counts(counts)
  metadata <- validate_metadata(metadata, counts)
  design <- validate_design(design, metadata)
  contrast <- validate_contrast(contrast, metadata)
  
  # Get pathways if database name provided
  if (is.character(pathways) && length(pathways) == 1) {
    message("Downloading ", pathways, " pathways...")
    pathways <- get_pathways(
      database = pathways,
      species = species,
      min_size = min_pathway_size,
      max_size = max_pathway_size
    )
  } else {
    pathways <- validate_gene_sets(pathways, min_pathway_size, max_pathway_size)
  }
  
  # Transform counts
  message("Transforming counts using ", transform_method, "...")
  if (transform_method == "log2") {
    counts_transformed <- log2(counts + 1)
  } else {
    counts_transformed <- transform_counts(counts, method = transform_method)
  }
  
  # Calculate pathway scores
  message("Calculating pathway scores using ", scoring_method, "...")
  pathway_scores <- score_pathway_activity(
    expression_matrix = counts_transformed,
    gene_sets = pathways,
    method = scoring_method,
    ...
  )
  
  # Ensure metadata order matches scores
  metadata <- metadata[colnames(pathway_scores), , drop = FALSE]
  
  # Perform differential analysis on pathway scores
  message("Testing for differential pathway activity using ", de_method, "...")
  
  if (de_method == "limma") {
    results <- test_pathways_limma(pathway_scores, metadata, design, contrast)
    
  } else if (de_method == "t.test") {
    results <- test_pathways_ttest(pathway_scores, metadata, contrast)
    
  } else if (de_method == "wilcox") {
    results <- test_pathways_wilcox(pathway_scores, metadata, contrast)
    
  } else {
    stop("Unsupported de_method: ", de_method)
  }
  
  # Add pathway information
  results$pathway <- rownames(results)
  results$n_genes <- sapply(pathways[rownames(results)], length)
  
  # Add genes present in the expression data for each pathway
  results$genes <- sapply(rownames(results), function(pathway_name) {
    pathway_genes <- pathways[[pathway_name]]
    genes_in_data <- intersect(pathway_genes, rownames(counts_transformed))
    paste(genes_in_data, collapse = ",")
  })
  
  # Reorder columns
  standard_cols <- c("pathway", "n_genes", "genes", "log2FoldChange", "pvalue", "padj")
  available_cols <- intersect(standard_cols, colnames(results))
  other_cols <- setdiff(colnames(results), standard_cols)
  results <- results[, c(available_cols, other_cols)]
  
  # Calculate summary statistics
  summary_stats <- list(
    n_pathways_tested = nrow(results),
    n_pathways_significant = sum(results$padj < alpha, na.rm = TRUE),
    n_upregulated = sum(results$padj < alpha & results$log2FoldChange > 0, na.rm = TRUE),
    n_downregulated = sum(results$padj < alpha & results$log2FoldChange < 0, na.rm = TRUE),
    median_n_genes = median(results$n_genes),
    scoring_method = scoring_method,
    de_method = de_method
  )
  
  # Create result object
  result <- list(
    results = results,
    pathway_scores = pathway_scores,
    pathways_used = pathways,
    counts_transformed = counts_transformed,
    metadata = metadata,
    design = design,
    contrast = contrast,
    method = list(scoring = scoring_method, de = de_method),
    params = list(
      min_pathway_size = min_pathway_size,
      max_pathway_size = max_pathway_size,
      transform_method = transform_method,
      alpha = alpha,
      level = "pathway"
    ),
    summary = summary_stats
  )
  
  class(result) <- c("functionalDE_result", "list")
  
  return(result)
}

#' Test pathway differential activity using limma
#' @keywords internal
test_pathways_limma <- function(pathway_scores, metadata, design, contrast) {
  check_package("limma", "limma-based pathway testing")
  
  # Create design matrix
  design_matrix <- model.matrix(design, data = metadata)
  
  # Fit linear model
  fit <- limma::lmFit(pathway_scores, design_matrix)
  fit <- limma::eBayes(fit)
  
  # Store original rownames to preserve them
  original_rownames <- rownames(pathway_scores)
  
  # Extract results
  if (is.null(contrast)) {
    # Use last coefficient
    results_df <- limma::topTable(fit, coef = ncol(design_matrix), 
                                  number = Inf, adjust.method = "BH")
  } else {
    # Create contrast
    var_name <- contrast[1]
    level1 <- contrast[2]
    level2 <- contrast[3]
    
    # Find contrast in design matrix
    contrast_name <- paste0(var_name, level1, " - ", var_name, level2)
    
    # Try to make contrast
    results_df <- tryCatch({
      contrast_matrix <- limma::makeContrasts(
        contrasts = contrast_name,
        levels = design_matrix
      )
      
      fit2 <- limma::contrasts.fit(fit, contrast_matrix)
      fit2 <- limma::eBayes(fit2)
      
      limma::topTable(fit2, number = Inf, adjust.method = "BH")
    }, error = function(e) {
      # Fallback to coefficient if contrast fails
      coef_idx <- grep(level1, colnames(design_matrix))
      if (length(coef_idx) == 0) coef_idx <- ncol(design_matrix)
      
      limma::topTable(fit, coef = coef_idx, 
                      number = Inf, adjust.method = "BH")
    })
  }
  
  # Ensure rownames are preserved correctly
  if (nrow(results_df) == length(original_rownames)) {
    rownames(results_df) <- original_rownames
  }
  
  # Standardize column names
  colnames(results_df)[colnames(results_df) == "logFC"] <- "log2FoldChange"
  colnames(results_df)[colnames(results_df) == "P.Value"] <- "pvalue"
  colnames(results_df)[colnames(results_df) == "adj.P.Val"] <- "padj"
  colnames(results_df)[colnames(results_df) == "t"] <- "statistic"
  
  return(results_df)
}

#' Test pathway differential activity using t-test
#' @keywords internal
test_pathways_ttest <- function(pathway_scores, metadata, contrast) {
  if (is.null(contrast)) {
    stop("contrast must be specified for t-test method")
  }
  
  var_name <- contrast[1]
  level1 <- contrast[2]
  level2 <- contrast[3]
  
  # Get sample groups
  group1_samples <- rownames(metadata)[metadata[[var_name]] == level1]
  group2_samples <- rownames(metadata)[metadata[[var_name]] == level2]
  
  # Initialize results
  n_pathways <- nrow(pathway_scores)
  results_df <- data.frame(
    log2FoldChange = numeric(n_pathways),
    pvalue = numeric(n_pathways),
    statistic = numeric(n_pathways),
    row.names = rownames(pathway_scores)
  )
  
  # Test each pathway
  for (i in seq_len(n_pathways)) {
    pathway_scores_i <- pathway_scores[i, ]
    
    group1_scores <- pathway_scores_i[group1_samples]
    group2_scores <- pathway_scores_i[group2_samples]
    
    # Perform t-test
    test_result <- t.test(group1_scores, group2_scores)
    
    results_df$log2FoldChange[i] <- mean(group1_scores) - mean(group2_scores)
    results_df$pvalue[i] <- test_result$p.value
    results_df$statistic[i] <- test_result$statistic
  }
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  return(results_df)
}

#' Test pathway differential activity using Wilcoxon test
#' @keywords internal
test_pathways_wilcox <- function(pathway_scores, metadata, contrast) {
  if (is.null(contrast)) {
    stop("contrast must be specified for wilcox method")
  }
  
  var_name <- contrast[1]
  level1 <- contrast[2]
  level2 <- contrast[3]
  
  # Get sample groups
  group1_samples <- rownames(metadata)[metadata[[var_name]] == level1]
  group2_samples <- rownames(metadata)[metadata[[var_name]] == level2]
  
  # Initialize results
  n_pathways <- nrow(pathway_scores)
  results_df <- data.frame(
    log2FoldChange = numeric(n_pathways),
    pvalue = numeric(n_pathways),
    statistic = numeric(n_pathways),
    row.names = rownames(pathway_scores)
  )
  
  # Test each pathway
  for (i in seq_len(n_pathways)) {
    pathway_scores_i <- pathway_scores[i, ]
    
    group1_scores <- pathway_scores_i[group1_samples]
    group2_scores <- pathway_scores_i[group2_samples]
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(group1_scores, group2_scores)
    
    results_df$log2FoldChange[i] <- median(group1_scores) - median(group2_scores)
    results_df$pvalue[i] <- test_result$p.value
    results_df$statistic[i] <- test_result$statistic
  }
  
  # Adjust p-values
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  return(results_df)
}