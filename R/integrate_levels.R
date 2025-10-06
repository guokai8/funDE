#' Integrate Results Across Multiple Functional Levels
#'
#' @description
#' Integrates and prioritizes findings across different functional levels
#' (genes, families, pathways) to provide comprehensive biological insights
#' and identify the most robust and significant patterns.
#'
#' @param results_list Named list of funDE results from different levels
#' @param integration_method Method for integration:
#'   - "rank_aggregation": Combine ranks across levels using RRA
#'   - "evidence_scoring": Score features based on evidence strength
#'   - "consensus": Find consensus significant features
#'   - "pathway_guided": Use pathway structure to guide integration
#' @param weights Named vector of weights for each level (optional)
#' @param alpha Significance threshold (default: 0.05)
#' @param min_levels Minimum levels a feature must be significant in (default: 1)
#' @param annotation Optional annotation data.frame for mapping between levels
#' @param pathway_db Pathway database for pathway-guided integration
#' @param ... Additional arguments for specific integration methods
#'
#' @return A list containing:
#'   \item{integrated_results}{Integrated prioritized results}
#'   \item{level_comparison}{Comparison across levels}
#'   \item{consensus_features}{Features significant across multiple levels}
#'   \item{method_specific}{Method-specific additional results}
#'   \item{summary}{Integration summary statistics}
#'
#' @details
#' **Integration Methods:**
#'
#' *Rank Aggregation:*
#' - Combines feature rankings across levels using Robust Rank Aggregation (RRA)
#' - Identifies features consistently highly ranked across levels
#' - Good for finding robust signals regardless of statistical method
#'
#' *Evidence Scoring:*
#' - Assigns evidence scores based on significance, effect size, and consistency
#' - Weights can be applied to prioritize certain levels
#' - Produces a unified ranking of biological importance
#'
#' *Consensus:*
#' - Identifies features significant in multiple levels
#' - Simple intersection-based approach
#' - Good for high-confidence findings
#'
#' *Pathway-Guided:*
#' - Uses pathway structure to inform integration
#' - Gene-level evidence informs pathway-level conclusions
#' - Family-level patterns guide gene interpretation
#'
#' **Biological Interpretation:**
#' - Features significant at multiple levels: High confidence
#' - Pathway-level only: Coordinated subtle changes
#' - Gene-level only: Strong individual effects
#' - Family-level: Functional redundancy or compensation
#'
#' @examples
#' # Run analyses at different levels
#' gene_results <- analyze_genes(...)
#' family_results <- analyze_families(...)
#' pathway_results <- analyze_pathways(...)
#'
#' results_list <- list(
#'   genes = gene_results,
#'   families = family_results,
#'   pathways = pathway_results
#' )
#'
#' # Evidence-based integration
#' integrated <- integrate_levels(
#'   results_list,
#'   integration_method = "evidence_scoring",
#'   weights = c(genes = 1, families = 1.5, pathways = 2)
#' )
#'
#' # Consensus approach
#' consensus <- integrate_levels(
#'   results_list,
#'   integration_method = "consensus",
#'   min_levels = 2
#' )
#'
#' # Pathway-guided integration
#' pathway_guided <- integrate_levels(
#'   results_list,
#'   integration_method = "pathway_guided",
#'   pathway_db = get_pathways("Hallmark")
#' )
#'
#' @seealso \code{\link{plot_multilevel}}, \code{\link{analyze_genes}}, \code{\link{analyze_pathways}}
#' @export
integrate_levels <- function(results_list,
                            integration_method = "evidence_scoring",
                            weights = NULL,
                            alpha = 0.05,
                            min_levels = 1,
                            annotation = NULL,
                            pathway_db = NULL,
                            ...) {
  
  # Validate inputs
  if (!is.list(results_list) || is.null(names(results_list))) {
    stop("results_list must be a named list")
  }
  
  if (length(results_list) < 2) {
    stop("At least 2 result sets required for integration")
  }
  
  # Set default weights
  if (is.null(weights)) {
    weights <- setNames(rep(1, length(results_list)), names(results_list))
  }
  
  # Perform integration based on method
  if (integration_method == "rank_aggregation") {
    integrated <- integrate_rank_aggregation(results_list, weights, alpha, ...)
    
  } else if (integration_method == "evidence_scoring") {
    integrated <- integrate_evidence_scoring(results_list, weights, alpha, ...)
    
  } else if (integration_method == "consensus") {
    integrated <- integrate_consensus(results_list, alpha, min_levels, ...)
    
  } else if (integration_method == "pathway_guided") {
    integrated <- integrate_pathway_guided(results_list, pathway_db, alpha, ...)
    
  } else {
    stop("Unsupported integration_method. Choose from: rank_aggregation, evidence_scoring, consensus, pathway_guided")
  }
  
  # Add common summary statistics
  integrated$summary <- create_integration_summary(results_list, integrated, alpha)
  integrated$method <- integration_method
  
  class(integrated) <- c("funDE_integration", "list")
  
  return(integrated)
}

#' Rank aggregation integration
#' @keywords internal
integrate_rank_aggregation <- function(results_list, weights, alpha, ...) {
  
  # This is a simplified implementation
  # In a real package, would use RobustRankAggreg or similar
  
  message("Note: This is a simplified rank aggregation implementation.")
  message("For production use, consider using the RobustRankAggreg package.")
  
  # Extract ranks for each level
  all_features <- unique(unlist(lapply(results_list, function(x) x$results$gene)))
  
  rank_matrix <- matrix(NA, nrow = length(all_features), ncol = length(results_list))
  rownames(rank_matrix) <- all_features
  colnames(rank_matrix) <- names(results_list)
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    level_ranks <- rank(results$padj, na.last = "keep")
    names(level_ranks) <- results$gene
    
    # Assign ranks
    common_features <- intersect(all_features, names(level_ranks))
    rank_matrix[common_features, level_name] <- level_ranks[common_features]
  }
  
  # Calculate weighted average ranks
  weighted_ranks <- numeric(length(all_features))
  names(weighted_ranks) <- all_features
  
  for (feature in all_features) {
    feature_ranks <- rank_matrix[feature, ]
    valid_ranks <- !is.na(feature_ranks)
    
    if (sum(valid_ranks) > 0) {
      level_weights <- weights[names(feature_ranks)[valid_ranks]]
      weighted_ranks[feature] <- weighted.mean(feature_ranks[valid_ranks], level_weights)
    } else {
      weighted_ranks[feature] <- Inf
    }
  }
  
  # Create integrated results
  integrated_results <- data.frame(
    feature = names(weighted_ranks),
    integrated_rank = weighted_ranks,
    n_levels = rowSums(!is.na(rank_matrix)),
    stringsAsFactors = FALSE
  ) %>%
    arrange(integrated_rank)
  
  # Add level-specific information
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    level_data <- results[, c("gene", "padj", "log2FoldChange")]
    colnames(level_data) <- c("feature", paste0(level_name, "_padj"), paste0(level_name, "_lfc"))
    
    integrated_results <- merge(integrated_results, level_data, by = "feature", all.x = TRUE)
  }
  
  return(list(
    integrated_results = integrated_results,
    rank_matrix = rank_matrix,
    weights_used = weights
  ))
}

#' Evidence scoring integration
#' @keywords internal
integrate_evidence_scoring <- function(results_list, weights, alpha, ...) {
  
  # Collect all features
  all_features <- unique(unlist(lapply(results_list, function(x) x$results$gene)))
  
  # Initialize evidence scores
  evidence_scores <- data.frame(
    feature = all_features,
    evidence_score = 0,
    n_significant = 0,
    mean_lfc = 0,
    min_padj = 1,
    stringsAsFactors = FALSE
  )
  
  # Calculate evidence scores
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    level_weight <- weights[level_name]
    
    for (i in seq_along(all_features)) {
      feature <- all_features[i]
      
      if (feature %in% results$gene) {
        feature_data <- results[results$gene == feature, ]
        
        # Significance component
        sig_score <- ifelse(feature_data$padj < alpha, level_weight, 0)
        
        # Effect size component
        effect_score <- abs(feature_data$log2FoldChange) * level_weight
        
        # Confidence component (inverse of p-value)
        conf_score <- -log10(feature_data$padj + 1e-300) * level_weight
        
        # Combined evidence score
        combined_score <- sig_score + 0.5 * effect_score + 0.3 * conf_score
        
        evidence_scores$evidence_score[i] <- evidence_scores$evidence_score[i] + combined_score
        
        if (feature_data$padj < alpha) {
          evidence_scores$n_significant[i] <- evidence_scores$n_significant[i] + 1
        }
        
        evidence_scores$mean_lfc[i] <- evidence_scores$mean_lfc[i] + feature_data$log2FoldChange
        evidence_scores$min_padj[i] <- min(evidence_scores$min_padj[i], feature_data$padj)
      }
    }
  }
  
  # Normalize mean LFC by number of levels
  n_levels_per_feature <- sapply(all_features, function(f) {
    sum(sapply(results_list, function(res) f %in% res$results$gene))
  })
  
  evidence_scores$mean_lfc <- evidence_scores$mean_lfc / n_levels_per_feature
  evidence_scores$n_levels <- n_levels_per_feature
  
  # Sort by evidence score
  evidence_scores <- evidence_scores[order(evidence_scores$evidence_score, decreasing = TRUE), ]
  
  return(list(
    integrated_results = evidence_scores,
    weights_used = weights
  ))
}

#' Consensus integration
#' @keywords internal
integrate_consensus <- function(results_list, alpha, min_levels, ...) {
  
  # Find significant features in each level
  sig_features_by_level <- list()
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    sig_features <- results$gene[results$padj < alpha]
    sig_features_by_level[[level_name]] <- sig_features
  }
  
  # Find consensus features
  all_features <- unique(unlist(sig_features_by_level))
  
  consensus_data <- data.frame(
    feature = all_features,
    n_levels_significant = 0,
    levels_significant = "",
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(all_features)) {
    feature <- all_features[i]
    sig_levels <- names(sig_features_by_level)[sapply(sig_features_by_level, function(x) feature %in% x)]
    
    consensus_data$n_levels_significant[i] <- length(sig_levels)
    consensus_data$levels_significant[i] <- paste(sig_levels, collapse = ",")
  }
  
  # Filter by minimum levels
  consensus_features <- consensus_data[consensus_data$n_levels_significant >= min_levels, ]
  consensus_features <- consensus_features[order(consensus_features$n_levels_significant, decreasing = TRUE), ]
  
  return(list(
    integrated_results = consensus_features,
    sig_features_by_level = sig_features_by_level,
    min_levels_used = min_levels
  ))
}

#' Pathway-guided integration
#' @keywords internal
integrate_pathway_guided <- function(results_list, pathway_db, alpha, ...) {
  
  if (is.null(pathway_db)) {
    stop("pathway_db required for pathway-guided integration")
  }
  
  # This is a simplified implementation
  message("Note: This is a simplified pathway-guided integration.")
  
  # Get gene-level results
  gene_results <- NULL
  if ("genes" %in% names(results_list)) {
    gene_results <- results_list$genes$results
  }
  
  # Get pathway-level results
  pathway_results <- NULL
  if ("pathways" %in% names(results_list)) {
    pathway_results <- results_list$pathways$results
  }
  
  if (is.null(gene_results) || is.null(pathway_results)) {
    stop("Both gene and pathway results required for pathway-guided integration")
  }
  
  # Create pathway-gene evidence matrix
  pathway_gene_evidence <- list()
  
  for (pathway_name in names(pathway_db)) {
    pathway_genes <- pathway_db[[pathway_name]]
    
    # Get evidence from gene level
    gene_evidence <- gene_results %>%
      filter(gene %in% pathway_genes, padj < alpha) %>%
      summarise(
        n_sig_genes = n(),
        mean_lfc = mean(log2FoldChange),
        min_padj = min(padj),
        .groups = "drop"
      )
    
    # Get evidence from pathway level
    pathway_evidence <- pathway_results %>%
      filter(pathway == pathway_name)
    
    if (nrow(pathway_evidence) > 0) {
      combined_evidence <- list(
        pathway = pathway_name,
        gene_level_support = gene_evidence,
        pathway_level_result = pathway_evidence,
        consistency = calculate_consistency(gene_evidence, pathway_evidence)
      )
      
      pathway_gene_evidence[[pathway_name]] <- combined_evidence
    }
  }
  
  return(list(
    integrated_results = pathway_gene_evidence,
    pathway_db_used = length(pathway_db)
  ))
}

#' Calculate consistency between gene and pathway evidence
#' @keywords internal
calculate_consistency <- function(gene_evidence, pathway_evidence) {
  
  if (nrow(gene_evidence) == 0 || nrow(pathway_evidence) == 0) {
    return(NA)
  }
  
  # Simple consistency metric based on direction
  gene_direction <- sign(gene_evidence$mean_lfc)
  pathway_direction <- sign(pathway_evidence$log2FoldChange)
  
  consistency <- ifelse(gene_direction == pathway_direction, "consistent", "inconsistent")
  
  return(consistency)
}

#' Create integration summary
#' @keywords internal
create_integration_summary <- function(results_list, integrated, alpha) {
  
  summary_stats <- list(
    n_levels = length(results_list),
    level_names = names(results_list),
    integration_method = integrated$method %||% "unknown"
  )
  
  # Add method-specific summaries
  if ("integrated_results" %in% names(integrated)) {
    summary_stats$n_integrated_features <- nrow(integrated$integrated_results)
  }
  
  # Level-specific summaries
  level_summaries <- list()
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    level_summaries[[level_name]] <- list(
      n_features = nrow(results),
      n_significant = sum(results$padj < alpha, na.rm = TRUE)
    )
  }
  
  summary_stats$level_summaries <- level_summaries
  
  return(summary_stats)
}

#' Print method for funDE_integration
#' @param x A funDE_integration object
#' @param ... Additional arguments (ignored)
#' @export
print.funDE_integration <- function(x, ...) {
  cat("funDE Multi-level Integration\n")
  cat("Method:", x$method, "\n")
  cat("Levels integrated:", paste(x$summary$level_names, collapse = ", "), "\n")
  
  if ("integrated_results" %in% names(x)) {
    cat("Integrated features:", nrow(x$integrated_results), "\n")
  }
}

#' Summary method for funDE_integration
#' @param object A funDE_integration object
#' @param ... Additional arguments (ignored)
#' @export
summary.funDE_integration <- function(object, ...) {
  print(object)
  
  cat("\nLevel-specific summary:\n")
  for (level_name in names(object$summary$level_summaries)) {
    level_sum <- object$summary$level_summaries[[level_name]]
    cat(sprintf("  %s: %d features, %d significant\n", 
                level_name, level_sum$n_features, level_sum$n_significant))
  }
  
  if ("integrated_results" %in% names(object)) {
    cat("\nTop integrated features:\n")
    if (object$method == "evidence_scoring") {
      top_features <- head(object$integrated_results[, c("feature", "evidence_score", "n_significant")], 10)
    } else if (object$method == "consensus") {
      top_features <- head(object$integrated_results[, c("feature", "n_levels_significant")], 10)
    } else {
      top_features <- head(object$integrated_results, 10)
    }
    print(top_features)
  }
}