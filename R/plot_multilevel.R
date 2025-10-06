#' Multi-level Visualization and Comparison
#'
#' @description
#' Creates visualizations comparing results across different functional levels
#' (genes, families, pathways, etc.) to show concordance and complementarity
#' of different analytical approaches.
#'
#' @param results_list Named list of funDE results from different levels
#' @param plot_type Type of comparison plot: 
#'   - "concordance": Correlation between fold changes
#'   - "upset": Overlap of significant features using UpSetR
#'   - "volcano_grid": Grid of volcano plots
#'   - "summary_bar": Bar plot of significant features per level
#'   - "rank_comparison": Rank correlation between levels
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_threshold Log fold change threshold (default: 0)
#' @param top_n Number of top features to include in some plots (default: 100)
#' @param color_palette Color palette for plots (default: NULL for default colors)
#' @param ... Additional arguments passed to specific plot functions
#'
#' @return A ggplot2 object or list of plots
#'
#' @details
#' **Plot Types:**
#'
#' *Concordance Plot:*
#' - Scatter plot comparing log fold changes between levels
#' - Shows correlation between different analytical approaches
#' - Points colored by significance status
#' - Useful for assessing consistency across levels
#'
#' *UpSet Plot:*
#' - Shows overlap of significant features across levels
#' - More informative than Venn diagrams for >3 sets
#' - Reveals unique and shared significant features
#' - Requires UpSetR package
#'
#' *Volcano Grid:*
#' - Grid of volcano plots for each level
#' - Easy comparison of effect sizes and significance
#' - Highlights level-specific patterns
#'
#' *Summary Bar Plot:*
#' - Bar plot showing number of significant features per level
#' - Split by up and down regulation
#' - Quick overview of detection power per level
#'
#' *Rank Comparison:*
#' - Compares ranking of features across levels
#' - Shows how top features from one level perform in others
#' - Useful for understanding complementarity
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
#' # Concordance plot
#' plot_multilevel(results_list, plot_type = "concordance")
#'
#' # UpSet plot of significant features
#' plot_multilevel(results_list, plot_type = "upset")
#'
#' # Summary comparison
#' plot_multilevel(results_list, plot_type = "summary_bar")
#'
#' @seealso \code{\link{analyze_genes}}, \code{\link{analyze_families}}, \code{\link{analyze_pathways}}
#' @export
plot_multilevel <- function(results_list,
                           plot_type = "concordance",
                           alpha = 0.05,
                           lfc_threshold = 0,
                           top_n = 100,
                           color_palette = NULL,
                           ...) {
  
  # Validate inputs
  if (!is.list(results_list) || is.null(names(results_list))) {
    stop("results_list must be a named list")
  }
  
  if (length(results_list) < 2) {
    stop("At least 2 result sets required for comparison")
  }
  
  # Check that all results have expected structure
  for (i in seq_along(results_list)) {
    if (!"results" %in% names(results_list[[i]])) {
      stop("Each result must have a 'results' component")
    }
  }
  
  # Create plots based on type
  if (plot_type == "concordance") {
    plot <- plot_concordance(results_list, alpha, lfc_threshold, color_palette, ...)
    
  } else if (plot_type == "upset") {
    plot <- plot_upset_multilevel(results_list, alpha, lfc_threshold, ...)
    
  } else if (plot_type == "volcano_grid") {
    plot <- plot_volcano_grid(results_list, alpha, lfc_threshold, ...)
    
  } else if (plot_type == "summary_bar") {
    plot <- plot_summary_bar(results_list, alpha, lfc_threshold, color_palette, ...)
    
  } else if (plot_type == "rank_comparison") {
    plot <- plot_rank_comparison(results_list, top_n, ...)
    
  } else {
    stop("Unsupported plot_type. Choose from: concordance, upset, volcano_grid, summary_bar, rank_comparison")
  }
  
  return(plot)
}

#' Plot concordance between levels
#' @keywords internal
plot_concordance <- function(results_list, alpha, lfc_threshold, color_palette, ...) {
  
  if (length(results_list) != 2) {
    stop("Concordance plot requires exactly 2 result sets")
  }
  
  level_names <- names(results_list)
  level1_name <- level_names[1]
  level2_name <- level_names[2]
  
  # Extract results
  res1 <- results_list[[1]]$results
  res2 <- results_list[[2]]$results
  
  # Find common features (this is simplified - in practice would need smarter matching)
  common_features <- intersect(res1$gene, res2$gene)
  
  if (length(common_features) == 0) {
    stop("No common features found between result sets")
  }
  
  # Create combined data
  combined_data <- data.frame(
    feature = common_features,
    lfc1 = res1$log2FoldChange[match(common_features, res1$gene)],
    lfc2 = res2$log2FoldChange[match(common_features, res2$gene)],
    padj1 = res1$padj[match(common_features, res1$gene)],
    padj2 = res2$padj[match(common_features, res2$gene)],
    stringsAsFactors = FALSE
  )
  
  # Create significance categories
  combined_data$significance <- ifelse(
    combined_data$padj1 < alpha & combined_data$padj2 < alpha,
    "Both significant",
    ifelse(
      combined_data$padj1 < alpha,
      paste("Only", level1_name),
      ifelse(
        combined_data$padj2 < alpha,
        paste("Only", level2_name),
        "Neither significant"
      )
    )
  )
  
  # Calculate correlation
  cor_val <- cor(combined_data$lfc1, combined_data$lfc2, use = "complete.obs")
  
  # Create plot
  p <- ggplot(combined_data, aes(x = lfc1, y = lfc2, color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(
      x = paste("Log2 Fold Change (", level1_name, ")", sep = ""),
      y = paste("Log2 Fold Change (", level2_name, ")", sep = ""),
      title = paste("Concordance between", level1_name, "and", level2_name),
      subtitle = paste("Correlation:", round(cor_val, 3)),
      color = "Significance"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Apply color palette if provided
  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  }
  
  return(p)
}

#' Plot UpSet diagram for multilevel overlap
#' @keywords internal
plot_upset_multilevel <- function(results_list, alpha, lfc_threshold, ...) {
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("UpSetR package required for upset plots. Install with: install.packages('UpSetR')")
  }
  
  # Get significant features for each level
  sig_lists <- list()
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    
    sig_features <- results %>%
      filter(padj < alpha, abs(log2FoldChange) > lfc_threshold) %>%
      pull(gene)
    
    sig_lists[[level_name]] <- sig_features
  }
  
  # Create UpSet plot
  upset_plot <- UpSetR::upset(
    UpSetR::fromList(sig_lists),
    order.by = "freq",
    nsets = length(sig_lists),
    ...
  )
  
  return(upset_plot)
}

#' Plot grid of volcano plots
#' @keywords internal
plot_volcano_grid <- function(results_list, alpha, lfc_threshold, ...) {
  
  # Create individual volcano plots
  volcano_plots <- list()
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    
    # Create significance categories
    results$significance <- ifelse(
      results$padj < alpha & abs(results$log2FoldChange) > lfc_threshold,
      ifelse(results$log2FoldChange > 0, "Upregulated", "Downregulated"),
      "Not significant"
    )
    
    # Create volcano plot
    p <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_manual(values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "gray"
      )) +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(alpha), linetype = "dashed", alpha = 0.5) +
      labs(
        x = "Log2 Fold Change",
        y = "-log10(adjusted p-value)",
        title = level_name,
        color = "Significance"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    
    volcano_plots[[level_name]] <- p
  }
  
  # Arrange in grid
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(volcano_plots, ncol = ceiling(sqrt(length(volcano_plots))))
  } else {
    warning("patchwork package recommended for grid layout")
    combined_plot <- volcano_plots
  }
  
  return(combined_plot)
}

#' Plot summary bar chart
#' @keywords internal
plot_summary_bar <- function(results_list, alpha, lfc_threshold, color_palette, ...) {
  
  # Calculate summary statistics
  summary_data <- data.frame(
    Level = character(),
    Direction = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    
    n_up <- sum(results$padj < alpha & results$log2FoldChange > lfc_threshold, na.rm = TRUE)
    n_down <- sum(results$padj < alpha & results$log2FoldChange < -lfc_threshold, na.rm = TRUE)
    
    summary_data <- rbind(summary_data, data.frame(
      Level = rep(level_name, 2),
      Direction = c("Upregulated", "Downregulated"),
      Count = c(n_up, n_down),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create bar plot
  p <- ggplot(summary_data, aes(x = Level, y = Count, fill = Direction)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.2, size = 3) +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    labs(
      x = "Analysis Level",
      y = "Number of Significant Features",
      title = "Significant Features by Analysis Level",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Plot rank comparison
#' @keywords internal
plot_rank_comparison <- function(results_list, top_n, ...) {
  
  if (length(results_list) != 2) {
    stop("Rank comparison requires exactly 2 result sets")
  }
  
  level_names <- names(results_list)
  
  # Get top features from first level
  top_features <- results_list[[1]]$results %>%
    arrange(padj) %>%
    head(top_n) %>%
    pull(gene)
  
  # Get ranks in both levels
  rank_data <- data.frame(
    feature = top_features,
    rank1 = match(top_features, results_list[[1]]$results$gene[order(results_list[[1]]$results$padj)]),
    rank2 = match(top_features, results_list[[2]]$results$gene[order(results_list[[2]]$results$padj)]),
    stringsAsFactors = FALSE
  )
  
  # Handle features not found in second level
  rank_data$rank2[is.na(rank_data$rank2)] <- nrow(results_list[[2]]$results) + 1
  
  # Create plot
  p <- ggplot(rank_data, aes(x = rank1, y = rank2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      x = paste("Rank in", level_names[1], "(log scale)"),
      y = paste("Rank in", level_names[2], "(log scale)"),
      title = paste("Rank Comparison:", level_names[1], "vs", level_names[2]),
      subtitle = paste("Top", top_n, "features from", level_names[1])
    ) +
    theme_minimal(base_size = 12)
  
  return(p)
}

#' Create Multi-level Summary Table
#'
#' @description
#' Creates a summary table comparing statistics across different analysis levels
#'
#' @param results_list Named list of funDE results
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_threshold Log fold change threshold (default: 0)
#'
#' @return Data.frame with summary statistics
#'
#' @examples
#' results_list <- list(genes = gene_results, pathways = pathway_results)
#' summary_table <- create_multilevel_summary(results_list)
#'
#' @export
create_multilevel_summary <- function(results_list, alpha = 0.05, lfc_threshold = 0) {
  
  summary_data <- data.frame(
    Level = character(),
    N_Features = integer(),
    N_Significant = integer(),
    N_Upregulated = integer(),
    N_Downregulated = integer(),
    Median_LFC = numeric(),
    Min_PValue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (level_name in names(results_list)) {
    results <- results_list[[level_name]]$results
    
    n_features <- nrow(results)
    n_sig <- sum(results$padj < alpha, na.rm = TRUE)
    n_up <- sum(results$padj < alpha & results$log2FoldChange > lfc_threshold, na.rm = TRUE)
    n_down <- sum(results$padj < alpha & results$log2FoldChange < -lfc_threshold, na.rm = TRUE)
    median_lfc <- median(abs(results$log2FoldChange[results$padj < alpha]), na.rm = TRUE)
    min_pval <- min(results$padj, na.rm = TRUE)
    
    summary_data <- rbind(summary_data, data.frame(
      Level = level_name,
      N_Features = n_features,
      N_Significant = n_sig,
      N_Upregulated = n_up,
      N_Downregulated = n_down,
      Median_LFC = median_lfc,
      Min_PValue = min_pval,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add percentage columns
  summary_data$Percent_Significant <- round(100 * summary_data$N_Significant / summary_data$N_Features, 1)
  
  return(summary_data)
}