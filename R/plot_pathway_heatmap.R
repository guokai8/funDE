#' Plot Pathway Activity Heatmap
#'
#' @description
#' Creates heatmaps of pathway activity scores with sample annotations,
#' clustering, and customizable features for exploring pathway-level
#' expression patterns.
#'
#' @param pathway_results Output from analyze_pathways() or pathway scores matrix
#' @param metadata Sample metadata for annotations (optional if included in pathway_results)
#' @param top_n Number of top significant pathways to include (default: 20)
#' @param cluster_rows Cluster pathways (default: TRUE)
#' @param cluster_cols Cluster samples (default: TRUE)
#' @param annotation_cols Column names in metadata to use for sample annotations
#' @param color_scheme Color scheme: "RdBu", "viridis", "RdYlBu" (default: "RdBu")
#' @param scale Scale data: "row", "column", "none" (default: "row")
#' @param show_pathway_names Show pathway names (default: TRUE)
#' @param fontsize_row Font size for pathway names (default: 8)
#' @param fontsize_col Font size for sample names (default: 8)
#' @param use_complexheatmap Use ComplexHeatmap instead of pheatmap (default: FALSE)
#' @param ... Additional arguments passed to heatmap function
#'
#' @return A heatmap plot object (grob for pheatmap, HeatmapList for ComplexHeatmap)
#'
#' @details
#' **Input Types:**
#' - pathway_results: Output from analyze_pathways() containing pathway scores
#' - Matrix: Direct pathway scores matrix (pathways x samples)
#'
#' **Clustering:**
#' - Hierarchical clustering reveals pathway modules and sample groups
#' - Distance metrics: correlation (default), euclidean, manhattan
#' - Linkage methods: complete (default), average, single, ward
#'
#' **Scaling:**
#' - Row scaling: Z-score normalization across samples for each pathway
#' - Column scaling: Z-score normalization across pathways for each sample
#' - No scaling: Raw pathway scores
#'
#' **Annotations:**
#' - Sample annotations from metadata
#' - Automatic color assignment for categorical variables
#' - Numeric variables shown as continuous color gradients
#'
#' **Color Schemes:**
#' - RdBu: Red-white-blue, good for centered data
#' - viridis: Perceptually uniform, colorblind-friendly
#' - RdYlBu: Red-yellow-blue, good for diverging data
#'
#' @examples
#' # Basic pathway heatmap
#' pathway_results <- analyze_pathways(...)
#' plot_pathway_heatmap(pathway_results, top_n = 30)
#'
#' # With sample annotations
#' plot_pathway_heatmap(
#'   pathway_results,
#'   annotation_cols = c("condition", "batch"),
#'   top_n = 25,
#'   color_scheme = "viridis"
#' )
#'
#' # Using ComplexHeatmap for advanced features
#' plot_pathway_heatmap(
#'   pathway_results,
#'   use_complexheatmap = TRUE,
#'   annotation_cols = "condition",
#'   cluster_rows = TRUE,
#'   cluster_cols = FALSE
#' )
#'
#' # Direct matrix input
#' pathway_scores <- score_pathway_activity(...)
#' plot_pathway_heatmap(
#'   pathway_scores,
#'   metadata = sample_metadata,
#'   annotation_cols = "treatment"
#' )
#'
#' @seealso \code{\link{analyze_pathways}}, \code{\link{score_pathway_activity}}
#' @export
plot_pathway_heatmap <- function(pathway_results,
                                metadata = NULL,
                                top_n = 20,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                annotation_cols = NULL,
                                color_scheme = "RdBu",
                                scale = "row",
                                show_pathway_names = TRUE,
                                fontsize_row = 8,
                                fontsize_col = 8,
                                use_complexheatmap = FALSE,
                                ...) {
  
  # Extract pathway scores matrix
  if (is.matrix(pathway_results)) {
    pathway_scores <- pathway_results
    if (is.null(metadata)) {
      metadata <- data.frame(sample = colnames(pathway_scores))
      rownames(metadata) <- colnames(pathway_scores)
    }
  } else if (is.list(pathway_results) && "pathway_scores" %in% names(pathway_results)) {
    pathway_scores <- pathway_results$pathway_scores
    if (is.null(metadata)) {
      metadata <- pathway_results$metadata
    }
  } else {
    stop("pathway_results must be a matrix or output from analyze_pathways()")
  }
  
  # Validate inputs
  if (!is.matrix(pathway_scores)) {
    pathway_scores <- as.matrix(pathway_scores)
  }
  
  if (is.null(metadata)) {
    metadata <- data.frame(sample = colnames(pathway_scores))
    rownames(metadata) <- colnames(pathway_scores)
  }
  
  # Select top pathways if pathway_results contains DE results
  if (is.list(pathway_results) && "results" %in% names(pathway_results)) {
    top_pathways <- pathway_results$results %>%
      arrange(padj) %>%
      head(top_n) %>%
      pull(pathway)
    
    # Filter to available pathways
    top_pathways <- intersect(top_pathways, rownames(pathway_scores))
  } else {
    # Select most variable pathways
    pathway_var <- apply(pathway_scores, 1, var, na.rm = TRUE)
    top_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:min(top_n, nrow(pathway_scores))]
  }
  
  # Subset pathway scores
  plot_data <- pathway_scores[top_pathways, , drop = FALSE]
  
  # Ensure metadata order matches
  metadata <- metadata[colnames(plot_data), , drop = FALSE]
  
  # Create heatmap based on backend choice
  if (use_complexheatmap) {
    heatmap_plot <- create_complex_heatmap(
      plot_data, metadata, annotation_cols, color_scheme, scale,
      cluster_rows, cluster_cols, show_pathway_names, fontsize_row, fontsize_col, ...
    )
  } else {
    heatmap_plot <- create_pheatmap(
      plot_data, metadata, annotation_cols, color_scheme, scale,
      cluster_rows, cluster_cols, show_pathway_names, fontsize_row, fontsize_col, ...
    )
  }
  
  return(heatmap_plot)
}

#' Create heatmap using pheatmap
#' @keywords internal
create_pheatmap <- function(plot_data, metadata, annotation_cols, color_scheme, scale,
                           cluster_rows, cluster_cols, show_pathway_names, 
                           fontsize_row, fontsize_col, ...) {
  
  check_package("pheatmap", "pheatmap visualization")
  
  # Prepare annotations
  annotation_df <- NULL
  if (!is.null(annotation_cols)) {
    available_cols <- intersect(annotation_cols, colnames(metadata))
    if (length(available_cols) > 0) {
      annotation_df <- metadata[, available_cols, drop = FALSE]
    }
  }
  
  # Set color palette
  if (color_scheme == "RdBu") {
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
  } else if (color_scheme == "viridis") {
    colors <- viridis::viridis(100)
  } else if (color_scheme == "RdYlBu") {
    colors <- colorRampPalette(c("red", "yellow", "blue"))(100)
  } else {
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
  }
  
  # Create heatmap
  heatmap_plot <- pheatmap::pheatmap(
    plot_data,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_col = annotation_df,
    color = colors,
    show_rownames = show_pathway_names,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    border_color = NA,
    ...
  )
  
  return(heatmap_plot)
}

#' Create heatmap using ComplexHeatmap
#' @keywords internal
create_complex_heatmap <- function(plot_data, metadata, annotation_cols, color_scheme, scale,
                                  cluster_rows, cluster_cols, show_pathway_names,
                                  fontsize_row, fontsize_col, ...) {
  
  check_package("ComplexHeatmap", "ComplexHeatmap visualization")
  
  # Prepare color function
  if (color_scheme == "RdBu") {
    col_fun <- circlize::colorRamp2(
      c(min(plot_data, na.rm = TRUE), 0, max(plot_data, na.rm = TRUE)),
      c("blue", "white", "red")
    )
  } else if (color_scheme == "viridis") {
    col_fun <- circlize::colorRamp2(
      seq(min(plot_data, na.rm = TRUE), max(plot_data, na.rm = TRUE), length = 3),
      viridis::viridis(3)
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(min(plot_data, na.rm = TRUE), 0, max(plot_data, na.rm = TRUE)),
      c("blue", "white", "red")
    )
  }
  
  # Prepare column annotations
  col_annotation <- NULL
  if (!is.null(annotation_cols)) {
    available_cols <- intersect(annotation_cols, colnames(metadata))
    if (length(available_cols) > 0) {
      col_annotation <- ComplexHeatmap::HeatmapAnnotation(
        df = metadata[, available_cols, drop = FALSE]
      )
    }
  }
  
  # Create main heatmap
  ht <- ComplexHeatmap::Heatmap(
    plot_data,
    name = "Pathway\nActivity",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_pathway_names,
    row_names_gp = grid::gpar(fontsize = fontsize_row),
    column_names_gp = grid::gpar(fontsize = fontsize_col),
    top_annotation = col_annotation,
    ...
  )
  
  return(ht)
}

#' Interactive Pathway Heatmap
#'
#' @description
#' Creates an interactive heatmap using plotly for exploring pathway activity
#'
#' @param pathway_results Output from analyze_pathways() or pathway scores matrix
#' @param metadata Sample metadata (optional)
#' @param top_n Number of top pathways to include (default: 20)
#' @param scale Scale data: "row", "column", "none" (default: "row")
#' @param color_scheme Color scheme (default: "RdBu")
#'
#' @return A plotly object
#'
#' @examples
#' pathway_results <- analyze_pathways(...)
#' interactive_heatmap <- plot_pathway_heatmap_interactive(pathway_results)
#'
#' @export
plot_pathway_heatmap_interactive <- function(pathway_results,
                                            metadata = NULL,
                                            top_n = 20,
                                            scale = "row",
                                            color_scheme = "RdBu") {
  
  check_package("plotly", "interactive heatmap")
  
  # Extract data (similar to main function)
  if (is.matrix(pathway_results)) {
    pathway_scores <- pathway_results
  } else {
    pathway_scores <- pathway_results$pathway_scores
    if (is.null(metadata)) {
      metadata <- pathway_results$metadata
    }
  }
  
  # Select top pathways
  if (is.list(pathway_results) && "results" %in% names(pathway_results)) {
    top_pathways <- pathway_results$results %>%
      arrange(padj) %>%
      head(top_n) %>%
      pull(pathway)
  } else {
    pathway_var <- apply(pathway_scores, 1, var, na.rm = TRUE)
    top_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:min(top_n, nrow(pathway_scores))]
  }
  
  plot_data <- pathway_scores[top_pathways, , drop = FALSE]
  
  # Scale data if requested
  if (scale == "row") {
    plot_data <- t(scale(t(plot_data)))
  } else if (scale == "column") {
    plot_data <- scale(plot_data)
  }
  
  # Create interactive heatmap
  p <- plotly::plot_ly(
    z = plot_data,
    x = colnames(plot_data),
    y = rownames(plot_data),
    type = "heatmap",
    colorscale = if (color_scheme == "RdBu") "RdBu" else "Viridis",
    hovertemplate = "Pathway: %{y}<br>Sample: %{x}<br>Score: %{z}<extra></extra>"
  ) %>%
    plotly::layout(
      title = "Pathway Activity Heatmap",
      xaxis = list(title = "Samples"),
      yaxis = list(title = "Pathways")
    )
  
  return(p)
}

#' Plot Pathway Activity Distribution
#'
#' @description
#' Creates density or violin plots showing distribution of pathway activity scores
#'
#' @param pathway_results Output from analyze_pathways()
#' @param pathways Character vector of pathway names to plot (default: top 6)
#' @param group_by Column name in metadata to group samples by
#' @param plot_type Type of plot: "violin", "boxplot", "density" (default: "violin")
#' @param ncol Number of columns for faceting (default: 3)
#'
#' @return A ggplot2 object
#'
#' @examples
#' pathway_results <- analyze_pathways(...)
#' plot_pathway_distribution(pathway_results, group_by = "condition")
#'
#' @export
plot_pathway_distribution <- function(pathway_results,
                                     pathways = NULL,
                                     group_by = NULL,
                                     plot_type = "violin",
                                     ncol = 3) {
  
  pathway_scores <- pathway_results$pathway_scores
  metadata <- pathway_results$metadata
  
  # Select pathways
  if (is.null(pathways)) {
    top_pathways <- pathway_results$results %>%
      arrange(padj) %>%
      head(6) %>%
      pull(pathway)
    pathways <- intersect(top_pathways, rownames(pathway_scores))
  }
  
  # Prepare data for plotting
  plot_data <- pathway_scores[pathways, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("pathway") %>%
    tidyr::pivot_longer(-pathway, names_to = "sample", values_to = "activity")
  
  # Add group information if specified
  if (!is.null(group_by) && group_by %in% colnames(metadata)) {
    sample_groups <- metadata[[group_by]]
    names(sample_groups) <- rownames(metadata)
    plot_data$group <- sample_groups[plot_data$sample]
  }
  
  # Create plot
  if (!is.null(group_by) && "group" %in% colnames(plot_data)) {
    p <- ggplot(plot_data, aes(x = group, y = activity, fill = group))
  } else {
    p <- ggplot(plot_data, aes(x = "", y = activity))
  }
  
  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  } else if (plot_type == "boxplot") {
    p <- p + geom_boxplot(alpha = 0.7)
  } else if (plot_type == "density") {
    p <- ggplot(plot_data, aes(x = activity, fill = group)) +
      geom_density(alpha = 0.7)
  }
  
  p <- p +
    facet_wrap(~pathway, scales = "free", ncol = ncol) +
    labs(
      x = if (!is.null(group_by)) group_by else "",
      y = "Pathway Activity Score",
      title = "Pathway Activity Distribution",
      fill = group_by
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10)
    )
  
  return(p)
}