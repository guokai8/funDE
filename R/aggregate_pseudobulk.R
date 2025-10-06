#' Aggregate Single-cell Data to Pseudobulk
#'
#' @description
#' Aggregates single-cell RNA-seq data to pseudobulk samples for differential
#' expression analysis. This function sums counts across cells within each 
#' sample and cell type combination, creating pseudobulk profiles suitable
#' for standard DE analysis tools.
#'
#' @param counts A sparse matrix of counts (genes x cells) or SingleCellExperiment object
#' @param metadata A data.frame containing cell metadata with rownames matching colnames of counts
#' @param sample_col Column name in metadata containing sample IDs
#' @param celltype_col Column name in metadata containing cell type annotations (optional)
#' @param min_cells Minimum cells required per pseudobulk sample (default: 10)
#' @param min_features Minimum features (genes) required per pseudobulk sample (default: 1000)
#' @param aggregate_by How to group cells: "sample" (default), "sample_celltype", or custom grouping
#' @param assay_name Name of assay to use if input is SingleCellExperiment (default: "counts")
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{counts}{Pseudobulk count matrix (genes x pseudobulk_samples)}
#'   \item{metadata}{Metadata for pseudobulk samples}
#'   \item{cell_counts}{Number of cells per pseudobulk sample}
#'   \item{original_cells}{Number of cells in original data}
#'   \item{aggregation_method}{Method used for aggregation}
#'
#' @details
#' **Aggregation Methods:**
#'
#' *By Sample:*
#' - Aggregate all cells from each sample
#' - Good for detecting sample-level differences
#' - Ignores cell type heterogeneity
#'
#' *By Sample and Cell Type:*
#' - Aggregate cells within each sample-celltype combination
#' - Enables cell type-specific DE analysis
#' - Requires cell type annotations
#' - More granular but requires more cells
#'
#' *Custom Grouping:*
#' - User defines grouping variable
#' - Maximum flexibility
#' - Column name in metadata specifying groups
#'
#' **Quality Filters:**
#' - Pseudobulk samples with < min_cells are removed
#' - Genes detected in < min_features are filtered
#' - Helps ensure reliable DE testing
#'
#' **Biological Considerations:**
#' - Pseudobulk approximates bulk RNA-seq from the same samples
#' - Loses single-cell resolution but gains statistical power
#' - Suitable for detecting overall expression changes
#' - May miss cell type-specific effects if aggregating across types
#'
#' @examples
#' # Basic aggregation by sample
#' pseudobulk <- aggregate_pseudobulk(
#'   counts = sc_counts,
#'   metadata = sc_metadata,
#'   sample_col = "sample_id",
#'   aggregate_by = "sample"
#' )
#'
#' # Aggregate by sample and cell type
#' pseudobulk_ct <- aggregate_pseudobulk(
#'   counts = sc_counts,
#'   metadata = sc_metadata,
#'   sample_col = "sample_id",
#'   celltype_col = "cell_type",
#'   aggregate_by = "sample_celltype",
#'   min_cells = 20
#' )
#'
#' # Use aggregated data for DE analysis
#' de_results <- analyze_genes(
#'   counts = pseudobulk$counts,
#'   metadata = pseudobulk$metadata,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control")
#' )
#'
#' # With SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   pseudobulk_sce <- aggregate_pseudobulk(
#'     counts = sce_object,
#'     sample_col = "sample_id",
#'     celltype_col = "cell_type",
#'     aggregate_by = "sample_celltype"
#'   )
#' }
#'
#' @seealso \code{\link{analyze_genes}}, \code{\link{analyze_pathways}}
#' @export
aggregate_pseudobulk <- function(counts,
                                 metadata = NULL,
                                 sample_col,
                                 celltype_col = NULL,
                                 min_cells = 10,
                                 min_features = 1000,
                                 aggregate_by = "sample",
                                 assay_name = "counts",
                                 verbose = TRUE) {
  
  # Handle SingleCellExperiment input
  if (is(counts, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package required for SCE input")
    }
    
    if (is.null(metadata)) {
      metadata <- as.data.frame(SummarizedExperiment::colData(counts))
    }
    counts <- SummarizedExperiment::assay(counts, assay_name)
  }
  
  # Handle Seurat input
  if (is(counts, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required for Seurat input")
    }
    
    if (is.null(metadata)) {
      metadata <- counts@meta.data
    }
    counts <- Seurat::GetAssayData(counts, slot = "counts")
  }
  
  # Validate inputs
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("counts must be a matrix, sparse matrix, or single-cell object")
  }
  
  if (is.null(metadata)) {
    stop("metadata must be provided or extractable from single-cell object")
  }
  
  metadata <- validate_metadata(metadata, counts)
  
  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  
  if (aggregate_by == "sample_celltype" && is.null(celltype_col)) {
    stop("celltype_col must be specified when aggregate_by = 'sample_celltype'")
  }
  
  if (!is.null(celltype_col) && !celltype_col %in% colnames(metadata)) {
    stop("celltype_col '", celltype_col, "' not found in metadata")
  }
  
  original_ncells <- ncol(counts)
  
  if (verbose) {
    message("Aggregating ", original_ncells, " cells to pseudobulk using method: ", aggregate_by)
  }
  
  # Create aggregation groups
  if (aggregate_by == "sample") {
    metadata$group <- metadata[[sample_col]]
    
  } else if (aggregate_by == "sample_celltype") {
    metadata$group <- paste(metadata[[sample_col]], metadata[[celltype_col]], sep = "_")
    
  } else {
    # Custom grouping
    if (!aggregate_by %in% colnames(metadata)) {
      stop("Custom aggregate_by column '", aggregate_by, "' not found in metadata")
    }
    metadata$group <- metadata[[aggregate_by]]
  }
  
  # Count cells per group
  cell_counts <- table(metadata$group)
  
  # Filter groups with insufficient cells
  valid_groups <- names(cell_counts)[cell_counts >= min_cells]
  
  if (length(valid_groups) == 0) {
    stop("No groups have >= ", min_cells, " cells")
  }
  
  if (verbose) {
    n_filtered <- sum(cell_counts < min_cells)
    if (n_filtered > 0) {
      message("Filtered ", n_filtered, " groups with < ", min_cells, " cells")
    }
    message("Proceeding with ", length(valid_groups), " pseudobulk samples")
  }
  
  # Filter cells to valid groups
  valid_cells <- metadata$group %in% valid_groups
  counts_filtered <- counts[, valid_cells, drop = FALSE]
  metadata_filtered <- metadata[valid_cells, ]
  
  # Aggregate counts by group
  groups <- metadata_filtered$group
  unique_groups <- unique(groups)
  
  # Create pseudobulk count matrix
  pseudobulk_counts <- matrix(0, 
                             nrow = nrow(counts_filtered), 
                             ncol = length(unique_groups))
  rownames(pseudobulk_counts) <- rownames(counts_filtered)
  colnames(pseudobulk_counts) <- unique_groups
  
  # Sum counts for each group
  for (group in unique_groups) {
    group_cells <- groups == group
    
    if (sum(group_cells) == 1) {
      pseudobulk_counts[, group] <- counts_filtered[, group_cells]
    } else {
      pseudobulk_counts[, group] <- Matrix::rowSums(counts_filtered[, group_cells, drop = FALSE])
    }
  }
  
  # Convert to regular matrix if sparse
  if (inherits(pseudobulk_counts, "Matrix")) {
    pseudobulk_counts <- as.matrix(pseudobulk_counts)
  }
  
  # Filter lowly detected genes
  genes_detected <- rowSums(pseudobulk_counts > 0)
  keep_genes <- genes_detected >= min_features / 1000 * ncol(pseudobulk_counts)
  
  if (sum(keep_genes) == 0) {
    warning("No genes pass the min_features filter")
    keep_genes <- genes_detected > 0
  }
  
  pseudobulk_counts <- pseudobulk_counts[keep_genes, , drop = FALSE]
  
  if (verbose) {
    message("Kept ", nrow(pseudobulk_counts), " genes after filtering")
  }
  
  # Create pseudobulk metadata
  pseudobulk_metadata <- metadata_filtered %>%
    group_by(group) %>%
    summarise(
      sample = first(!!sym(sample_col)),
      n_cells = n(),
      .groups = "drop"
    )
  
  # Add cell type information if available
  if (!is.null(celltype_col)) {
    celltype_info <- metadata_filtered %>%
      group_by(group) %>%
      summarise(
        cell_type = first(!!sym(celltype_col)),
        .groups = "drop"
      )
    
    pseudobulk_metadata <- merge(pseudobulk_metadata, celltype_info, by = "group")
  }
  
  # Add any other metadata columns (taking first value per group)
  other_cols <- setdiff(colnames(metadata_filtered), 
                       c("group", sample_col, celltype_col))
  
  if (length(other_cols) > 0) {
    other_metadata <- metadata_filtered %>%
      group_by(group) %>%
      summarise(across(all_of(other_cols), first), .groups = "drop")
    
    pseudobulk_metadata <- merge(pseudobulk_metadata, other_metadata, by = "group")
  }
  
  # Set rownames for metadata
  rownames(pseudobulk_metadata) <- pseudobulk_metadata$group
  pseudobulk_metadata <- pseudobulk_metadata[colnames(pseudobulk_counts), ]
  
  # Create result object
  result <- list(
    counts = pseudobulk_counts,
    metadata = pseudobulk_metadata,
    cell_counts = cell_counts[valid_groups],
    original_cells = original_ncells,
    aggregation_method = aggregate_by,
    parameters = list(
      sample_col = sample_col,
      celltype_col = celltype_col,
      min_cells = min_cells,
      min_features = min_features,
      n_pseudobulk_samples = ncol(pseudobulk_counts),
      n_genes_kept = nrow(pseudobulk_counts)
    )
  )
  
  class(result) <- c("funDE_pseudobulk", "list")
  
  if (verbose) {
    message("Pseudobulk aggregation complete:")
    message("  ", original_ncells, " cells -> ", ncol(pseudobulk_counts), " pseudobulk samples")
    message("  ", nrow(pseudobulk_counts), " genes retained")
  }
  
  return(result)
}

#' Print method for funDE_pseudobulk
#' @param x A funDE_pseudobulk object
#' @param ... Additional arguments (ignored)
#' @export
print.funDE_pseudobulk <- function(x, ...) {
  cat("funDE Pseudobulk Aggregation\n")
  cat("Method:", x$aggregation_method, "\n")
  cat("Original cells:", x$original_cells, "\n")
  cat("Pseudobulk samples:", ncol(x$counts), "\n")
  cat("Genes:", nrow(x$counts), "\n")
  cat("Cell count range:", min(x$cell_counts), "-", max(x$cell_counts), "\n")
}

#' Summary method for funDE_pseudobulk
#' @param object A funDE_pseudobulk object
#' @param ... Additional arguments (ignored)
#' @export
summary.funDE_pseudobulk <- function(object, ...) {
  print(object)
  
  cat("\nPseudobulk sample details:\n")
  print(object$metadata[, c("sample", "n_cells")])
  
  cat("\nCell count distribution:\n")
  print(summary(object$cell_counts))
}

#' Quality Control for Pseudobulk Data
#'
#' @description
#' Performs quality control checks on pseudobulk data to identify potential issues
#'
#' @param pseudobulk_result Output from aggregate_pseudobulk()
#' @param plot Logical, create QC plots (default: TRUE)
#'
#' @return List with QC metrics and plots
#'
#' @examples
#' pseudobulk <- aggregate_pseudobulk(...)
#' qc_results <- qc_pseudobulk(pseudobulk)
#'
#' @export
qc_pseudobulk <- function(pseudobulk_result, plot = TRUE) {
  
  counts <- pseudobulk_result$counts
  metadata <- pseudobulk_result$metadata
  
  # Calculate QC metrics
  qc_metrics <- data.frame(
    sample = rownames(metadata),
    n_cells = metadata$n_cells,
    total_counts = colSums(counts),
    n_detected_genes = colSums(counts > 0),
    median_counts = apply(counts, 2, median),
    stringsAsFactors = FALSE
  )
  
  # Add sample information
  qc_metrics <- merge(qc_metrics, metadata, by.x = "sample", by.y = "row.names")
  
  # Identify potential outliers
  qc_metrics$total_counts_outlier <- abs(scale(qc_metrics$total_counts)) > 2
  qc_metrics$n_cells_outlier <- abs(scale(qc_metrics$n_cells)) > 2
  
  # Summary statistics
  qc_summary <- list(
    n_samples = nrow(qc_metrics),
    total_cells = sum(qc_metrics$n_cells),
    median_cells_per_sample = median(qc_metrics$n_cells),
    total_count_range = range(qc_metrics$total_counts),
    n_outlier_samples = sum(qc_metrics$total_counts_outlier | qc_metrics$n_cells_outlier)
  )
  
  result <- list(
    metrics = qc_metrics,
    summary = qc_summary
  )
  
  return(result)
}