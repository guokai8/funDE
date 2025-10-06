#' Validate count matrix
#' @param counts Count matrix 
#' @param allow_missing Allow missing values
#' @return Validated count matrix
#' @keywords internal
validate_counts <- function(counts, allow_missing = FALSE) {
  if (is.null(counts)) {
    stop("counts cannot be NULL")
  }
  
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts must be a matrix or data.frame")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  
  # Check for numeric data
  if (!is.numeric(counts)) {
    stop("counts must contain numeric values")
  }
  
  # Check for negative values
  if (any(counts < 0, na.rm = TRUE)) {
    stop("counts cannot contain negative values")
  }
  
  # Check for missing values
  if (!allow_missing && any(is.na(counts))) {
    stop("counts cannot contain missing values")
  }
  
  # Check dimensions
  if (any(dim(counts) == 0)) {
    stop("counts cannot have zero dimensions")
  }
  
  # Ensure row and column names
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("Gene_", seq_len(nrow(counts)))
    warning("Row names not provided, using Gene_1, Gene_2, ...")
  }
  
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("Sample_", seq_len(ncol(counts)))
    warning("Column names not provided, using Sample_1, Sample_2, ...")
  }
  
  return(counts)
}

#' Validate metadata
#' @param metadata Sample metadata data.frame
#' @param counts Count matrix for sample matching
#' @return Validated metadata
#' @keywords internal
validate_metadata <- function(metadata, counts) {
  if (is.null(metadata)) {
    stop("metadata cannot be NULL")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  if (nrow(metadata) == 0) {
    stop("metadata cannot be empty")
  }
  
  # Check sample names match
  if (!is.null(counts)) {
    sample_names <- colnames(counts)
    
    if (is.null(rownames(metadata))) {
      if (nrow(metadata) != length(sample_names)) {
        stop("Number of rows in metadata must match number of samples in counts")
      }
      rownames(metadata) <- sample_names
    } else {
      missing_samples <- setdiff(sample_names, rownames(metadata))
      if (length(missing_samples) > 0) {
        stop("Missing samples in metadata: ", paste(missing_samples, collapse = ", "))
      }
      
      # Reorder metadata to match counts
      metadata <- metadata[sample_names, , drop = FALSE]
    }
  }
  
  return(metadata)
}

#' Validate design formula
#' @param design Design formula
#' @param metadata Sample metadata
#' @return Validated design
#' @keywords internal
validate_design <- function(design, metadata) {
  if (is.null(design)) {
    stop("design cannot be NULL")
  }
  
  if (!inherits(design, "formula")) {
    stop("design must be a formula")
  }
  
  # Check variables exist in metadata
  design_vars <- all.vars(design)
  missing_vars <- setdiff(design_vars, colnames(metadata))
  
  if (length(missing_vars) > 0) {
    stop("Design variables not found in metadata: ", paste(missing_vars, collapse = ", "))
  }
  
  return(design)
}

#' Validate contrast specification
#' @param contrast Contrast specification
#' @param metadata Sample metadata
#' @return Validated contrast
#' @keywords internal
validate_contrast <- function(contrast, metadata) {
  if (is.null(contrast)) {
    return(NULL)
  }
  
  if (!is.character(contrast) || length(contrast) != 3) {
    stop("contrast must be a character vector of length 3: c('variable', 'level1', 'level2')")
  }
  
  var_name <- contrast[1]
  level1 <- contrast[2]
  level2 <- contrast[3]
  
  if (!var_name %in% colnames(metadata)) {
    stop("Contrast variable '", var_name, "' not found in metadata")
  }
  
  available_levels <- unique(metadata[[var_name]])
  missing_levels <- setdiff(c(level1, level2), available_levels)
  
  if (length(missing_levels) > 0) {
    stop("Contrast levels not found in metadata: ", paste(missing_levels, collapse = ", "),
         "\nAvailable levels: ", paste(available_levels, collapse = ", "))
  }
  
  return(contrast)
}

#' Validate gene sets
#' @param gene_sets Named list of gene sets
#' @param min_size Minimum set size
#' @param max_size Maximum set size
#' @return Validated gene sets
#' @keywords internal
validate_gene_sets <- function(gene_sets, min_size = 1, max_size = Inf) {
  if (!is.list(gene_sets)) {
    stop("gene_sets must be a list")
  }
  
  if (is.null(names(gene_sets))) {
    stop("gene_sets must be a named list")
  }
  
  # Filter by size
  set_sizes <- sapply(gene_sets, length)
  valid_sets <- set_sizes >= min_size & set_sizes <= max_size
  
  if (sum(valid_sets) == 0) {
    stop("No gene sets meet size criteria (min: ", min_size, ", max: ", max_size, ")")
  }
  
  gene_sets <- gene_sets[valid_sets]
  
  return(gene_sets)
}

#' Check if package is available
#' @param package Package name
#' @param function_name Function name for error message
#' @return TRUE if available, stops if not
#' @keywords internal
check_package <- function(package, function_name = "This function") {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(function_name, " requires the '", package, "' package. ",
         "Install it with: install.packages('", package, "')")
  }
  return(TRUE)
}