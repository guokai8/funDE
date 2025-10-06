#' Convert gene IDs between different formats
#' @param genes Character vector of gene identifiers
#' @param from Source identifier type: "symbol", "ensembl", "entrez"
#' @param to Target identifier type: "symbol", "ensembl", "entrez"
#' @param species Character, "human" or "mouse"
#' @param remove_unmapped Logical, remove genes that couldn't be mapped
#' @return Named character vector with converted IDs
#' @export
convert_gene_ids <- function(genes, from = "symbol", to = "ensembl", 
                            species = "human", remove_unmapped = TRUE) {
  
  # This is a simplified implementation
  # In a real package, this would use biomaRt or org.Hs.eg.db
  
  warning("convert_gene_ids is a simplified implementation. ",
          "For production use, implement proper ID conversion using biomaRt or annotation packages.")
  
  # Return input with warning for now
  converted <- setNames(genes, genes)
  
  if (remove_unmapped) {
    # Remove any NAs or empty strings
    converted <- converted[!is.na(converted) & converted != ""]
  }
  
  return(converted)
}

#' Safely extract elements from list with defaults
#' @param x List
#' @param name Element name
#' @param default Default value if element doesn't exist
#' @return Element value or default
#' @keywords internal
safe_extract <- function(x, name, default = NULL) {
  if (name %in% names(x)) {
    return(x[[name]])
  } else {
    return(default)
  }
}

#' Check if all required columns exist in data.frame
#' @param df Data.frame to check
#' @param required_cols Character vector of required column names
#' @param df_name Name of data.frame for error messages
#' @return TRUE if all columns exist, stops otherwise
#' @keywords internal
check_required_columns <- function(df, required_cols, df_name = "data.frame") {
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", df_name, ": ", 
         paste(missing_cols, collapse = ", "))
  }
  
  return(TRUE)
}

#' Format p-values for display
#' @param pvals Numeric vector of p-values
#' @param threshold Threshold below which to show as "< threshold"
#' @param digits Number of significant digits
#' @return Character vector of formatted p-values
#' @keywords internal
format_pvalues <- function(pvals, threshold = 0.001, digits = 3) {
  formatted <- ifelse(
    pvals < threshold,
    paste0("< ", threshold),
    format(pvals, digits = digits, scientific = TRUE)
  )
  
  # Handle NAs
  formatted[is.na(pvals)] <- "NA"
  
  return(formatted)
}

#' Format log fold changes for display
#' @param lfc Numeric vector of log fold changes
#' @param digits Number of decimal places
#' @return Character vector of formatted LFCs
#' @keywords internal
format_logfc <- function(lfc, digits = 2) {
  formatted <- sprintf(paste0("%.", digits, "f"), lfc)
  formatted[is.na(lfc)] <- "NA"
  return(formatted)
}

#' Create a summary table from results
#' @param results_df Data.frame with DE results
#' @param alpha Significance threshold
#' @param lfc_threshold Log fold change threshold
#' @return Data.frame with summary statistics
#' @keywords internal
create_summary_table <- function(results_df, alpha = 0.05, lfc_threshold = 0) {
  summary_df <- data.frame(
    Category = c(
      "Total tested",
      "Significant (padj < alpha)",
      "Upregulated",
      "Downregulated",
      "Large effect (|LFC| > threshold)"
    ),
    Count = c(
      nrow(results_df),
      sum(results_df$padj < alpha, na.rm = TRUE),
      sum(results_df$padj < alpha & results_df$log2FoldChange > lfc_threshold, na.rm = TRUE),
      sum(results_df$padj < alpha & results_df$log2FoldChange < -lfc_threshold, na.rm = TRUE),
      sum(results_df$padj < alpha & abs(results_df$log2FoldChange) > lfc_threshold, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  summary_df$Percentage <- round(100 * summary_df$Count / summary_df$Count[1], 1)
  
  return(summary_df)
}

#' Calculate quantiles of a numeric vector
#' @param x Numeric vector
#' @param probs Quantile probabilities
#' @param na.rm Remove NA values
#' @return Named numeric vector of quantiles
#' @keywords internal
calculate_quantiles <- function(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE) {
  if (all(is.na(x))) {
    result <- rep(NA, length(probs))
    names(result) <- paste0("Q", probs * 100)
    return(result)
  }
  
  quantiles <- quantile(x, probs = probs, na.rm = na.rm)
  names(quantiles) <- paste0("Q", probs * 100)
  
  return(quantiles)
}

#' Find overlapping elements between lists
#' @param list1 First list
#' @param list2 Second list
#' @param method Method for overlap: "intersect", "jaccard", "dice"
#' @return List with overlap information
#' @keywords internal
calculate_overlap <- function(list1, list2, method = "intersect") {
  overlap <- intersect(list1, list2)
  
  result <- list(
    overlap = overlap,
    n_overlap = length(overlap),
    n_list1 = length(list1),
    n_list2 = length(list2)
  )
  
  if (method == "jaccard") {
    union_size <- length(union(list1, list2))
    result$jaccard <- length(overlap) / union_size
    
  } else if (method == "dice") {
    result$dice <- 2 * length(overlap) / (length(list1) + length(list2))
  }
  
  return(result)
}

#' Create a named vector from two vectors
#' @param names Character vector of names
#' @param values Vector of values
#' @param remove_na Remove NA values
#' @return Named vector
#' @keywords internal
create_named_vector <- function(names, values, remove_na = TRUE) {
  if (length(names) != length(values)) {
    stop("names and values must have same length")
  }
  
  result <- setNames(values, names)
  
  if (remove_na) {
    result <- result[!is.na(result)]
  }
  
  return(result)
}

#' Check if object has specific class
#' @param x Object to check
#' @param class_name Name of class
#' @return Logical
#' @keywords internal
has_class <- function(x, class_name) {
  return(class_name %in% class(x))
}

#' Get top N elements from named vector
#' @param x Named numeric vector
#' @param n Number of elements to return
#' @param decreasing Sort in decreasing order
#' @return Named vector with top N elements
#' @keywords internal
get_top_n <- function(x, n = 10, decreasing = TRUE) {
  if (length(x) == 0) {
    return(x)
  }
  
  sorted_x <- sort(x, decreasing = decreasing, na.last = TRUE)
  n <- min(n, length(sorted_x))
  
  return(sorted_x[1:n])
}

#' Create file paths safely
#' @param ... Path components
#' @param create_dirs Create directories if they don't exist
#' @return Character path
#' @keywords internal
safe_file_path <- function(..., create_dirs = FALSE) {
  path <- file.path(...)
  
  if (create_dirs) {
    dir_path <- dirname(path)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
  
  return(path)
}

#' Round numbers intelligently
#' @param x Numeric vector
#' @param digits Number of digits
#' @return Rounded numeric vector
#' @keywords internal
smart_round <- function(x, digits = 3) {
  if (all(is.na(x))) {
    return(x)
  }
  
  # Use more digits for very small numbers
  small_numbers <- abs(x) < 0.001 & abs(x) > 0
  x[small_numbers] <- signif(x[small_numbers], digits)
  x[!small_numbers] <- round(x[!small_numbers], digits)
  
  return(x)
}

#' Get memory usage in human readable format
#' @param object Object to measure
#' @return Character string with memory usage
#' @keywords internal
get_memory_usage <- function(object) {
  size_bytes <- as.numeric(object.size(object))
  
  if (size_bytes < 1024) {
    return(paste(size_bytes, "bytes"))
  } else if (size_bytes < 1024^2) {
    return(paste(round(size_bytes / 1024, 2), "KB"))
  } else if (size_bytes < 1024^3) {
    return(paste(round(size_bytes / 1024^2, 2), "MB"))
  } else {
    return(paste(round(size_bytes / 1024^3, 2), "GB"))
  }
}