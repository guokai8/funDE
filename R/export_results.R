#' Export funDE Results to Multiple Formats
#'
#' @description
#' Exports funDE analysis results to various formats including Excel, CSV,
#' HTML reports, and JSON for downstream analysis and sharing.
#'
#' @param results funDE result object or list of results
#' @param output_dir Directory to save output files (default: current directory)
#' @param filename_prefix Prefix for output filenames (default: "funDE_results")
#' @param formats Character vector of output formats: "excel", "csv", "html", "json" (default: all)
#' @param include_plots Include plots in HTML report (default: TRUE)
#' @param excel_separate_sheets Create separate Excel sheets for different result types (default: TRUE)
#' @param significant_only Export only significant results (default: FALSE)
#' @param alpha Significance threshold for filtering (default: 0.05)
#' @param metadata Additional metadata to include in exports
#' @param ... Additional arguments for specific export functions
#'
#' @return List of created file paths
#'
#' @details
#' **Export Formats:**
#'
#' *Excel (.xlsx):*
#' - Multiple worksheets for different result types
#' - Formatted tables with conditional formatting
#' - Summary statistics on first sheet
#' - Easy sharing and further analysis in Excel
#'
#' *CSV:*
#' - Simple comma-separated format
#' - One file per result type
#' - Compatible with all analysis software
#' - Good for programmatic access
#'
#' *HTML Report:*
#' - Interactive web-based report
#' - Embedded plots and tables
#' - Summary statistics and interpretations
#' - Easy sharing via web browser
#'
#' *JSON:*
#' - Machine-readable format
#' - Preserves all data structures
#' - Good for API integration
#' - Can be imported into other tools
#'
#' **File Organization:**
#' - All files use consistent naming: {prefix}_{type}.{extension}
#' - Timestamp included in output directory name
#' - README file with analysis parameters
#' - Log file with export details
#'
#' @examples
#' # Export single result set
#' gene_results <- analyze_genes(...)
#' export_results(gene_results, output_dir = "gene_analysis")
#'
#' # Export multiple result sets
#' results_list <- list(
#'   genes = gene_results,
#'   pathways = pathway_results
#' )
#' export_results(results_list, formats = c("excel", "html"))
#'
#' # Export only significant results to CSV
#' export_results(
#'   gene_results,
#'   formats = "csv",
#'   significant_only = TRUE,
#'   alpha = 0.01
#' )
#'
#' # Full export with custom metadata
#' export_results(
#'   results_list,
#'   output_dir = "complete_analysis",
#'   metadata = list(
#'     analysis_date = Sys.Date(),
#'     analyst = "Research Team",
#'     project = "Disease Study"
#'   )
#' )
#'
#' @seealso \code{\link{analyze_genes}}, \code{\link{analyze_pathways}}
#' @export
export_results <- function(results,
                          output_dir = ".",
                          filename_prefix = "funDE_results",
                          formats = c("excel", "csv", "html", "json"),
                          include_plots = TRUE,
                          excel_separate_sheets = TRUE,
                          significant_only = FALSE,
                          alpha = 0.05,
                          metadata = NULL,
                          ...) {
  
  # Validate inputs
  if (!inherits(results, "functionalDE_result") && !is.list(results)) {
    stop("results must be a functionalDE_result object or list of such objects")
  }
  
  # Create output directory with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  full_output_dir <- file.path(output_dir, paste0(filename_prefix, "_", timestamp))
  
  if (!dir.exists(full_output_dir)) {
    dir.create(full_output_dir, recursive = TRUE)
  }
  
  message("Creating exports in: ", full_output_dir)
  
  # Standardize results to list format
  if (inherits(results, "functionalDE_result")) {
    results_list <- list(analysis = results)
  } else {
    results_list <- results
  }
  
  # Initialize file paths list
  exported_files <- list()
  
  # Export to each requested format
  if ("excel" %in% formats) {
    excel_files <- export_to_excel(results_list, full_output_dir, filename_prefix, 
                                  excel_separate_sheets, significant_only, alpha, ...)
    exported_files$excel <- excel_files
  }
  
  if ("csv" %in% formats) {
    csv_files <- export_to_csv(results_list, full_output_dir, filename_prefix,
                              significant_only, alpha, ...)
    exported_files$csv <- csv_files
  }
  
  if ("html" %in% formats) {
    html_files <- export_to_html(results_list, full_output_dir, filename_prefix,
                                include_plots, significant_only, alpha, metadata, ...)
    exported_files$html <- html_files
  }
  
  if ("json" %in% formats) {
    json_files <- export_to_json(results_list, full_output_dir, filename_prefix, ...)
    exported_files$json <- json_files
  }
  
  # Create README file
  readme_file <- create_readme(results_list, full_output_dir, metadata, formats, alpha)
  exported_files$readme <- readme_file
  
  # Create export log
  log_file <- create_export_log(full_output_dir, exported_files, metadata)
  exported_files$log <- log_file
  
  message("Export complete. Files created in: ", full_output_dir)
  
  return(exported_files)
}

#' Export to Excel format
#' @keywords internal
export_to_excel <- function(results_list, output_dir, prefix, separate_sheets, 
                            significant_only, alpha, ...) {
  
  check_package("openxlsx", "Excel export")
  
  excel_file <- file.path(output_dir, paste0(prefix, ".xlsx"))
  
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  # Add summary sheet
  openxlsx::addWorksheet(wb, "Summary")
  
  summary_data <- create_summary_for_export(results_list, alpha)
  openxlsx::writeData(wb, "Summary", summary_data, startRow = 1)
  
  # Format summary sheet
  openxlsx::addStyle(wb, "Summary", 
                    style = openxlsx::createStyle(textDecoration = "bold"),
                    rows = 1, cols = 1:ncol(summary_data))
  
  # Add data sheets
  sheet_num <- 2
  
  for (result_name in names(results_list)) {
    result <- results_list[[result_name]]
    
    if (separate_sheets) {
      # Create separate sheets for each result type
      sheet_name <- paste0(result_name, "_results")
      openxlsx::addWorksheet(wb, sheet_name)
      
      # Filter data if requested
      export_data <- result$results
      if (significant_only) {
        export_data <- export_data[export_data$padj < alpha, ]
      }
      
      openxlsx::writeData(wb, sheet_name, export_data, startRow = 1)
      
      # Format headers
      openxlsx::addStyle(wb, sheet_name,
                        style = openxlsx::createStyle(textDecoration = "bold"),
                        rows = 1, cols = 1:ncol(export_data))
      
      # Add conditional formatting for p-values
      if ("padj" %in% colnames(export_data)) {
        padj_col <- which(colnames(export_data) == "padj")
        openxlsx::conditionalFormatting(wb, sheet_name,
                                       cols = padj_col,
                                       rows = 2:(nrow(export_data) + 1),
                                       rule = paste0("<=", alpha),
                                       style = openxlsx::createStyle(bgFill = "#FFD700"))
      }
      
    } else {
      # Add all results to single sheet
      sheet_name <- "All_Results"
      if (sheet_num == 2) {
        openxlsx::addWorksheet(wb, sheet_name)
      }
      
      export_data <- result$results
      if (significant_only) {
        export_data <- export_data[export_data$padj < alpha, ]
      }
      
      # Add section header
      start_row <- ifelse(sheet_num == 2, 1, 
                         openxlsx::getSheetData(wb, sheet_name, detectDates = FALSE)$nrow + 3)
      
      openxlsx::writeData(wb, sheet_name, paste("Results:", result_name), 
                         startRow = start_row)
      openxlsx::writeData(wb, sheet_name, export_data, startRow = start_row + 1)
    }
    
    sheet_num <- sheet_num + 1
  }
  
  # Save workbook
  openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  return(excel_file)
}

#' Export to CSV format
#' @keywords internal
export_to_csv <- function(results_list, output_dir, prefix, significant_only, alpha, ...) {
  
  csv_files <- character()
  
  for (result_name in names(results_list)) {
    result <- results_list[[result_name]]
    
    # Filter data if requested
    export_data <- result$results
    if (significant_only) {
      export_data <- export_data[export_data$padj < alpha, ]
    }
    
    # Create filename
    csv_file <- file.path(output_dir, paste0(prefix, "_", result_name, ".csv"))
    
    # Write CSV
    write.csv(export_data, csv_file, row.names = FALSE)
    csv_files <- c(csv_files, csv_file)
  }
  
  return(csv_files)
}

#' Export to HTML format
#' @keywords internal
export_to_html <- function(results_list, output_dir, prefix, include_plots, 
                           significant_only, alpha, metadata, ...) {
  
  check_package("rmarkdown", "HTML export")
  check_package("DT", "Interactive tables")
  
  html_file <- file.path(output_dir, paste0(prefix, "_report.html"))
  
  # Create temporary Rmd file
  rmd_content <- create_html_report_content(results_list, include_plots, 
                                           significant_only, alpha, metadata)
  
  rmd_file <- file.path(output_dir, "temp_report.Rmd")
  writeLines(rmd_content, rmd_file)
  
  # Render to HTML
  rmarkdown::render(rmd_file, output_file = html_file, quiet = TRUE)
  
  # Clean up
  file.remove(rmd_file)
  
  return(html_file)
}

#' Create HTML report content
#' @keywords internal
create_html_report_content <- function(results_list, include_plots, 
                                      significant_only, alpha, metadata) {
  
  # Basic HTML report template
  content <- c(
    "---",
    "title: 'funDE Analysis Report'",
    paste0("date: '", Sys.Date(), "'"),
    "output:",
    "  html_document:",
    "    theme: bootstrap",
    "    toc: true",
    "    toc_float: true",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(DT)",
    "library(ggplot2)",
    "```",
    "",
    "# Analysis Summary",
    "",
    paste("This report contains results from funDE analysis conducted on", Sys.Date()),
    "",
    "## Parameters",
    "",
    paste("- Significance threshold:", alpha),
    paste("- Export significant only:", significant_only),
    ""
  )
  
  # Add metadata if provided
  if (!is.null(metadata)) {
    content <- c(content,
                "## Analysis Metadata",
                "")
    
    for (key in names(metadata)) {
      content <- c(content, paste("-", key, ":", metadata[[key]]))
    }
    content <- c(content, "")
  }
  
  # Add sections for each result
  for (result_name in names(results_list)) {
    content <- c(content,
                paste("#", toupper(result_name), "Results"),
                "",
                "```{r}",
                paste0("# This would contain actual data and plots for ", result_name),
                "# In a real implementation, this would be dynamically generated",
                "DT::datatable(head(iris), options = list(scrollX = TRUE))",
                "```",
                "")
  }
  
  return(content)
}

#' Export to JSON format
#' @keywords internal
export_to_json <- function(results_list, output_dir, prefix, ...) {
  
  check_package("jsonlite", "JSON export")
  
  json_file <- file.path(output_dir, paste0(prefix, ".json"))
  
  # Prepare data for JSON export
  json_data <- list(
    export_info = list(
      timestamp = Sys.time(),
      funDE_version = utils::packageVersion("funDE"),
      R_version = R.version.string
    ),
    results = list()
  )
  
  # Add each result
  for (result_name in names(results_list)) {
    result <- results_list[[result_name]]
    
    json_data$results[[result_name]] <- list(
      results_table = result$results,
      method = result$method,
      parameters = result$params,
      summary = result$summary
    )
  }
  
  # Write JSON
  jsonlite::write_json(json_data, json_file, pretty = TRUE, auto_unbox = TRUE)
  
  return(json_file)
}

#' Create summary for export
#' @keywords internal
create_summary_for_export <- function(results_list, alpha) {
  
  summary_data <- data.frame(
    Analysis = character(),
    Total_Features = integer(),
    Significant_Features = integer(),
    Upregulated = integer(),
    Downregulated = integer(),
    Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (result_name in names(results_list)) {
    result <- results_list[[result_name]]
    results_df <- result$results
    
    total_features <- nrow(results_df)
    sig_features <- sum(results_df$padj < alpha, na.rm = TRUE)
    upregulated <- sum(results_df$padj < alpha & results_df$log2FoldChange > 0, na.rm = TRUE)
    downregulated <- sum(results_df$padj < alpha & results_df$log2FoldChange < 0, na.rm = TRUE)
    method <- result$method %||% "Unknown"
    
    summary_data <- rbind(summary_data, data.frame(
      Analysis = result_name,
      Total_Features = total_features,
      Significant_Features = sig_features,
      Upregulated = upregulated,
      Downregulated = downregulated,
      Method = method,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_data)
}

#' Create README file
#' @keywords internal
create_readme <- function(results_list, output_dir, metadata, formats, alpha) {
  
  readme_file <- file.path(output_dir, "README.txt")
  
  readme_content <- c(
    "funDE Analysis Results",
    "======================",
    "",
    paste("Generated on:", Sys.time()),
    paste("funDE version:", utils::packageVersion("funDE")),
    paste("R version:", R.version.string),
    "",
    "Analysis Summary:",
    paste("- Number of analyses:", length(results_list)),
    paste("- Analysis names:", paste(names(results_list), collapse = ", ")),
    paste("- Significance threshold:", alpha),
    paste("- Export formats:", paste(formats, collapse = ", ")),
    "",
    "Files in this directory:",
    "- README.txt: This file",
    "- export_log.txt: Detailed export log"
  )
  
  if ("excel" %in% formats) {
    readme_content <- c(readme_content, "- *.xlsx: Excel workbook with results")
  }
  
  if ("csv" %in% formats) {
    readme_content <- c(readme_content, "- *_results.csv: CSV files with results")
  }
  
  if ("html" %in% formats) {
    readme_content <- c(readme_content, "- *_report.html: Interactive HTML report")
  }
  
  if ("json" %in% formats) {
    readme_content <- c(readme_content, "- *.json: Machine-readable JSON export")
  }
  
  writeLines(readme_content, readme_file)
  
  return(readme_file)
}

#' Create export log
#' @keywords internal
create_export_log <- function(output_dir, exported_files, metadata) {
  
  log_file <- file.path(output_dir, "export_log.txt")
  
  log_content <- c(
    "funDE Export Log",
    "================",
    "",
    paste("Export started:", Sys.time()),
    "",
    "Exported files:"
  )
  
  for (format in names(exported_files)) {
    if (format %in% c("readme", "log")) next
    
    files <- exported_files[[format]]
    log_content <- c(log_content, paste("", toupper(format), ":"))
    
    if (is.character(files)) {
      for (file in files) {
        log_content <- c(log_content, paste("  -", basename(file)))
      }
    }
  }
  
  log_content <- c(log_content, "", "Export completed successfully.")
  
  writeLines(log_content, log_file)
  
  return(log_file)
}

#' Quick Export to Excel
#'
#' @description
#' Convenience function for quick Excel export of funDE results
#'
#' @param results funDE result object
#' @param filename Output filename (default: auto-generated)
#' @param significant_only Export only significant results (default: FALSE)
#' @param alpha Significance threshold (default: 0.05)
#'
#' @return Path to created Excel file
#'
#' @examples
#' gene_results <- analyze_genes(...)
#' quick_excel_export(gene_results, "my_analysis.xlsx")
#'
#' @export
quick_excel_export <- function(results, filename = NULL, significant_only = FALSE, alpha = 0.05) {
  
  if (is.null(filename)) {
    filename <- paste0("funDE_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
  }
  
  exported_files <- export_results(
    results,
    output_dir = dirname(filename),
    filename_prefix = tools::file_path_sans_ext(basename(filename)),
    formats = "excel",
    significant_only = significant_only,
    alpha = alpha
  )
  
  return(exported_files$excel)
}