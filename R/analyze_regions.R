#' Differential Expression Analysis at Regional Level
#'
#' @description
#' Performs differential analysis on transcripts and isoforms using DESeq2
#' or edgeR at the transcript level. For more advanced exon-level and DTU
#' analysis, consider using specialized packages externally.
#'
#' @param counts A matrix of counts at transcript level (transcripts x samples)
#' @param metadata A data.frame containing sample information
#' @param design A formula specifying the design matrix
#' @param contrast Character vector specifying the contrast
#' @param annotation Data.frame containing transcript annotations with columns:
#'   transcript_id, gene_id, gene_name
#' @param method Method for DE: "DESeq2" (default) or "edgeR"
#' @param filter_low Logical, filter low-expressed transcripts (default: TRUE)
#' @param min_count Minimum count threshold (default: 10)
#' @param alpha Significance threshold (default: 0.05)
#' @param ... Additional arguments passed to the DE method
#'
#' @return A list of class "functionalDE_result" containing:
#'   \item{results}{Data.frame with transcript-level DE statistics}
#'   \item{gene_summary}{Data.frame summarizing results per gene}
#'   \item{annotation}{Transcript annotation used}
#'   \item{method}{Method used}
#'
#' @details
#' This simplified function performs differential expression at the transcript
#' level. For advanced analyses like:
#' - Exon-level differential usage: Use DEXSeq package directly
#' - Differential transcript usage (DTU): Use DRIMSeq or satuRn directly
#' - Splice junction analysis: Use leafcutter or other specialized tools
#'
#' **Biological Interpretation:**
#' - Multiple transcripts from same gene with different patterns suggests isoform switching
#' - Compare transcript-level results with gene-level to identify DTU
#' - Use annotation to understand functional consequences (domains, UTRs, etc.)
#'
#' @examples
#' # Transcript-level analysis
#' transcript_results <- analyze_regions(
#'   counts = transcript_counts,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   annotation = transcript_annotation,
#'   method = "DESeq2"
#' )
#'
#' # Identify potential isoform switches
#' # (genes with multiple significant transcripts with opposite directions)
#' gene_summary <- transcript_results$gene_summary
#' potential_switches <- gene_summary %>%
#'   filter(n_sig_transcripts >= 2, 
#'          has_opposite_direction == TRUE)
#'
#' @export
analyze_regions <- function(counts,
                           metadata,
                           design,
                           contrast,
                           annotation = NULL,
                           method = "DESeq2",
                           filter_low = TRUE,
                           min_count = 10,
                           alpha = 0.05,
                           ...) {
  
  message("Note: This function provides transcript-level DE analysis.")
  message("For advanced DTU/exon analysis, consider using DRIMSeq/DEXSeq directly.")
  
  # Validate inputs
  counts <- validate_counts(counts)
  metadata <- validate_metadata(metadata, counts)
  
  # Perform standard DE analysis at transcript level
  results <- analyze_genes(
    counts = counts,
    metadata = metadata,
    design = design,
    contrast = contrast,
    method = method,
    filter_low_counts = filter_low,
    min_count = min_count,
    alpha = alpha,
    ...
  )
  
  # Add annotation if provided
  if (!is.null(annotation)) {
    # Validate annotation
    required_cols <- c("transcript_id", "gene_id")
    if (!all(required_cols %in% colnames(annotation))) {
      warning("annotation should contain columns: transcript_id, gene_id")
    } else {
      # Merge with results
      results$results <- results$results %>%
        left_join(annotation, by = c("gene" = "transcript_id"))
    }
  }
  
  # Create gene-level summary if annotation is available
  if (!is.null(annotation) && "gene_id" %in% colnames(annotation)) {
    
    # Ensure results have gene_id
    if (!"gene_id" %in% colnames(results$results)) {
      warning("Could not create gene summary - gene_id not found in results")
      gene_summary <- NULL
    } else {
      gene_summary <- results$results %>%
        filter(!is.na(gene_id)) %>%
        group_by(gene_id) %>%
        summarise(
          n_transcripts = n(),
          n_sig_transcripts = sum(padj < alpha, na.rm = TRUE),
          has_opposite_direction = length(unique(sign(log2FoldChange[padj < alpha]))) > 1,
          max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
          min_padj = min(padj, na.rm = TRUE),
          most_sig_transcript = gene[which.min(padj)],
          mean_expression = mean(baseMean, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(min_padj)
    }
  } else {
    gene_summary <- NULL
  }
  
  # Add region-specific information to results
  results$gene_summary <- gene_summary
  results$annotation <- annotation
  results$params$level <- "transcript"
  results$params$min_count <- min_count
  
  # Update summary with transcript-specific metrics
  if (!is.null(gene_summary)) {
    region_summary <- list(
      n_genes_with_transcripts = nrow(gene_summary),
      n_genes_with_sig_transcripts = sum(gene_summary$n_sig_transcripts > 0),
      n_potential_isoform_switches = sum(gene_summary$has_opposite_direction, na.rm = TRUE),
      median_transcripts_per_gene = median(gene_summary$n_transcripts),
      max_transcripts_per_gene = max(gene_summary$n_transcripts)
    )
    
    results$summary <- c(results$summary, region_summary)
  }
  
  return(results)
}

#' Identify Potential Isoform Switches
#'
#' @description
#' Identifies genes with potential isoform switching based on transcript-level
#' differential expression results.
#'
#' @param transcript_results Output from analyze_regions()
#' @param min_transcripts Minimum transcripts per gene to consider (default: 2)
#' @param alpha Significance threshold (default: 0.05)
#' @param min_lfc Minimum log fold change threshold (default: 0.5)
#'
#' @return Data.frame with potential isoform switches
#'
#' @examples
#' transcript_results <- analyze_regions(...)
#' switches <- identify_isoform_switches(transcript_results)
#'
#' @export
identify_isoform_switches <- function(transcript_results,
                                     min_transcripts = 2,
                                     alpha = 0.05,
                                     min_lfc = 0.5) {
  
  if (is.null(transcript_results$gene_summary)) {
    stop("Gene summary not available. Ensure annotation was provided to analyze_regions()")
  }
  
  # Filter genes with potential switches
  potential_switches <- transcript_results$gene_summary %>%
    filter(
      n_transcripts >= min_transcripts,
      n_sig_transcripts >= 2,
      has_opposite_direction == TRUE,
      max_abs_lfc >= min_lfc
    )
  
  # Add detailed transcript information for each switch
  if (nrow(potential_switches) > 0) {
    switch_details <- list()
    
    for (i in seq_len(nrow(potential_switches))) {
      gene_id <- potential_switches$gene_id[i]
      
      gene_transcripts <- transcript_results$results %>%
        filter(gene_id == !!gene_id, padj < alpha) %>%
        arrange(padj) %>%
        select(gene, log2FoldChange, padj, baseMean)
      
      switch_details[[gene_id]] <- gene_transcripts
    }
    
    potential_switches$transcript_details <- switch_details
  }
  
  return(potential_switches)
}

#' Calculate Transcript Usage Proportions
#'
#' @description
#' Calculates the proportion of each transcript relative to total gene expression
#'
#' @param counts Transcript count matrix
#' @param annotation Data.frame with transcript_id and gene_id columns
#'
#' @return Matrix of transcript usage proportions
#'
#' @examples
#' proportions <- calculate_transcript_proportions(transcript_counts, annotation)
#'
#' @export
calculate_transcript_proportions <- function(counts, annotation) {
  
  if (is.null(annotation)) {
    stop("annotation required with transcript_id and gene_id columns")
  }
  
  # Match transcripts to annotation
  common_transcripts <- intersect(rownames(counts), annotation$transcript_id)
  
  if (length(common_transcripts) == 0) {
    stop("No transcripts found in both counts and annotation")
  }
  
  counts_subset <- counts[common_transcripts, , drop = FALSE]
  annotation_subset <- annotation[annotation$transcript_id %in% common_transcripts, ]
  
  # Calculate gene-level expression
  gene_expression <- aggregate_genes(
    counts = counts_subset,
    grouping = setNames(annotation_subset$gene_id, annotation_subset$transcript_id),
    method = "sum"
  )
  
  # Calculate proportions
  proportions <- matrix(0, nrow = nrow(counts_subset), ncol = ncol(counts_subset))
  rownames(proportions) <- rownames(counts_subset)
  colnames(proportions) <- colnames(counts_subset)
  
  for (transcript in rownames(counts_subset)) {
    gene_id <- annotation_subset$gene_id[annotation_subset$transcript_id == transcript]
    
    if (length(gene_id) > 0 && gene_id %in% rownames(gene_expression)) {
      gene_counts <- gene_expression[gene_id, ]
      transcript_counts <- counts_subset[transcript, ]
      
      # Avoid division by zero
      proportions[transcript, ] <- ifelse(gene_counts > 0, 
                                         transcript_counts / gene_counts, 
                                         0)
    }
  }
  
  return(proportions)
}

#' Compare Transcript and Gene-level Results
#'
#' @description
#' Compares results from transcript-level and gene-level analyses to identify
#' cases where transcript-level analysis provides additional insights
#'
#' @param transcript_results Output from analyze_regions()
#' @param gene_results Output from analyze_genes()
#' @param alpha Significance threshold (default: 0.05)
#'
#' @return List with comparison results
#'
#' @examples
#' transcript_results <- analyze_regions(...)
#' gene_results <- analyze_genes(...)
#' comparison <- compare_transcript_gene_results(transcript_results, gene_results)
#'
#' @export
compare_transcript_gene_results <- function(transcript_results,
                                           gene_results,
                                           alpha = 0.05) {
  
  if (is.null(transcript_results$gene_summary)) {
    stop("Gene summary not available in transcript results")
  }
  
  # Get significant genes and transcripts
  sig_genes <- gene_results$results %>%
    filter(padj < alpha) %>%
    pull(gene)
  
  sig_transcripts_by_gene <- transcript_results$gene_summary %>%
    filter(n_sig_transcripts > 0)
  
  # Find genes significant at transcript but not gene level
  transcript_specific <- sig_transcripts_by_gene %>%
    filter(!gene_id %in% sig_genes)
  
  # Find genes with isoform switches
  isoform_switches <- sig_transcripts_by_gene %>%
    filter(has_opposite_direction == TRUE)
  
  # Calculate overlap statistics
  overlap_stats <- list(
    n_genes_tested = nrow(gene_results$results),
    n_genes_sig = length(sig_genes),
    n_genes_with_sig_transcripts = nrow(sig_transcripts_by_gene),
    n_transcript_specific = nrow(transcript_specific),
    n_isoform_switches = nrow(isoform_switches),
    genes_both_significant = length(intersect(sig_genes, sig_transcripts_by_gene$gene_id))
  )
  
  result <- list(
    overlap_stats = overlap_stats,
    transcript_specific_genes = transcript_specific,
    isoform_switches = isoform_switches,
    comparison_summary = data.frame(
      Analysis = c("Gene-level", "Transcript-level", "Both", "Transcript-specific"),
      N_Significant = c(
        overlap_stats$n_genes_sig,
        overlap_stats$n_genes_with_sig_transcripts,
        overlap_stats$genes_both_significant,
        overlap_stats$n_transcript_specific
      )
    )
  )
  
  return(result)
}