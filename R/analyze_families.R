#' Gene Family-level Differential Expression Analysis
#'
#' @description
#' Performs differential expression analysis at the gene family level by
#' aggregating related genes and testing for coordinated changes. This can
#' detect functional changes that might be missed at individual gene level.
#'
#' @param counts A matrix of raw counts (genes x samples)
#' @param metadata A data.frame containing sample information
#' @param design A formula specifying the design matrix
#' @param contrast Character vector of length 3: c("variable", "level1", "level2")
#' @param families Character database name or named list or data.frame of gene families
#' @param species Character, "human" or "mouse" (used when families is a database name)
#' @param aggregation_method Method for aggregating genes: "mean" (default), "median", "sum"
#' @param min_family_size Minimum genes required in family (default: 2)
#' @param max_family_size Maximum genes allowed in family (default: 100)
#' @param method Method for DE analysis: "DESeq2" (default), "edgeR", "limma"
#' @param filter_low_counts Logical, filter lowly expressed families (default: TRUE)
#' @param min_count Minimum count threshold for filtering (default: 10)
#' @param alpha Significance threshold (default: 0.05)
#' @param ... Additional arguments passed to DE method
#'
#' @return A list of class "functionalDE_result" containing:
#'   \item{results}{Data.frame with family-level DE statistics}
#'   \item{family_counts}{Aggregated count matrix at family level}
#'   \item{families_used}{Gene family mappings used}
#'   \item{gene_mapping}{Mapping from genes to families}
#'   \item{counts_raw}{Original gene-level counts}
#'   \item{metadata}{Sample metadata}
#'   \item{method}{Method used}
#'   \item{params}{Parameters used}
#'   \item{summary}{Summary statistics}
#'
#' @details
#' **Gene Family Sources:**
#'
#' *HGNC (Human Gene Nomenclature Committee):*
#' - Authoritative gene family classifications
#' - Based on sequence similarity and function
#' - Available for human genes
#'
#' *Custom families:*
#' - User-provided gene groupings
#' - Can be functional groups, protein complexes, etc.
#' - Provided as named list or data.frame
#'
#' **Aggregation Methods:**
#'
#' *Mean:*
#' - Average expression across family members
#' - Good for detecting overall family activity
#' - Reduces noise from individual genes
#'
#' *Median:*
#' - More robust to outliers
#' - Less sensitive to highly variable family members
#'
#' *Sum:*
#' - Total family expression
#' - Appropriate for count data
#' - Maintains count distributions
#'
#' **Biological Interpretation:**
#' - Family-level changes suggest coordinated regulation
#' - Can detect compensation effects within families
#' - May reveal functional changes missed at gene level
#' - Consider both family-level and gene-level results
#'
#' @examples
#' # Basic family analysis with HGNC families
#' family_results <- analyze_families(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   families = "HGNC",
#'   species = "human"
#' )
#'
#' # View results
#' head(family_results$results)
#' summary(family_results)
#'
#' # With custom families
#' custom_families <- list(
#'   "Ribosomal_proteins" = c("RPS1", "RPS2", "RPL1", "RPL2"),
#'   "Heat_shock_proteins" = c("HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B"),
#'   "Immunoglobulins" = c("IGHG1", "IGHG2", "IGHA1", "IGHA2")
#' )
#'
#' family_results <- analyze_families(
#'   counts = count_matrix,
#'   metadata = sample_info,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   families = custom_families,
#'   aggregation_method = "sum"
#' )
#'
#' @seealso \code{\link{get_gene_families}}, \code{\link{analyze_genes}}, \code{\link{analyze_pathways}}
#' @export
analyze_families <- function(counts,
                            metadata,
                            design,
                            contrast = NULL,
                            families = "HGNC",
                            species = "human",
                            aggregation_method = "mean",
                            min_family_size = 2,
                            max_family_size = 100,
                            method = "DESeq2",
                            filter_low_counts = TRUE,
                            min_count = 10,
                            alpha = 0.05,
                            ...) {
  
  # Validate inputs
  counts <- validate_counts(counts)
  metadata <- validate_metadata(metadata, counts)
  design <- validate_design(design, metadata)
  contrast <- validate_contrast(contrast, metadata)
  
  # Get gene families if database name provided
  if (is.character(families) && length(families) == 1) {
    message("Downloading ", families, " gene families...")
    families <- get_gene_families(
      database = families,
      species = species,
      min_size = min_family_size,
      max_size = max_family_size
    )
  }
  
  # Convert families to gene mapping
  gene_mapping <- create_gene_mapping(families)
  
  # Filter to genes in count matrix
  available_genes <- rownames(counts)
  gene_mapping <- gene_mapping[gene_mapping$gene %in% available_genes, ]
  
  if (nrow(gene_mapping) == 0) {
    stop("No genes from families found in count matrix")
  }
  
  # Create family mapping (gene -> family)
  family_mapping <- setNames(gene_mapping$family, gene_mapping$gene)
  
  # Aggregate genes to families
  message("Aggregating ", length(unique(gene_mapping$family)), 
          " families using ", aggregation_method, " method...")
  
  family_counts <- aggregate_genes(
    counts = counts,
    grouping = family_mapping,
    method = aggregation_method
  )
  
  # Perform differential expression analysis on family counts
  message("Analyzing families using ", method, "...")
  
  family_results <- analyze_genes(
    counts = family_counts,
    metadata = metadata,
    design = design,
    contrast = contrast,
    method = method,
    filter_low_counts = filter_low_counts,
    min_count = min_count,
    alpha = alpha,
    ...
  )
  
  # Add family information to results
  family_results$results$family <- family_results$results$gene
  family_results$results$n_genes <- sapply(families[family_results$results$family], length)
  
  # Add family-specific summary
  family_summary <- list(
    n_families_total = length(families),
    n_families_analyzed = nrow(family_results$results),
    n_families_significant = sum(family_results$results$padj < alpha, na.rm = TRUE),
    median_family_size = median(family_results$results$n_genes),
    aggregation_method = aggregation_method
  )
  
  # Update result object
  family_results$family_counts <- family_counts
  family_results$families_used <- families
  family_results$gene_mapping <- gene_mapping
  family_results$counts_raw <- counts
  family_results$params$level <- "family"
  family_results$params$aggregation_method <- aggregation_method
  family_results$params$min_family_size <- min_family_size
  family_results$params$max_family_size <- max_family_size
  family_results$summary <- c(family_results$summary, family_summary)
  
  return(family_results)
}

#' Convert family definitions to gene mapping
#' @keywords internal
create_gene_mapping <- function(families) {
  if (is.list(families)) {
    # Named list format
    gene_mapping <- data.frame(
      gene = character(),
      family = character(),
      stringsAsFactors = FALSE
    )
    
    for (family_name in names(families)) {
      family_genes <- families[[family_name]]
      gene_mapping <- rbind(gene_mapping, data.frame(
        gene = family_genes,
        family = family_name,
        stringsAsFactors = FALSE
      ))
    }
    
  } else if (is.data.frame(families)) {
    # Data.frame format - expect columns: gene, family
    if (!all(c("gene", "family") %in% colnames(families))) {
      stop("families data.frame must have columns 'gene' and 'family'")
    }
    gene_mapping <- families[, c("gene", "family")]
    
  } else {
    stop("families must be a named list or data.frame")
  }
  
  # Remove duplicates and missing values
  gene_mapping <- gene_mapping[!is.na(gene_mapping$gene) & 
                              !is.na(gene_mapping$family), ]
  gene_mapping <- gene_mapping[!duplicated(gene_mapping), ]
  
  return(gene_mapping)
}

#' Get Gene Family Annotations
#'
#' @description
#' Retrieves gene family annotations from HGNC or other sources.
#' Returns gene family mappings for use with family-level analysis.
#'
#' @param database Character specifying database: "HGNC" (default)
#' @param species Character, "human" or "mouse"
#' @param min_size Minimum family size (default: 2)
#' @param max_size Maximum family size (default: 100)
#' @param include_groups Character vector of family types to include (default: all)
#'
#' @return Named list where names are family names and values are gene vectors
#'
#' @details
#' **HGNC Gene Families:**
#' - Based on sequence similarity and functional relationships
#' - Includes protein families, gene clusters, pseudogenes
#' - Regularly updated by HGNC
#' - Available family types include:
#'   - "protein_coding": Protein-coding gene families
#'   - "pseudogene": Pseudogene families
#'   - "lncRNA": Long non-coding RNA families
#'   - "miRNA": microRNA families
#'
#' **Custom Families:**
#' - Users can provide custom family definitions
#' - Useful for functional groups, complexes, or pathways
#' - Should be provided as named list or data.frame
#'
#' @examples
#' # Get HGNC families for human
#' hgnc_families <- get_gene_families("HGNC", species = "human")
#' length(hgnc_families)
#'
#' # Get only protein-coding families
#' protein_families <- get_gene_families(
#'   "HGNC", 
#'   species = "human",
#'   include_groups = "protein_coding"
#' )
#'
#' @export
get_gene_families <- function(database = "HGNC",
                             species = "human",
                             min_size = 2,
                             max_size = 100,
                             include_groups = NULL) {
  
  if (database == "HGNC") {
    families <- get_hgnc_families(species, include_groups)
  } else {
    stop("Currently only HGNC database is supported")
  }
  
  # Filter by size
  family_sizes <- sapply(families, length)
  families <- families[family_sizes >= min_size & family_sizes <= max_size]
  
  message("Retrieved ", length(families), " gene families")
  message("Size range: ", min(sapply(families, length)), " - ", 
          max(sapply(families, length)), " genes")
  
  # Add attributes
  attr(families, "database") <- database
  attr(families, "species") <- species
  attr(families, "date_retrieved") <- Sys.Date()
  attr(families, "n_families") <- length(families)
  
  return(families)
}

#' Get HGNC gene families
#' @keywords internal
get_hgnc_families <- function(species, include_groups) {
  # This is a simplified implementation
  # In a real package, this would fetch from HGNC REST API or local database
  
  if (species != "human") {
    stop("HGNC families are only available for human. For mouse, consider using ortholog mapping.")
  }
  
  # Example families (in real implementation, this would be comprehensive)
  families <- list(
    "Ribosomal_proteins_cytosolic_40S" = c(
      "RPS1", "RPS2", "RPS3", "RPS4X", "RPS4Y1", "RPS5", "RPS6", "RPS7",
      "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15",
      "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23",
      "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS29"
    ),
    "Ribosomal_proteins_cytosolic_60S" = c(
      "RPL1", "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8",
      "RPL9", "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A",
      "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21",
      "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A",
      "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35",
      "RPL35A", "RPL36", "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39",
      "RPL40", "RPL41"
    ),
    "Heat_shock_proteins_70kDa" = c(
      "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5",
      "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14"
    ),
    "Heat_shock_proteins_90kDa" = c(
      "HSP90AA1", "HSP90AB1", "HSP90B1", "TRAP1"
    ),
    "Immunoglobulin_heavy_chains" = c(
      "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
      "IGHM", "IGHD1-1", "IGHD1-7", "IGHD1-14", "IGHD1-20", "IGHD1-26",
      "IGHD2-2", "IGHD2-8", "IGHD2-15", "IGHD2-21", "IGHD3-3", "IGHD3-9",
      "IGHD3-10", "IGHD3-16", "IGHD3-22", "IGHD4-4", "IGHD4-11", "IGHD4-17",
      "IGHD4-23", "IGHD5-5", "IGHD5-12", "IGHD5-18", "IGHD5-24",
      "IGHD6-6", "IGHD6-13", "IGHD6-19", "IGHD6-25", "IGHD7-27"
    ),
    "Immunoglobulin_kappa_chains" = c(
      "IGKC", "IGKJ1", "IGKJ2", "IGKJ3", "IGKJ4", "IGKJ5",
      "IGKV1-5", "IGKV1-6", "IGKV1-8", "IGKV1-9", "IGKV1-12",
      "IGKV1-13", "IGKV1-16", "IGKV1-17", "IGKV1-27", "IGKV1-33",
      "IGKV1-37", "IGKV1-39", "IGKV1D-8", "IGKV1D-12", "IGKV1D-13",
      "IGKV1D-16", "IGKV1D-17", "IGKV1D-22", "IGKV1D-27", "IGKV1D-32",
      "IGKV1D-33", "IGKV1D-37", "IGKV1D-39", "IGKV1D-42", "IGKV1D-43"
    ),
    "Histones_H1" = c(
      "HIST1H1A", "HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E",
      "HIST1H1T", "H1F0", "H1FX"
    ),
    "Histones_H2A" = c(
      "HIST1H2AA", "HIST1H2AB", "HIST1H2AC", "HIST1H2AD", "HIST1H2AE",
      "HIST1H2AG", "HIST1H2AI", "HIST1H2AJ", "HIST1H2AK", "HIST1H2AL",
      "HIST1H2AM", "HIST2H2AA3", "HIST2H2AA4", "HIST2H2AC", "H2AFX",
      "H2AFY", "H2AFY2", "H2AFZ", "MACROH2A1", "MACROH2A2"
    ),
    "Histones_H2B" = c(
      "HIST1H2BA", "HIST1H2BB", "HIST1H2BC", "HIST1H2BD", "HIST1H2BE",
      "HIST1H2BF", "HIST1H2BG", "HIST1H2BH", "HIST1H2BI", "HIST1H2BJ",
      "HIST1H2BK", "HIST1H2BL", "HIST1H2BM", "HIST1H2BN", "HIST1H2BO",
      "HIST2H2BE", "HIST2H2BF"
    ),
    "Histones_H3" = c(
      "HIST1H3A", "HIST1H3B", "HIST1H3C", "HIST1H3D", "HIST1H3E",
      "HIST1H3F", "HIST1H3G", "HIST1H3H", "HIST1H3I", "HIST1H3J",
      "HIST2H3A", "HIST2H3C", "HIST2H3D", "H3F3A", "H3F3B", "H3F3C"
    ),
    "Histones_H4" = c(
      "HIST1H4A", "HIST1H4B", "HIST1H4C", "HIST1H4D", "HIST1H4E",
      "HIST1H4F", "HIST1H4G", "HIST1H4H", "HIST1H4I", "HIST1H4J",
      "HIST1H4K", "HIST1H4L", "HIST2H4A", "HIST2H4B", "HIST4H4"
    ),
    "Metallothioneins" = c(
      "MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H", "MT1JP", "MT1L",
      "MT1M", "MT1X", "MT2A", "MT3", "MT4"
    )
  )
  
  # Filter by groups if specified
  if (!is.null(include_groups)) {
    warning("Group filtering not implemented in this simplified version")
  }
  
  return(families)
}