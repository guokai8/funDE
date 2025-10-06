#' Example Count Data
#'
#' @description
#' A simulated RNA-seq count matrix for demonstrating funDE functionality.
#' Contains 2000 genes across 24 samples with two experimental conditions.
#'
#' @format A matrix with 2000 rows (genes) and 24 columns (samples):
#' \describe{
#'   \item{rows}{Gene symbols (e.g., "GENE_0001", "GENE_0002", ...)}
#'   \item{columns}{Sample IDs (e.g., "Sample_01", "Sample_02", ...)}
#'   \item{values}{Raw RNA-seq read counts (integers)}
#' }
#'
#' @details
#' This dataset was simulated to include:
#' - 2000 genes total
#' - 12 samples per condition (treated vs control)
#' - ~200 truly differentially expressed genes
#' - Realistic count distributions and overdispersion
#' - Gene names that can be mapped to pathway databases
#'
#' The data follows a negative binomial distribution typical of RNA-seq data
#' and includes both highly and lowly expressed genes to represent realistic
#' experimental conditions.
#'
#' @source Simulated data created for funDE package demonstration
#'
#' @examples
#' data(example_counts)
#' dim(example_counts)
#' head(example_counts[, 1:6])
#'
#' # Basic analysis
#' data(example_metadata)
#' gene_results <- analyze_genes(
#'   counts = example_counts,
#'   metadata = example_metadata,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control")
#' )
#'
"example_counts"

#' Example Sample Metadata
#'
#' @description
#' Sample metadata corresponding to the example_counts dataset.
#' Contains experimental design information for 24 samples.
#'
#' @format A data frame with 24 rows and 4 columns:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifiers matching colnames of example_counts}
#'   \item{condition}{Factor. Experimental condition: "treated" or "control"}
#'   \item{batch}{Factor. Sequencing batch: "batch1", "batch2", or "batch3"}
#'   \item{replicate}{Integer. Biological replicate number within condition (1-12)}
#' }
#'
#' @details
#' This metadata represents a balanced experimental design with:
#' - 12 treated samples
#' - 12 control samples  
#' - 3 sequencing batches (8 samples each)
#' - Balanced representation across batches and conditions
#'
#' The rownames match the column names of example_counts to enable
#' proper sample matching in funDE analyses.
#'
#' @source Simulated metadata created for funDE package demonstration
#'
#' @examples
#' data(example_metadata)
#' str(example_metadata)
#' table(example_metadata$condition)
#' table(example_metadata$batch, example_metadata$condition)
#'
#' # Use with example_counts
#' data(example_counts)
#' all(colnames(example_counts) == rownames(example_metadata))
#'
"example_metadata"

#' HGNC Gene Family Mappings
#'
#' @description
#' A subset of Human Gene Nomenclature Committee (HGNC) gene family
#' classifications for use with family-level analysis.
#'
#' @format A named list with gene family mappings:
#' \describe{
#'   \item{names}{Gene family names (e.g., "Ribosomal_proteins_cytosolic_40S")}
#'   \item{values}{Character vectors of gene symbols belonging to each family}
#' }
#'
#' @details
#' This dataset contains curated gene family assignments from HGNC including:
#' - Ribosomal protein families (40S and 60S subunits)
#' - Heat shock protein families (HSP70, HSP90)
#' - Immunoglobulin families (heavy and light chains)
#' - Histone protein families (H1, H2A, H2B, H3, H4)
#' - Metallothionein family
#'
#' These families represent well-characterized functional groups that
#' demonstrate coordinated regulation in many biological contexts.
#'
#' @source Human Gene Nomenclature Committee (HGNC) via genenames.org
#'
#' @examples
#' data(hgnc_families)
#' length(hgnc_families)
#' names(hgnc_families)[1:5]
#' head(hgnc_families$Ribosomal_proteins_cytosolic_40S)
#'
#' # Use for family analysis
#' data(example_counts, example_metadata)
#' family_results <- analyze_families(
#'   counts = example_counts,
#'   metadata = example_metadata,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   families = hgnc_families
#' )
#'
"hgnc_families"

#' Hallmark Pathway Gene Sets
#'
#' @description
#' MSigDB Hallmark gene sets representing well-defined biological states
#' and processes. This collection represents coherent gene sets for
#' pathway-level analysis.
#'
#' @format A named list with pathway mappings:
#' \describe{
#'   \item{names}{Hallmark pathway names (e.g., "HALLMARK_GLYCOLYSIS")}
#'   \item{values}{Character vectors of gene symbols in each pathway}
#' }
#'
#' @details
#' The Hallmark collection contains 50 gene sets representing:
#' - Metabolic processes (glycolysis, oxidative phosphorylation)
#' - Signaling pathways (mTORC1, p53, TNF-α)
#' - Cellular processes (DNA repair, cell cycle)
#' - Immune responses (inflammatory response, complement)
#' - Stress responses (hypoxia, reactive oxygen species)
#'
#' These pathways are designed to be:
#' - Coherent and interpretable
#' - Minimally redundant
#' - Well-suited for GSEA and pathway analysis
#'
#' @source MSigDB v7.5.1 (Broad Institute)
#' @references Liberzon et al. (2015) Cell Systems 1:417-425
#'
#' @examples
#' data(hallmark_pathways)
#' length(hallmark_pathways)
#' names(hallmark_pathways)[1:5]
#' head(hallmark_pathways$HALLMARK_GLYCOLYSIS)
#'
#' # Use for pathway analysis
#' data(example_counts, example_metadata)
#' pathway_results <- analyze_pathways(
#'   counts = example_counts,
#'   metadata = example_metadata,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   pathways = hallmark_pathways,
#'   scoring_method = "GSVA"
#' )
#'
"hallmark_pathways"