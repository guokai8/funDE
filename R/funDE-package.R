#' funDE: Functional Differential Expression Analysis
#'
#' @description
#' funDE provides differential expression analysis at multiple functional
#' levels: genes, gene families, and pathways. It supports both bulk and 
#' single-cell RNA-seq data with a unified, easy-to-use interface.
#'
#' @section Main Analysis Functions:
#' 
#' **Core DE Analysis:**
#' \itemize{
#'   \item \code{\link{analyze_genes}}: Gene-level differential expression
#'   \item \code{\link{analyze_families}}: Gene family-level analysis
#'   \item \code{\link{analyze_pathways}}: Pathway-level analysis
#'   \item \code{\link{analyze_regions}}: Transcript/isoform-level analysis
#' }
#'
#' @section Data Preparation:
#' \itemize{
#'   \item \code{\link{aggregate_pseudobulk}}: Single-cell to pseudobulk aggregation
#'   \item \code{\link{score_pathway_activity}}: Calculate pathway scores
#' }
#'
#' @section Enrichment & Annotation:
#' \itemize{
#'   \item \code{\link{test_enrichment}}: Gene set enrichment testing
#'   \item \code{\link{get_gene_families}}: Retrieve gene family annotations
#'   \item \code{\link{get_pathways}}: Access pathway databases
#' }
#'
#' @section Visualization:
#' \itemize{
#'   \item \code{\link{plot_multilevel}}: Compare across analytical levels
#'   \item \code{\link{plot_pathway_heatmap}}: Pathway activity heatmaps
#'   \item \code{\link{plot_enrichment}}: Enrichment bar/dot plots
#' }
#'
#' @section Integration & Utilities:
#' \itemize{
#'   \item \code{\link{integrate_levels}}: Combine multi-level results
#'   \item \code{\link{compare_methods}}: Compare pathway scoring methods
#'   \item \code{\link{export_results}}: Export to multiple formats
#' }
#'
#' @section Getting Started:
#' 
#' See \code{vignette("introduction", package = "funDE")} for a quick start guide.
#' 
#' **Basic workflow:**
#' \preformatted{
#' # Gene-level analysis
#' gene_res <- analyze_genes(counts, metadata, design = ~ condition,
#'                          contrast = c("condition", "treat", "ctrl"))
#' 
#' # Pathway analysis  
#' path_res <- analyze_pathways(counts, metadata, design = ~ condition,
#'                             contrast = c("condition", "treat", "ctrl"),
#'                             pathways = "Hallmark")
#' 
#' # Enrichment
#' enrichment <- test_enrichment(
#'   genes = rownames(gene_res$results)[gene_res$results$padj < 0.05],
#'   universe = rownames(gene_res$results),
#'   gene_sets = get_pathways("KEGG")
#' )
#' }
#'
#' @docType package
#' @name funDE_package
#' @aliases funDE
#'
#' @author Your Name \email{your.email@@example.com}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats model.matrix p.adjust cor sd median quantile
#' @importFrom stats fisher.test phyper dhyper var
#' @importFrom utils head tail write.csv write.table
#' @importFrom methods is as new
#' @import ggplot2
#' @import dplyr
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom GSVA gsva
#' @importFrom fgsea fgsea
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom Matrix Matrix colSums rowSums
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_wrap str_replace_all
## usethis namespace: end
NULL