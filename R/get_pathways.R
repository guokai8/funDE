#' Retrieve Pathway Gene Sets
#'
#' @description
#' Downloads pathway annotations from major databases using msigdbr and 
#' custom sources. Returns gene sets for use with pathway analysis functions.
#'
#' @param database Character specifying pathway database:
#'   - "Hallmark": MSigDB Hallmark gene sets (50 well-defined processes)
#'   - "KEGG": KEGG pathways via msigdbr
#'   - "KEGG_latest": Latest KEGG pathways via KEGG REST API (keggLink)
#'   - "Reactome": Reactome pathways via msigdbr
#'   - "GO_BP": Gene Ontology Biological Process
#'   - "GO_MF": Gene Ontology Molecular Function  
#'   - "GO_CC": Gene Ontology Cellular Component
#'   - "WikiPathways": WikiPathways via msigdbr
#'   - "BioCarta": BioCarta pathways via msigdbr
#'   - "PID": Pathway Interaction Database via msigdbr
#' @param species Character, "human" or "mouse"
#' @param min_size Minimum pathway size (default: 5)
#' @param max_size Maximum pathway size (default: 500)
#' @param id_type Type of gene identifiers: "symbol" (default), "ensembl", "entrez"
#' @param update Logical, download latest version (default: FALSE, uses cached)
#'
#' @return Named list where:
#'   - Names are pathway IDs/names
#'   - Values are character vectors of gene identifiers
#'
#' @details
#' This function uses the msigdbr package to access MSigDB collections.
#' For databases not in msigdbr, it provides custom download functionality.
#'
#' **Database Characteristics:**
#'
#' *Hallmark (H):*
#' - 50 gene sets representing well-defined biological states
#' - Best for initial exploration
#' - Coherent and interpretable
#'
#' *KEGG:*
#' - Metabolic and signaling pathways
#' - Well-curated, widely used
#' - ~300 pathways
#'
#' *Reactome:*
#' - Most detailed pathway database
#' - Hierarchical structure
#' - ~2000 pathways
#'
#' *GO:*
#' - Broadest coverage
#' - Can be redundant
#' - Consider filtering by evidence codes
#'
#' @examples
#' # Get Hallmark pathways (recommended starting point)
#' hallmark <- get_pathways("Hallmark", species = "human")
#' length(hallmark)  # 50 pathways
#'
#' # Get KEGG pathways
#' kegg <- get_pathways("KEGG", species = "mouse")
#'
#' # Get GO Biological Process with size filtering
#' go_bp <- get_pathways(
#'   "GO_BP", 
#'   species = "human",
#'   min_size = 10,
#'   max_size = 200
#' )
#'
#' @export
get_pathways <- function(database = "Hallmark",
                        species = "human",
                        min_size = 5,
                        max_size = 500,
                        id_type = "symbol",
                        update = FALSE) {
  
  # Check if msigdbr is available
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required but not installed. ",
         "Install it with: install.packages('msigdbr')")
  }
  
  # Map species to msigdbr format
  species_map <- c(
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus"
  )
  
  if (!species %in% names(species_map)) {
    stop("species must be 'human' or 'mouse'")
  }
  
  msigdb_species <- species_map[species]
  
  # Map database to msigdbr categories
  db_map <- list(
    "Hallmark" = list(category = "H", subcategory = NULL),
    "KEGG" = list(category = "C2", subcategory = "CP:KEGG"),
    "Reactome" = list(category = "C2", subcategory = "CP:REACTOME"),
    "WikiPathways" = list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
    "BioCarta" = list(category = "C2", subcategory = "CP:BIOCARTA"),
    "PID" = list(category = "C2", subcategory = "CP:PID"),
    "GO_BP" = list(category = "C5", subcategory = "GO:BP"),
    "GO_MF" = list(category = "C5", subcategory = "GO:MF"),
    "GO_CC" = list(category = "C5", subcategory = "GO:CC")
  )
  
  # Special handling for KEGG_latest using KEGG REST API
  if (database == "KEGG_latest") {
    return(get_kegg_latest(species, min_size, max_size, id_type))
  }
  
  if (!database %in% names(db_map)) {
    stop("Unsupported database. Choose from: ", 
         paste(c(names(db_map), "KEGG_latest"), collapse = ", "))
  }
  
  db_info <- db_map[[database]]
  
  # Download from msigdbr
  message("Downloading ", database, " pathways for ", species, "...")
  
  if (is.null(db_info$subcategory)) {
    msigdb_data <- msigdbr::msigdbr(
      species = msigdb_species,
      category = db_info$category
    )
  } else {
    msigdb_data <- msigdbr::msigdbr(
      species = msigdb_species,
      category = db_info$category,
      subcategory = db_info$subcategory
    )
  }
  
  # Convert to gene list format based on id_type
  gene_col <- switch(id_type,
                    "symbol" = "gene_symbol",
                    "ensembl" = "ensembl_gene",
                    "entrez" = "entrez_gene",
                    stop("id_type must be 'symbol', 'ensembl', or 'entrez'"))
  
  # Convert to named list
  pathways <- msigdb_data %>%
    select(gs_name, !!sym(gene_col)) %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(!!sym(gene_col))), .groups = "drop") %>%
    deframe()
  
  # Convert list column to character vectors
  pathways <- lapply(pathways, function(x) as.character(unlist(x)))
  
  # Filter by size
  pathway_sizes <- sapply(pathways, length)
  pathways <- pathways[pathway_sizes >= min_size & pathway_sizes <= max_size]
  
  message("Retrieved ", length(pathways), " pathways")
  message("Size range: ", min(sapply(pathways, length)), " - ", 
          max(sapply(pathways, length)), " genes")
  
  # Add attributes
  attr(pathways, "database") <- database
  attr(pathways, "species") <- species
  attr(pathways, "id_type") <- id_type
  attr(pathways, "date_retrieved") <- Sys.Date()
  attr(pathways, "n_pathways") <- length(pathways)
  
  return(pathways)
}


#' List Available Pathway Databases
#'
#' @description
#' Shows all available pathway databases and their characteristics
#'
#' @export
list_pathway_databases <- function() {
  databases <- data.frame(
    Database = c("Hallmark", "KEGG", "KEGG_latest", "Reactome", "WikiPathways", 
                "BioCarta", "PID", "GO_BP", "GO_MF", "GO_CC"),
    N_Pathways = c("50", "~300", "~300", "~2000", "~600", 
                   "~300", "~200", "~7500", "~1000", "~800"),
    Focus = c("General biological states", 
             "Metabolism & signaling",
             "Latest KEGG pathways",
             "Detailed processes",
             "Community-curated",
             "Signaling pathways",
             "Pathway interactions",
             "Biological processes",
             "Molecular functions",
             "Cellular components"),
    Recommended_For = c("Initial exploration",
                       "Metabolism studies",
                       "Latest KEGG data",
                       "Detailed mechanistic",
                       "Disease pathways",
                       "Signaling (legacy)",
                       "Signaling (legacy)",
                       "Broad coverage",
                       "Enzymatic functions",
                       "Localization"),
    stringsAsFactors = FALSE
  )
  
  print(databases, row.names = FALSE)
  invisible(databases)
}

#' Get latest KEGG pathways using KEGG REST API
#' @keywords internal
get_kegg_latest <- function(species = "human", min_size = 5, max_size = 500, id_type = "symbol") {
  
  # Check if KEGGREST is available
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("Package 'KEGGREST' is required for KEGG_latest. ",
         "Install it with: BiocManager::install('KEGGREST')")
  }
  
  message("Downloading latest KEGG pathways for ", species, " using KEGG REST API...")
  
  # Map species to KEGG organism codes
  kegg_species_map <- c(
    "human" = "hsa",
    "mouse" = "mmu"
  )
  
  if (!species %in% names(kegg_species_map)) {
    stop("KEGG_latest only supports human and mouse")
  }
  
  kegg_org <- kegg_species_map[species]
  
  # Get all pathway IDs for the organism
  pathway_list <- KEGGREST::keggList("pathway", kegg_org)
  pathway_ids <- names(pathway_list)
  
  message("Found ", length(pathway_ids), " KEGG pathways")
  
  # Get genes for each pathway
  pathway_genes <- list()
  
  for (i in seq_along(pathway_ids)) {
    pathway_id <- pathway_ids[i]
    pathway_name <- pathway_list[pathway_id]
    
    if (i %% 10 == 0) {
      message("Processing pathway ", i, " of ", length(pathway_ids))
    }
    
    tryCatch({
      # Get genes linked to this pathway
      genes <- KEGGREST::keggLink("genes", pathway_id)
      
      if (length(genes) > 0) {
        # Extract gene symbols/IDs
        if (id_type == "symbol") {
          # Get gene symbols
          gene_info <- KEGGREST::keggList(genes)
          gene_symbols <- sub("^[^;]*; ", "", gene_info)  # Extract gene symbols
          gene_symbols <- sub(" \\[.*", "", gene_symbols)  # Remove descriptions
          gene_symbols <- unique(gene_symbols)
        } else {
          # Use KEGG gene IDs directly
          gene_symbols <- sub(paste0(kegg_org, ":"), "", genes)
        }
        
        # Filter by size
        if (length(gene_symbols) >= min_size && length(gene_symbols) <= max_size) {
          clean_name <- sub("^path:", "", pathway_id)
          clean_name <- sub(paste0("^", kegg_org), "", clean_name)
          pathway_genes[[clean_name]] <- gene_symbols
        }
      }
    }, error = function(e) {
      # Skip pathways that fail
      message("Skipping pathway ", pathway_id, " due to error")
    })
    
    # Add small delay to be respectful to KEGG servers
    Sys.sleep(0.1)
  }
  
  message("Successfully retrieved ", length(pathway_genes), " pathways")
  
  return(pathway_genes)
}