test_that("test_enrichment works correctly", {
  data(hallmark_pathways, package = "funDE")
  
  # Create test data
  test_genes <- c("ALDOA", "ENO1", "GAPDH", "HK1", "LDHA", "PGK1", "PKM", "TPI1")
  universe_genes <- unique(unlist(hallmark_pathways))
  
  # Test basic functionality
  expect_no_error({
    enrichment_result <- test_enrichment(
      genes = test_genes,
      universe = universe_genes,
      gene_sets = hallmark_pathways,
      method = "hypergeometric"
    )
  })
  
  # Test result structure
  enrichment_result <- test_enrichment(
    genes = test_genes,
    universe = universe_genes,
    gene_sets = hallmark_pathways,
    method = "hypergeometric"
  )
  
  expect_s3_class(enrichment_result, "functionalDE_enrichment")
  expect_true(is.data.frame(enrichment_result))
  expect_true("gene_set" %in% colnames(enrichment_result))
  expect_true("pvalue" %in% colnames(enrichment_result))
  expect_true("padj" %in% colnames(enrichment_result))
  expect_true("overlap" %in% colnames(enrichment_result))
  expect_true("enrichment" %in% colnames(enrichment_result))
})

test_that("test_enrichment validates inputs correctly", {
  data(hallmark_pathways, package = "funDE")
  
  # Test invalid genes
  expect_error(
    test_enrichment(
      genes = 123,
      universe = c("A", "B", "C"),
      gene_sets = hallmark_pathways
    )
  )
  
  # Test invalid universe
  expect_error(
    test_enrichment(
      genes = c("A", "B"),
      universe = 123,
      gene_sets = hallmark_pathways
    )
  )
  
  # Test invalid gene_sets
  expect_error(
    test_enrichment(
      genes = c("A", "B"),
      universe = c("A", "B", "C"),
      gene_sets = "not_a_list"
    )
  )
})

test_that("test_enrichment works with different methods", {
  data(hallmark_pathways, package = "funDE")
  
  test_genes <- c("ALDOA", "ENO1", "GAPDH", "HK1", "LDHA")
  universe_genes <- unique(unlist(hallmark_pathways))
  
  # Test Fisher's exact test
  expect_no_error({
    enrichment_fisher <- test_enrichment(
      genes = test_genes,
      universe = universe_genes,
      gene_sets = hallmark_pathways,
      method = "fisher"
    )
  })
  
  # Test hypergeometric test
  expect_no_error({
    enrichment_hyper <- test_enrichment(
      genes = test_genes,
      universe = universe_genes,
      gene_sets = hallmark_pathways,
      method = "hypergeometric"
    )
  })
})