test_that("analyze_genes works with example data", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  
  # Test basic functionality
  expect_no_error({
    result <- analyze_genes(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      method = "DESeq2"
    )
  })
  
  # Test result structure
  result <- analyze_genes(
    counts = example_counts,
    metadata = example_metadata,
    design = ~ condition,
    contrast = c("condition", "treated", "control"),
    method = "DESeq2"
  )
  
  expect_s3_class(result, "functionalDE_result")
  expect_true("results" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true(is.data.frame(result$results))
  expect_true(nrow(result$results) > 0)
})

test_that("analyze_genes validates inputs correctly", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  
  # Test invalid counts
  expect_error(
    analyze_genes(
      counts = "not_a_matrix",
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control")
    )
  )
  
  # Test mismatched samples
  bad_metadata <- example_metadata[1:10, ]
  expect_error(
    analyze_genes(
      counts = example_counts,
      metadata = bad_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control")
    )
  )
  
  # Test invalid contrast
  expect_error(
    analyze_genes(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "wrong_level", "control")
    )
  )
})

test_that("analyze_genes works with different methods", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  
  # Test edgeR method
  expect_no_error({
    result_edger <- analyze_genes(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      method = "edgeR"
    )
  })
  
  # Test limma method
  expect_no_error({
    result_limma <- analyze_genes(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      method = "limma"
    )
  })
})