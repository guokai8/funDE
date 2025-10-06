test_that("analyze_pathways works with example data", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  data(hallmark_pathways, package = "funDE")
  
  # Test basic functionality
  expect_no_error({
    result <- analyze_pathways(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      pathways = hallmark_pathways,
      scoring_method = "mean"
    )
  })
  
  # Test result structure
  result <- analyze_pathways(
    counts = example_counts,
    metadata = example_metadata,
    design = ~ condition,
    contrast = c("condition", "treated", "control"),
    pathways = hallmark_pathways,
    scoring_method = "mean"
  )
  
  expect_s3_class(result, "functionalDE_result")
  expect_true("results" %in% names(result))
  expect_true("pathway_scores" %in% names(result))
  expect_true("pathways_used" %in% names(result))
  expect_true(is.data.frame(result$results))
  expect_true(is.matrix(result$pathway_scores))
})

test_that("analyze_pathways validates inputs correctly", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  
  # Test invalid pathways
  expect_error(
    analyze_pathways(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      pathways = "not_a_list"
    )
  )
  
  # Test empty pathways
  empty_pathways <- list()
  expect_error(
    analyze_pathways(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      pathways = empty_pathways
    )
  )
})

test_that("analyze_pathways works with different scoring methods", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  data(hallmark_pathways, package = "funDE")
  
  # Test mean scoring
  expect_no_error({
    result_mean <- analyze_pathways(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      pathways = hallmark_pathways,
      scoring_method = "mean"
    )
  })
  
  # Test median scoring
  expect_no_error({
    result_median <- analyze_pathways(
      counts = example_counts,
      metadata = example_metadata,
      design = ~ condition,
      contrast = c("condition", "treated", "control"),
      pathways = hallmark_pathways,
      scoring_method = "median"
    )
  })
})