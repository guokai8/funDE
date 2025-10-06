test_that("validation functions work correctly", {
  data(example_counts, package = "funDE")
  data(example_metadata, package = "funDE")
  
  # Test validate_counts
  expect_no_error(validate_counts(example_counts))
  
  # Test with data.frame input
  counts_df <- as.data.frame(example_counts)
  validated <- validate_counts(counts_df)
  expect_true(is.matrix(validated))
  
  # Test invalid counts
  expect_error(validate_counts("not_a_matrix"))
  expect_error(validate_counts(matrix(c(-1, 2, 3, 4), nrow = 2)))
  
  # Test validate_metadata
  expect_no_error(validate_metadata(example_metadata, example_counts))
  
  # Test validate_design
  design <- ~ condition
  expect_no_error(validate_design(design, example_metadata))
  
  # Test invalid design
  bad_design <- ~ nonexistent_column
  expect_error(validate_design(bad_design, example_metadata))
  
  # Test validate_contrast
  contrast <- c("condition", "treated", "control")
  expect_no_error(validate_contrast(contrast, example_metadata))
  
  # Test invalid contrast
  bad_contrast <- c("condition", "wrong_level", "control")
  expect_error(validate_contrast(bad_contrast, example_metadata))
})

test_that("statistical utility functions work correctly", {
  data(example_counts, package = "funDE")
  
  # Test filter_low_counts
  filtered <- filter_low_counts(example_counts, min_count = 10, min_samples = 2)
  expect_true(is.matrix(filtered))
  expect_true(nrow(filtered) <= nrow(example_counts))
  
  # Test calculate_size_factors
  size_factors <- calculate_size_factors(example_counts[1:100, ])
  expect_true(is.numeric(size_factors))
  expect_equal(length(size_factors), ncol(example_counts))
  expect_true(all(size_factors > 0))
  
  # Test normalize_counts
  normalized <- normalize_counts(example_counts[1:100, ])
  expect_true(is.matrix(normalized))
  expect_equal(dim(normalized), dim(example_counts[1:100, ]))
})

test_that("general utility functions work correctly", {
  # Test format_pvalues
  pvals <- c(0.0001, 0.01, 0.05, 0.1, NA)
  formatted <- format_pvalues(pvals)
  expect_true(is.character(formatted))
  expect_equal(length(formatted), length(pvals))
  
  # Test format_logfc
  lfc <- c(-2.5, -0.5, 0, 0.5, 2.5, NA)
  formatted_lfc <- format_logfc(lfc)
  expect_true(is.character(formatted_lfc))
  expect_equal(length(formatted_lfc), length(lfc))
  
  # Test get_top_n
  test_vector <- c(a = 5, b = 2, c = 8, d = 1, e = 6)
  top_3 <- get_top_n(test_vector, n = 3)
  expect_equal(length(top_3), 3)
  expect_equal(names(top_3)[1], "c")  # highest value
})