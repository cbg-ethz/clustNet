# Generate a test dataset
set.seed(123)
n_samples <- 100
n_variables <- 6
myData <- matrix(sample(0:1, n_samples * n_variables, replace = TRUE), ncol = n_variables)

# Test basic functionality
test_that("get_clusters returns a list with the expected elements", {
  result <- get_clusters(myData, k_clust = 2)
  expect_type(result, "list")
  expect_true("clustermembership" %in% names(result))
  expect_true("assignprogress" %in% names(result))
  expect_true("DAGs" %in% names(result))
  expect_true("probs" %in% names(result))
})

# Test basic functionality
test_that("get_clusters returns a list with the expected elements", {
  result <- get_clusters(myData, k_clust = 2, EMseeds = 1)
  expect_type(result, "list")
  expect_true("clustermembership" %in% names(result))
  expect_true("assignprogress" %in% names(result))
  expect_true("DAGs" %in% names(result))
  expect_true("probs" %in% names(result))
})


# Test that the returned object has the expected dimensions
test_that("get_clusters returns objects with the expected dimensions", {
  result <- get_clusters(myData, k_clust = 2)
  expect_equal(length(result$clustermembership), n_samples)
  expect_equal(length(result$DAGs), 2)
  expect_equal(dim(result$probs), c(n_samples, 2))
})

# Test that the returned object has the expected dimensions
test_that("get_clusters returns objects with the expected dimensions", {
  result <- get_clusters(myData, k_clust = 2, EMseeds = 1)
  expect_equal(length(result$clustermembership), n_samples)
  expect_equal(length(result$DAGs), 2)
  expect_equal(dim(result$probs), c(n_samples, 2))
})

# Test that the function can handle missing inputs
test_that("get_clusters handles missing inputs correctly", {
  expect_error(get_clusters(), "Need a categorical matrix as input to cluster.")
})

# Test that the function can handle non-integer inputs
test_that("get_clusters handles non-integer inputs correctly", {
  expect_error(get_clusters(matrix(c(0, 0.5, 1, 0), ncol = 2)), "All categorical variables need to be specified as integers. Binary variables can be 0 or 1.")
})
