######################################################
#                                                    #
# Lab 2 assignment: Instructor solution              #
#                                                    #
# Solutions for Tasks 3-5.                           #
# Work from the Lab2/assignment folder.              #
#                                                    #
######################################################

if (!exists("validate_pairwise_inputs", mode = "function") ||
    !exists("validate_threshold", mode = "function")) {
  helpers_path <- file.path("R", "assignment_helpers.r")

  if (!file.exists(helpers_path)) {
    stop(
      paste0(
        "Cannot find '", helpers_path, "'. ",
        "Set the working directory to the Lab2/assignment folder and try again."
      ),
      call. = FALSE
    )
  }

  source(helpers_path)
}


# Task 3: Thresholded pairwise distances with nested loops

thresholded_distance_loop <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  distances <- matrix(
    0,
    nrow = length(x),
    ncol = length(y)
  )

  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      distance <- abs(x[i] - y[j])

      if (distance < threshold) {
        distances[i, j] <- 0
      } else {
        distances[i, j] <- distance
      }
    }
  }

  distances
}


# Task 4: Thresholded pairwise distances with vectorized code using outer()

thresholded_distance_vectorized <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  distances <- abs(outer(x, y, FUN = "-"))
  distances[distances < threshold] <- 0

  distances
}


# Task 4 alternative: Hand-coded vectorized version without outer()

thresholded_distance_vectorized_handcoded <- function(
    x,
    y,
    threshold = 0.5
) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  x_matrix <- matrix(
    x,
    nrow = length(x),
    ncol = length(y)
  )

  distances <- abs(
    sweep(x_matrix, MARGIN = 2L, STATS = y, FUN = "-")
  )

  distances[distances < threshold] <- 0

  distances
}


# Task 5: Thresholded vectorized distances with sparse storage
#
# The final object is sparse, but Task 4 still creates a dense intermediate.

thresholded_distance_sparse <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "The `Matrix` package is required for the sparse implementation.",
      call. = FALSE
    )
  }

  distances <- thresholded_distance_vectorized(x, y, threshold)

  Matrix::Matrix(distances, sparse = TRUE)
}
