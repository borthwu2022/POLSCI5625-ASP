######################################################
#                                                    #
# Lab 2 assignment: Student starter                  #
#                                                    #
# R/lab2_assignment.r                                #
# Complete Tasks 3-5.                                #
#                                                    #
######################################################

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


# Task 3: Thresholded pairwise distances with nested loops

thresholded_distance_loop <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  # TODO: preallocate a length(x) by length(y) dense matrix.
  # TODO: calculate abs(x[i] - y[j]) inside two nested loops.
  # TODO: store zero if the distance is strictly less than threshold.
  # TODO: otherwise preserve the distance.
}


# Task 4: Thresholded pairwise distances with vectorized code

thresholded_distance_vectorized <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  # TODO: calculate all pairwise distances without nested loops.
  # TODO: apply the strict threshold.
}


# Task 5: Thresholded vectorized distances with sparse storage

thresholded_distance_sparse <- function(x, y, threshold = 0.5) {
  validate_pairwise_inputs(x, y)
  validate_threshold(threshold)

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install the Matrix package before using this function.",
         call. = FALSE)
  }

  # TODO: call thresholded_distance_vectorized().
  # TODO: convert the result with Matrix::Matrix(..., sparse = TRUE).
}
