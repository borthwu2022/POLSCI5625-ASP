######################################################
#                                                    #
# Lab 2 assignment: Tasks 3-5 tests                  #
#                                                    #
# tests/lab2_assignment_test.r                       #
#                                                    #
######################################################

solution_path <- file.path("R", "lab2_assignment_solution.r")

if (!file.exists(solution_path)) {
  stop(
    paste0(
      "Cannot find '", solution_path, "'. ",
      "Set the working directory to the Lab2/assignment folder and try again."
    ),
    call. = FALSE
  )
}

source(solution_path)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop(
    "The Matrix package is required. Install it with install.packages('Matrix').",
    call. = FALSE
  )
}

x <- c(0, 1)
y <- c(0, 0.5, 2)
threshold <- 0.5

expected <- matrix(
  c(0, 0.5, 2, 1, 0.5, 1),
  nrow = 2L,
  byrow = TRUE
)

d_loop <- thresholded_distance_loop(x, y, threshold)
d_outer <- thresholded_distance_vectorized(x, y, threshold)
d_handcoded <- thresholded_distance_vectorized_handcoded(x, y, threshold)
d_sparse <- thresholded_distance_sparse(x, y, threshold)

stopifnot(
  identical(dim(d_loop), c(2L, 3L)),
  isTRUE(all.equal(d_loop, expected, tolerance = 1e-12)),
  isTRUE(all.equal(d_loop, d_outer, tolerance = 1e-12)),
  isTRUE(all.equal(d_loop, d_handcoded, tolerance = 1e-12)),
  isTRUE(all.equal(d_loop, as.matrix(d_sparse), tolerance = 1e-12)),
  d_outer[1, 2] == threshold,
  Matrix::nnzero(d_sparse) == sum(d_outer != 0)
)

cat("PASS: all thresholded-distance implementations agree.\n")

invalid_threshold_error <- tryCatch(
  {
    thresholded_distance_vectorized(x, y, threshold = -1)
    NA_character_
  },
  error = function(e) conditionMessage(e)
)

stopifnot(
  !is.na(invalid_threshold_error),
  grepl("non-negative finite number", invalid_threshold_error, fixed = TRUE)
)

cat("PASS: an invalid threshold returns an informative error.\n")
cat("All Tasks 3-5 assignment tests passed.\n")
