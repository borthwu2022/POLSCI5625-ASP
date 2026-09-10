######################################################
#                                                    #
# Lab 2 assignment: Tasks 3-5 benchmark              #
#                                                    #
# R/run_assignment.r                                 #
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

data_path <- file.path("data-raw", "lab2_vectors.rds")

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Cannot find '", data_path, "'. ",
      "Run source(file.path('R', 'gendata.r')) first."
    ),
    call. = FALSE
  )
}

lab2_vectors <- readRDS(data_path)
x <- lab2_vectors$x
y <- lab2_vectors$y
threshold <- 0.5

if (length(x) != 3000L || length(y) != 5000L) {
  stop("The assignment vectors must have lengths 3000 and 5000.",
       call. = FALSE)
}

# Check correctness before benchmarking.
gc()
d_loop <- thresholded_distance_loop(x, y, threshold)
d_dense <- thresholded_distance_vectorized(x, y, threshold)
d_handcoded <- thresholded_distance_vectorized_handcoded(x, y, threshold)
d_sparse <- thresholded_distance_sparse(x, y, threshold)

stopifnot(
  identical(dim(d_loop), c(3000L, 5000L)),
  isTRUE(all.equal(d_loop, d_dense, tolerance = 1e-12)),
  isTRUE(all.equal(d_dense, d_handcoded, tolerance = 1e-12)),
  isTRUE(all.equal(d_dense, as.matrix(d_sparse), tolerance = 1e-12)),
  all(d_dense == 0 | d_dense >= threshold),
  Matrix::nnzero(d_sparse) == sum(d_dense != 0)
)

stored_result_bytes <- c(
  as.numeric(object.size(d_loop)),
  as.numeric(object.size(d_dense)),
  as.numeric(object.size(d_sparse))
)

nonzero_entries <- c(
  sum(d_loop != 0),
  sum(d_dense != 0),
  Matrix::nnzero(d_sparse)
)

rm(d_loop, d_dense, d_handcoded, d_sparse)
gc()

loop_time <- elapsed_seconds(
  thresholded_distance_loop,
  x,
  y,
  threshold
)

gc()
vectorized_time <- elapsed_seconds(
  thresholded_distance_vectorized,
  x,
  y,
  threshold
)

gc()
sparse_time <- elapsed_seconds(
  thresholded_distance_sparse,
  x,
  y,
  threshold
)

benchmark <- data.frame(
  experiment = "3000 x 5000 thresholded distances",
  method = c("nested loop", "vectorized dense", "vectorized sparse"),
  elapsed_seconds = c(loop_time, vectorized_time, sparse_time),
  stored_result_bytes = stored_result_bytes,
  nonzero_entries = nonzero_entries,
  stringsAsFactors = FALSE
)

benchmark$time_gain_over_loop <-
  benchmark$elapsed_seconds[1] / benchmark$elapsed_seconds

benchmark$time_reduction_percent <-
  100 * (1 - benchmark$elapsed_seconds / benchmark$elapsed_seconds[1])

benchmark$memory_reduction_from_dense_percent <-
  100 * (1 - benchmark$stored_result_bytes /
           benchmark$stored_result_bytes[2])

dir.create("results", showWarnings = FALSE)

write.csv(
  benchmark,
  file = file.path("results", "large_benchmark.csv"),
  row.names = FALSE
)

cat("Tasks 3-5 timing and stored-result memory:\n")
print(benchmark)

cat(
  "\nTime gain, Task 3 to Task 4:",
  benchmark$time_gain_over_loop[2],
  "x\n"
)

cat(
  "Stored-result memory reduction, Task 4 to Task 5:",
  benchmark$memory_reduction_from_dense_percent[3],
  "%\n"
)
