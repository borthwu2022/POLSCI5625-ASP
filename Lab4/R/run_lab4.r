######################################################
#                                                    #
# Lab 4 in-class exercise: Tasks 1-3                 #
#                                                    #
# R/run_lab4.r                                       #
#                                                    #
# Check and benchmark nested loops, vectorization,   #
# and sequential and parallel foreach.               #
#                                                    #
######################################################

functions_path <- file.path("R", "lab4.r")
data_path <- file.path("data-raw", "counties.csv")

source(functions_path)

library(doParallel)
library(microbenchmark)

n_counties <- 3109L
n_workers <- 2L

run_lab4 <- function() {
  counties <- read.csv(data_path)
  mat <- as.matrix(
    counties[seq_len(n_counties), c("latitude", "longitude")]
  )

  # Set up the parallel backend before timing.
  cl <- parallel::makeCluster(n_workers, type = "PSOCK")
  on.exit({
    parallel::stopCluster(cl)
    registerDoSEQ()
  }, add = TRUE)

  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
  registerDoParallel(cl)

  # Check that all four versions return the same distances.
  reference <- unname(as.matrix(dist(mat)))

  stopifnot(
    isTRUE(all.equal(distance_loop(mat), reference)),
    isTRUE(all.equal(distance_vectorized(mat), reference)),
    isTRUE(all.equal(distance_foreach(mat), reference)),
    isTRUE(all.equal(distance_foreach(mat, parallel = TRUE), reference))
  )

  # Compare elapsed time with five repetitions.
  benchmark <- microbenchmark(
    nested = distance_loop(mat),
    vectorized = distance_vectorized(mat),
    foreach_sequential = distance_foreach(mat),
    foreach_parallel = distance_foreach(mat, parallel = TRUE),
    times = 5L
  )

  measured <- summary(benchmark, unit = "s")

  timing <- data.frame(
    n = n_counties,
    method = as.character(measured$expr),
    median_seconds = measured$median
  )

  elapsed <- setNames(timing$median_seconds, timing$method)

  timing$gain_over_nested <-
    elapsed["nested"] / timing$median_seconds

  # Save the timing comparison.
  dir.create("results", showWarnings = FALSE)

  write.csv(
    timing,
    file = file.path("results", "timing.csv"),
    row.names = FALSE
  )

  cat("Tasks 1-3 timing:\n")
  print(timing, row.names = FALSE)

  # Calculate vectorization gain, parallel speedup, and efficiency.
  vectorization_gain <- elapsed["nested"] / elapsed["vectorized"]
  parallel_speedup <-
    elapsed["foreach_sequential"] / elapsed["foreach_parallel"]
  parallel_efficiency <- parallel_speedup / n_workers

  cat("\nVectorization gain:", vectorization_gain, "\n")
  cat("Parallel speedup:", parallel_speedup, "\n")
  cat("Parallel efficiency:", parallel_efficiency, "\n")
}

run_lab4()
