######################################################
#                                                    #
# Lab 4 in-class exercise: Tasks 1-3 tests           #
#                                                    #
# tests/lab4_test.r                                  #
#                                                    #
######################################################

functions_path <- file.path("R", "lab4.r")

source(functions_path)

library(doParallel)

run_checks <- function() {
  cl <- parallel::makeCluster(2L, type = "PSOCK")
  on.exit({
    parallel::stopCluster(cl)
    registerDoSEQ()
  }, add = TRUE)

  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
  registerDoParallel(cl)

  mat <- rbind(c(0, 0), c(3, 4), c(0, 0))

  expected <- matrix(
    c(0, 5, 0, 5, 0, 5, 0, 5, 0),
    nrow = 3L,
    ncol = 3L
  )

  d_loop <- distance_loop(mat)
  d_vectorized <- distance_vectorized(mat)
  d_sequential <- distance_foreach(mat)
  d_parallel <- distance_foreach(mat, parallel = TRUE)

  stopifnot(
    isTRUE(all.equal(d_loop, expected)),
    isTRUE(all.equal(d_vectorized, expected)),
    isTRUE(all.equal(d_sequential, expected)),
    isTRUE(all.equal(d_parallel, expected))
  )

  cat("PASS: nested loops, vectorization, and serial/parallel foreach agree.\n")
}

run_checks()
