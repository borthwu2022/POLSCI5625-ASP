######################################################
#                                                    #
# Lab 6: Gaussian mixture EM                         #
#                                                    #
# R/gendata.r                                        #
#                                                    #
# Generate three simulated Gaussian-mixture datasets.#
# Use exact class counts and a common variance of 1. #
#                                                    #
######################################################

sample_sizes <- c(100L, 1000L, 10000L)
seeds <- rep(123L, 3L)
means <- c(-3, 0, 4)
variance <- 1

dir.create("data-raw", showWarnings = FALSE)

for (i in seq_along(sample_sizes)) {
  set.seed(seeds[i])

  n <- sample_sizes[i]
  true_class <- rep(1:3, times = c(n / 4, n / 2, n / 4))

  # rnorm() takes a standard deviation rather than a variance.
  x <- rnorm(n, mean = means[true_class], sd = sqrt(variance))
  rows <- sample.int(n)

  lab6_data <- data.frame(
    id = seq_len(n),
    x = x[rows],
    true_class = true_class[rows]
  )

  data_path <- file.path("data-raw", paste0("mixture_n", n, ".csv"))

  write.csv(
    lab6_data,
    file = data_path,
    row.names = FALSE
  )

  cat("Wrote", nrow(lab6_data), "observations to", data_path, "\n")
}
