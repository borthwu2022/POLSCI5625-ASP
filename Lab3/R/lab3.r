######################################################
#                                                    #
# Lab 3: Bootstrap inference                         #
#                                                    #
# R/lab3.r                                           #
#                                                    #
# Bootstrap the mean of X and the OLS slope using:   #
#   1. an explicit loop; and                         #
#   2. boot::boot().                                 #
#                                                    #
######################################################

validate_lab3_data <- function(data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  required <- c("X", "Y")
  missing_columns <- setdiff(required, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "`data` is missing required column(s): ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (nrow(data) < 2L) {
    stop("`data` must contain at least two observations.", call. = FALSE)
  }

  valid_values <-
    is.numeric(data$X) &&
    is.numeric(data$Y) &&
    all(is.finite(data$X)) &&
    all(is.finite(data$Y))

  if (!valid_values) {
    stop("Columns `X` and `Y` must contain only finite numeric values.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_bootstrap_options <- function(B, conf) {
  valid_B <-
    length(B) == 1L &&
    is.numeric(B) &&
    is.finite(B) &&
    B >= 2L &&
    B == floor(B)

  if (!valid_B) {
    stop("`B` must be a whole number of at least 2.", call. = FALSE)
  }

  valid_conf <-
    length(conf) == 1L &&
    is.numeric(conf) &&
    is.finite(conf) &&
    conf > 0 &&
    conf < 1

  if (!valid_conf) {
    stop("`conf` must be a number strictly between 0 and 1.", call. = FALSE)
  }

  invisible(TRUE)
}

# These statistics use the same (data, indices) interface expected by boot().

mean_x_statistic <- function(data, indices) {
  mean(data$X[indices])
}

slope_statistic <- function(data, indices) {
  resampled_data <- data[indices, , drop = FALSE]
  unname(coef(lm(Y ~ X, data = resampled_data))[["X"]])
}

# Nonparametric case bootstrap written as an explicit loop.

bootstrap_loop <- function(data, statistic, B = 1000L, seed = 123L) {
  validate_lab3_data(data)
  validate_bootstrap_options(B, conf = 0.95)

  if (!is.function(statistic)) {
    stop("`statistic` must be a function.", call. = FALSE)
  }

  set.seed(seed)

  n <- nrow(data)
  bootstrap_estimates <- numeric(B)

  for (b in seq_len(B)) {
    indices <- sample.int(n, size = n, replace = TRUE)
    bootstrap_estimates[b] <- statistic(data, indices)
  }

  if (!all(is.finite(bootstrap_estimates))) {
    stop("A bootstrap replication returned a non-finite estimate.", call. = FALSE)
  }

  structure(
    list(
      t0 = statistic(data, seq_len(n)),
      t = bootstrap_estimates,
      R = as.integer(B),
      seed = seed
    ),
    class = "lab3_loop_boot"
  )
}

# Summarize loop replicates using the formulas in the Week 3 slides.
# The bias-corrected interval is also called the basic or reverse-percentile
# interval in the boot package.

summarize_loop_bootstrap <- function(result, statistic_name, conf = 0.95) {
  validate_bootstrap_options(result$R, conf)

  alpha <- 1 - conf
  percentile_interval <- unname(
    quantile(
      result$t,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )
  )

  basic_interval <- c(
    2 * result$t0 - percentile_interval[2],
    2 * result$t0 - percentile_interval[1]
  )

  bootstrap_mean <- mean(result$t)

  data.frame(
    statistic = statistic_name,
    method = "explicit loop",
    estimate = result$t0,
    bootstrap_mean = bootstrap_mean,
    bias = bootstrap_mean - result$t0,
    bias_corrected_estimate = 2 * result$t0 - bootstrap_mean,
    std_error = sd(result$t),
    percentile_lower = percentile_interval[1],
    percentile_upper = percentile_interval[2],
    bias_corrected_lower = basic_interval[1],
    bias_corrected_upper = basic_interval[2],
    row.names = NULL,
    check.names = FALSE
  )
}

# Run the same nonparametric case bootstrap with the boot package.

bootstrap_with_boot <- function(data, statistic, B = 1000L, seed = 123L) {
  validate_lab3_data(data)
  validate_bootstrap_options(B, conf = 0.95)

  if (!is.function(statistic)) {
    stop("`statistic` must be a function.", call. = FALSE)
  }

  if (!requireNamespace("boot", quietly = TRUE)) {
    stop(
      "Package `boot` is required. Install it with install.packages('boot').",
      call. = FALSE
    )
  }

  set.seed(seed)

  boot::boot(
    data = data,
    statistic = statistic,
    R = as.integer(B),
    sim = "ordinary",
    stype = "i"
  )
}

summarize_boot_package <- function(result, statistic_name, conf = 0.95) {
  validate_bootstrap_options(result$R, conf)

  intervals <- boot::boot.ci(
    result,
    conf = conf,
    type = c("basic", "perc")
  )

  bootstrap_estimates <- as.numeric(result$t[, 1L])
  estimate <- as.numeric(result$t0[[1L]])
  bootstrap_mean <- mean(bootstrap_estimates)

  data.frame(
    statistic = statistic_name,
    method = "boot::boot",
    estimate = estimate,
    bootstrap_mean = bootstrap_mean,
    bias = bootstrap_mean - estimate,
    bias_corrected_estimate = 2 * estimate - bootstrap_mean,
    std_error = sd(bootstrap_estimates),
    percentile_lower = intervals$percent[1L, 4L],
    percentile_upper = intervals$percent[1L, 5L],
    bias_corrected_lower = intervals$basic[1L, 4L],
    bias_corrected_upper = intervals$basic[1L, 5L],
    row.names = NULL,
    check.names = FALSE
  )
}
