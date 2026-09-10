######################################################
#                                                    #
# Lab 3: Bootstrap inference tests                   #
#                                                    #
# tests/lab3_test.r                                  #
#                                                    #
######################################################

functions_path <- file.path("R", "lab3.r")
data_path <- file.path("data-raw", "lab3.csv")

if (!file.exists(functions_path) || !file.exists(data_path)) {
  stop(
    paste0(
      "Run this script with Lab3 as the working directory, and run ",
      "R/gendata.r before the tests."
    ),
    call. = FALSE
  )
}

source(functions_path)
lab3_data <- read.csv(data_path)

stopifnot(
  nrow(lab3_data) == 500L,
  identical(names(lab3_data), c("id", "X", "Y")),
  all(is.finite(lab3_data$X)),
  all(is.finite(lab3_data$Y))
)

cat("PASS: the simulated data contain 500 complete observations.\n")

# Use a smaller number of replications for quick tests. The analysis script
# uses the required B = 1000.
test_B <- 199L

mean_loop_one <- bootstrap_loop(
  lab3_data,
  mean_x_statistic,
  B = test_B,
  seed = 123L
)

mean_loop_two <- bootstrap_loop(
  lab3_data,
  mean_x_statistic,
  B = test_B,
  seed = 123L
)

stopifnot(
  identical(mean_loop_one$t, mean_loop_two$t),
  length(mean_loop_one$t) == test_B,
  isTRUE(all.equal(mean_loop_one$t0, mean(lab3_data$X), tolerance = 1e-12))
)

cat("PASS: the loop bootstrap is reproducible and returns B replicates.\n")

mean_loop_summary <- summarize_loop_bootstrap(
  mean_loop_one,
  statistic_name = "mean(X)"
)

expected_percentile <- unname(
  quantile(mean_loop_one$t, c(0.025, 0.975), type = 7)
)
expected_basic <- c(
  2 * mean_loop_one$t0 - expected_percentile[2],
  2 * mean_loop_one$t0 - expected_percentile[1]
)

stopifnot(
  isTRUE(all.equal(
    mean_loop_summary$std_error,
    sd(mean_loop_one$t),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unlist(mean_loop_summary[c("percentile_lower", "percentile_upper")]),
    expected_percentile,
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    unlist(mean_loop_summary[c("bias_corrected_lower", "bias_corrected_upper")]),
    expected_basic,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))
)

cat("PASS: loop standard error and both intervals match their formulas.\n")

slope_loop <- bootstrap_loop(
  lab3_data,
  slope_statistic,
  B = test_B,
  seed = 123L
)

observed_slope <- unname(coef(lm(Y ~ X, data = lab3_data))[["X"]])

stopifnot(
  isTRUE(all.equal(slope_loop$t0, observed_slope, tolerance = 1e-12)),
  all(is.finite(slope_loop$t))
)

cat("PASS: the loop bootstrap resamples complete cases for the OLS slope.\n")

mean_boot <- bootstrap_with_boot(
  lab3_data,
  mean_x_statistic,
  B = test_B,
  seed = 123L
)

mean_boot_summary <- summarize_boot_package(
  mean_boot,
  statistic_name = "mean(X)"
)

direct_ci <- boot::boot.ci(mean_boot, type = c("basic", "perc"))

stopifnot(
  inherits(mean_boot, "boot"),
  nrow(mean_boot$t) == test_B,
  isTRUE(all.equal(mean_boot$t0[[1L]], mean(lab3_data$X), tolerance = 1e-12)),
  isTRUE(all.equal(
    mean_boot_summary$percentile_lower,
    unname(direct_ci$percent[1L, 4L]),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    mean_boot_summary$bias_corrected_lower,
    unname(direct_ci$basic[1L, 4L]),
    tolerance = 1e-12
  ))
)

cat("PASS: boot::boot() and boot::boot.ci() return the requested results.\n")

invalid_B_error <- tryCatch(
  {
    bootstrap_loop(lab3_data, mean_x_statistic, B = 1L)
    NA_character_
  },
  error = function(e) conditionMessage(e)
)

stopifnot(
  !is.na(invalid_B_error),
  grepl("at least 2", invalid_B_error, fixed = TRUE)
)

cat("PASS: invalid replication counts return an informative error.\n")
cat("All Lab 3 tests passed.\n")
