######################################################
#                                                    #
# Lab 3: Bootstrap inference                         #
#                                                    #
# R/run_lab3.r                                       #
#                                                    #
# Run B = 1000 bootstrap replications for the mean   #
# of X and the OLS slope, using a loop and boot().   #
#                                                    #
######################################################

functions_path <- file.path("R", "lab3.r")
data_path <- file.path("data-raw", "lab3.csv")

if (!file.exists(functions_path)) {
  stop(
    paste0(
      "Cannot find '", functions_path, "'. ",
      "Set the working directory to the Lab3 folder and try again."
    ),
    call. = FALSE
  )
}

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Cannot find '", data_path, "'. ",
      "Run source(file.path('R', 'gendata.r')) first."
    ),
    call. = FALSE
  )
}

source(functions_path)
lab3_data <- read.csv(data_path)

B <- 1000L
confidence_level <- 0.95
bootstrap_seed <- 123L

# Task 1: Bootstrap the sample mean of X with an explicit loop.
mean_loop <- bootstrap_loop(
  data = lab3_data,
  statistic = mean_x_statistic,
  B = B,
  seed = bootstrap_seed
)

# Task 2: Bootstrap the OLS slope on X with an explicit loop.
slope_loop <- bootstrap_loop(
  data = lab3_data,
  statistic = slope_statistic,
  B = B,
  seed = bootstrap_seed
)

# Task 3: Repeat both bootstraps with boot::boot().
mean_boot <- bootstrap_with_boot(
  data = lab3_data,
  statistic = mean_x_statistic,
  B = B,
  seed = bootstrap_seed
)

slope_boot <- bootstrap_with_boot(
  data = lab3_data,
  statistic = slope_statistic,
  B = B,
  seed = bootstrap_seed
)

summary_table <- rbind(
  summarize_loop_bootstrap(mean_loop, "mean(X)", confidence_level),
  summarize_loop_bootstrap(slope_loop, "slope on X", confidence_level),
  summarize_boot_package(mean_boot, "mean(X)", confidence_level),
  summarize_boot_package(slope_boot, "slope on X", confidence_level)
)

replicate_table <- data.frame(
  replication = seq_len(B),
  mean_loop = mean_loop$t,
  slope_loop = slope_loop$t,
  mean_boot = as.numeric(mean_boot$t[, 1L]),
  slope_boot = as.numeric(slope_boot$t[, 1L])
)

dir.create("results", showWarnings = FALSE)

write.csv(
  summary_table,
  file = file.path("results", "bootstrap_summary.csv"),
  row.names = FALSE
)

write.csv(
  replicate_table,
  file = file.path("results", "bootstrap_replicates.csv"),
  row.names = FALSE
)

cat("Lab 3 bootstrap summary (", B, " replications):\n", sep = "")
print(summary_table, digits = 5, row.names = FALSE)
cat("\nWrote results/bootstrap_summary.csv and results/bootstrap_replicates.csv.\n")
