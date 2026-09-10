######################################################
#                                                    #
# Lab 3: Bootstrap inference                         #
#                                                    #
# R/gendata.r                                        #
#                                                    #
# Generate the simulated regression data.            #
#                                                    #
######################################################

set.seed(123)

N <- 500L
beta_0 <- 1.5
beta_1 <- -2

# The notation N(0, 2) is interpreted as mean 0 and variance 2.
X <- rnorm(N, mean = 0, sd = sqrt(2))
error <- rnorm(N, mean = 0, sd = 1)
Y <- beta_0 + beta_1 * X + error

lab3_data <- data.frame(
  id = seq_len(N),
  X = X,
  Y = Y
)

dir.create("data-raw", showWarnings = FALSE)

write.csv(
  lab3_data,
  file = file.path("data-raw", "lab3.csv"),
  row.names = FALSE
)

cat("Wrote", nrow(lab3_data), "observations to data-raw/lab3.csv.\n")
