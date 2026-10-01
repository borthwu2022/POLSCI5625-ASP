######################################################
#                                                    #
# Lab 6: Gaussian mixture EM                         #
#                                                    #
# tests/lab6_test.r                                  #
#                                                    #
# Check simulated data and the hand-coded EM updates. #
# Test likelihood histories and final estimates.     #
#                                                    #
######################################################

source(file.path("R", "lab6.r"))

# Check the saved datasets against the specified counts and seeds.
sample_sizes <- c(100L, 1000L, 10000L)
data_paths <- file.path("data-raw", paste0("mixture_n", sample_sizes, ".csv"))

if (!all(file.exists(data_paths))) {
  stop("Missing simulated data. Run source(file.path('R', 'gendata.r')) first.",
       call. = FALSE)
}

for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  data <- read.csv(data_paths[i])

  set.seed(123L)
  true_class <- rep(1:3, times = c(n / 4, n / 2, n / 4))
  expected_x <- rnorm(n, mean = c(-3, 0, 4)[true_class], sd = 1)
  rows <- sample.int(n)

  stopifnot(
    nrow(data) == n,
    identical(data$id, seq_len(n)),
    identical(as.integer(table(data$true_class)), as.integer(c(n / 4, n / 2, n / 4))),
    identical(data$true_class, true_class[rows]),
    isTRUE(all.equal(data$x, expected_x[rows], tolerance = 1e-12))
  )
}


cat("PASS: saved data have exact class counts and reproduce the specified seeds.\n")

# Independently evaluate Bayes' rule and the mixture density on a small sample.
x <- c(-2, 0, 3)
p <- c(0.2, 0.5, 0.3)
mu <- c(-3, 0, 4)
direct <- sapply(1:3, function(k) p[k] * dnorm(x, mu[k], sqrt(2)))
e <- e_step(x, p, mu, 2)

stopifnot(
  isTRUE(all.equal(e$responsibilities, direct / rowSums(direct))),
  isTRUE(all.equal(e$loglik, sum(log(rowSums(direct)))))
)
# Ordinary density calculation underflows here; log-sum-exp stays finite.
extreme <- e_step(c(-1000, 1000), p, mu, 1)

stopifnot(
  is.finite(extreme$loglik),
  all(abs(rowSums(extreme$responsibilities) - 1) < 1e-12)
)

cat("PASS: Bayes responsibilities, observed likelihood, and log-scale stability.\n")

# Hand-calculated M-step: pi=(3/8,1/4,3/8), mu=(-1/3,0,1/3), var=11/12.
w <- rbind(c(0.5, 0.25, 0.25), c(0.25, 0.25, 0.5))
m <- m_step(c(-1, 1), w)

stopifnot(
  isTRUE(all.equal(m$proportions, c(3 / 8, 1 / 4, 3 / 8))),
  isTRUE(all.equal(m$means, c(-1 / 3, 0, 1 / 3))),
  isTRUE(all.equal(m$variance, 11 / 12))
)

cat("PASS: pooled variance uses updated means and denominator n.\n")

# Check all 15 requested fits and their returned likelihood/responsibilities.
# Test the five initializations specified in R/run_lab6.r.
starts <- data.frame(
  start = 1:5,
  pi1 = c(1 / 3, 0.2, 1 / 3, 0.5, 0.2),
  pi2 = c(1 / 3, 0.3, 1 / 3, 0.3, 0.5),
  pi3 = c(1 / 3, 0.5, 1 / 3, 0.2, 0.3),
  mu1 = c(-4, -2, -0.2, -5, 5),
  mu2 = c(-0.5, 1, 0, -3, 1),
  mu3 = c(5, 3, 0.2, -1, -2),
  variance = c(2, 4, 8, 2, 3)
)
for (i in 1:3) {
  x <- read.csv(data_paths[i])$x
  for (s in 1:5) {
    p <- as.numeric(starts[s, c("pi1", "pi2", "pi3")])
    mu <- as.numeric(starts[s, c("mu1", "mu2", "mu3")])
    fit <- em_mixture(x, p, mu, starts$variance[s])
    final <- e_step(x, fit$proportions, fit$means, fit$variance)
    stopifnot(
      fit$converged, all(diff(fit$history$loglik) >= -1e-8),
      nrow(fit$history) == fit$iter + 1L, fit$history$iteration[1] == 0,
      isTRUE(all.equal(fit$history$loglik[1],
                       e_step(x, p, mu, starts$variance[s])$loglik)),
      isTRUE(all.equal(tail(fit$history$loglik, 1), final$loglik)),
      isTRUE(all.equal(fit$responsibilities, final$responsibilities)),
      all(diff(fit$means) >= 0), fit$variance > 0,
      abs(sum(fit$proportions) - 1) < 1e-12
    )
  }
}

cat("PASS: all 15 fits converge with nondecreasing observed likelihoods.\n")

# Permuting labels must leave the likelihood and aligned estimates unchanged.
x <- read.csv(data_paths[2])$x
p <- c(0.2, 0.5, 0.3)
mu <- c(-4, 0.5, 5)
a <- em_mixture(x, p, mu, 2)
b <- em_mixture(x, rev(p), rev(mu), 2)

stopifnot(
  isTRUE(all.equal(a$means, b$means, tolerance = 1e-7)),
  isTRUE(all.equal(a$proportions, b$proportions, tolerance = 1e-7)),
  isTRUE(all.equal(a$loglik, b$loglik))
)

# Independent local optimizer: verify that a converged EM fit is stationary.
# Softmax weights use class 3 as reference; exp() makes the variance positive.
objective <- function(theta) {
  logits <- c(theta[1:2], 0)
  weights <- exp(logits - max(logits))
  weights <- weights / sum(weights)
  density <- sapply(1:3, function(k) {
    weights[k] * dnorm(x, theta[k + 2L], sqrt(exp(theta[6])))
  })
  -sum(log(rowSums(density)))
}
initial <- c(log(a$proportions[1:2] / a$proportions[3]), a$means, log(a$variance))
opt <- optim(
  initial,
  objective,
  method = "BFGS",
  control = list(reltol = 1e-12, maxit = 1000L)
)

stopifnot(
  opt$convergence == 0L, abs(-opt$value - a$loglik) < 1e-4
)

cat("PASS: label invariance and agreement with independent local optimization.\n")

# Iteration limits and invalid/degenerate inputs are visible to the caller.
warned <- FALSE
short <- withCallingHandlers(
  em_mixture(x, p, mu, 2, iter.max = 1L),
  warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)

stopifnot(
  warned, !short$converged, short$iter == 1L,
  inherits(try(e_step(x, p, mu, 0), silent = TRUE), "try-error"),
  inherits(try(m_step(c(-1, 1), cbind(c(1, 1), 0, 0)),
               silent = TRUE), "try-error")
)

cat("PASS: iteration limit, invalid variance, and empty components.\n")
cat("All Lab 6 tests passed.\n")
