######################################################
#                                                    #
# Lab 6: Gaussian mixture EM                         #
#                                                    #
# R/lab6.r                                           #
#                                                    #
# In-class reference implementations for:            #
#   Task 1. E-step and observed-data log-likelihood  #
#   Task 2. M-step with a shared variance            #
#   Task 3. repeat EM until convergence              #
#                                                    #
######################################################

# x is a numeric vector of observations.
# Each component has its own proportion and mean, but shares one variance.


# Task 1. E-step and observed-data log-likelihood

e_step <- function(x, proportions, means, variance) {
  stopifnot(
    is.numeric(x), length(x) > 0L, all(is.finite(x)),
    length(proportions) == 3L, all(is.finite(proportions)),
    all(proportions > 0), abs(sum(proportions) - 1) < 1e-10,
    length(means) == 3L, all(is.finite(means)),
    length(variance) == 1L, is.finite(variance), variance > 0
  )

  log_joint <- matrix(0, nrow = length(x), ncol = 3L)

  for (k in 1:3) {
    log_joint[, k] <- log(proportions[k]) +
                      dnorm(x, mean = means[k], sd = sqrt(variance), log = TRUE)
  }

  row_max <- apply(log_joint, 1, max)

  if (any(!is.finite(row_max))) {
    stop("Non-finite log densities: choose less extreme starting values.",
         call. = FALSE)
  }

  # Subtract the row maximum before exponentiating to avoid underflow.
  # R recycles the n-vector down each column of the n by 3 matrix.
  shifted <- exp(log_joint - row_max)
  denominator <- rowSums(shifted)
  log_density <- row_max + log(denominator)

  list(
    responsibilities = shifted / denominator,
    loglik = sum(log_density)
  )
}


# Task 2. Update proportions, means, and the common variance

m_step <- function(x, responsibilities) {
  stopifnot(
    is.matrix(responsibilities),
    identical(dim(responsibilities), c(length(x), 3L)),
    all(is.finite(responsibilities)), all(responsibilities >= 0),
    all(abs(rowSums(responsibilities) - 1) < 1e-10)
  )

  n <- length(x)
  effective_n <- colSums(responsibilities)

  if (any(effective_n <= 0)) {
    stop("A component has zero effective membership; choose another start.",
         call. = FALSE)
  }

  proportions <- effective_n / n
  means <- colSums(responsibilities * x) / effective_n
  residuals <- outer(x, means, "-")

  # Keep responsibilities fixed. Pool across components using the new means.
  # The maximum-likelihood variance uses n, not n - 1, as the denominator.
  variance <- sum(responsibilities * residuals^2) / n

  if (!is.finite(variance) || variance <= 0) {
    stop("Degenerate shared variance; choose another start or inspect the data.",
         call. = FALSE)
  }

  list(
    proportions = proportions,
    means = means,
    variance = variance
  )
}


# Task 3. Repeat until the observed-data log-likelihood converges

em_mixture <- function(x, proportions, means, variance,
                       tol = 1e-10, iter.max = 2000L,
                       monotone_tol = 1e-8) {
  stopifnot(
    length(tol) == 1L, is.finite(tol), tol > 0,
    length(iter.max) == 1L, is.finite(iter.max), iter.max >= 1L,
    iter.max == floor(iter.max), length(monotone_tol) == 1L,
    is.finite(monotone_tol), monotone_tol >= 0
  )

  current <- e_step(x, proportions, means, variance)
  loglik <- numeric(iter.max + 1L)
  loglik[1] <- current$loglik  # Iteration 0 is the starting likelihood.
  converged <- FALSE

  for (iteration in seq_len(iter.max)) {
    updated <- m_step(x, current$responsibilities)
    next_step <- e_step(
      x,
      proportions = updated$proportions,
      means = updated$means,
      variance = updated$variance
    )

    loglik[iteration + 1L] <- next_step$loglik
    gain <- next_step$loglik - current$loglik

    if (gain < -monotone_tol) {
      stop("Observed-data log-likelihood decreased beyond rounding tolerance.",
           call. = FALSE)
    }

    converged <- abs(gain) <= tol * (1 + abs(current$loglik))
    current <- next_step

    if (converged) {
      break
    }
  }

  if (!converged) {
    warning("Iteration limit reached; increase iter.max.", call. = FALSE)
  }

  # Component labels are arbitrary. Align all outputs by increasing mean.
  component_order <- order(updated$means)

  list(
    proportions = updated$proportions[component_order],
    means = updated$means[component_order],
    variance = updated$variance,
    responsibilities = current$responsibilities[, component_order, drop = FALSE],
    loglik = current$loglik,
    iter = iteration,
    converged = converged,
    history = data.frame(
      iteration = 0:iteration,
      loglik = loglik[seq_len(iteration + 1L)]
    )
  )
}
