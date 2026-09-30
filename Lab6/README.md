# Lab 6: Expectation-Maximization for Gaussian Mixtures

## Overview

Implement EM for a three-component, univariate Gaussian mixture. Estimate three
class proportions, three means, and one shared variance. Track the observed-data
log-likelihood at every iteration, check that it does not decrease, and compare
five starting values at each of three sample sizes.

This lab follows Week 6's complete-data likelihood, class probabilities,
parameter updates, log-scale calculations, monotone convergence, and label
switching (displayed slides 2-5, 8-11, and 20). It extends Lab 5's hard K-means
assignments to soft probabilities. The slides' component-specific covariance
update is adapted here to **one pooled variance**.

## Learning objectives

1. Simulate normal data with exact component counts.
2. Calculate soft class probabilities using Bayes' rule and log-sum-exp.
3. Update proportions, means, and a common variance.
4. Store and plot the observed-data log-likelihood throughout EM.
5. Compare starting values after aligning component labels.

## Project structure

```text
Lab6/
|-- data-raw/
|   |-- mixture_n100.csv
|   |-- mixture_n1000.csv
|   |-- mixture_n10000.csv
|-- R/
|   |-- gendata.r
|   |-- lab6.r
|   |-- run_lab6.r
|-- results/
|   |-- starting_values.csv
|   |-- start_comparison.csv
|   |-- best_estimates.csv
|   |-- loglik_history.csv
|   |-- loglik_history.pdf
|   |-- mixture_fits.pdf
|-- tests/
|   |-- lab6_test.r
|-- README.md
```

## Requirements

R only; no additional packages. Paths are relative to `Lab6`.

## Simulation and model

For each $n \in \{100, 1000, 10000\}$, generate:

| Component | Observations | Generating distribution | True proportion |
|---|---:|---|---:|
| 1 | $n/4$ | $N(-3, 1)$ | $0.25$ |
| 2 | $n/2$ | $N(0, 1)$ | $0.50$ |
| 3 | $n/4$ | $N(4, 1)$ | $0.25$ |

**The generating variance is $\sigma^2 = 1$; it is unknown during estimation.** Estimate a
single shared variance rather than fixing it to 1 or estimating three separate
variances. `dnorm()` and `rnorm()` take a standard deviation, so use
`sd = sqrt(variance)`.

Use exactly the specified counts, then shuffle observations. Keep `true_class`
only to verify the simulation; fit using `x` alone. Class proportions are free
parameters even though simulation counts were fixed. We maximize the usual
mixture likelihood, treating memberships as latent and not imposing the known
simulation counts on the fitting procedure.

Let $Z_i \in \{1,2,3\}$ denote the unobserved component membership. The model is

$$
P(Z_i=k)=\pi_k,
\qquad
X_i \mid Z_i=k \sim N(\mu_k,\sigma^2),
\qquad
\sum_{k=1}^{3}\pi_k=1.
$$

Write the parameters as $\theta=(\boldsymbol{\pi},\boldsymbol{\mu},\sigma^2)$,
with $\pi_k>0$ and $\sigma^2>0$. The normal density is

$$
\phi(x;\mu,\sigma^2)
=
\frac{1}{\sqrt{2\pi\sigma^2}}
\exp\left[-\frac{(x-\mu)^2}{2\sigma^2}\right].
$$

The observed-data log-likelihood is

$$
\ell(\theta)
=
\sum_{i=1}^{n}
\log\left[
\sum_{k=1}^{3}\pi_k\,\phi(x_i;\mu_k,\sigma^2)
\right].
$$

The logarithm is outside the component sum. This is different from the
expected complete-data log-likelihood $Q$ that the M-step maximizes with fixed
responsibilities.

## In-class exercise

Write your functions before consulting the reference implementations. The
scripts follow the same separation of data generation, functions, and analysis
as the previous labs:

| Script | Role |
|---|---|
| `R/gendata.r` | Setup: generate and save the three datasets |
| `R/lab6.r` | Tasks 1-3: define `e_step()`, `m_step()`, and `em_mixture()` |
| `R/run_lab6.r` | Task 4: define the `starts` data frame, read the saved data, fit all 15 models, and save comparisons and plots |

The five starting values are defined directly in `R/run_lab6.r`.

### Setup: Simulate the three datasets

Write the data-generation script in `R/gendata.r`, following the separate
`gendata.r` workflow in Labs 1 and 3. Use `rnorm()` for the three component
counts and `sample.int()` to shuffle rows. Check the counts with `table()`.
Use seed 123 for each $n \in \{100,1000,10000\}$, and save the three CSVs in
`data-raw/`. The generator resets the same seed for each sample size, so these
datasets should not be treated as independent simulation replications.

The analysis runner reads these saved CSVs and reuses each dataset for all five
starts. Run data generation before the tests or analysis; the runner reports
missing files and directs you to `R/gendata.r`.

### Task 1: E-step and observed-data log-likelihood

Write `e_step(x, proportions, means, variance)`. At iteration $t$, calculate
the responsibilities using the current parameter estimates:

$$
\gamma_{ik}^{(t)}
=
P(Z_i=k\mid x_i,\theta^{(t)})
=
\frac{\pi_k^{(t)}\phi(x_i;\mu_k^{(t)},\sigma^{2(t)})}
{\sum_{h=1}^{3}\pi_h^{(t)}\phi(x_i;\mu_h^{(t)},\sigma^{2(t)})}.
$$

Implement this calculation on the log scale. Suppressing the iteration
superscript for readability, compute

$$
\begin{aligned}
a_{ik} &= \log\pi_k + \log\phi(x_i;\mu_k,\sigma^2), \\
m_i &= \max_{1\leq k\leq3} a_{ik}, \\
L_i &= m_i + \log\left[\sum_{k=1}^{3}\exp(a_{ik}-m_i)\right], \\
\gamma_{ik} &= \exp(a_{ik}-L_i), \\
\ell(\theta) &= \sum_{i=1}^{n}L_i.
\end{aligned}
$$

In R, `dnorm(x, mean, sd, log = TRUE)` returns the log normal density directly.
Use `sd = sqrt(variance)`.

Return an $n \times 3$ responsibility matrix and the observed-data
log-likelihood. Check that

$$
\gamma_{ik}\geq0,
\qquad
\sum_{k=1}^{3}\gamma_{ik}=1.
$$

Subtracting the row maximum avoids underflow when densities are very small.
Unlike K-means, each observation contributes fractionally to multiple classes.

### Task 2: M-step with one shared variance

Write `m_step(x, responsibilities)`. Hold the E-step responsibilities fixed
and maximize

$$
Q(\theta\mid\theta^{(t)})
=
\sum_{i=1}^{n}\sum_{k=1}^{3}
\gamma_{ik}^{(t)}
\left[\log\pi_k+\log\phi(x_i;\mu_k,\sigma^2)\right].
$$

The resulting updates are

$$
\begin{aligned}
N_k^{(t)} &= \sum_{i=1}^{n}\gamma_{ik}^{(t)}, \\
\pi_k^{(t+1)} &= \frac{N_k^{(t)}}{n}, \\
\mu_k^{(t+1)} &= \frac{\sum_{i=1}^{n}\gamma_{ik}^{(t)}x_i}{N_k^{(t)}}, \\
\sigma^{2(t+1)} &= \frac{1}{n}
\sum_{i=1}^{n}\sum_{k=1}^{3}
\gamma_{ik}^{(t)}\left(x_i-\mu_k^{(t+1)}\right)^2.
\end{aligned}
$$

Use the **updated means** in the variance formula. Pool the weighted residual
sum of squares across all components and divide by $n$, not $n-1$ or $n-3$.
This is a maximum-likelihood update. Do not recompute responsibilities between
the mean and variance updates. Stop informatively if an effective class count
or the variance becomes zero; this implementation does not silently restart or
clip the variance.

### Task 3: Repeat and check monotonicity

Write `em_mixture()` to alternate the two steps. Store the initial likelihood
at iteration 0, then evaluate and store it after every full M-step using the
updated parameters. Return final estimates, final responsibilities, the history,
iteration count, and convergence flag.

Let $\ell^{(t)}=\ell(\theta^{(t)})$. Use at most 2000 iterations and stop when

$$
\left|\ell^{(t+1)}-\ell^{(t)}\right|
\leq
10^{-10}\left(1+\left|\ell^{(t)}\right|\right).
$$

Exact EM updates should satisfy

$$
\ell^{(t+1)}\geq\ell^{(t)}.
$$

Check all differences separately from the stopping rule, allowing a numerical
tolerance of $10^{-8}$:

```r
all(diff(fit$history$loglik) >= -1e-8)
plot(fit$history$iteration, fit$history$loglik, type = "l",
     xlab = "Completed EM iteration", ylab = "Observed-data log-likelihood")
```

The $10^{-8}$ allowance is for floating-point rounding in the total likelihood.
A larger decrease is an error, not evidence of convergence. Preserve the actual
likelihood history; do not replace it with a cumulative maximum. An iteration
limit must produce a warning and `converged = FALSE`.

### Task 4: Five starting values and aligned comparisons

Define these five full initializations in `R/run_lab6.r` as the `starts` data
frame, and use them for each $n$:

| Start | Initial proportions | Initial means | Initial variance |
|---|---|---|---:|
| 1 | $(1/3, 1/3, 1/3)$ | $(-4, -0.5, 5)$ | 2 |
| 2 | $(0.2, 0.3, 0.5)$ | $(-2, 1, 3)$ | 4 |
| 3 | $(1/3, 1/3, 1/3)$ | $(-0.2, 0, 0.2)$ | 8 |
| 4 | $(0.5, 0.3, 0.2)$ | $(-5, -3, -1)$ | 2 |
| 5 | $(0.2, 0.5, 0.3)$ | $(5, 1, -2)$ | 3 |

These include close means, means concentrated on one side, and reversed labels.
The starts are deterministic; only data generation uses randomness. Fit all
15 combinations. Within each dataset, compare final likelihoods, iteration
counts, convergence flags, proportions, means, and shared variance.

Component numbers are arbitrary. Order each final fit by increasing mean,
reordering its proportions and responsibility columns with the same permutation,
before comparing parameters. Compare aligned estimates with

$$
\boldsymbol{\pi}=(0.25,0.50,0.25),
\qquad
\boldsymbol{\mu}=(-3,0,4),
\qquad
\sigma^2=1.
$$

Choose the largest final likelihood among the five converged fits for each $n$. Tiny likelihood gaps can reflect the stopping tolerance.
Agreement across starts does not establish a global maximum. For start
$s\in\{1,\ldots,5\}$, the runner calculates the likelihood gap within each dataset:

$$
\operatorname{gap}_s
=
\max_{1\leq r\leq5}\ell(\widehat\theta_r)
-
\ell(\widehat\theta_s).
$$

A gap of zero identifies a fit with the largest final likelihood among the five
starts.

## How to run

From the course repository root in R:

```r
setwd("Lab6")
source(file.path("R", "gendata.r"))
source(file.path("tests", "lab6_test.r"))
source(file.path("R", "run_lab6.r"))
```

Or, with `Lab6` as the shell working directory:

```sh
Rscript R/gendata.r
Rscript tests/lab6_test.r
Rscript R/run_lab6.r
```

The data-generation script saves the three datasets. The runner reads them
and saves four result CSVs and two PDFs:

- `starting_values.csv`: the five complete initial parameter sets.
- `start_comparison.csv`: all 15 aligned fits, their minimum likelihood gains,
  monotonicity flags, and likelihood gaps from the best start within each $n$.
- `best_estimates.csv`: one selected fit per $n$; exact likelihood ties use the
  first start.
- `loglik_history.csv`: sample size, start, iteration, and observed-data log-likelihood,
  including iteration 0.
- `loglik_history.pdf`: three pages, one per $n$, with a separate panel per start.
  Panel axes vary so every trace is visible.
- `mixture_fits.pdf`: histograms with the true and best fitted mixture densities.

## Reference results

With seed 123 for each sample size, all 15 fits converge and every recorded likelihood
increment is nonnegative. After aligning labels, the five starts give essentially
the same estimates within each $n$; the largest within-dataset likelihood gap is
less than $10^{-6}$. The best-of-five estimates, rounded to four decimals, are:

| $n$ | Proportions (ordered by mean) | Means | Shared variance |
|---:|---|---|---:|
| 100 | $(0.2383, 0.5055, 0.2562)$ | $(-3.1373, 0.0099, 4.2264)$ | 0.7754 |
| 1000 | $(0.2627, 0.4896, 0.2477)$ | $(-2.9681, 0.0914, 4.0414)$ | 0.9411 |
| 10000 | $(0.2521, 0.4999, 0.2480)$ | $(-2.9781, 0.0152, 3.9955)$ | 0.9974 |

Initialization affects speed substantially:

| $n$ | Start 1 | Start 2 | Start 3 | Start 4 | Start 5 |
|---:|---:|---:|---:|---:|---:|
| 100 | 15 | 32 | 267 | 96 | 22 |
| 1000 | 22 | 42 | 328 | 151 | 30 |
| 10000 | 25 | 47 | 336 | 154 | 34 |

Start 3, with nearby initial means and a large variance, is slowest for every $n$.
Its long plateaus show why apparent visual flatness alone is a poor stopping
rule. All starts agree here, but other data or starts can reach different local
solutions. Numerical stopping decisions can vary slightly across R platforms.

## Discuss the results

1. Why can an observation have positive membership in all three classes?
2. Why is the shared-variance denominator $n$? How would the update differ with
   three separate variances?
3. Why should each exact EM iteration avoid decreasing the observed likelihood?
4. Do the five starts reach the same likelihood and aligned estimates? Which
   initialization takes the most iterations?
5. Are estimates closer to the generating values at larger $n$? A single dataset
   per $n$ illustrates sampling variation; it does not establish a general trend
   or estimate bias and RMSE. Repeated simulations would be needed for that.
6. Why can estimates differ from truth even with exact class counts? Why are
   raw likelihood totals across different $n$ not a measure of which fit is better?
7. Does a flat likelihood trace prove a global maximum? What happens if all
   three initial means are exactly equal?

## Reproducibility and checks

Tests verify exact simulation counts, seeded reproducibility, Bayes' rule,
log-scale stability, a hand-calculated pooled variance, all 15 likelihood traces,
label invariance, iteration-limit reporting, and agreement with an independent
local `optim()` check. The optimizer check assesses a stationary fit, not global
optimality. Data and results are generated with base R.

## Author

Xiangyu Song

## Date

2026-09-29
