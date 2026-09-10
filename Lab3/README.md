# Lab 3: Bootstrap Inference for a Mean and Regression Coefficient

## Overview

This lab applies the nonparametric bootstrap from Week 3 to two smooth
statistics: the sample mean of `X` and the OLS slope on `X`. Students first
write the bootstrap as an explicit loop and then repeat the analysis with
`boot::boot()`.

The simulated data follow

\[
X_i \sim N(0, 2), \qquad
Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i,
\]

where `N(0, 2)` uses 2 as the variance, so `sd = sqrt(2)`, and

- \(\varepsilon_i \sim N(0, 1)\);
- \(\beta_0 = 1.5\);
- \(\beta_1 = -2\);
- \(N = 500\); and
- \(B = 1000\) bootstrap replications.

Each bootstrap resample draws 500 rows with replacement. Resampling complete
rows preserves the relationship between `X` and `Y` for the regression.

## Learning objectives

After completing this lab, students should be able to:

1. generate reproducible data from a linear regression model;
2. implement a nonparametric bootstrap with a preallocated loop;
3. bootstrap a sample mean and an OLS regression coefficient;
4. calculate a bootstrap standard error and bias estimate;
5. construct percentile and bias-corrected percentile intervals; and
6. reproduce the calculation with `boot::boot()` and `boot::boot.ci()`.

## Project structure

```text
Lab3/
|-- data-raw/
|   |-- lab3.csv
|-- R/
|   |-- gendata.r
|   |-- lab3.r
|   |-- run_lab3.r
|-- results/
|   |-- bootstrap_replicates.csv
|   |-- bootstrap_summary.csv
|-- tests/
|   |-- lab3_test.r
|-- README.md
```

## Requirements

- R 4.0 or later.
- The recommended package `boot`.

Install `boot` once if it is not available:

```r
install.packages("boot")
```

## Bootstrap quantities

For a statistic \(\hat\theta\) and bootstrap replicates
\(\hat\theta^{*(1)}, \ldots, \hat\theta^{*(B)}\), the lab calculates:

- bootstrap standard error: the sample standard deviation of the replicates;
- estimated bias: \(\bar\theta^* - \hat\theta\);
- bias-corrected point estimate: \(2\hat\theta - \bar\theta^*\);
- percentile interval: \([q^*_{.025}, q^*_{.975}]\); and
- bias-corrected percentile interval:
  \([2\hat\theta-q^*_{.975}, 2\hat\theta-q^*_{.025}]\).

The last interval is commonly called the **basic** or **reverse-percentile**
interval. It corresponds to `type = "basic"` in `boot::boot.ci()`. It is not
the bias-corrected and accelerated (BCa) interval.

## Tasks

### Task 1: Generate the data

Run `R/gendata.r`. Check that the result contains 500 observations and that
the simulated relationship between `X` and `Y` has a negative slope near -2.

### Task 2: Bootstrap the mean with a loop

Use a loop over `b = 1, ..., B`. For each replication:

1. sample 500 row indices with replacement;
2. calculate the mean of the resampled `X` values; and
3. store the result in a preallocated numeric vector.

Use the resulting distribution to calculate the standard error, bias,
bias-corrected estimate, percentile interval, and bias-corrected interval.

### Task 3: Bootstrap the regression coefficient with a loop

Repeat Task 2 for the coefficient on `X` from `lm(Y ~ X)`. Resample rows, not
`X` and `Y` separately. Independent resampling would break the observed
relationship between the predictor and outcome.

### Task 4: Use the `boot` package

Write statistic functions with the `(data, indices)` interface. Run
`boot::boot()` for the mean and slope, again with `R = 1000`, and obtain the
percentile and basic intervals with:

```r
boot::boot.ci(result, type = c("basic", "perc"))
```

### Task 5: Compare the methods

Compare the loop and package results. They need not be numerically identical
because Monte Carlo samples can differ, but the point estimates, standard
errors, and intervals should be close. Increasing `B` reduces Monte Carlo
error; it does not remove bootstrap approximation error due to finite `N`.

## How to run the materials

Start R with `Lab3` as the working directory. If the course repository is the
current working directory, run:

```r
setwd("Lab3")
```

Generate the data:

```r
source(file.path("R", "gendata.r"))
```

Run the checks:

```r
source(file.path("tests", "lab3_test.r"))
```

Run the complete analysis:

```r
source(file.path("R", "run_lab3.r"))
```

The analysis writes a four-row summary to
`results/bootstrap_summary.csv` and all four sets of bootstrap replicates to
`results/bootstrap_replicates.csv`.



## Reproducibility

The data generator uses seed `123`, and the bootstrap analysis uses seed
`123`. All paths are project-relative and assume that `Lab3` is the working
directory.

## Author

Xiangyu Song

## Date

2026-09-08
