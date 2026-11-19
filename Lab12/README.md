# Lab 12: Ideal Points with MCMC and EM

## Overview

Repeat the House roll-call analysis from Lab 11 using `pscl::ideal()`, then
estimate one-dimensional ideal points with `emIRT::binIRT()`. Add posterior
means, posterior medians, and EM estimates to a legislator table. Use scatter
plots to compare the posterior mean with the posterior median for each
legislator and to compare estimates across methods.

This lab builds on Lab 11's roll-call data, party-based starting values,
normalization, and Bayesian estimation, and returns to the EM algorithm
introduced in Lab 6.

## Learning objectives

1. Repeat the ideal-point estimation code from Lab 11.
2. Extract retained MCMC draws and calculate a posterior median per legislator.
3. Prepare the same roll-call sample for estimation with `emIRT`.
4. Distinguish EM point estimates from MCMC posterior summaries.
5. Compare posterior means and medians with a scatter plot.
6. Align location, scale, and sign before comparing the two fitted axes.

## Project structure

```text
Lab12/
|-- R/
|   |-- lab12.r       # in-class estimation and scatter plots
|-- README.md
```

## Requirements and setup

Use R with `Rvoteview`, `data.table`, `pscl`, and `emIRT`, plus internet access
for the Voteview search and download. Install packages once if needed:

```r
install.packages(c("data.table", "pscl", "emIRT", "remotes"))
remotes::install_github("voteview/Rvoteview")
```

If the Lab 11 packages are already installed, only `emIRT` is new. The plots
use base R.

From the course repository root in R:

```r
setwd("Lab12")
```

## In-class exercise

Open `R/lab12.r` and run the code section by section. Tasks 1-3 repeat the
Lab 11 workflow; the remaining tasks extend it.

### Task 1: Search for and download roll-call votes

```r
library("Rvoteview")
library('data.table')
library(pscl)
library(emIRT)

res <- voteview_search(chamber = "House", congress = 119)
rc <- voteview_download(res$id)
```

As in Lab 11, `res` contains the search results, `res$id` supplies the vote
identifiers, and `rc` is a `rollcall` object. Rows of `rc$votes` represent
legislators; columns represent roll calls. `rc$legis.data` contains metadata.

1. What does one entry of `rc$votes` represent?
2. Why should both estimators use the same roll calls and legislators?

### Task 2: Construct starting values

```r
n <- dim(rc$legis.data)[1]
x0 <- rep(0,n)
x0[rc$legis.data$party_code=="100"] <- -1
x0[rc$legis.data$party_code=="200"] <- 1
```

Democrats start at -1, Republicans at 1, and other party codes at zero.
These values initialize estimation; they do not fix legislators' positions
or specify party-specific priors.

1. Why does `x0` need one value per legislator?
2. Must the fitted positions equal their starting values?

### Task 3: Repeat Bayesian ideal-point estimation

```r
set.seed(5625)
id1 <- ideal(rc,
             d=1,
             startvals=list(x=x0),
             normalize=TRUE,
             maxiter=10000,
             burnin=100,
             thin=100)

id1$n
id1$m
id1$xbar
id1$codes
```

The estimation arguments are the same as in Lab 11. `d=1` selects one latent
dimension, `normalize=TRUE` identifies location and scale, and `ideal()`
uses MCMC to estimate legislator and roll-call parameters jointly. The seed
makes the stochastic run repeatable with the same data and software.

1. What do `maxiter`, `burnin`, and `thin` control?
2. Why does normalization leave the substantive direction of the axis open?
3. How are `id1$xbar` and `x0` different?

### Task 4: Attach posterior means and medians

```r
rc_fit <- extractRollCallObject(id1)
stopifnot(identical(rownames(rc_fit$votes), rownames(id1$xbar)))
legis119 <- data.table(rc_fit$legis.data)
legis119$ideal <- as.numeric(id1$xbar[, 1])

dim(id1$x)
keep <- as.numeric(dimnames(id1$x)[[1]]) > id1$call$burnin
x_draws <- id1$x[keep, , 1, drop = FALSE]
legis119$ideal_median <- apply(x_draws, 2, median)
```

`extractRollCallObject()` recovers the sample used by `ideal()`, which drops
unanimous roll calls by default. The row-name check verifies that the table
and fitted estimates refer to the same legislators in the same order.

The dimensions of `id1$x` are saved iterations, legislators, and latent
dimensions. `keep` retains iterations strictly after the burn-in boundary,
matching the rule used for `id1$xbar`. With these settings, the stored draw
at iteration 100 is excluded, leaving 99 draws at iterations 200-10,000.
`drop = FALSE` preserves the array dimensions. `apply(..., 2, median)`
summarizes the draws separately for each legislator.

The script also checks that means calculated from these retained draws
reproduce `id1$xbar`. The `ideal` column keeps the name used in Lab 11.

1. Why do we summarize over draws rather than over legislators?
2. What would `median(id1$xbar)` calculate instead?
3. When might a legislator's posterior mean and median differ?

### Task 5: Estimate ideal points with emIRT

```r
rc_em <- convertRC(rc_fit)
rc_em$votes[is.na(rc_em$votes) | rc_em$votes == 9] <- 0
stopifnot(all(rc_em$votes %in% c(-1, 0, 1)))

p <- makePriors(rc_em$n, rc_em$m, 1)
s <- getStarts(rc_em$n, rc_em$m, 1)
fit_index <- match(rownames(rc_fit$votes), rownames(rc$votes))
stopifnot(!anyNA(fit_index))
s$x <- matrix(x0[fit_index], ncol = 1)

id2 <- binIRT(.rc = rc_em,
              .starts = s,
              .priors = p,
              .D = 1,
              .control = list(threads = 1,
                              verbose = TRUE,
                              maxit = 1000,
                              thresh = 1e-6,
                              asEM = TRUE))

id2$runtime
legis119$ideal_em <- as.numeric(id2$means$x[, 1])
```

`convertRC()` prepares the vote coding for `binIRT()`: Yea is 1, Nay is -1,
and missing is 0. The extra recoding treats any remaining not-in-legislature
code (9) or `NA` as missing. A missing vote must not be treated as a Nay.

`makePriors()` creates normal priors; `getStarts()` creates starting-value
matrices for the legislators and roll-call parameters. Replacing `s$x`
uses the same party-based positions as the MCMC fit, matched by legislator.

| Argument | Meaning |
|---|---|
| `.rc = rc_em` | Use the recoded analysis sample |
| `.starts = s` | Supply starting-value matrices |
| `.priors = p` | Supply prior parameters |
| `.D = 1` | Estimate one latent dimension |
| `threads = 1` | Use one computational thread |
| `maxit = 1000` | Allow up to 1,000 EM iterations |
| `thresh = 1e-6` | Set the stopping tolerance |
| `asEM = TRUE` | Use EM estimation |

Inspect `id2$runtime$conv`: 1 means the stopping threshold was reached;
0 means the iteration limit was reached. The script warns about the latter;
increase `maxit` and rerun Task 5 before interpreting that comparison.

`id2$means$x` contains one EM point estimate per legislator. Despite the
component name `means`, this EM run does not produce posterior draws from
which to calculate a separate posterior median. The mean-median comparison
in the next task uses the MCMC fit.

1. How does EM estimation differ from MCMC sampling?
2. Why must we convert the vote codes before calling `binIRT()`?
3. Why is reaching an iteration limit different from meeting a stopping rule?

### Task 6: Plot posterior means against posterior medians

```r
lim <- range(c(legis119$ideal, legis119$ideal_median))
plot(legis119$ideal, legis119$ideal_median,
     xlab = "Posterior mean ideal point",
     ylab = "Posterior median ideal point",
     main = "Posterior means and medians",
     xlim = lim, ylim = lim, asp = 1,
     pch = 19, col = party_col)
abline(a = 0, b = 1, lty = 2)
```

The script creates `party_col` immediately before this plot: blue for
Democrats, red for Republicans, and gray for other party codes. It also adds
a legend. Each point represents one legislator. The dashed line is the
45-degree equality line, not a fitted regression line. Both axes use the
same limits and units.

The script reports the correlation and summarizes the median-minus-mean
differences. Close agreement is consistent with roughly symmetric marginal
posteriors, but does not establish symmetry or MCMC convergence. Differences
can also reflect Monte Carlo error from the small number of retained draws.

1. What does a point above the dashed line tell us?
2. Can correlation be high even when some points are away from equality?
3. Why is agreement between means and medians not a convergence diagnostic?

### Task 7: Compare MCMC and EM estimates

```r
legis119$ideal_std <- as.numeric(scale(legis119$ideal))
legis119$ideal_em_std <- as.numeric(scale(legis119$ideal_em))
if (cor(legis119$ideal_std, legis119$ideal_em_std) < 0) {
  legis119$ideal_em_std <- -legis119$ideal_em_std
}
```

Standardization centers each set of estimates and divides by its standard
deviation across legislators. The sign check puts EM in the same plotting
direction as MCMC. It does not establish which end is liberal or conservative.
The original estimates remain in `ideal` and `ideal_em`.

Run the second scatter plot in the script. Its horizontal axis is the
standardized MCMC posterior mean; its vertical axis is the standardized,
sign-aligned EM estimate. The dashed line again represents equality.

1. Why might unaligned estimates have a negative correlation?
2. Why does comparison require attention to location and scale as well as sign?
3. Why can EM point estimates differ from MCMC posterior means even after alignment?

## How to run

From the course repository root in R:

```r
setwd("Lab12")
source(file.path("R", "lab12.r"), echo = TRUE)
```

For in-class work, run each section interactively and inspect the first
scatter plot before creating the second. In RStudio, the Plots pane's back
arrow returns to the previous plot. The download and model fits may take time.

## Reproducibility and interpretation

The script queries the live Voteview database, so the downloaded votes can
change. The seed controls MCMC randomness for a fixed input; it does not
freeze the data. The script leaves the fitted objects, `legis119`, and plots
available in R and does not explicitly save files.

The MCMC settings repeat the Lab 11 demonstration. Longer runs, multiple
chains, and convergence checks are needed for substantive inference. EM's
stopping flag also does not establish a global optimum. Compare the two
methods using the same analysis sample, and distinguish agreement in point
estimates from agreement about uncertainty.

For the package interfaces and output definitions, see the
[emIRT reference manual](https://kosukeimai.r-universe.dev/emIRT/doc/manual.html)
and the [pscl reference manual](https://cran.r-project.org/web/packages/pscl/pscl.pdf).

## Author

Xiangyu Song

## Date

2026-11-17
