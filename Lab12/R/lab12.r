######################################################
#                                                    #
# Lab 12: Ideal points with MCMC and EM               #
#                                                    #
# R/lab12.r                                          #
#                                                    #
# In-class implementation:                           #
#   Task 1. search and download roll-call votes       #
#   Task 2. construct party-based starting values    #
#   Task 3. repeat Bayesian ideal-point estimation   #
#   Task 4. attach posterior means and medians       #
#   Task 5. estimate ideal points with emIRT         #
#   Task 6. plot posterior means against medians     #
#   Task 7. compare MCMC and EM estimates            #
#                                                    #
######################################################

# Use the packages from Lab 11 and add emIRT for EM estimation.
library("Rvoteview")
library('data.table')
library(pscl)
library(emIRT)


# Task 1. Search for and download roll-call votes.
# Repeat the House/Congress 119 search from Lab 11.
res <- voteview_search(chamber = "House", congress = 119)
rc <- voteview_download(res$id)
# In rc$votes, rows are legislators and columns are roll-call votes.
# rc$legis.data contains information about those legislators.
# Check: what does one entry in rc$votes represent?


# Task 2. Construct party-based starting values.
n <- dim(rc$legis.data)[1]
x0 <- rep(0,n)
x0[rc$legis.data$party_code=="100"] <- -1
x0[rc$legis.data$party_code=="200"] <- 1
# Party code 100 = Democrat; party code 200 = Republican.
# These are initial values, not fixed positions or party-specific priors.


# Task 3. Repeat the one-dimensional model from Lab 11.
# Set a seed so the MCMC run can be repeated with the same downloaded data.
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
# normalize=TRUE identifies location and scale, but not substantive direction.
# The short burn-in and thinning repeat the in-class demonstration settings;
# finishing the run does not establish MCMC convergence.


# Task 4. Attach posterior means and medians.
# ideal() drops unanimous roll calls by default. Recover its analysis sample
# so both methods use the same legislators, roll calls, and missing votes.
rc_fit <- extractRollCallObject(id1)
stopifnot(identical(rownames(rc_fit$votes), rownames(id1$xbar)))
legis119 <- data.table(rc_fit$legis.data)
legis119$ideal <- as.numeric(id1$xbar[, 1])

# id1$x has three dimensions: saved iterations, legislators, and dimensions.
dim(id1$x)
# Match the strict burn-in rule used to compute id1$xbar.
keep <- as.numeric(dimnames(id1$x)[[1]]) > id1$call$burnin
x_draws <- id1$x[keep, , 1, drop = FALSE]
# Margin 2 keeps legislators; median() summarizes each one's retained draws.
legis119$ideal_median <- apply(x_draws, 2, median)
# Confirm that averaging those same draws reproduces the supplied means.
stopifnot(isTRUE(all.equal(as.numeric(apply(x_draws, 2, mean)),
                          legis119$ideal)))
head(legis119[, .(ideal, ideal_median)])
# Check: why is median(id1$xbar) a different quantity?


# Task 5. Estimate ideal points with emIRT.
# convertRC() recodes the fitted rollcall object for binIRT().
# emIRT uses 1 = Yea, -1 = Nay, and 0 = missing, not Voteview's raw codes.
rc_em <- convertRC(rc_fit)
# Treat any remaining not-in-legislature code (9) or NA as missing too.
rc_em$votes[is.na(rc_em$votes) | rc_em$votes == 9] <- 0
stopifnot(all(rc_em$votes %in% c(-1, 0, 1)))

# Create priors and correctly sized starting-value matrices.
p <- makePriors(rc_em$n, rc_em$m, 1)
s <- getStarts(rc_em$n, rc_em$m, 1)
# Match the original starting positions to the fitted legislators by row name.
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
# conv = 1 means the stopping threshold was met; 0 means maxit was reached.
if (id2$runtime$conv != 1) {
  warning("EM reached maxit. Increase maxit and rerun Task 5 before interpreting the comparison.")
}
# binIRT() returns one point estimate per input row, in the same order.
stopifnot(nrow(id2$means$x) == nrow(legis119))
legis119$ideal_em <- as.numeric(id2$means$x[, 1])
# Despite the name 'means', this EM output is not a set of posterior draws.
# Check: how does EM estimation differ from the MCMC run in Task 3?


# Task 6. Plot posterior means against posterior medians.
# One point represents one legislator. Use party colors to read the plot.
party_col <- rep("gray50", nrow(legis119))
party_col[legis119$party_code == "100"] <- "steelblue"
party_col[legis119$party_code == "200"] <- "firebrick"

lim <- range(c(legis119$ideal, legis119$ideal_median))
plot(legis119$ideal, legis119$ideal_median,
     xlab = "Posterior mean ideal point",
     ylab = "Posterior median ideal point",
     main = "Posterior means and medians",
     xlim = lim, ylim = lim, asp = 1,
     pch = 19, col = party_col)
abline(a = 0, b = 1, lty = 2)
legend("topleft", legend = c("Democrat", "Republican", "Other"),
       col = c("steelblue", "firebrick", "gray50"), pch = 19, bty = "n")

cor(legis119$ideal, legis119$ideal_median)
summary(legis119$ideal_median - legis119$ideal)
# The dashed line is equality: a point on it has the same mean and median.
# Check: does agreement between these summaries prove convergence?


# Task 7. Compare MCMC and EM estimates.
# The two fitted axes can differ in location, scale, and sign.
# Standardize both sets of point estimates across legislators for this plot.
legis119$ideal_std <- as.numeric(scale(legis119$ideal))
legis119$ideal_em_std <- as.numeric(scale(legis119$ideal_em))
# Reverse EM's plotting direction if it is opposite to the MCMC direction.
if (cor(legis119$ideal_std, legis119$ideal_em_std) < 0) {
  legis119$ideal_em_std <- -legis119$ideal_em_std
}

lim <- range(c(legis119$ideal_std, legis119$ideal_em_std))
plot(legis119$ideal_std, legis119$ideal_em_std,
     xlab = "Standardized MCMC posterior mean",
     ylab = "Standardized EM estimate (aligned sign)",
     main = "MCMC and EM ideal points",
     xlim = lim, ylim = lim, asp = 1,
     pch = 19, col = party_col)
abline(a = 0, b = 1, lty = 2)
legend("topleft", legend = c("Democrat", "Republican", "Other"),
       col = c("steelblue", "firebrick", "gray50"), pch = 19, bty = "n")

cor(legis119$ideal_std, legis119$ideal_em_std)
# Check: why can the estimates still differ after standardization?
# Inspect legis119 for posterior means, posterior medians, and EM estimates.
