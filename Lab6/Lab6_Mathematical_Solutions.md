---
author: Xiangyu Song
tags: 
created: 2026-09-30 23:13
modified: 2026-09-30 23:38:41
fontsize: 11pt
header-includes:
    - \usepackage[notes]{note}
    - \linespread{1}
---

# Lab 6: Mathematical Solutions

## Model and notation

Let $x_1,\ldots,x_n$ be independent observations from a three-component
Gaussian mixture model. Introduce the unobserved component label
$Z_i\in\{1,2,3\}$. The model is

```math
\Pr(Z_i=k)=\pi_k,
\qquad
X_i\mid Z_i=k\sim N(\mu_k,\sigma^2),
\qquad k=1,2,3,
```

where

```math
\pi_k>0,\qquad
\sum_{k=1}^3\pi_k=1,\qquad
\sigma^2>0.
```

The three components have different mixing proportions and means, but they
share the same variance. Write the full parameter vector as

```math
\theta=(\boldsymbol{\pi},\boldsymbol{\mu},\sigma^2),
```

where

```math
\boldsymbol{\pi}=(\pi_1,\pi_2,\pi_3)
\quad\text{and}\quad
\boldsymbol{\mu}=(\mu_1,\mu_2,\mu_3).
```

The Gaussian density is

```math
\phi(x;\mu,\sigma^2)
=
\frac{1}{\sqrt{2\pi\sigma^2}}
\exp\left\{-\frac{(x-\mu)^2}{2\sigma^2}\right\}.
```

Define the component indicators

```math
Z_{ik}=\mathbf{1}\{Z_i=k\}.
```

For each observation, exactly one of the three indicators equals one:

```math
\sum_{k=1}^3 Z_{ik}=1.
```

## Complete-data likelihood

If the component labels were observed, observation $i$ would contribute

```math
\prod_{k=1}^3
\left[\pi_k\phi(x_i;\mu_k,\sigma^2)\right]^{Z_{ik}}
```

to the likelihood. Therefore, the complete-data likelihood is

```math
L_c(\theta)
=
\prod_{i=1}^n\prod_{k=1}^3
\left[\pi_k\phi(x_i;\mu_k,\sigma^2)\right]^{Z_{ik}}.
```

Taking logarithms gives the complete-data log-likelihood:

```math
\ell_c(\theta)
=
\sum_{i=1}^n\sum_{k=1}^3
Z_{ik}
\left[
\log\pi_k
-\frac{1}{2}\log(2\pi\sigma^2)
-\frac{(x_i-\mu_k)^2}{2\sigma^2}
\right].
```

The indicators $Z_{ik}$ are not observed. The EM algorithm handles this
missing information by replacing each $Z_{ik}$ with its conditional
expectation under the current parameter values.

## Task 1: E-step and observed-data log-likelihood

At iteration $t$, define the responsibility

```math
\gamma_{ik}^{(t)}
=
\mathrm{E}\left(Z_{ik}\mid x_i,\theta^{(t)}\right)
=
\Pr\left(Z_i=k\mid x_i,\theta^{(t)}\right).
```

The responsibility is the posterior probability that observation $i$ belongs
to component $k$.

By Bayes' rule,

```math
\Pr(Z_i=k\mid x_i,\theta^{(t)})
=
\frac{
\Pr(Z_i=k\mid\theta^{(t)})
f(x_i\mid Z_i=k,\theta^{(t)})
}{
\sum_{j=1}^3
\Pr(Z_i=j\mid\theta^{(t)})
f(x_i\mid Z_i=j,\theta^{(t)})
}.
```

Substituting the mixture proportions and Gaussian densities yields the
E-step formula

```math
\gamma_{ik}^{(t)}
=
\frac{
\pi_k^{(t)}
\phi\left(x_i;\mu_k^{(t)},(\sigma^2)^{(t)}\right)
}{
\sum_{j=1}^3
\pi_j^{(t)}
\phi\left(x_i;\mu_j^{(t)},(\sigma^2)^{(t)}\right)
}.
```

For every $i$, the responsibilities satisfy

```math
0\leq\gamma_{ik}^{(t)}\leq 1
\quad\text{and}\quad
\sum_{k=1}^3\gamma_{ik}^{(t)}=1.
```

Thus, the responsibilities are soft component assignments. An observation can
receive positive probability for more than one component.

### Observed-data likelihood

Because the label $Z_i$ is unobserved, the marginal density of $x_i$ is found
by summing over all possible components:

```math
f(x_i\mid\theta)
=
\sum_{k=1}^3
\pi_k\phi(x_i;\mu_k,\sigma^2).
```

The observed-data likelihood is therefore

```math
L(\theta)
=
\prod_{i=1}^n
\left[
\sum_{k=1}^3
\pi_k\phi(x_i;\mu_k,\sigma^2)
\right],
```

and the observed-data log-likelihood is

```math
\ell(\theta)
=
\sum_{i=1}^n
\log\left[
\sum_{k=1}^3
\pi_k\phi(x_i;\mu_k,\sigma^2)
\right].
```

This log-likelihood is used to monitor convergence. In exact arithmetic, a
correct EM iteration does not decrease it:

```math
\ell\left(\theta^{(t+1)}\right)
\geq
\ell\left(\theta^{(t)}\right).
```

For numerical stability, the inner sum can be evaluated on the log scale. If

```math
a_{ik}
=
\log\pi_k
\log\phi(x_i;\mu_k,\sigma^2),
```

then

```math
\log\sum_{k=1}^3 e^{a_{ik}}
=
m_i+\log\sum_{k=1}^3 e^{a_{ik}-m_i},
\qquad
m_i=\max_k a_{ik}.
```

Subtracting $m_i$ before exponentiation reduces the risk of numerical
underflow.

## Task 2: M-step derivation

The EM auxiliary function is the conditional expectation of the complete-data
log-likelihood:

```math
Q(\theta\mid\theta^{(t)})
=
\mathrm{E}_{Z\mid X,\theta^{(t)}}
\left[\ell_c(\theta)\right].
```

Since

```math
\mathrm{E}
\left(Z_{ik}\mid x_i,\theta^{(t)}\right)
=
\gamma_{ik}^{(t)},
```

we replace $Z_{ik}$ by $\gamma_{ik}^{(t)}$:

```math
Q(\theta\mid\theta^{(t)})
=
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}
\left[
\log\pi_k
-\frac{1}{2}\log(2\pi\sigma^2)
-\frac{(x_i-\mu_k)^2}{2\sigma^2}
\right].
```

Define the effective membership count of component $k$ by

```math
N_k^{(t)}
=
\sum_{i=1}^n\gamma_{ik}^{(t)}.
```

Because the responsibilities sum to one for every observation,

```math
\sum_{k=1}^3N_k^{(t)}
=
\sum_{i=1}^n\sum_{k=1}^3\gamma_{ik}^{(t)}
=n.
```

### Updating the mixing proportions

The part of $Q$ that depends on $\boldsymbol{\pi}$ is

```math
Q_\pi
=
\sum_{k=1}^3N_k^{(t)}\log\pi_k.
```

We maximize this expression subject to
$\sum_{k=1}^3\pi_k=1$. Introduce a Lagrange multiplier $\lambda$:

```math
\mathcal{L}(\boldsymbol{\pi},\lambda)
=
\sum_{k=1}^3N_k^{(t)}\log\pi_k
+\lambda\left(1-\sum_{k=1}^3\pi_k\right).
```

The first-order condition for $\pi_k$ is

```math
\frac{\partial\mathcal{L}}{\partial\pi_k}
=
\frac{N_k^{(t)}}{\pi_k}-\lambda=0.
```

Hence,

```math
\pi_k=\frac{N_k^{(t)}}{\lambda}.
```

Summing over $k$ and applying the constraint gives

```math
1
=
\sum_{k=1}^3\pi_k
=
\frac{1}{\lambda}\sum_{k=1}^3N_k^{(t)}
=
\frac{n}{\lambda}.
```

Therefore $\lambda=n$, and the update is

```math
\pi_k^{(t+1)}
=
\frac{N_k^{(t)}}{n}
=
\frac{1}{n}\sum_{i=1}^n\gamma_{ik}^{(t)}.
```

Thus, the updated mixture proportion is the effective fraction of
observations assigned to component $k$.

### Updating the component means

For a fixed component $k$, the terms involving $\mu_k$ are

```math
Q_{\mu_k}
=
-\frac{1}{2\sigma^2}
\sum_{i=1}^n
\gamma_{ik}^{(t)}(x_i-\mu_k)^2.
```

Differentiate with respect to $\mu_k$:

```math
\frac{\partial Q_{\mu_k}}{\partial\mu_k}
=
\frac{1}{\sigma^2}
\sum_{i=1}^n
\gamma_{ik}^{(t)}(x_i-\mu_k).
```

Setting this derivative equal to zero gives

```math
\sum_{i=1}^n\gamma_{ik}^{(t)}x_i
-\mu_k\sum_{i=1}^n\gamma_{ik}^{(t)}
=0.
```

Solving for $\mu_k$ produces

```math
\mu_k^{(t+1)}
=
\frac{
\sum_{i=1}^n\gamma_{ik}^{(t)}x_i
}{
\sum_{i=1}^n\gamma_{ik}^{(t)}
}
=
\frac{
\sum_{i=1}^n\gamma_{ik}^{(t)}x_i
}{
N_k^{(t)}
}.
```

The updated mean is a responsibility-weighted average of the observations.

### Updating the shared variance

Let $v=\sigma^2$. After inserting the updated means, define the total
responsibility-weighted residual sum of squares

```math
S^{(t+1)}
=
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}
\left(x_i-\mu_k^{(t+1)}\right)^2.
```

The terms of $Q$ that depend on $v$ are

```math
Q_v
=
-\frac{1}{2}
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}\log v
-\frac{1}{2v}S^{(t+1)}
+C,
```

where $C$ does not depend on $v$. Since

```math
\sum_{i=1}^n\sum_{k=1}^3\gamma_{ik}^{(t)}=n,
```

this simplifies to

```math
Q_v
=
-\frac{n}{2}\log v
-\frac{S^{(t+1)}}{2v}
+C.
```

Differentiate:

```math
\frac{\partial Q_v}{\partial v}
=
-\frac{n}{2v}
+\frac{S^{(t+1)}}{2v^2}.
```

Setting the derivative equal to zero gives

```math
-\frac{n}{2v}
+\frac{S^{(t+1)}}{2v^2}
=0.
```

Multiplying by $2v^2$ yields

```math
-nv+S^{(t+1)}=0,
```

so

```math
(\sigma^2)^{(t+1)}
=
\frac{1}{n}
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}
\left(x_i-\mu_k^{(t+1)}\right)^2.
```

The denominator is $n$, not $3n$, because the three responsibilities for each
observation sum to one. It is also not $n-1$: this is a maximum-likelihood
update rather than an unbiased sample-variance estimator.

### Complete M-step

The three M-step updates are

```math
\begin{aligned}
\pi_k^{(t+1)}
&=
\frac{1}{n}\sum_{i=1}^n\gamma_{ik}^{(t)},\\
\mu_k^{(t+1)}
&=
\frac{
\sum_{i=1}^n\gamma_{ik}^{(t)}x_i
}{
\sum_{i=1}^n\gamma_{ik}^{(t)}
},\\
(\sigma^2)^{(t+1)}
&=
\frac{1}{n}
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}
\left(x_i-\mu_k^{(t+1)}\right)^2.
\end{aligned}
```

## Task 3: EM algorithm and convergence

### Initialization

Choose starting values

```math
\theta^{(0)}
=
\left(
\boldsymbol{\pi}^{(0)},
\boldsymbol{\mu}^{(0)},
(\sigma^2)^{(0)}
\right),
```

with positive proportions that sum to one and a positive variance.

### Iterative procedure

For $t=0,1,2,\ldots$, repeat the following steps.

1. **E-step.** Compute all posterior responsibilities:

   $$
   \gamma_{ik}^{(t)}
   =
   \frac{
   \pi_k^{(t)}
   \phi\left(x_i;\mu_k^{(t)},(\sigma^2)^{(t)}\right)
   }{
   \sum_{j=1}^3
   \pi_j^{(t)}
   \phi\left(x_i;\mu_j^{(t)},(\sigma^2)^{(t)}\right)
   }.
   $$

2. **M-step.** Use the responsibilities to update
   $\pi_k^{(t+1)}$, $\mu_k^{(t+1)}$, and $(\sigma^2)^{(t+1)}$ with
   the formulas derived above.

3. **Evaluation.** Compute the observed-data log-likelihood
   $\ell(\theta^{(t+1)})$.

4. **Stopping rule.** Stop when the log-likelihood change is sufficiently
   small, for example,

   $$
   \left|
   \ell(\theta^{(t+1)})
   -\ell(\theta^{(t)})
   \right|
   <\varepsilon.
   $$

A relative stopping rule is often more meaningful when the magnitude of the
log-likelihood is large:

```math
\frac{
\left|\ell(\theta^{(t+1)})-\ell(\theta^{(t)})\right|
}{
1+\left|\ell(\theta^{(t)})\right|
}
<\varepsilon.
```

It is useful to impose a maximum number of iterations as a safeguard.

### Why the log-likelihood does not decrease

At iteration $t$, EM constructs a lower bound on the observed-data
log-likelihood that touches the log-likelihood at $\theta^{(t)}$. The E-step
chooses the conditional distribution of the missing labels under
$\theta^{(t)}$. The M-step maximizes the resulting function
$Q(\theta\mid\theta^{(t)})$. Consequently,

```math
\ell(\theta^{(t+1)})
\geq
\ell(\theta^{(t)}).
```

Small numerical decreases can occur because of finite-precision arithmetic,
but a substantial decrease usually indicates an implementation error.

### What convergence means

The monotonicity property does not guarantee convergence to the global
maximum. A finite mixture likelihood can have multiple stationary points, and
EM can converge to different solutions from different starting values.
Therefore, one should run the algorithm from several initializations and
retain the converged fit with the largest observed-data log-likelihood.

## Task 4: Comparing starting values and fitted solutions

Suppose the algorithm is run from several sets of initial values. For each run,
record

- the starting proportions, means, and variance;
- the converged proportions, means, and variance;
- the final observed-data log-likelihood;
- the number of iterations; and
- whether the convergence criterion was satisfied.

### Label switching

Mixture-component labels have no intrinsic meaning. If a permutation $p$ of
$\{1,2,3\}$ is applied to the component-specific parameters, then

```math
\sum_{k=1}^3
\pi_k\phi(x;\mu_k,\sigma^2)
=
\sum_{k=1}^3
\pi_{p(k)}\phi(x;\mu_{p(k)},\sigma^2).
```

Therefore, two runs can describe exactly the same fitted density while listing
the components in different orders. Before comparing parameter estimates
across runs, align the labels. A simple convention in one dimension is

```math
\mu_1<\mu_2<\mu_3,
```

with each $\pi_k$ reordered together with its corresponding $\mu_k$.

### Comparing objective values

Let $\ell_r$ be the final log-likelihood for run $r$, and define

```math
\ell_{\max}=\max_r\ell_r.
```

The log-likelihood gap for run $r$ is

```math
\Delta_r=\ell_{\max}-\ell_r.
```

Interpretation:

- If $\Delta_r$ is approximately zero, the run reached the best solution found.
- If two runs have equal log-likelihoods after relabeling, they are effectively
  the same fitted solution.
- If $\Delta_r$ is clearly positive, the run likely converged to an inferior
  local optimum or stationary point.

The best reported fit should be the converged run with the largest final
observed-data log-likelihood, not necessarily the run that converged in the
fewest iterations.

### Interpreting the fitted parameters

After ordering the component means, interpret

```math
\widehat{\pi}_k
```

as the estimated population fraction in component $k$,

```math
\widehat{\mu}_k
```

as the estimated center of that component, and

```math
\widehat{\sigma}
=
\sqrt{\widehat{\sigma^2}}
```

as the common within-component standard deviation.

The effective fitted sample size of component $k$ is

```math
\widehat{N}_k
=
\sum_{i=1}^n\widehat{\gamma}_{ik}
=
n\widehat{\pi}_k.
```

A very small $\widehat{N}_k$ indicates that the component is supported by few
observations. Such a component can be unstable and especially sensitive to
initialization.

## Final summary

For this three-component Gaussian mixture with a shared variance, one EM
iteration consists of the following updates:

```math
\gamma_{ik}^{(t)}
=
\frac{
\pi_k^{(t)}
\phi\left(x_i;\mu_k^{(t)},(\sigma^2)^{(t)}\right)
}{
\sum_{j=1}^3
\pi_j^{(t)}
\phi\left(x_i;\mu_j^{(t)},(\sigma^2)^{(t)}\right)
},
```

```math
\pi_k^{(t+1)}
=
\frac{1}{n}\sum_{i=1}^n\gamma_{ik}^{(t)},
```

```math
\mu_k^{(t+1)}
=
\frac{
\sum_{i=1}^n\gamma_{ik}^{(t)}x_i
}{
\sum_{i=1}^n\gamma_{ik}^{(t)}
},
```

and

```math
(\sigma^2)^{(t+1)}
=
\frac{1}{n}
\sum_{i=1}^n\sum_{k=1}^3
\gamma_{ik}^{(t)}
\left(x_i-\mu_k^{(t+1)}\right)^2.
```

Repeat these steps until the observed-data log-likelihood stabilizes. Use
multiple starting values, align component labels before comparing runs, and
select the converged solution with the largest final log-likelihood.
