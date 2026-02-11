# Mathematical Foundations: IRT Models and Dirichlet Process Mixture Priors

## 1. Overview

This vignette provides a rigorous mathematical treatment of the theory
underlying the DPMirt package. It is intended for statisticians and
methodologically-inclined researchers who wish to understand the full
probabilistic framework connecting Item Response Theory (IRT)
measurement models, Bayesian hierarchical priors, Dirichlet Process
Mixture (DPM) extensions, identification strategies, and posterior
summary methods.

We cover the following topics in depth:

1.  The Bayesian Rasch model and its hierarchical structure
2.  The 2PL and 3PL extensions
3.  Identification indeterminacy and post-hoc rescaling
4.  Dirichlet Process Mixture priors for the latent trait distribution
5.  The concentration parameter $`\alpha`$ and its hyperprior
6.  Posterior summary theory: PM, CB, and GR estimators
7.  DP density reconstruction from MCMC output

Throughout, we carefully distinguish between **established results**
from the IRT, Bayesian nonparametric, and decision-theoretic literatures
and **novel contributions** of this work. The DPMirt package synthesizes
these components into a unified computational framework; the novelty
lies primarily in the integration and in extending the CB and GR
estimators to the semiparametric IRT setting.

> **Note:** This vignette contains no live MCMC computation. All code
> blocks show model specifications and mathematical formulas for
> reference. For applied examples with actual model fitting, see the
> companion vignettes `models-and-workflow` and `posterior-summaries`.

## 2. The Bayesian Rasch Model

### 2.1 Measurement Model

The Rasch model (Rasch, 1960) specifies the probability that person
$`p`$$`(p = 1, \ldots, N)`$ endorses item $`i`$$`(i = 1, \ldots, I)`$ as
a function of two parameters: a person ability $`\theta_p`$ and an item
difficulty $`\beta_i`$.

The measurement model is:
``` math
y_{ip} \sim \text{Bernoulli}(\pi_{ip}),
```
where the item response function (IRF) maps the latent linear predictor
to the probability scale via the logistic link:
``` math
\text{logit}(\pi_{ip}) = \theta_p - \beta_i.
```

Equivalently, the probability of a correct response is:
``` math
\pi_{ip} = P(y_{ip} = 1 \mid \theta_p, \beta_i) =
\frac{\exp(\theta_p - \beta_i)}{1 + \exp(\theta_p - \beta_i)}.
```

The Rasch model embodies the principle of *specific objectivity*: the
difference between any two persons’ abilities can be estimated
independently of which items are administered, and conversely for item
difficulties.

In NIMBLE, this measurement model is expressed as:

``` r

for (j in 1:N) {
  for (i in 1:I) {
    y[j, i] ~ dbern(pi[j, i])
    logit(pi[j, i]) <- eta[j] - beta[i]
  }
}
```

### 2.2 Sufficient Statistics

A fundamental property of the Rasch model is that the total score
$`r_p = \sum_{i=1}^{I} y_{ip}`$ is a **sufficient statistic** for
$`\theta_p`$. The joint likelihood factorizes as:
``` math
P(\mathbf{y}_p \mid \theta_p, \boldsymbol{\beta}) =
\frac{\exp\left(r_p \theta_p - \sum_i y_{ip} \beta_i\right)}
     {\prod_{i=1}^{I} \left(1 + \exp(\theta_p - \beta_i)\right)},
```
which, by the Neyman-Fisher factorization theorem, establishes that the
posterior for $`\theta_p`$ depends on the data only through $`r_p`$.

> **Key insight:** The Rasch model is the *only* standard unidimensional
> IRT model with this sufficiency property. Two persons with the same
> total score have identical likelihoods for $`\theta`$.

### 2.3 Hierarchical Framework

The full Bayesian Rasch model has three levels:

**Level 1 (Measurement model):**
``` math
p(\mathbf{Y} \mid \boldsymbol{\theta}, \boldsymbol{\beta}) =
\prod_{p=1}^{N} \prod_{i=1}^{I}
\pi_{ip}^{y_{ip}} (1 - \pi_{ip})^{1 - y_{ip}},
```
where $`\pi_{ip}`$ is defined by the Rasch IRF above.

**Level 2 (Population model):** Person abilities are drawn from
$`\theta_p \stackrel{iid}{\sim} G`$ and item difficulties from
$`\beta_i \stackrel{iid}{\sim} \text{N}(\mu_\beta, \sigma^2_\beta)`$.
The choice of $`G`$ is where the two approaches diverge: under the
parametric prior $`G = \text{N}(\mu_\theta, \sigma^2_\theta)`$, and
under the DPM prior $`G \sim \text{DP}(\alpha, G_0)`$ (Section 5).

**Level 3 (Hyperpriors):** For the Normal prior,
$`\mu_\theta \sim \text{N}(0, 3)`$ and
$`\sigma^2_\theta \sim \text{Inv-Gamma}(2.01, 1.01)`$. The DPM
hyperprior structure is detailed in Section 5.

### 2.4 Posterior Mean and Shrinkage

Under a Normal population prior
$`G = \text{N}(\mu_\theta, \sigma^2_\theta)`$, the posterior mean
(expected a posteriori, or EAP) estimate of $`\theta_p`$ exhibits
**shrinkage** toward the population mean.

The EAP can be approximately characterized as a weighted combination of
the likelihood-based estimate and the prior mean:
``` math
\hat{\theta}_p^{\text{EAP}} \approx w_p \hat{\theta}_p^{\text{MLE}} +
(1 - w_p) \mu_\theta,
```
where $`\hat{\theta}_p^{\text{MLE}}`$ is the maximum likelihood estimate
and $`w_p`$ is the **shrinkage weight**:
``` math
w_p = \frac{\sigma^2_\theta}{\sigma^2_\theta + \text{se}(\hat{\theta}_p)^2}.
```

Here $`\text{se}(\hat{\theta}_p)^2 = 1 / I(\hat{\theta}_p)`$ is the
squared standard error. Key consequences: persons with extreme scores
are pulled more strongly toward $`\mu_\theta`$ (larger standard errors);
the weight $`w_p \in (0, 1)`$ is person-specific; and more items mean
less shrinkage.

> **Key insight:** Shrinkage produces **underdispersion** of the
> posterior means relative to the true population distribution. The
> variance of $`\{\hat{\theta}_p^{\text{EAP}}\}_{p=1}^N`$ is strictly
> less than the population variance $`\sigma^2_\theta`$. This
> underdispersion is the primary motivation for the CB and GR estimators
> (Section 7).

### 2.5 Test Reliability

The average shrinkage weight across persons provides a measure of test
reliability analogous to classical reliability coefficients:
``` math
\bar{w} = \frac{\sigma^2_\theta}{\sigma^2_\theta + \overline{\text{MSEM}}},
```
where $`\overline{\text{MSEM}}`$ is the mean squared error of
measurement averaged across persons.

The **Fisher information** for person $`p`$ under the Rasch model is:
``` math
I(\theta_p) = \sum_{i=1}^{I} P_i(\theta_p)(1 - P_i(\theta_p)),
```
where $`P_i(\theta_p) = \text{logistic}(\theta_p - \beta_i)`$ is the
probability of correct response. The information is maximized when
$`\theta_p = \beta_i`$ (i.e., when person ability matches item
difficulty) and equals $`I/4`$ when all items have the same difficulty
as the person’s ability.

This connects to the Rasch separation index:
``` math
G = \frac{\text{SD}_{\text{adjusted}}}{\text{RMSE}},
```
and the associated strata count $`H = (4G + 1)/3`$, which estimates the
number of statistically distinguishable ability levels the test can
discriminate.

## 3. The 2PL and 3PL Extensions

### 3.1 Two-Parameter Logistic Model (2PL IRT)

The 2PL model (Birnbaum, 1968) extends the Rasch model by introducing
item-specific discrimination parameters $`\lambda_i > 0`$:
``` math
\text{logit}(\pi_{ip}) = \lambda_i (\theta_p - \beta_i).
```

The discrimination parameter $`\lambda_i`$ controls the steepness of the
item characteristic curve (ICC). Items with higher $`\lambda_i`$ are
better at distinguishing between persons near the difficulty level
$`\beta_i`$.

In NIMBLE:

``` r

for (i in 1:I) {
  for (j in 1:N) {
    y[j, i] ~ dbern(pi[j, i])
    logit(pi[j, i]) <- lambda[i] * (eta[j] - beta[i])
  }
}
for (i in 1:I) {
  log(lambda[i]) ~ dnorm(0.5, var = 0.5)
  beta[i] ~ dnorm(0, var = sigma2_beta)
}
```

Note that $`\lambda_i`$ is constrained to be positive via the log-normal
prior: $`\log(\lambda_i) \sim \text{N}(0.5, 0.5)`$. This prior places
the median discrimination at $`\exp(0.5) \approx 1.65`$ with reasonable
spread over practically relevant values.

> **Key difference from Rasch:** The total score $`r_p`$ is no longer
> sufficient for $`\theta_p`$ in the 2PL model. The specific pattern of
> correct and incorrect responses matters because items contribute
> differentially to the likelihood.

### 3.2 Slope-Intercept (SI) Parameterization

An alternative parameterization of the 2PL replaces
$`(\lambda_i, \beta_i)`$ with $`(\lambda_i, \gamma_i)`$, where
$`\gamma_i`$ is the intercept:
``` math
\text{logit}(\pi_{ip}) = \gamma_i + \lambda_i \theta_p.
```

The two parameterizations are related by:
``` math
\gamma_i = -\lambda_i \beta_i \quad \Longleftrightarrow \quad
\beta_i = -\gamma_i / \lambda_i.
```

In NIMBLE:

``` r

for (i in 1:I) {
  for (j in 1:N) {
    y[j, i] ~ dbern(pi[j, i])
    logit(pi[j, i]) <- lambda[i] * eta[j] + gamma[i]
  }
}
```

The SI parameterization is sometimes preferred for computational reasons
because the linear predictor is additive in $`\theta_p`$, which can
improve mixing in certain MCMC samplers. DPMirt supports both
parameterizations.

### 3.3 Three-Parameter Logistic Model (3PL)

The 3PL model (Birnbaum, 1968) adds a lower asymptote parameter
$`\delta_i \in [0, 1]`$ representing the probability that a person with
very low ability still endorses the item (the “guessing” parameter):
``` math
\pi_{ip} = \delta_i + (1 - \delta_i) \cdot
\text{logistic}(\lambda_i (\theta_p - \beta_i)).
```

The item response probability approaches $`\delta_i`$ (rather than 0) as
$`\theta_p \to -\infty`$. The upper asymptote remains at 1.

In NIMBLE:

``` r

for (i in 1:I) {
  for (j in 1:N) {
    y[j, i] ~ dbern(pi[j, i])
    pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
    logit(linearReg[j, i]) <- lambda[i] * (eta[j] - beta[i])
  }
}
for (i in 1:I) {
  delta[i] ~ dbeta(4, 12)
}
```

The prior $`\delta_i \sim \text{Beta}(4, 12)`$ has mean $`4/16 = 0.25`$
and concentrates most mass between 0.05 and 0.45, reflecting the prior
belief that guessing parameters are typically modest. This is a standard
informative prior for the lower asymptote (De Ayala, 2022).

## 4. Identification and Rescaling

### 4.1 The Identification Problem

IRT models contain **identification indeterminacy**: certain
transformations of the parameters leave the likelihood invariant.
Without constraints, the posterior is improper along these directions.

**Rasch model (location indeterminacy):** For any constant $`c`$, the
transformation $`\theta_p \mapsto \theta_p + c`$,
$`\beta_i \mapsto \beta_i + c`$ leaves $`\theta_p - \beta_i`$ unchanged
and hence the likelihood invariant. Only one degree of freedom
(location) is indeterminate.

**2PL and 3PL models (location + scale indeterminacy):** For constants
$`c`$ and $`d > 0`$, the transformation
$`\theta_p \mapsto (\theta_p + c) / d`$,
$`\beta_i \mapsto (\beta_i + c) / d`$,
$`\lambda_i \mapsto \lambda_i \cdot d`$ leaves
$`\lambda_i(\theta_p - \beta_i)`$ invariant. Two degrees of freedom
(location and scale) are indeterminate.

The choice of identification strategy affects both the interpretation of
parameters and the efficiency of MCMC sampling.

### 4.2 Three Identification Strategies

DPMirt supports three identification strategies, selected via the
`identification` argument to
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md):

**Strategy 1: `constrained_ability`** – Fix
$`\theta \sim \text{N}(0, 1)`$. Resolves both location and scale
indeterminacy but constrains the latent distribution and is incompatible
with the DPM prior.

**Strategy 2: `constrained_item`** – Center item parameters during MCMC.
For Rasch:
``` math
\beta_i^* = \beta_i^{\text{tmp}} - \frac{1}{I}\sum_{j=1}^{I} \beta_j^{\text{tmp}}.
```

For 2PL/3PL, both $`\beta`$ (or $`\gamma`$) and $`\log(\lambda)`$ are
centered:
``` math
\beta_i^* = \beta_i^{\text{tmp}} - \bar{\beta}^{\text{tmp}}, \quad
\log(\lambda_i^*) = \log(\lambda_i^{\text{tmp}}) - \overline{\log \lambda}^{\text{tmp}}.
```

This approach identifies the model within the MCMC without constraining
the ability distribution, making it compatible with both Normal and DPM
priors.

**Strategy 3: `unconstrained`** – Place no identification constraints
during MCMC. Instead, apply post-hoc rescaling to the posterior samples
to achieve identification. This approach, advocated by Paganin et al.
(2023), has the advantage of not interfering with the sampler’s
geometry.

### 4.3 Post-hoc Rescaling Formulas

For **unconstrained** models, rescaling is applied
iteration-by-iteration to the MCMC output. Let superscript $`(s)`$
denote MCMC iteration $`s`$.

**Rasch model (location shift only):**
``` math
\beta_i^{*(s)} = \beta_i^{(s)} - \bar{\beta}^{(s)}, \quad
\theta_p^{*(s)} = \theta_p^{(s)} - \bar{\beta}^{(s)},
```
where $`\bar{\beta}^{(s)} = \frac{1}{I}\sum_{i=1}^{I} \beta_i^{(s)}`$ is
the mean item difficulty at iteration $`s`$.

**2PL/3PL IRT parameterization (location + scale):**
``` math
c^{(s)} = \bar{\beta}^{(s)}, \quad
d^{(s)} = \left(\prod_{i=1}^{I} \lambda_i^{(s)}\right)^{-1/I},
```
``` math
\beta_i^{*(s)} = \frac{\beta_i^{(s)} - c^{(s)}}{d^{(s)}}, \quad
\lambda_i^{*(s)} = \lambda_i^{(s)} \cdot d^{(s)}, \quad
\theta_p^{*(s)} = \frac{\theta_p^{(s)} - c^{(s)}}{d^{(s)}}.
```

The scale factor $`d^{(s)}`$ is the inverse geometric mean of the
discriminations, ensuring that the geometric mean of the rescaled
$`\lambda^*`$ equals 1.

**2PL/3PL SI parameterization:**
``` math
c^{(s)} = \frac{\sum_i \gamma_i^{(s)}}{\sum_i \lambda_i^{(s)}}, \quad
d^{(s)} = \left(\prod_{i=1}^{I} \lambda_i^{(s)}\right)^{-1/I},
```
``` math
\gamma_i^{*(s)} = \gamma_i^{(s)} - \lambda_i^{(s)} \cdot c^{(s)}, \quad
\lambda_i^{*(s)} = \lambda_i^{(s)} \cdot d^{(s)}, \quad
\theta_p^{*(s)} = \frac{\theta_p^{(s)} + c^{(s)}}{d^{(s)}}.
```

Note the sign difference for $`\theta`$ in the SI case: because
$`\gamma_i = -\lambda_i \beta_i`$, the location shift has the opposite
sign compared to the IRT parameterization.

> **Key insight:** The 3PL guessing parameter $`\delta_i`$ is
> scale-invariant and does not require rescaling. It enters the model
> multiplicatively outside the logistic function and is unaffected by
> affine transformations of the ability scale.

### 4.4 Efficiency of Unconstrained Sampling

Paganin et al. (2023) found that unconstrained sampling with post-hoc
rescaling is generally the most efficient strategy, because constraining
parameters during sampling can distort the posterior geometry. DPMirt
uses unconstrained identification as the default for 2PL and 3PL models,
and constrained_item as the default for Rasch (where the difference is
smaller).

## 5. Dirichlet Process Mixture Priors

### 5.1 The Dirichlet Process

The Dirichlet Process (DP; Ferguson, 1973) is a distribution over
probability distributions. A random measure
$`G \sim \text{DP}(\alpha, G_0)`$ satisfies: for every finite measurable
partition $`(A_1, \ldots, A_m)`$,
``` math
(G(A_1), \ldots, G(A_m)) \sim
\text{Dirichlet}(\alpha G_0(A_1), \ldots, \alpha G_0(A_m)).
```

Three properties are essential: (1) **Almost-sure discreteness** – draws
from the DP are discrete with probability one, inducing clustering. (2)
**Centering at $`G_0`$** – $`\mathbb{E}[G(A)] = G_0(A)`$ for any
measurable $`A`$. (3) **Concentration control** – as
$`\alpha \to \infty`$, $`G \to G_0`$; as $`\alpha \to 0`$, $`G`$
concentrates on a single atom.

### 5.2 The Chinese Restaurant Process

The Chinese Restaurant Process (CRP; Blackwell & MacQueen, 1973)
characterizes the partition structure induced by the DP. Customers
(persons) arrive sequentially at a restaurant with infinitely many
tables (clusters). Customer $`p`$ joins existing table $`k`$ with
probability $`n_k / (\alpha + p - 1)`$ or starts a new table with
probability $`\alpha / (\alpha + p - 1)`$, where $`n_k`$ is the current
table size.

> **Key insight:** The CRP exhibits a “rich get richer” property: tables
> with more customers are proportionally more likely to attract new
> customers. This produces a power-law distribution of cluster sizes and
> explains why DP mixtures naturally create a few large clusters and
> many small ones.

The CRP partition is **exchangeable** (independent of arrival order),
enabling Gibbs-based MCMC. The expected number of clusters is:
``` math
\mathbb{E}[K_N \mid \alpha] = \sum_{p=1}^{N} \frac{\alpha}{\alpha + p - 1}
\approx \alpha \log(N) \quad \text{for large } N.
```

This logarithmic growth is a reasonable modeling assumption for latent
trait distributions.

### 5.3 DPM-IRT Model Specification

The DPM-IRT model replaces the parametric Normal prior with a DP
Mixture. The measurement model is unchanged; the DPM prior specifies:
``` math
\theta_p \mid z_p, \tilde{\boldsymbol{\mu}}, \tilde{\boldsymbol{\sigma}}^2
\sim \text{N}(\tilde{\mu}_{z_p}, \tilde{\sigma}^2_{z_p}),
```
where $`z_p \in \{1, 2, \ldots\}`$ is the cluster assignment for person
$`p`$.

**CRP prior for cluster assignments:**
``` math
\mathbf{z} = (z_1, \ldots, z_N) \sim \text{CRP}(\alpha, N).
```

**Base measure for cluster parameters:**
``` math
G_0 = \text{N}(0, \sigma^2_\mu) \times \text{Inv-Gamma}(\nu_1, \nu_2),
```
meaning that for each cluster $`m`$:
``` math
\tilde{\mu}_m \sim \text{N}(0, \sigma^2_\mu), \quad
\tilde{\sigma}^2_m \sim \text{Inv-Gamma}(\nu_1, \nu_2).
```

**Default hyperparameters** (following Paganin et al., 2023):
``` math
\sigma^2_\mu = 2, \quad \nu_1 = 2.01, \quad \nu_2 = 1.01.
```

These values place the cluster variance prior just above the boundary
for finite mean ($`\nu_1 > 2`$) with
$`\mathbb{E}[\tilde{\sigma}^2] \approx 1`$.

In NIMBLE, the DPM prior is implemented via the CRP representation:

``` r

# CRP cluster assignments
zi[1:N] ~ dCRP(alpha, size = N)

# Concentration parameter hyperprior
alpha ~ dgamma(a, b)

# Person abilities: Normal kernel with cluster-specific parameters
for (j in 1:N) {
  eta[j] ~ dnorm(mu_j[j], var = s2_j[j])
  mu_j[j]  <- muTilde[zi[j]]
  s2_j[j]  <- s2Tilde[zi[j]]
}

# Cluster parameters drawn from base measure G0
for (m in 1:M) {
  muTilde[m]  ~ dnorm(0, var = s2_mu)
  s2Tilde[m]  ~ dinvgamma(nu1, nu2)
}
```

### 5.4 Stick-Breaking vs. CRP

An alternative to the CRP is Sethuraman’s (1994) **stick-breaking**
construction, which represents $`G`$ as an infinite discrete measure:
``` math
G = \sum_{h=1}^{\infty} w_h \delta_{\theta^*_h}, \quad
w_h = v_h \prod_{\ell < h}(1 - v_\ell), \quad
v_h \stackrel{iid}{\sim} \text{Beta}(1, \alpha),
```
where $`\{\theta^*_h\} \stackrel{iid}{\sim} G_0`$. DPMirt uses the **CRP
representation** instead, because it is natively supported by NIMBLE’s
`dCRP` distribution and offers direct access to cluster assignments
$`z_p`$ without requiring explicit truncation of the number of
components.

### 5.5 Practical Truncation

Although the CRP is defined over infinitely many potential clusters,
NIMBLE requires a finite upper bound $`M`$ on pre-allocated cluster
parameter vectors. The default in DPMirt is $`M = 50`$. Since
$`\mathbb{E}[K_N \mid \alpha] \approx \alpha \log(N)`$ and typical
applications have $`\alpha \in [0.5, 5]`$ and $`N \leq 2{,}000`$, the
default $`M = 50`$ is conservative. If the MCMC sampler creates clusters
approaching $`M`$, the user should increase it via the `M` argument to
[`dpmirt_spec()`](https://joonho112.github.io/DPMirt/reference/dpmirt_spec.md).

## 6. The Concentration Parameter ($`\alpha`$)

### 6.1 Role of $`\alpha`$

The concentration parameter $`\alpha`$ controls (1) the expected number
of clusters $`\mathbb{E}[K_N \mid \alpha]`$, (2) the degree of departure
from $`G_0`$, and (3) the resolution of density estimation. For large
$`N`$:
``` math
\mathbb{E}[K_N \mid \alpha] = \sum_{p=1}^{N} \frac{\alpha}{\alpha + p - 1}
\approx \alpha \log(N).
```

Illustrative values for $`N = 500`$:

| $`\alpha`$ | $`\mathbb{E}[K_N]`$ |   Interpretation    |
|:----------:|:-------------------:|:-------------------:|
|    0.5     |    $`\approx 3`$    | Few broad clusters  |
|    1.0     |    $`\approx 6`$    | Moderate clustering |
|    3.0     |   $`\approx 19`$    |    Many clusters    |
|    10.0    |   $`\approx 62`$    |   Near-continuous   |

### 6.2 Gamma Hyperprior

DPMirt places $`\alpha \sim \text{Gamma}(a, b)`$ (shape-rate
parameterization, $`\mathbb{E}[\alpha] = a/b`$). The default is
$`\alpha \sim \text{Gamma}(1, 3)`$ following Paganin et al. (2023), with
$`\mathbb{E}[\alpha] = 1/3`$. This conservative prior favors
parsimonious models with few clusters, allowing data to drive the
discovery of multiple subpopulations.

### 6.3 Principled Elicitation via DPprior

The DPprior package (Lee, 2026) translates beliefs about $`\mu_K =
\mathbb{E}[K_N]`$ into Gamma parameters $`(a, b)`$ via moment-matching.
DPMirt integrates through
[`dpmirt_alpha_prior()`](https://joonho112.github.io/DPMirt/reference/dpmirt_alpha_prior.md):

``` r

alpha_ab <- dpmirt_alpha_prior(N = 500, mu_K = 5, confidence = "medium")
spec <- dpmirt_spec(data, model = "rasch", prior = "dpm",
                    alpha_prior = alpha_ab)
```

If DPprior is not installed, DPMirt falls back to Gamma(1, 3).

> **Key insight:** The $`\alpha`$ hyperprior meaningfully impacts
> inference for moderate sample sizes. Confirmatory applications benefit
> from principled elicitation encoding domain knowledge about population
> heterogeneity.

## 7. Posterior Summary Theory

A central theme of the DPMirt package is that the choice of posterior
summary method should be guided by the inferential goal. Different loss
functions lead to different optimal estimators, and no single estimator
is universally best.

### 7.1 Loss Functions and Optimality

Consider $`N`$ persons with true abilities
$`\theta_1, \ldots, \theta_N`$ and estimates
$`\hat{\theta}_1, \ldots, \hat{\theta}_N`$. We define three loss
functions corresponding to three inferential goals.

**Goal 1: Individual estimation accuracy.**
``` math
\text{MSEL} = \frac{1}{N} \sum_{p=1}^{N}
(\hat{\theta}_p - \theta_p)^2.
```

The posterior mean (PM) minimizes MSEL in expectation under the
posterior. This is the classical Bayes estimator under squared error
loss.

**Goal 2: Ranking quality.**
``` math
\text{MSELP} = \frac{1}{N} \sum_{p=1}^{N}
\left(\frac{\hat{R}_p}{N} - \frac{R_p}{N}\right)^2,
```
where $`\hat{R}_p = \text{rank}(\hat{\theta}_p)`$ and
$`R_p = \text{rank}(\theta_p)`$ are the estimated and true ranks,
respectively. This loss measures how well the estimator preserves the
ordering of persons.

**Goal 3: Distribution recovery.**
``` math
\text{KS} = \sup_{t} \left|\hat{G}_N(t) - G_N(t)\right|,
```
where
$`\hat{G}_N(t) = \frac{1}{N}\sum_p \mathbf{1}(\hat{\theta}_p \leq t)`$
and $`G_N(t) = \frac{1}{N}\sum_p \mathbf{1}(\theta_p \leq t)`$ are the
empirical distribution functions of the estimates and true values. This
loss measures how well the set of estimates reproduces the shape of the
true ability distribution.

> **Key insight:** The PM estimator is optimal for Goal 1 but performs
> poorly on Goals 2 and 3 due to shrinkage-induced underdispersion. The
> CB estimator corrects the underdispersion (Goal 3), and the GR
> estimator simultaneously addresses ranking and distribution recovery
> (Goals 2 and 3).

### 7.2 Posterior Mean (PM)

The posterior mean is the most familiar Bayesian point estimate:
``` math
\hat{\theta}_p^{\text{PM}} = \mathbb{E}[\theta_p \mid \mathbf{Y}] =
\frac{1}{S} \sum_{s=1}^{S} \theta_p^{(s)},
```
where $`\theta_p^{(s)}`$ is the value of $`\theta_p`$ at MCMC iteration
$`s`$ and $`S`$ is the total number of post-burnin samples.

The PM minimizes the posterior expected squared error loss for each
person individually. However, the set of posterior means
$`\{\hat{\theta}_p^{\text{PM}}\}_{p=1}^N`$ is **underdispersed**: their
empirical variance is strictly less than the posterior expectation of
the population variance. This underdispersion distorts the shape of the
estimated ability distribution and can bias percentile-based inferences.

### 7.3 Constrained Bayes (CB)

The Constrained Bayes estimator, introduced by Ghosh (1992) and
developed further by Louis (1984) and Ghosh and Kim (2002), addresses
the underdispersion of the posterior mean by imposing two constraints:

1.  **Match the marginal mean:**
    $`\frac{1}{N}\sum_p \hat{\theta}_p^{\text{CB}}
    = \bar{\eta}`$ (same as PM).
2.  **Match the marginal variance:**
    $`\text{Var}(\hat{\theta}_p^{\text{CB}})
    = V_\eta + \bar{\lambda}`$ (match the posterior expected population
    variance).

These constraints yield the following closed-form estimator:
``` math
\hat{\theta}_p^{\text{CB}} = \bar{\eta} +
(\eta_p - \bar{\eta}) \sqrt{1 + \frac{\bar{\lambda}}{V_\eta}},
```
where:

- $`\eta_p = \hat{\theta}_p^{\text{PM}}`$ is the posterior mean for
  person $`p`$,
- $`\bar{\eta} = \frac{1}{N}\sum_p \eta_p`$ is the grand mean of
  posterior means,
- $`\bar{\lambda} = \frac{1}{N}\sum_p \text{Var}(\theta_p \mid \mathbf{Y})`$
  is the mean posterior variance,
- $`V_\eta = \text{Var}(\eta_1, \ldots, \eta_N)`$ is the sample variance
  of the posterior means (computed with R’s
  [`var()`](https://rdrr.io/r/stats/cor.html), i.e., $`N-1`$
  denominator).

The CB estimator inflates the deviations of each posterior mean from the
grand mean by the factor $`\sqrt{1 + \bar{\lambda}/V_\eta}`$. This
factor is always greater than 1 (since $`\bar{\lambda} > 0`$), so the CB
estimates are more dispersed than the posterior means.

> **Interpretation of the CB scaling factor:** The ratio
> $`\bar{\lambda}/V_\eta`$ measures the “severity of shrinkage.” When
> individual posterior variances $`\bar{\lambda}`$ are large relative to
> the between-person variation $`V_\eta`$ (i.e., the test is
> unreliable), the scaling factor is large and CB makes substantial
> corrections. When the test is highly reliable
> ($`\bar{\lambda} \ll V_\eta`$), CB and PM estimates are nearly
> identical.

### 7.4 Triple-Goal (GR)

The Triple-Goal estimator, introduced by Shen and Louis (1998),
simultaneously addresses all three goals: individual estimation,
ranking, and distribution recovery. It is the most computationally
intensive of the three methods but provides the most balanced trade-off.

The GR algorithm proceeds in four steps:

**Step 1: Posterior mean ranks.** For each MCMC iteration $`s`$, rank
all $`N`$ persons by their $`\theta^{(s)}`$ values. Average these ranks
across iterations:
``` math
\bar{R}_p = \frac{1}{S} \sum_{s=1}^{S} R_p^{(s)}, \quad
\text{where } R_p^{(s)} = \text{rank}(\theta_p^{(s)})
\text{ among } \theta_1^{(s)}, \ldots, \theta_N^{(s)}.
```

Equivalently, $`\bar{R}_p`$ can be computed as:
``` math
\bar{R}_p = \sum_{q=1}^{N} P(\theta_q \leq \theta_p \mid \mathbf{Y}),
```
which is the expected rank of person $`p`$ under the posterior.

**Step 2: Integer ranks.** Convert the (generally non-integer) posterior
mean ranks to integer ranks:
``` math
\hat{R}_p = \text{rank}(\bar{R}_p) \in \{1, 2, \ldots, N\}.
```

Ties in $`\bar{R}_p`$ are broken randomly.

**Step 3: ISEL empirical distribution function.** Compute the
**integrated squared error loss** (ISEL) optimal estimator of the
population EDF:
``` math
\hat{G}_N(t) = \frac{1}{N} \sum_{p=1}^{N} P(\theta_p \leq t \mid \mathbf{Y}).
```

This is the posterior expectation of the true EDF and is implemented by
pooling all posterior samples across all persons.

**Step 4: Quantile inversion.** Assign each person the value from
$`\hat{G}_N`$ corresponding to their integer rank:
``` math
\hat{\theta}_p^{\text{GR}} = \hat{G}_N^{-1}\left(\frac{2\hat{R}_p - 1}{2N}\right).
```

The fraction $`(2\hat{R}_p - 1)/(2N)`$ maps rank $`\hat{R}_p`$ to the
midpoint of the corresponding probability bin, preventing boundary
artifacts.

**Computational complexity.** Step 1 requires $`O(S \times N \log N)`$
operations. For $`S = 5{,}000`$ and $`N = 500`$, this runs in seconds.
In DPMirt, the GR estimator is implemented in `.triple_goal()`, adapted
from the HETOP package:

``` r

# Step 1: Posterior mean ranks
rbar <- apply(
  t(apply(s, 1, rank, ties.method = "average")),
  2, mean
)

# Step 2: Integer ranks
rhat <- rank(rbar, ties.method = "random")

# Steps 3-4: ISEL EDF + quantile inversion
theta_gr <- quantile(
  c(s),                            # Pool all posterior samples
  probs = (2 * rhat - 1) / (2 * K),
  type = quantile_type
)
```

### 7.5 Comparison of Estimators

The following table summarizes the properties of the three estimators:

| Property                               |    PM     |    CB     |        GR        |
|:---------------------------------------|:---------:|:---------:|:----------------:|
| Optimal for individual MSE (Goal 1)    |    Yes    |    No     |        No        |
| Preserves population variance (Goal 3) |    No     |    Yes    |       Yes        |
| Preserves ranks (Goal 2)               |  Approx.  |  Approx.  |       Yes        |
| Closed-form from posterior samples     |    Yes    |    Yes    |       Yes        |
| Computational cost                     | $`O(SN)`$ | $`O(SN)`$ | $`O(SN \log N)`$ |

> **Practical guidance:** For educational testing applications where
> percentile-based interpretations are important (e.g., reporting a
> student’s standing relative to the population), GR estimates are
> preferred. For applications focused on individual point estimates
> (e.g., adaptive testing), PM estimates are optimal. CB provides a
> useful middle ground with the computational simplicity of PM.

## 8. DP Density Reconstruction

### 8.1 Mixture Density per Iteration

At each MCMC iteration $`s`$, the DP posterior induces a mixture
density:
``` math
f_s(x) = \sum_{k=1}^{K_s} w_{s,k} \,
\phi(x; \tilde{\mu}_{s,k}, \tilde{\sigma}^2_{s,k}),
```
where $`K_s`$ is the number of occupied clusters, $`w_{s,k}`$ is the
cluster weight, and $`\phi(\cdot; \mu, \sigma^2)`$ denotes the Normal
density.

### 8.2 Posterior Predictive Density and Credible Bands

The posterior predictive density averages over MCMC iterations:
``` math
\hat{f}(x) = \frac{1}{S} \sum_{s=1}^{S} f_s(x),
```
integrating over all sources of uncertainty. Pointwise credible bands
are obtained by taking quantiles of $`\{f_s(x)\}_{s=1}^S`$ at each grid
point:
``` math
\hat{f}^{\text{lower}}(x) = Q_{0.025}(\{f_s(x)\}), \quad
\hat{f}^{\text{upper}}(x) = Q_{0.975}(\{f_s(x)\}).
```

> **Caveat:** These are *pointwise* credible bands, not simultaneous
> bands. The probability that the true density lies within the band at
> every point simultaneously is less than the nominal level.

### 8.3 Implementation via NIMBLE

DPMirt implements density reconstruction through the
[`dpmirt_dp_density()`](https://joonho112.github.io/DPMirt/reference/dpmirt_dp_density.md)
function, following Paganin et al.’s (2023) approach: extract DP
posterior samples, call NIMBLE’s
[`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
to obtain stick-breaking weights and atoms, evaluate $`f_s(x)`$ on a
grid, and apply rescaling.

For unconstrained Rasch models, the rescaling is a location shift
(Jacobian = 1):
``` math
f_s^{\text{rescaled}}(x) = f_s(x + \bar{\beta}^{(s)}).
```

For 2PL/3PL models, both location and scale adjustments are required:
``` math
f_s^{\text{rescaled}}(x) = d^{(s)} \cdot f_s(d^{(s)} x + c^{(s)}),
```
where $`c^{(s)}`$ is the location shift and $`d^{(s)}`$ is the scale
factor.

## 9. Attribution Table

The following table clarifies the provenance of each major component
used in the DPMirt package.

| Component | Source | Status |
|:---|:---|:---|
| Rasch model | Rasch (1960); De Ayala (2022) | Established |
| 2PL / 3PL models | Birnbaum (1968); De Ayala (2022) | Established |
| Dirichlet Process | Ferguson (1973) | Established |
| Chinese Restaurant Process | Blackwell & MacQueen (1973) | Established |
| Stick-breaking construction | Sethuraman (1994) | Established |
| DPM-IRT models + NIMBLE code | Paganin et al. (2023) | Established (adapted) |
| Post-hoc rescaling | Paganin et al. (2023) | Established (adapted) |
| Posterior Mean (PM) | Standard Bayes (Berger, 1985) | Established |
| Constrained Bayes (CB) | Ghosh (1992); Louis (1984) | Established |
| Triple-Goal (GR) | Shen & Louis (1998) | Established |
| Combined model-estimator framework | Lee & Wind (APM) | Novel integration |
| Reliability-targeted simulation | Lee (IRTsimrel) | Novel contribution |
| Alpha elicitation (DPprior) | Lee (2026, DPprior) | Novel contribution |

## 10. References

Antoniak, C. E. (1974). Mixtures of Dirichlet processes with
applications to Bayesian nonparametric problems. *The Annals of
Statistics*, 2(6), 1152–1174.

Berger, J. O. (1985). *Statistical Decision Theory and Bayesian
Analysis* (2nd ed.). Springer.

Birnbaum, A. (1968). Some latent trait models and their use in inferring
an examinee’s ability. In F. M. Lord & M. R. Novick (Eds.), *Statistical
Theories of Mental Test Scores* (pp. 397–479). Addison-Wesley.

Blackwell, D., & MacQueen, J. B. (1973). Ferguson distributions via
Polya urn schemes. *The Annals of Statistics*, 1(2), 353–355.

De Ayala, R. J. (2022). *The Theory and Practice of Item Response
Theory* (2nd ed.). Guilford Press.

Ferguson, T. S. (1973). A Bayesian analysis of some nonparametric
problems. *The Annals of Statistics*, 1(2), 209–230.

Fox, J.-P. (2010). *Bayesian Item Response Modeling: Theory and
Applications*. Springer.

Ghosh, M. (1992). Constrained Bayes estimation with applications.
*Journal of the American Statistical Association*, 87(418), 533–540.

Ghosh, M., & Kim, M.-H. (2002). The Bayes and constrained Bayes
estimators under balanced loss. *Statistics & Probability Letters*,
59(2), 175–183.

Lee, J. & Wind, S. Targeting toward inferential goals in Bayesian Rasch
models for estimating person-specific latent traits. *OSF Preprint*.
<https://doi.org/10.31219/osf.io/qrw4n>

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. arXiv preprint arXiv:2512.16012.
<https://arxiv.org/abs/2512.16012>

Lee, J. (2026). Design-conditional prior elicitation for Dirichlet
Process mixtures: A unified framework for cluster counts and weight
control. arXiv preprint arXiv:2602.06301.
<https://arxiv.org/abs/2602.06301>

Louis, T. A. (1984). Estimating a population of parameter values using
Bayes and empirical Bayes methods. *Journal of the American Statistical
Association*, 79(386), 393–398.

Paganin, S., Paciorek, C. J., Wehrhahn, C., Rodríguez, A., Rabe-Hesketh,
S., & de Valpine, P. (2023). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics*, 48(2), 147–188.
<https://doi.org/10.3102/10769986221136105>

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Danish Institute for Educational Research.

Sethuraman, J. (1994). A constructive definition of Dirichlet priors.
*Statistica Sinica*, 4(2), 639–650.

Shen, W., & Louis, T. A. (1998). Triple-goal estimates in two-stage
hierarchical models. *Journal of the Royal Statistical Society: Series B
(Statistical Methodology)*, 60(2), 455–471.

Wright, B. D., & Masters, G. N. (1982). *Rating Scale Analysis: Rasch
Measurement*. MESA Press.

------------------------------------------------------------------------

## What’s Next?

This vignette covered the mathematical foundations. For practical
guidance on using these models and methods, see:

- **[`vignette("models-and-workflow", package = "DPMirt")`](https://joonho112.github.io/DPMirt/articles/models-and-workflow.md)**
  – Step-by-step guide to fitting Rasch, 2PL, and 3PL models with both
  Normal and DPM priors. Covers the compile-once, sample-many workflow.

- **[`vignette("posterior-summaries", package = "DPMirt")`](https://joonho112.github.io/DPMirt/articles/posterior-summaries.md)**
  – Practical comparison of PM, CB, and GR estimators with real data
  examples. Demonstrates when each estimator is preferred.

- **[`vignette("nimble-internals", package = "DPMirt")`](https://joonho112.github.io/DPMirt/articles/nimble-internals.md)**
  – Under-the-hood details of NIMBLE integration, custom samplers, and
  the
  [`getSamplesDPmeasure()`](https://rdrr.io/pkg/nimble/man/getSamplesDPmeasure.html)
  workflow for DP density reconstruction.

------------------------------------------------------------------------

*For questions about this vignette or the DPMirt package, please visit
the [GitHub repository](https://github.com/joonho112/DPMirt).*
