# Beta-binomial regression in Shiba

!!! Under development

	This document is currently under development and may contain inaccuracies or incomplete information. Please refer to the latest version on GitHub for updates.

This document summarizes the statistical theory and the current implementation of beta-binomial regression used in Shiba.

The implementation targets event-wise differential splicing testing between two groups, using read counts directly.

## Why beta-binomial

PSI is a ratio derived from count data, and per-sample uncertainty depends strongly on read depth.

If we model PSI alone, depth information is partially lost. Beta-binomial modeling keeps both:

- success counts $k_i$ (inclusion-side reads)
- total counts $n_i$ (effective total reads)

for sample $i$.

This allows overdispersion beyond a simple binomial model while preserving count-scale uncertainty.

## Data structure for one event

For one splicing event and one sample $i$:

- $x_i \in \{0,1\}$: group indicator (reference $=0$, alternative $=1$)
- $k_i$: success reads
- $n_i$: total reads, with $0 \le k_i \le n_i$

The model is fitted per event using all valid samples from both groups.

### Event-wise definition of success and total counts

Shiba computes event-specific $k_i$ and $n_i$ from junction counts as follows.

Let $a,b,c$ denote event-specific junction counts in each sample.

1. FIVE / THREE:

$$
k = a, \quad n = a + b
$$

2. AFE / ALE:

$$
k = \sum a_j, \quad n = \sum a_j + \sum b_j
$$

3. MXE:

$$
k = a_1 + a_2, \quad n = a_1 + a_2 + b_1 + b_2
$$

4. MSE:

$$
k = \sum a_j, \quad n = \sum a_j + c
$$

where $a_j$ are inclusion junctions and $c$ is exclusion junction count.

5. SE:

$$
k = a + b, \quad n = a + b + 2c
$$

This matches the PSI structure

$$
\psi = \frac{(a+b)/2}{(a+b)/2 + c} = \frac{a+b}{a+b+2c} = \frac{k}{n}.
$$

6. RI:

$$
k = s + e, \quad n = s + e + 2a
$$

where $s,e$ are splice-junction counts at intron start/end, and $a$ is intron-retention junction count.

This also preserves

$$
\psi = \frac{(s+e)/2}{(s+e)/2 + a} = \frac{s+e}{s+e+2a} = \frac{k}{n}.
$$

## Statistical model

For sample $i$:

$$
K_i \mid p_i \sim \mathrm{Binomial}(n_i, p_i), \quad p_i \sim \mathrm{Beta}(\alpha_i, \beta_i).
$$

Marginally:

$$
K_i \sim \mathrm{BetaBinomial}(n_i, \alpha_i, \beta_i).
$$

The PMF is

$$
\Pr(K_i = k_i)
= \binom{n_i}{k_i}
\frac{B(k_i + \alpha_i,\; n_i-k_i+\beta_i)}{B(\alpha_i,\beta_i)}.
$$

### Mean and overdispersion parameterization

Shiba uses a mean-dispersion parameterization:

$$
\mu_i = \mathbb{E}[p_i], \quad \rho \in (0,1).
$$

Regression part:

$$
\mathrm{logit}(\mu_i) = \beta_0 + \beta_1 x_i.
$$

Dispersion reparameterization:

$$
\rho = \sigma(\theta) = \frac{1}{1+e^{-\theta}},
$$

$$
\phi = \frac{1-\rho}{\rho},
$$

$$
\alpha_i = \mu_i\phi, \quad \beta_i = (1-\mu_i)\phi.
$$

When $\rho \to 0$, the model approaches binomial-like behavior (low extra-binomial variation).

### What is $\rho$? (core intuition)

In this model, $\rho$ controls **extra-binomial variability**.

If we had a plain binomial model, the success probability is fixed within sample $i$:

$$
K_i \sim \mathrm{Binomial}(n_i, \mu_i), \quad
\mathrm{Var}(K_i) = n_i\mu_i(1-\mu_i).
$$

In beta-binomial, the probability itself fluctuates across biological/technical contexts:

$$
p_i \sim \mathrm{Beta}(\alpha_i,\beta_i), \quad K_i\mid p_i \sim \mathrm{Binomial}(n_i,p_i).
$$

After marginalization, variance is inflated to

$$
\mathrm{Var}(K_i)
= n_i\mu_i(1-\mu_i)\{1 + (n_i-1)\rho\}.
$$

So $\rho$ is the inflation driver:

1. $\rho=0$: no inflation, exactly binomial variance.
2. $\rho>0$: overdispersion; counts are more spread than binomial.
3. Larger $n_i$ amplifies the impact via $(n_i-1)\rho$.

### Equivalent interpretation as intra-class correlation

$\rho$ can also be interpreted as the correlation between two Bernoulli trials within the same sample/event under the beta-binomial hierarchy.

Roughly speaking:

- small $\rho$: reads are close to conditionally independent given $\mu_i$
- large $\rho$: reads co-move more strongly (latent heterogeneity is strong)

This is why beta-binomial is often preferable for RNA splicing counts: technical and biological heterogeneity creates positive correlation and overdispersion that binomial cannot capture.

### Relationship between $\rho$ and $\phi$

Shiba uses

$$
\phi = \frac{1-\rho}{\rho} \quad \Longleftrightarrow \quad \rho = \frac{1}{\phi+1}.
$$

Hence:

1. Large $\phi$ means concentrated Beta distribution around $\mu_i$ and small overdispersion ($\rho \approx 0$).
2. Small $\phi$ means broad Beta distribution and strong overdispersion (large $\rho$).

### Why optimize $\theta$ instead of $\rho$ directly

Implementation uses

$$
\rho = \sigma(\theta) = \frac{1}{1+e^{-\theta}},
$$

so unconstrained $\theta\in\mathbb{R}$ maps smoothly into valid $\rho\in(0,1)$.
This avoids invalid updates during numerical optimization.

## Likelihood used in implementation

For one event, with samples $i=1,\dots,m$, the log-likelihood is

$$
\ell(\Theta)
= \sum_{i=1}^{m}
\left[
\log\binom{n_i}{k_i}
+ \log B(k_i+\alpha_i, n_i-k_i+\beta_i)
- \log B(\alpha_i,\beta_i)
\right].
$$

In code, this is computed via gamma/beta log functions:

$$
\log\binom{n}{k}
= \log\Gamma(n+1)-\log\Gamma(k+1)-\log\Gamma(n-k+1).
$$

Two models are fitted by numerical optimization:

1. Null model $H_0$: $\beta_1 = 0$
2. Full model $H_1$: $\beta_1$ free

Negative log-likelihoods are minimized with L-BFGS-B.

## Hypothesis test

Primary hypothesis:

$$
H_0: \beta_1 = 0
\quad\text{vs}\quad
H_1: \beta_1 \ne 0.
$$

Likelihood-ratio statistic:

$$
\Lambda = 2\{\ell(\widehat\Theta_1) - \ell(\widehat\Theta_0)\}.
$$

Implementation uses equivalent form with minimized negative log-likelihoods:

$$
\Lambda = 2\{\mathrm{NLL}_{0} - \mathrm{NLL}_{1}\}, \quad \Lambda \leftarrow \max(\Lambda,0).
$$

P-value per event:

$$
p_{\beta} = \Pr\left(\chi^2_{(1)} \ge \Lambda\right).
$$

Then BH correction is applied across valid events to create $q_{\beta}$.

## Practical safeguards in Shiba implementation

The implementation includes several robustification steps.

1. Sample count requirement:

- if either group has fewer than 2 valid samples, result is `NaN`.

2. Early exits for near-constant events:

- if variance of $k_i/n_i$ is extremely small, return $\Lambda=0$.
- if group means are almost identical, return $\Lambda=0$.

3. Parameter clipping and bounds:

- $\mu$ clipped to $(10^{-10}, 1-10^{-10})$
- $\rho$ clipped to $(10^{-8}, 1-10^{-8})$
- optimization bounds:
  - null: $\beta_0\in[-20,20],\;\theta\in[-10,10]$
  - full: $\beta_0\in[-20,20],\;\beta_1\in[-20,20],\;\theta\in[-10,10]$

4. Numerical validity checks:

- if optimizer returns non-finite objective values, event result is `NaN`.

5. P-value underflow handling:

- when $p_{\beta}$ underflows to zero, it is replaced with `np.finfo(float).tiny`.

## Relationship to PSI output

Important distinction:

- Statistical inference now uses direct count inputs $(k,n)$.
- PSI is still reported as an interpretable summary metric.

So PSI remains useful for biological interpretation, while hypothesis testing is count-aware and depth-aware.

## End-to-end flow in Shiba

For each event type:

1. Compute group-level PSI and detect candidate differential events.
2. For candidate events, compute sample-level PSI and count columns (`_success`, `_total_reads`).
3. Run beta-binomial LRT to get `p_beta`.
4. Apply BH correction to get `q_beta`.
