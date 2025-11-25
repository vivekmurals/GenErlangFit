# FAQ

## Fitting Algorithms and Assumptions

This package computes the maximum likelihood estimates (MLE) of an
Erlang distribution (K, λ) and/or MLE of an Erlang-Exponential
distribution (K, λ_(Erlang), λ_(Exponential)) from observed data.  
\## Erlang Fit Algorithm

The algorithm calculates the maximum likelihood estimates (MLE) of an
Erlang distribution from observed data. Since the Erlang shape parameter
K must be an integer, direct optimization over K using common
optimization libraries is non-trivial. To address this, we first
minimize the negative log-likelihood of the observed data under a Gamma
distribution with respect to its continuous shape parameter K, using
numerical optimization. This provides a candidate value of K which can
be adjusted to nearest integers to identify the maximum likelihood
estimates for the Erlang distribution. The exact details of this are
explained below.

This search starts from an initial guess for K specified by
`InitialguessK`. The cost function to be minimized with respect to K is
the negative log likelihood of the observed data given a Gamma
distribution:

``` math

- \log \left( \prod_{i=1}^N \frac{\lambda^K x_i^{K - 1} e^{-\lambda x_i}}{(K - 1)!} \right)
```

In the function above, $`x_i`$ refers to the observed data, N refers to
the size of the observed data while K and λ are the shape and scale
parameters of the Gamma distribution. The corresponding scale parameter,
λ used in the evaluation of the negative log likelihood is calculated
from the initial guess K via the following relationship, where μ is the
mean of the observed data:

``` math

\frac{K}{\lambda} = \mu
```

This relationship exists at the MLE of a Gamma distribution (and by
extension an Erlang). Specifically, this relationship arises from the
Gamma MLE conditions: differentiating the log-likelihood with respect to
λ and setting that to 0 yields λ = K / μ. We leverage this relationship
having to only search over one variable K.

After finding the optimal continuous K, the floor ⌊K⌋ and ceiling ⌈K⌉ of
this value are evaluated as candidate integer K values for finding the
MLE of the Erlang Distribution. For each candidate K, the respective λ
is computed from λ = K / μ which holds at the MLE, and the
log-likelihood of the observed data given an Erlang distribution is
computed. The pair of integer K\* and λ\* that gives the higher
log-likelihood is selected as the MLE for the Erlang Distribution.

In the special case where ⌈K⌉ = 0, the candidate K = 0 is not admissible
in an Erlang Distribution. Therefore, K = 1 and λ = 1 / μ will be the
MLE for the distribution. This also corresponds to an exponential
distribution with mean and variance λ.

Assumptions: The likelihood of the observed data given a Gamma
distribution forms a 2D surface over continuous variables K and λ. At
the MLE, these parameters satisfy K / λ = μ, reducing our search space
to a 1D line through this surface. Once the optimal continuous (K, λ)
for the Gamma is found, the floor and ceiling of K are selected as the
two corresponding Erlang candidates.

This method assumes: (1) the maximum likelihood found for the continuous
Gamma distribution along this line is the true global maximum, and (2)
the likelihood function along the segment between floor ⌊K⌋ and ceiling
⌈K⌉ is monotonic and continuous, allowing us to identify the MLE of an
Erlang distribution by checking only these two integer neighbors.

### Erlang-Exp Fit Algorithm

#### What is an Erlang-Exponential Distribution?

The Erlang–Exponential distribution describes the sum of 2 independent
random variables: an erlang random variable with shape parameter K and
rate parameter, $`\lambda_{\text{Erlang}}`$ and an exponential random
variable with rate parameter $`\lambda_{\text{Exponential}}`$.
Essentially, this distribution generalizes the Erlang distribution by
adding an extra exponential stage, which can be either faster or slower
than the Erlang stages. This allows for better approximation of
empirical distributions with varying skewness and tail behavior.

#### Fit Algorithm

By default, this algorithm performs an adaptive search across a range of
K-values centered around the input value of K to best identify K – along
with the corresponding $`\lambda_{\text{Erlang}}`$ and
$`\lambda_{\text{Exponential}}`$ – that maximizes the likelihood of the
observed data. If ‘FixedK’, true is specified, the algorithm instead
computes the MLEs ($`\lambda_{\text{Erlang}}`$ and
$`\lambda_{\text{Exponential}}`$) for a fixed value of K (specified
input value of K).

This section details how the algorithm computes the MLE for a fixed
input argument K. This is the backbone of the overall procedure, since
the adaptive search described earlier simply repeats this fixed-K
estimation across multiple candidate K values.

The algorithm calculates the MLE ($`\lambda_{\text{Erlang}}`$ &
$`\lambda_{\text{Exponential}}`$) of an Erlang-Exponential distribution
from observed data for a given input shape parameter K through
optimizing a cost function that minimizes the negative log likelihood of
the observed data given an Erlang-Exponential distribution. This
optimization starts from an initial guess for
$`\lambda_{\text{Erlang}}`$ specified by `InitialguessErLam`.

The negative log likelihood itself is computed using a custom function
which outputs both the pointwise PDF of the Erlang–Exponential
distribution and the log-likelihood of the observed dataset under the
model.

Note: This custom function performs direct numerical integration of the
convolution of an Erlang distribution with an Exponential distribution.
The expression for this convolution is as follows:

``` math

\frac{\lambda_{Er}^k \lambda_{Exp}}{(k - 1)!} e^{-\lambda_{Exp} z} \int_0^z x^{k - 1} e^{-(\lambda_{Er} - \lambda_{Exp}) x} \, dx
```

In computing the negative log likelihood, the corresponding
$`\lambda_{\text{Exponential}}`$ is calculated via the following
relationship, which exists at the MLE of this function:

``` math

\frac{K}{\lambda_{Erlang}} + \frac{1}{\lambda_{Exponential}} = \mu
```

The pair of $`\lambda_{\text{Erlang}}`$ &
$`\lambda_{\text{Exponential}}`$ that maximizes the log-likelihood is
selected as the MLE for the Erlang-Exponential Distribution from the
observed data.

------------------------------------------------------------------------

## Goodness of Fit Evaluation

### Bootstrap Hypothesis Testing

Assessment of the goodness-of-fit of the estimated Erlang or Erlang-Exp
model is evaluated through parametric bootstrap hypothesis testing.

The hypotheses are:

- **$`H_0`$:** The observed data could have come from the fitted
  distribution.
- **$`H_1`$:** The observed data is not consistent with having come from
  the fitted distribution.

In assessing this, a test statistic, (T\*) is first computed between the
observed data and the best-fit model under the null hypothesis. Next, n
(default n = 10,000) synthetic datasets are generated from the best-fit
model, each being the same size as the empirical dataset. The same
statistics are then computed for each synthetic dataset against the
best-fit model (without refitting) which yields an empirical null
distribution of the statistic under the fitted model. The observed
statistic (T\*) is then compared to this null distribution to obtain a
bootstrap p-value (p\*), defined as the proportion of bootstrap
statistics at least as extreme as the observed one. This p\* is further
mapped to a binary decision variable q\* where q\* = 0 indicates
rejection of the null hypothesis at a significance level α (specified by
‘Alpha’), and q\* = 1 indicates otherwise.

The specific choice of the test statistic is determined by the
user-specified `pvaloption` (where `pvaloption` = `KS`, `AD`, `CvM`
correspond to distance based test statistics such as the
Kolmogorov–Smirnov (KS), Anderson–Darling (AD), and the Cramér–von Mises
(CvM) statistic, respectively).

As an example to illustrate this, when the default `pvaloption` = `KS`
is specified, the Kolmogorov–Smirnov (KS) test statistic is used, which
calculates the maximum distance between empirical data CDF and the
best-fit model CDF. This statistic (T\*) is first computed for the
observed data CDF against the best-fit model CDF. This same statistic is
then computed for each of the n synthetic datasets against best-fit
model CDF, producing an empirical null distribution under the assumption
that the null hypothesis is true. The p-value (p\*) is calculated as the
proportion of bootstrap statistics greater than or equal to the observed
statistic (T\*), i.e., quantifies the probability of observing a
statistic at least as extreme as the empirical one under the null (H0).
At a chosen significance level α, the null is rejected if p\* \< α and
retained otherwise. In implementation, this decision is encoded using a
binary indicator q where q = 0 denotes rejection of the null (model
fails) while q = 1 denotes non-rejection of the null (model passes)

#### Additional note for the Erlang-Exponential Fit

Since there is no closed-form CDF for the Erlang–Exponential
distribution, the cumulative distribution is obtained through numerical
integration of the convolution expression.

This function first approximates the Erlang–Exponential PDF on a uniform
grid of points 0:interval:max(observed data) (with default interval =
0.01). These approximate PDF values are then cumulatively integrated
across the grid to obtain the CDF. Interpolation is used to estimate CDF
values at the observed data points. This approximate CDF is subsequently
used to evaluate the distance-based test statistics, (i.e., maximum
distance between empirical data CDF and the fitted model CDF).

Note that this procedure involves two layers of numerical
integration—first to approximate the PDF and second to obtain the CDF—so
small numerical errors can accumulate, and results should be interpreted
with this in mind. The default interval is set to 0.01, which was chosen
as a compromise between numerical accuracy and computational efficiency:
it is sufficiently fine to ensure smooth numerical integration of the
PDF, while remaining coarse enough to avoid unnecessary computational
overhead.

------------------------------------------------------------------------

## Smallest statistically acceptable K Model (Smallest K)

### Erlang Smallest K

When the `SmallestK` option is passed as an input, the algorithm
performs a search over decreasing integer values of K to identify the
smallest K that still yields an adequate model fit. Starting from the
best fit K\*, the procedure repeatedly lowers K (with λ updated to
preserve the MLE relationship) and re-evaluates the goodness-of-fit.
This continues until the goodness-of-fit no longer passes at the
significance level α (i.e., p\* \< α), at which point the last accepted
K is recorded as the smallest admissible K. The corresponding λ,
p-value, binary decision q, and log-likelihood are then reported under
the results.

### Erlang-Exp Smallest K

This works the same way as above when the `SmallestK` option is passed
as an input but internally, it requires a starting value of K. If an
additional input `SmallestKValue` is provided, this is used directly and
sets the starting value of K for the iterative search. This can help
reduce computational time if users suspect that the smallest acceptable
K is much lower than the best-fitting K. Otherwise the algorithm runs
the Erlang-Fit function to identify the smallest admissible K and uses
this as the starting value for the Erlang–Exponential model.
