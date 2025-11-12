# Generalized gamma distribution

Density, distribution function, hazards, quantile function and random
generation for the generalized gamma distribution, using the
parameterisation originating from Prentice (1974). Also known as the
(generalized) log-gamma distribution.

## Usage

``` r
dgengamma(x, mu = 0, sigma = 1, Q, log = FALSE)

pgengamma(q, mu = 0, sigma = 1, Q, lower.tail = TRUE, log.p = FALSE)

Hgengamma(x, mu = 0, sigma = 1, Q)

hgengamma(x, mu = 0, sigma = 1, Q)

qgengamma(p, mu = 0, sigma = 1, Q, lower.tail = TRUE, log.p = FALSE)

rgengamma(n, mu = 0, sigma = 1, Q)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  Vector of “location” parameters.

- sigma:

  Vector of “scale” parameters. Note the inconsistent meanings of the
  term “scale” - this parameter is analogous to the (log-scale) standard
  deviation of the log-normal distribution, “sdlog” in
  [`dlnorm`](https://rdrr.io/r/stats/Lognormal.html), rather than the
  “scale” parameter of the gamma distribution
  [`dgamma`](https://rdrr.io/r/stats/GammaDist.html). Constrained to be
  positive.

- Q:

  Vector of shape parameters.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X \> x)\\.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dgengamma` gives the density, `pgengamma` gives the distribution
function, `qgengamma` gives the quantile function, `rgengamma` generates
random deviates, `Hgengamma` retuns the cumulative hazard and
`hgengamma` the hazard.

## Details

If \\\gamma \sim Gamma(Q^{-2}, 1)\\ , and \\w = log(Q^2 \gamma) / Q\\,
then \\x = \exp(\mu + \sigma w)\\ follows the generalized gamma
distribution with probability density function

\$\$f(x \| \mu, \sigma, Q) = \frac{\|Q\|(Q^{-2})^{Q^{-2}}}{\sigma x
\Gamma(Q^{-2})} \exp(Q^{-2}(Qw - \exp(Qw)))\$\$

This parameterisation is preferred to the original parameterisation of
the generalized gamma by Stacy (1962) since it is more numerically
stable near to \\Q=0\\ (the log-normal distribution), and allows
\\Q\<=0\\. The original is available in this package as
[`dgengamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md),
for the sake of completion and compatibility with other software - this
is implicitly restricted to `Q`\>0 (or `k`\>0 in the original notation).
The parameters of `dgengamma` and
[`dgengamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md)
are related as follows.

`dgengamma.orig(x, shape=shape, scale=scale, k=k) = `

`dgengamma(x, mu=log(scale) + log(k)/shape, sigma=1/(shape*sqrt(k)), Q=1/sqrt(k))`

The generalized gamma distribution simplifies to the gamma, log-normal
and Weibull distributions with the following parameterisations:

|                                    |     |                                                       |
|------------------------------------|-----|-------------------------------------------------------|
| `dgengamma(x, mu, sigma, Q=0)`     | `=` | `dlnorm(x, mu, sigma)`                                |
| `dgengamma(x, mu, sigma, Q=1)`     | `=` | `dweibull(x, shape=1/sigma, scale=exp(mu))`           |
| `dgengamma(x, mu, sigma, Q=sigma)` | `=` | `dgamma(x, shape=1/sigma^2, rate=exp(-mu) / sigma^2)` |

The properties of the generalized gamma and its applications to survival
analysis are discussed in detail by Cox (2007).

The generalized F distribution
[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md)
extends the generalized gamma to four parameters.

## References

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

Farewell, V. T. and Prentice, R. L. (1977). A study of distributional
shape in life testing. Technometrics 19(1):69-75.

Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
distributions. Technometrics 22(3):409-419.

Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007). Parametric
survival analysis and taxonomy of hazard functions for the generalized
gamma distribution. Statistics in Medicine 26:4252-4374

Stacy, E. W. (1962). A generalization of the gamma distribution. Annals
of Mathematical Statistics 33:1187-92

## See also

[`GenGamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md),
[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md),
[`Lognormal`](https://rdrr.io/r/stats/Lognormal.html),
[`GammaDist`](https://rdrr.io/r/stats/GammaDist.html),
[`Weibull`](https://rdrr.io/r/stats/Weibull.html).

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
