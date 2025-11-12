# Mean and restricted mean survival functions

Mean and restricted mean survival time functions for distributions which
are built into flexsurv.

## Usage

``` r
mean_exp(rate = 1)

rmst_exp(t, rate = 1, start = 0)

mean_gamma(shape, rate = 1)

rmst_gamma(t, shape, rate = 1, start = 0)

rmst_genf(t, mu, sigma, Q, P, start = 0)

mean_genf(mu, sigma, Q, P)

rmst_genf.orig(t, mu, sigma, s1, s2, start = 0)

mean_genf.orig(mu, sigma, s1, s2)

rmst_gengamma(t, mu = 0, sigma = 1, Q, start = 0)

mean_gengamma(mu = 0, sigma = 1, Q)

rmst_gengamma.orig(t, shape, scale = 1, k, start = 0)

mean_gengamma.orig(shape, scale = 1, k)

rmst_gompertz(t, shape, rate = 1, start = 0)

mean_gompertz(shape, rate = 1)

mean_llogis(shape = 1, scale = 1)

rmst_llogis(t, shape = 1, scale = 1, start = 0)

mean_lnorm(meanlog = 0, sdlog = 1)

rmst_lnorm(t, meanlog = 0, sdlog = 1, start = 0)

mean_weibull(shape, scale = 1)

rmst_weibull(t, shape, scale = 1, start = 0)

rmst_weibullPH(t, shape, scale = 1, start = 0)

mean_weibullPH(shape, scale = 1)
```

## Arguments

- rate:

  Rate parameter (exponential and gamma)

- t:

  Vector of times to which restricted mean survival time is evaluated

- start:

  Optional left-truncation time or times. The returned restricted mean
  survival will be conditioned on survival up to this time.

- shape:

  Shape parameter (Weibull, gamma, log-logistic, generalized gamma
  \[orig\], generalized F \[orig\])

- mu:

  Mean on the log scale (generalized gamma, generalized F)

- sigma:

  Standard deviation on the log scale (generalized gamma, generalized F)

- Q:

  Vector of first shape parameters (generalized gamma, generalized F)

- P:

  Vector of second shape parameters (generalized F)

- s1:

  Vector of first F shape parameters (generalized F \[orig\])

- s2:

  vector of second F shape parameters (generalized F \[orig\])

- scale:

  Scale parameter (Weibull, log-logistic, generalized gamma \[orig\],
  generalized F \[orig\])

- k:

  vector of shape parameters (generalized gamma \[orig\]).

- meanlog:

  Mean on the log scale (log normal)

- sdlog:

  Standard deviation on the log scale (log normal)

## Value

mean survival (functions beginning 'mean') or restricted mean survival
(functions beginning 'rmst\_').

## Details

For the exponential, Weibull, log-logistic, lognormal, and gamma, mean
survival is provided analytically. Restricted mean survival for the
exponential distribution is also provided analytically. Mean and
restricted means for other distributions are calculated via numeric
integration.

## See also

[`dexp`](https://rdrr.io/r/stats/Exponential.html),[`dweibull`](https://rdrr.io/r/stats/Weibull.html),[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),[`dlnorm`](https://ehar.se/r/eha/reference/Lognormal.html),[`dgompertz`](http://chjackson.github.io/flexsurv-dev/reference/Gompertz.md),[`dgengamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),[`dgenf`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
