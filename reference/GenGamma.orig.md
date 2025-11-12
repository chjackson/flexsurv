# Generalized gamma distribution (original parameterisation)

Density, distribution function, hazards, quantile function and random
generation for the generalized gamma distribution, using the original
parameterisation from Stacy (1962).

## Usage

``` r
dgengamma.orig(x, shape, scale = 1, k, log = FALSE)

pgengamma.orig(q, shape, scale = 1, k, lower.tail = TRUE, log.p = FALSE)

Hgengamma.orig(x, shape, scale = 1, k)

hgengamma.orig(x, shape, scale = 1, k)

qgengamma.orig(p, shape, scale = 1, k, lower.tail = TRUE, log.p = FALSE)

rgengamma.orig(n, shape, scale = 1, k)
```

## Arguments

- x, q:

  vector of quantiles.

- shape:

  vector of “Weibull” shape parameters.

- scale:

  vector of scale parameters.

- k:

  vector of “Gamma” shape parameters.

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

`dgengamma.orig` gives the density, `pgengamma.orig` gives the
distribution function, `qgengamma.orig` gives the quantile function,
`rgengamma.orig` generates random deviates, `Hgengamma.orig` retuns the
cumulative hazard and `hgengamma.orig` the hazard.

## Details

If \\w \sim Gamma(k,1)\\, then \\x = \exp(w/shape + \log(scale))\\
follows the original generalised gamma distribution with the
parameterisation given here (Stacy 1962). Defining `shape`\\=b\>0\\,
`scale`\\=a\>0\\, \\x\\ has probability density

\$\$f(x \| a, b, k) = \frac{b}{\Gamma(k)} \frac{x^{bk - 1}}{a^{bk}}
\$\$\$\$ \exp(-(x/a)^b)\$\$

The original generalized gamma distribution simplifies to the gamma,
exponential and Weibull distributions with the following
parameterisations:

|                                          |     |                                                                         |
|------------------------------------------|-----|-------------------------------------------------------------------------|
| `dgengamma.orig(x, shape, scale, k=1)`   | `=` | [`dweibull`](https://rdrr.io/r/stats/Weibull.html)`(x, shape, scale)`   |
| `dgengamma.orig(x, shape=1, scale, k)`   | `=` | [`dgamma`](https://rdrr.io/r/stats/GammaDist.html)`(x, shape=k, scale)` |
| `dgengamma.orig(x, shape=1, scale, k=1)` | `=` | [`dexp`](https://rdrr.io/r/stats/Exponential.html)`(x, rate=1/scale)`   |

Also as k tends to infinity, it tends to the log normal (as in
[`dlnorm`](https://rdrr.io/r/stats/Lognormal.html)) with the following
parameters (Lawless, 1980):

`dlnorm(x, meanlog=log(scale) + log(k)/shape, sdlog=1/(shape*sqrt(k)))`

For more stable behaviour as the distribution tends to the log-normal,
an alternative parameterisation was developed by Prentice (1974). This
is given in
[`dgengamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),
and is now preferred for statistical modelling. It is also more
flexible, including a further new class of distributions with negative
shape `k`.

The generalized F distribution
[`GenF.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenF.orig.md),
and its similar alternative parameterisation
[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md),
extend the generalized gamma to four parameters.

## References

Stacy, E. W. (1962). A generalization of the gamma distribution. Annals
of Mathematical Statistics 33:1187-92.

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
distributions. Technometrics 22(3):409-419.

## See also

[`GenGamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),
[`GenF.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenF.orig.md),
[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md),
[`Lognormal`](https://rdrr.io/r/stats/Lognormal.html),
[`GammaDist`](https://rdrr.io/r/stats/GammaDist.html),
[`Weibull`](https://rdrr.io/r/stats/Weibull.html).

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
