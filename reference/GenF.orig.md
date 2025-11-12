# Generalized F distribution (original parameterisation)

Density, distribution function, quantile function and random generation
for the generalized F distribution, using the less flexible original
parameterisation described by Prentice (1975).

## Usage

``` r
dgenf.orig(x, mu = 0, sigma = 1, s1, s2, log = FALSE)

pgenf.orig(q, mu = 0, sigma = 1, s1, s2, lower.tail = TRUE, log.p = FALSE)

Hgenf.orig(x, mu = 0, sigma = 1, s1, s2)

hgenf.orig(x, mu = 0, sigma = 1, s1, s2)

qgenf.orig(p, mu = 0, sigma = 1, s1, s2, lower.tail = TRUE, log.p = FALSE)

rgenf.orig(n, mu = 0, sigma = 1, s1, s2)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  Vector of location parameters.

- sigma:

  Vector of scale parameters.

- s1:

  Vector of first F shape parameters.

- s2:

  vector of second F shape parameters.

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

`dgenf.orig` gives the density, `pgenf.orig` gives the distribution
function, `qgenf.orig` gives the quantile function, `rgenf.orig`
generates random deviates, `Hgenf.orig` retuns the cumulative hazard and
`hgenf.orig` the hazard.

## Details

If \\y \sim F(2s_1, 2s_2)\\, and \\w = \log(y)\\ then \\x =
\exp(w\sigma + \mu)\\ has the original generalized F distribution with
location parameter \\\mu\\, scale parameter \\\sigma\>0\\ and shape
parameters \\s_1\>0,s_2\>0\\. The probability density function of \\x\\
is

\$\$f(x \| \mu, \sigma, s_1, s_2) = \frac{(s_1/s_2)^{s_1} e^{s_1
w}}{\sigma x (1 + s_1 e^w/s_2) ^ {(s_1 + s_2)} B(s_1, s_2)}\$\$

where \\w = (\log(x) - \mu)/\sigma\\, and \\B(s_1,s_2) =
\Gamma(s_1)\Gamma(s_2)/\Gamma(s_1+s_2)\\ is the beta function.

As \\s_2 \rightarrow \infty\\, the distribution of \\x\\ tends towards
an original generalized gamma distribution with the following
parameters:

[`dgengamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md)`(x, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1)`

See
[`GenGamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md)
for how this includes several other common distributions as special
cases.

The alternative parameterisation of the generalized F distribution,
originating from Prentice (1975) and given in this package as
[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md), is
preferred for statistical modelling, since it is more stable as \\s_1\\
tends to infinity, and includes a further new class of distributions
with negative first shape parameter. The original is provided here for
the sake of completion and compatibility.

## References

R. L. Prentice (1975). Discrimination among some parametric models.
Biometrika 62(3):607-614.

## See also

[`GenF`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md),
[`GenGamma.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.orig.md),
[`GenGamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
