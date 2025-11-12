# Generalized F distribution

Density, distribution function, hazards, quantile function and random
generation for the generalized F distribution, using the
reparameterisation by Prentice (1975).

## Usage

``` r
dgenf(x, mu = 0, sigma = 1, Q, P, log = FALSE)

pgenf(q, mu = 0, sigma = 1, Q, P, lower.tail = TRUE, log.p = FALSE)

Hgenf(x, mu = 0, sigma = 1, Q, P)

hgenf(x, mu = 0, sigma = 1, Q, P)

qgenf(p, mu = 0, sigma = 1, Q, P, lower.tail = TRUE, log.p = FALSE)

rgenf(n, mu = 0, sigma = 1, Q, P)
```

## Arguments

- x, q:

  Vector of quantiles.

- mu:

  Vector of location parameters.

- sigma:

  Vector of scale parameters.

- Q:

  Vector of first shape parameters.

- P:

  Vector of second shape parameters.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X \> x)\\.

- p:

  Vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dgenf` gives the density, `pgenf` gives the distribution function,
`qgenf` gives the quantile function, `rgenf` generates random deviates,
`Hgenf` retuns the cumulative hazard and `hgenf` the hazard.

## Details

If \\y \sim F(2s_1, 2s_2)\\, and \\w = \\\\ \log(y)\\ then \\x =
\exp(w\sigma + \mu)\\ has the original generalized F distribution with
location parameter \\\mu\\, scale parameter \\\sigma\>0\\ and shape
parameters \\s_1,s_2\\.

In this more stable version described by Prentice (1975), \\s_1,s_2\\
are replaced by shape parameters \\Q,P\\, with \\P\>0\\, and

\$\$s_1 = 2(Q^2 + 2P + Q\delta)^{-1}, \quad s_2 = 2(Q^2 + 2P -
Q\delta)^{-1}\$\$ equivalently \$\$Q = \left(\frac{1}{s_1} -
\frac{1}{s_2}\right)\left(\frac{1}{s_1} + \frac{1}{s_2}\right)^{-1/2},
\quad P = \frac{2}{s_1 + s_2} \$\$

Define \\\delta = (Q^2 + 2P)^{1/2}\\, and \\w = (\log(x) - \mu)\delta
/\sigma\\, then the probability density function of \\x\\ is \$\$ f(x) =
\frac{\delta (s_1/s_2)^{s_1} e^{s_1 w}}{\sigma x (1 + s_1 e^w/s_2) ^
{(s_1 + s_2)} B(s_1, s_2)} \$\$ The original parameterisation is
available in this package as
[`dgenf.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenF.orig.md),
for the sake of completion / compatibility. With the above definitions,

`dgenf(x, mu=mu, sigma=sigma, Q=Q, P=P) = dgenf.orig(x, mu=mu, sigma=sigma/delta, s1=s1, s2=s2)`

The generalized F distribution with `P=0` is equivalent to the
generalized gamma distribution
[`dgengamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),
so that `dgenf(x, mu, sigma, Q, P=0)` equals
`dgengamma(x, mu, sigma, Q)`. The generalized gamma reduces further to
several common distributions, as described in the
[`GenGamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md)
help page.

The generalized F distribution includes the log-logistic distribution
(see [`eha::dllogis`](https://ehar.se/r/eha/reference/Loglogistic.html))
as a further special case:

`dgenf(x, mu=mu, sigma=sigma, Q=0, P=1) = `[`eha::dllogis`](https://ehar.se/r/eha/reference/Loglogistic.html)`(x, shape=sqrt(2)/sigma, scale=exp(mu))`

The range of hazard trajectories available under this distribution are
discussed in detail by Cox (2008). Jackson et al. (2010) give an
application to modelling oral cancer survival for use in a health
economic evaluation of screening.

## Note

The parameters `Q` and `P` are usually called \\q\\ and \\p\\ in the
literature - they were made upper-case in these R functions to avoid
clashing with the conventional arguments `q` in the probability function
and `p` in the quantile function.

## References

R. L. Prentice (1975). Discrimination among some parametric models.
Biometrika 62(3):607-614.

Cox, C. (2008). The generalized \\F\\ distribution: An umbrella for
parametric survival analysis. Statistics in Medicine 27:4301-4312.

Jackson, C. H. and Sharples, L. D. and Thompson, S. G. (2010). Survival
models in health economic evaluations: balancing fit and parsimony to
improve prediction. International Journal of Biostatistics 6(1):Article
34.

## See also

[`GenF.orig`](http://chjackson.github.io/flexsurv-dev/reference/GenF.orig.md),
[`GenGamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
