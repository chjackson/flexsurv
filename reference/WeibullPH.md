# Weibull distribution in proportional hazards parameterisation

Density, distribution function, hazards, quantile function and random
generation for the Weibull distribution in its proportional hazards
parameterisation.

## Usage

``` r
dweibullPH(x, shape, scale = 1, log = FALSE)

pweibullPH(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)

qweibullPH(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)

hweibullPH(x, shape, scale = 1, log = FALSE)

HweibullPH(x, shape, scale = 1, log = FALSE)

rweibullPH(n, shape, scale = 1)
```

## Arguments

- x, q:

  Vector of quantiles.

- shape:

  Vector of shape parameters.

- scale:

  Vector of scale parameters.

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

`dweibullPH` gives the density, `pweibullPH` gives the distribution
function, `qweibullPH` gives the quantile function, `rweibullPH`
generates random deviates, `HweibullPH` retuns the cumulative hazard and
`hweibullPH` the hazard.

## Details

The Weibull distribution in proportional hazards parameterisation with
\`shape' parameter a and \`scale' parameter m has density given by

\$\$f(x) = a m x^{a-1} exp(- m x^a) \$\$

cumulative distribution function \\F(x) = 1 - exp( -m x^a )\\, survivor
function \\S(x) = exp( -m x^a )\\, cumulative hazard \\m x^a\\ and
hazard \\a m x^{a-1}\\.

[`dweibull`](https://rdrr.io/r/stats/Weibull.html) in base R has the
alternative 'accelerated failure time' (AFT) parameterisation with shape
a and scale b. The shape parameter \\a\\ is the same in both versions.
The scale parameters are related as \\b = m^{-1/a}\\, equivalently m =
b^-a.

In survival modelling, covariates are typically included through a
linear model on the log scale parameter. Thus, in the proportional
hazards model, the coefficients in such a model on \\m\\ are interpreted
as log hazard ratios.

In the AFT model, covariates on \\b\\ are interpreted as time
acceleration factors. For example, doubling the value of a covariate
with coefficient \\beta=log(2)\\ would give half the expected survival
time. These coefficients are related to the log hazard ratios \\\gamma\\
as \\\beta = -\gamma / a\\.

## See also

[`dweibull`](https://rdrr.io/r/stats/Weibull.html)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
