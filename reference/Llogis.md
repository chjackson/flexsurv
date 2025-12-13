# The log-logistic distribution

Density, distribution function, hazards, quantile function and random
generation for the log-logistic distribution.

## Usage

``` r
dllogis(x, shape = 1, scale = 1, log = FALSE)

pllogis(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)

qllogis(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)

rllogis(n, shape = 1, scale = 1)

hllogis(x, shape = 1, scale = 1, log = FALSE)

Hllogis(x, shape = 1, scale = 1, log = FALSE)
```

## Arguments

- x, q:

  vector of quantiles.

- shape, scale:

  vector of shape and scale parameters.

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

`dllogis` gives the density, `pllogis` gives the distribution function,
`qllogis` gives the quantile function, `hllogis` gives the hazard
function, `Hllogis` gives the cumulative hazard function, and `rllogis`
generates random deviates.

## Details

The log-logistic distribution with `shape` parameter \\a\>0\\ and
`scale` parameter \\b\>0\\ has probability density function

\$\$f(x \| a, b) = (a/b) (x/b)^{a-1} / (1 + (x/b)^a)^2\$\$

and hazard

\$\$h(x \| a, b) = (a/b) (x/b)^{a-1} / (1 + (x/b)^a)\$\$

for \\x\>0\\. The hazard is decreasing for shape \\a\leq 1\\, and
unimodal for \\a \> 1\\.

The probability distribution function is \$\$F(x \| a, b) = 1 - 1 / (1 +
(x/b)^a)\$\$

If \\a \> 1\\, the mean is \\b c / sin(c)\\, and if \\a \> 2\\ the
variance is \\b^2 \* (2\*c/sin(2\*c) - c^2/sin(c)^2)\\, where \\c =
\pi/a\\, otherwise these are undefined.

## Note

Various different parameterisations of this distribution are used. In
the one used here, the interpretation of the parameters is the same as
in the standard Weibull distribution
([`dweibull`](https://rdrr.io/r/stats/Weibull.html)). Like the Weibull,
the survivor function is a transformation of \\(x/b)^a\\ from the
non-negative real line to \[0,1\], but with a different link function.
Covariates on \\b\\ represent time acceleration factors, or ratios of
expected survival.

The same parameterisation is also used in
[`eha::dllogis`](https://rdrr.io/pkg/eha/man/Loglogistic.html) in the
eha package.

## References

Stata Press (2007) Stata release 10 manual: Survival analysis and
epidemiological tables.

## See also

[`dweibull`](https://rdrr.io/r/stats/Weibull.html)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
