# The Gompertz distribution

Density, distribution function, hazards, quantile function and random
generation for the Gompertz distribution with unrestricted shape.

## Usage

``` r
dgompertz(x, shape, rate = 1, log = FALSE)

pgompertz(q, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)

qgompertz(p, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)

rgompertz(n, shape = 1, rate = 1)

hgompertz(x, shape, rate = 1, log = FALSE)

Hgompertz(x, shape, rate = 1, log = FALSE)
```

## Arguments

- x, q:

  vector of quantiles.

- shape, rate:

  vector of shape and rate parameters.

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

`dgompertz` gives the density, `pgompertz` gives the distribution
function, `qgompertz` gives the quantile function, `hgompertz` gives the
hazard function, `Hgompertz` gives the cumulative hazard function, and
`rgompertz` generates random deviates.

## Details

The Gompertz distribution with `shape` parameter \\a\\ and `rate`
parameter \\b\\ has probability density function

\$\$f(x \| a, b) = be^{ax}\exp(-b/a (e^{ax} - 1))\$\$

and hazard

\$\$h(x \| a, b) = b e^{ax}\$\$

The hazard is increasing for shape \\a\>0\\ and decreasing for \\a\<0\\.
For \\a=0\\ the Gompertz is equivalent to the exponential distribution
with constant hazard and rate \\b\\.

The probability distribution function is \$\$F(x \| a, b) = 1 -
\exp(-b/a (e^{ax} - 1))\$\$

Thus if \\a\\ is negative, letting \\x\\ tend to infinity shows that
there is a non-zero probability \\\exp(b/a)\\ of living forever. On
these occasions `qgompertz` and `rgompertz` will return `Inf`.

## Note

Some implementations of the Gompertz restrict \\a\\ to be strictly
positive, which ensures that the probability of survival decreases to
zero as \\x\\ increases to infinity. The more flexible implementation
given here is consistent with `streg` in Stata.

The functions
[`eha::dgompertz`](https://rdrr.io/pkg/eha/man/Gompertz.html) and
similar available in the package eha label the parameters the other way
round, so that what is called the `shape` there is called the `rate`
here, and what is called `1 / scale` there is called the `shape` here.
The terminology here is consistent with the exponential
[`dexp`](https://rdrr.io/r/stats/Exponential.html) and Weibull
[`dweibull`](https://rdrr.io/r/stats/Weibull.html) distributions in R.

## References

Stata Press (2007) Stata release 10 manual: Survival analysis and
epidemiological tables.

## See also

[`dexp`](https://rdrr.io/r/stats/Exponential.html)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
