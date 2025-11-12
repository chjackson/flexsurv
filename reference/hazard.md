# Hazard and cumulative hazard functions

Hazard and cumulative hazard functions for distributions which are built
into flexsurv, and whose distribution functions are in base R.

## Usage

``` r
hexp(x, rate = 1, log = FALSE)

Hexp(x, rate = 1, log = FALSE)

hgamma(x, shape, rate = 1, log = FALSE)

Hgamma(x, shape, rate = 1, log = FALSE)

hlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)

Hlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)

hweibull(x, shape, scale = 1, log = FALSE)

Hweibull(x, shape, scale = 1, log = FALSE)
```

## Arguments

- x:

  Vector of quantiles

- rate:

  Rate parameter (exponential and gamma)

- log:

  Compute log hazard or log cumulative hazard

- shape:

  Shape parameter (Weibull and gamma)

- meanlog:

  Mean on the log scale (log normal)

- sdlog:

  Standard deviation on the log scale (log normal)

- scale:

  Scale parameter (Weibull)

## Value

Hazard (functions beginning 'h') or cumulative hazard (functions
beginning 'H').

## Details

For the exponential and the Weibull these are available analytically,
and so are programmed here in numerically stable and efficient forms.

For the gamma and log-normal, these are simply computed as minus the log
of the survivor function (cumulative hazard) or the ratio of the density
and survivor function (hazard), so are not expected to be robust to
extreme values or quick to compute.

## See also

[`dexp`](https://rdrr.io/r/stats/Exponential.html),[`dweibull`](https://rdrr.io/r/stats/Weibull.html),[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),[`dlnorm`](https://ehar.se/r/eha/reference/Lognormal.html),[`dgompertz`](http://chjackson.github.io/flexsurv-dev/reference/Gompertz.md),[`dgengamma`](http://chjackson.github.io/flexsurv-dev/reference/GenGamma.md),[`dgenf`](http://chjackson.github.io/flexsurv-dev/reference/GenF.md)

## Author

Christopher Jackson \<chris.jackson@mrc-bsu.cam.ac.uk\>
