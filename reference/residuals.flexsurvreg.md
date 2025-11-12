# Calculate residuals for flexible survival models

Calculates residuals for
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
or
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
model fits.

## Usage

``` r
# S3 method for class 'flexsurvreg'
residuals(object, type = "response", ...)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- type:

  Character string for the type of residual desired. Currently only
  `"response"` and `"coxsnell"` are supported. More residual types may
  become available in future versions.

- ...:

  Not currently used.

## Value

Numeric vector with the same length as `nobs(object)`.

## Details

Residuals of `type = "response"` are calculated as the naive difference
between the observed survival and the covariate-specific predicted mean
survival from
[`predict.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/predict.flexsurvreg.md),
ignoring whether the event time is observed or censored.

`type="coxsnell"` returns the Cox-Snell residual, defined as the
estimated cumulative hazard at each data point. To check the fit of the
A more fully featured utility for this is provided in the function
[`coxsnell_flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/coxsnell_flexsurvreg.md).

## See also

[`predict.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/predict.flexsurvreg.md)

## Examples

``` r
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
residuals(fitg, type="response")
#>  [1]  -187.21560   -88.78121  -255.36822  -873.56963 -1255.83883  -541.96330
#>  [7]  -483.01162  -258.61757   -25.66308  -541.58483  -324.11343  -977.14761
#> [13]    20.81491  -167.52130 -3640.36131 -2316.90107 -3552.79251 -1681.40682
#> [19]  -105.56084 -1679.31116   476.11922    64.41307 -2839.10079  -193.52725
#> [25]  -126.81187  -462.83150


```
