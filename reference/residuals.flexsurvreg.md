# Calculate residuals for flexible survival models

Calculates residuals for
[`flexsurvreg`](https://chjackson.github.io/flexsurv/reference/flexsurvreg.md)
or
[`flexsurvspline`](https://chjackson.github.io/flexsurv/reference/flexsurvspline.md)
model fits.

## Usage

``` r
# S3 method for class 'flexsurvreg'
residuals(object, type = "response", ...)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](https://chjackson.github.io/flexsurv/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](https://chjackson.github.io/flexsurv/reference/flexsurvspline.md),
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
[`predict.flexsurvreg`](https://chjackson.github.io/flexsurv/reference/predict.flexsurvreg.md),
ignoring whether the event time is observed or censored.

`type="coxsnell"` returns the Cox-Snell residual, defined as the
estimated cumulative hazard at each data point. To check the fit of the
A more fully featured utility for this is provided in the function
[`coxsnell_flexsurvreg`](https://chjackson.github.io/flexsurv/reference/coxsnell_flexsurvreg.md).

## See also

[`predict.flexsurvreg`](https://chjackson.github.io/flexsurv/reference/predict.flexsurvreg.md)

## Examples

``` r

fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
residuals(fitg, type="response")
#>  [1]  -187.21563   -88.78124  -255.36823  -873.56932 -1255.83832  -541.96312
#>  [7]  -483.01145  -258.61749   -25.66306  -541.58460  -324.11326  -977.14709
#> [13]    20.81500  -167.52113 -3640.35905 -2316.89969 -3552.79014 -1681.40569
#> [19]  -105.56055 -1679.30996   476.11931    64.41304 -2839.09941  -193.52722
#> [25]  -126.81185  -462.83138


```
