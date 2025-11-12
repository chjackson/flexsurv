# Tidy a flexsurv model object

Tidy summarizes information about the components of the model into a
tidy data frame.

## Usage

``` r
# S3 method for class 'flexsurvreg'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  pars = "all",
  transform = "none",
  ...
)
```

## Arguments

- x:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- conf.int:

  Logical. Should confidence intervals be returned? Default is `FALSE`.

- conf.level:

  The confidence level to use for the confidence interval if
  `conf.int = TRUE`. Default is `0.95`.

- pars:

  One of `"all"`, `"baseline"`, or `"coefs"` for all parameters,
  baseline distribution parameters, or covariate effects (i.e.
  regression betas), respectively. Default is `"all"`.

- transform:

  Character vector of transformations to apply to requested `pars`.
  Default is `"none"`, which returns `pars` as-is.

  Users can specify one or both types of transformations:

  - `"baseline.real"` which transforms the baseline distribution
    parameters to the real number line used for estimation.

  - `"coefs.exp"` which exponentiates the covariate effects.

  See `Details` for a more complete explanation.

- ...:

  Not currently used.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing the columns: `term`, `estimate`, `std.error`, `statistic`,
`p.value`, `conf.low`, and `conf.high`, by default.

`statistic` and `p.value` are only provided for covariate effects (`NA`
for baseline distribution parameters). These are computed as Wald-type
test statistics with p-values from a standard normal distribution.

## Details

`flexsurvreg` models estimate two types of coefficients, baseline
distribution parameters, and covariate effects which act on the baseline
distribution. By design, `flexsurvreg` returns distribution parameters
on the same scale as is found in the relevant `d/p/q/r` functions.
Covariate effects are returned on the log-scale, which represents either
log-time ratios (accelerated failure time models) or log-hazard ratios
for proportional hazard models. By default,
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) will return
baseline distribution parameters on their natural scale and covariate
effects on the log-scale.

To transform the baseline distribution parameters to the real-value
number line (the scale used for estimation), pass the character argument
`"baseline.real"` to `transform`. To get time ratios or hazard ratios,
pass `"coefs.exp"` to `transform`. These transformations may be done
together by submitting both arguments as a character vector.

## Examples

``` r
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
tidy(fitg)
#> # A tibble: 4 × 5
#>   term  estimate std.error statistic   p.value
#>   <chr>    <dbl>     <dbl>     <dbl>     <dbl>
#> 1 mu     11.7       1.66       NA    NA       
#> 2 sigma   0.751     0.244      NA    NA       
#> 3 Q       0.295     0.912      NA    NA       
#> 4 age    -0.0875    0.0250     -3.50  0.000467
tidy(fitg, pars = "coefs", transform = "coefs.exp")
#> # A tibble: 1 × 5
#>   term  estimate std.error statistic  p.value
#>   <chr>    <dbl>     <dbl>     <dbl>    <dbl>
#> 1 age      0.916      1.03     -3.50 0.000467
```
