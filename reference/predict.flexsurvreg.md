# Predictions from flexible survival models

Predict outcomes from flexible survival models at the covariate values
specified in `newdata`.

## Usage

``` r
# S3 method for class 'flexsurvreg'
predict(
  object,
  newdata,
  type = "response",
  times,
  start = 0,
  conf.int = FALSE,
  conf.level = 0.95,
  se.fit = FALSE,
  p = c(0.1, 0.9),
  B = 1000,
  ...
)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- newdata:

  Data frame containing covariate values at which to produce fitted
  values. There must be a column for every covariate in the model
  formula used to fit `object`, and one row for every combination of
  covariate values at which to obtain the fitted predictions.

  If `newdata` is omitted, then the original data used to fit the model
  are used, as extracted by `model.frame(object)`. However this will
  currently not work if the model formula contains functions, e.g.
  `~ factor(x)`. The names of the model frame must correspond to
  variables in the original data.

- type:

  Character vector for the type of predictions desired.

  - `"response"` for mean survival time (the default). `"mean"` is an
    acceptable synonym

  - `"quantile"` for quantiles of the survival distribution as specified
    by `p`

  - `"rmst"` for restricted mean survival time

  - `"survival"` for survival probabilities

  - `"cumhaz"` for cumulative hazards

  - `"hazard"` for hazards

  - `"link"` for fitted values of the location parameter, analogous to
    the linear predictor in generalized linear models (`type = "lp"` and
    `type = "linear"` are acceptable synonyms). This is on the natural
    scale of the parameter, not the log scale.

- times:

  Vector of time horizons at which to compute fitted values. Only
  applies when `type` is `"survival"`, `"cumhaz"`, `"hazard"`, or
  `"rmst"`. Will be silently ignored for all other types.

  If not specified, predictions for `"survival"`, `"cumhaz"`, and
  `"hazard"` will be made at each observed event time in
  `model.frame(object)`.

  For `"rmst"`, when `times` is not specified predictions will be made
  at the maximum observed event time from the data used to fit `object`.
  Specifying `times = Inf` is valid, and will return mean survival
  (equal to `type = "response"`).

- start:

  Optional left-truncation time or times. The returned survival, hazard,
  or cumulative hazard will be conditioned on survival up to this time.
  `start` must be length 1 or the same length as `times`. Predicted
  times returned with `type` `"rmst"` or `"quantile"` will be times
  since time zero, not times since the `start` time.

- conf.int:

  Logical. Should confidence intervals be returned? Default is `FALSE`.

- conf.level:

  Width of symmetric confidence intervals, relative to 1.

- se.fit:

  Logical. Should standard errors of fitted values be returned? Default
  is `FALSE`.

- p:

  Vector of quantiles at which to return fitted values when
  `type = "quantile"`. Default is `c(0.1, 0.9)`.

- B:

  Number of simulations from the normal asymptotic distribution of the
  estimates used to calculate confidence intervals or standard errors.
  Decrease for greater speed at the expense of accuracy, or set `B=0` to
  turn off calculation of CIs and SEs.

- ...:

  Not currently used.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
same number of rows as `newdata` and in the same order. If multiple
predictions are requested, a
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing a single list-column of data frames.

For the list-column of data frames - the dimensions of each data frame
will be identical. Rows are added for each value of `times` or `p`
requested.

This function is a wrapper around
[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md),
designed to help flexsurv to integrate with the "tidymodels" ecosystem,
in particular the censored package.
[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md)
returns the same results but in a more conventional format.

## See also

[`summary.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/summary.flexsurvreg.md),
[`residuals.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/residuals.flexsurvreg.md)

## Author

Matthew T. Warkentin (<https://github.com/mattwarkentin>)

## Examples

``` r
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")

## Simplest prediction: mean or median, for covariates defined by original dataset
predict(fitg)
#> # A tibble: 26 × 1
#>    .pred_time
#>         <dbl>
#>  1       246.
#>  2       204.
#>  3       411.
#>  4      1295.
#>  5      1687.
#>  6       990.
#>  7       947.
#>  8       734.
#>  9       503.
#> 10      1105.
#> # ℹ 16 more rows
predict(fitg, type = "quantile", p = 0.5)
#> # A tibble: 26 × 2
#>    .quantile .pred_quantile
#>        <dbl>          <dbl>
#>  1       0.5           194.
#>  2       0.5           161.
#>  3       0.5           325.
#>  4       0.5          1022.
#>  5       0.5          1331.
#>  6       0.5           781.
#>  7       0.5           747.
#>  8       0.5           579.
#>  9       0.5           397.
#> 10       0.5           872.
#> # ℹ 16 more rows

## Simple prediction for user-defined covariate values
predict(fitg, newdata = data.frame(age = c(40, 50, 60)))
#> # A tibble: 3 × 1
#>   .pred_time
#>        <dbl>
#> 1      4169.
#> 2      1738.
#> 3       724.
predict(fitg, type = "quantile", p = 0.5, newdata = data.frame(age = c(40,50,60)))
#> # A tibble: 3 × 2
#>   .quantile .pred_quantile
#>       <dbl>          <dbl>
#> 1       0.5          3291.
#> 2       0.5          1372.
#> 3       0.5           572.

## Predict multiple quantiles and unnest
require(tidyr)
#> Loading required package: tidyr
pr <- predict(fitg, type = "survival", times = c(600, 800))
tidyr::unnest(pr, .pred)
#> # A tibble: 52 × 2
#>    .eval_time .pred_survival
#>         <dbl>          <dbl>
#>  1        600        0.0548 
#>  2        800        0.0202 
#>  3        600        0.0292 
#>  4        800        0.00926
#>  5        600        0.200  
#>  6        800        0.104  
#>  7        600        0.751  
#>  8        800        0.624  
#>  9        600        0.841  
#> 10        800        0.742  
#> # ℹ 42 more rows
```
