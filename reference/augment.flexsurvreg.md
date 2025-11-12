# Augment data with information from a flexsurv model object

Augment accepts a model object and a dataset and adds information about
each observation in the dataset. Most commonly, this includes predicted
values in the `.fitted` column, residuals in the `.resid` column, and
standard errors for the fitted values in a `.se.fit` column. New columns
always begin with a . prefix to avoid overwriting columns in the
original dataset.

## Usage

``` r
# S3 method for class 'flexsurvreg'
augment(
  x,
  data = NULL,
  newdata = NULL,
  type.predict = "response",
  type.residuals = "response",
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

- data:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  containing the original data that was used to produce the object `x`.

- newdata:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  containing all the original predictors used to create `x`. Defaults to
  `NULL`, indicating that nothing has been passed to `newdata`. If
  `newdata` is specified, the `data` argument will be ignored.

- type.predict:

  Character indicating type of prediction to use. Passed to the `type`
  argument of the [`predict`](https://rdrr.io/r/stats/predict.html)
  generic. Allowed arguments vary with model class, so be sure to read
  the `predict.my_class` documentation.

- type.residuals:

  Character indicating type of residuals to use. Passed to the type
  argument of [`residuals`](https://rdrr.io/r/stats/residuals.html)
  generic. Allowed arguments vary with model class, so be sure to read
  the `residuals.my_class` documentation.

- ...:

  Additional arguments. Not currently used.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing `data` or `newdata` and possible additional columns:

- `.fitted` Fitted values of model

- `.se.fit` Standard errors of fitted values

- `.resid` Residuals (not present if `newdata` specified)

## Details

If neither of `data` or `newdata` are specified, then `model.frame(x)`
will be used. It is worth noting that `model.frame(x)` will include a
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object and not the
original time-to-event variables used when fitting the `flexsurvreg`
object. If the original data is desired, specify `data`.

## Examples

``` r
fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "exp")
augment(fit, data = ovarian)
#> # A tibble: 26 × 9
#>    futime fustat   age resid.ds    rx ecog.ps .pred_time .std_error  .resid
#>     <dbl>  <dbl> <dbl>    <dbl> <dbl>   <dbl>      <dbl>      <dbl>   <dbl>
#>  1     59      1  72.3        2     1       1       213.       113.  -154. 
#>  2    115      1  74.5        2     1       1       165.       103.   -49.8
#>  3    156      1  66.5        2     1       2       427.       149.  -271. 
#>  4    421      0  53.4        2     2       1      2016.       879. -1595. 
#>  5    431      1  50.3        2     1       1      2886.      1575. -2455. 
#>  6    448      0  56.4        1     1       2      1402.       495.  -954. 
#>  7    464      1  56.9        2     2       2      1320.       452.  -856. 
#>  8    475      1  59.9        2     2       2       934.       281.  -459. 
#>  9    477      0  64.2        2     1       1       560.       174.   -82.8
#> 10    563      1  55.2        1     2       2      1626.       623. -1063. 
#> # ℹ 16 more rows
```
