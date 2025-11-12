# Glance at a flexsurv model object

Glance accepts a model object and returns a tibble with exactly one row
of model summaries.

## Usage

``` r
# S3 method for class 'flexsurvreg'
glance(x, ...)
```

## Arguments

- x:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- ...:

  Not currently used.

## Value

A one-row [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing columns:

- `N` Number of observations used in fitting

- `events` Number of events

- `censored` Number of censored events

- `trisk` Total length of time-at-risk (i.e. follow-up)

- `df` Degrees of freedom (i.e. number of estimated parameters)

- `logLik` Log-likelihood

- `AIC` Akaike's "An Information Criteria"

- `BIC` Bayesian Information Criteria

## Examples

``` r
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
glance(fitg)
#> # A tibble: 1 × 8
#>       N events censored trisk    df logLik   AIC   BIC
#>   <int>  <int>    <int> <dbl> <int>  <dbl> <dbl> <dbl>
#> 1    26     12       14 15588     4  -89.7  187.  192.
```
