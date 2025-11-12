# Log likelihood from a flexsurvreg model

Log likelihood from a flexsurvreg model

## Usage

``` r
# S3 method for class 'flexsurvreg'
logLik(object, ...)
```

## Arguments

- object:

  A fitted model object of class
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
  e.g. as returned by `flexsurvreg` or `flexsurvspline`.

- ...:

  Other arguments (currently unused).

## Value

Log-likelihood (numeric) with additional attributes `df` (degrees of
freedom, or number of parameters that were estimated), and number of
observations `nobs` (including observed events and censored
observations).
