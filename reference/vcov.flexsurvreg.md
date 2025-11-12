# Variance-covariance matrix from a flexsurvreg model

Variance-covariance matrix from a flexsurvreg model

## Usage

``` r
# S3 method for class 'flexsurvreg'
vcov(object, ...)
```

## Arguments

- object:

  A fitted model object of class
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
  e.g. as returned by `flexsurvreg` or `flexsurvspline`.

- ...:

  Other arguments (currently unused).

## Value

Variance-covariance matrix of the estimated parameters, on the scale
that they were estimated on (for positive parameters this is the log
scale).
