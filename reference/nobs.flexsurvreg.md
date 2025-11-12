# Number of observations contributing to a fitted flexible survival model

Number of observations contributing to a fitted flexible survival model

## Usage

``` r
# S3 method for class 'flexsurvreg'
nobs(object, cens = TRUE, ...)
```

## Arguments

- object:

  Output from
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md),
  representing a fitted survival model object.

- cens:

  Include censored observations in the number. `TRUE` by default. If
  `FALSE` then the number of observed events is returned. See
  [`BIC.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/BIC.flexsurvreg.md)
  for a discussion of the issues with defining the sample size for
  censored data.

- ...:

  Further arguments passed to or from other methods. Currently unused.

## Value

This returns the `mod$N` component of the fitted model object `mod`. See
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
for full documentation of all components.

## Details

By default, this matches the behaviour of the `nobs` method for
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) objects,
including both censored and uncensored observations.

If a weighted `flexsurvreg` analysis was done, then this function
returns the sum of the weights.

## See also

[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
[`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Author

C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
