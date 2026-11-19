# Model frame from a flexsurvmix object

Returns a list of data frames, with each component containing the data
that were used for the model fit for that mixture component.

## Usage

``` r
# S3 method for class 'flexsurvmix'
model.frame(formula, ...)
```

## Arguments

- formula:

  Fitted model object from
  [`flexsurvmix`](https://chjackson.github.io/flexsurv/reference/flexsurvmix.md).

- ...:

  Further arguments (currently unused).

## Value

A list of data frames
