# Tidy a standsurv object.

This function is used internally by `standsurv` and tidy data.frames are
automatically returned by the function.

## Usage

``` r
# S3 method for class 'standsurv'
tidy(x, ...)
```

## Arguments

- x:

  A standsurv object.

- ...:

  Not currently used.

## Value

Returns additional tidy data.frames (tibbles) stored as attributes named
standpred_at and standpred_contrast.
