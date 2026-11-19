# Plot nonparametric estimates of survival from right-truncated data.

`plot.survrtrunc` creates a new plot, while `lines.survrtrunc` adds
lines to an exising plot.

## Usage

``` r
# S3 method for class 'survrtrunc'
plot(x, ...)

# S3 method for class 'survrtrunc'
lines(x, ...)
```

## Arguments

- x:

  Object of class `"survrtrunc"` as returned by
  [`survrtrunc`](https://chjackson.github.io/flexsurv/reference/survrtrunc.md).

- ...:

  Other arguments to be passed to
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html)
  or
  [`lines.survfit`](https://rdrr.io/pkg/survival/man/lines.survfit.html).
