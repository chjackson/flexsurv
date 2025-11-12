# Constructor for a mixture multi-state model based on flexsurvmix

Constructor for a mixture multi-state model based on flexsurvmix

## Usage

``` r
fmixmsm(...)
```

## Arguments

- ...:

  Named arguments. Each argument should be a fitted model as returned by
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).
  The name of each argument names the starting state for that model.

## Value

A list of
[`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md)
objects, with the following attribute(s):

`pathways` A list of all potential pathways until absorption, for models
without cycles. For models with cycles this will have an element
`has_cycle=TRUE`, plus the pathways discovered before the function found
the cycle and gave up.
