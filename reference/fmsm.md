# Construct a multi-state model from a set of parametric survival models

Construct a multi-state model from a set of parametric survival models

## Usage

``` r
fmsm(..., trans)
```

## Arguments

- ...:

  Objects returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  representing fitted survival models.

- trans:

  A matrix of integers specifying which models correspond to which
  transitions. The \\r,s\\ entry is `i` if the \\i\\th argument
  specified in `...` is the model for the state \\r\\ to state \\s\\
  transition. The entry should be `NA` if the transition is disallowed.

## Value

A list containing the objects given in `...`, and with attributes
`"trans"` and `"statenames"` defining the transition structure matrix
and state names, and with list components named to describe the
transitions they correspond to. If any of the arguments in `...` are
named, then these are used to define the transition names, otherwise
default names are chosen based on the state names.
