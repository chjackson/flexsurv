# Akaike's information criterion from a flexible parametric multistate model

Defined as the sum of the AICs of the transition-specific models.

## Usage

``` r
# S3 method for class 'fmsm'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  Object returned by
  [`fmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmsm.md)
  representing a multistate model.

- ...:

  Further arguments (currently unused).

- k:

  Penalty applied to number of parameters (defaults to the standard 2).

## Value

The sum of the AICs of the transition-specific models.
