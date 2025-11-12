# Simulate times to competing events from a mixture multi-state model

Simulate times to competing events from a mixture multi-state model

## Usage

``` r
simt_flexsurvmix(x, newdata = NULL, n)
```

## Arguments

- x:

  Fitted model object returned from
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- newdata:

  Data frame or list of covariate values. If omitted for a model with
  covariates, a default is used, defined by all combinations of factors
  if the only covariates in the model are factors, or all covariate
  values of zero if there are any non-factor covariates in the model.

- n:

  Number of simulations

## Value

Data frame with `n*m` rows and a column for each competing event, where
`m` is the number of alternative covariate values, that is the number of
rows of `newdata`. The simulated time represents the time to that event
conditionally on that event being the one that occurs. This function
doesn't simulate which event occurs.
