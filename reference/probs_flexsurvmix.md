# Probabilities of competing events from a flexsurvmix model

Probabilities of competing events from a flexsurvmix model

## Usage

``` r
probs_flexsurvmix(x, newdata = NULL, B = NULL)
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

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

## Value

A data frame containing the probability that each of the competing
events will occur next, by event and by any covariate values specified
in `newdata`.
