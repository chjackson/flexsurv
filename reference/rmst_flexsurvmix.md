# Restricted mean times to events from a flexsurvmix model

This returns the restricted mean of each event-specific parametric
time-to-event distribution in the mixture model, which is the mean time
to event conditionally on that event being the one that happens, and
conditionally on the event time being less than some time horizon `tot`.

## Usage

``` r
rmst_flexsurvmix(x, newdata = NULL, tot = Inf, B = NULL)
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

- tot:

  Time horizon to compute the restricted mean until.

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

## Value

Restricted mean times to next event conditionally on each alternative
event, given the specified covariate values.
