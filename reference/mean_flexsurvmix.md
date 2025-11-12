# Mean times to events from a flexsurvmix model

This returns the mean of each event-specific parametric time-to-event
distribution in the mixture model, which is the mean time to event
conditionally on that event being the one that happens.

## Usage

``` r
mean_flexsurvmix(x, newdata = NULL, B = NULL)
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

Mean times to next event conditionally on each alternative event, given
the specified covariate values.
