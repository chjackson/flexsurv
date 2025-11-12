# Fitted densities for times to events in a flexsurvmix model

This returns an estimate of the probability density for the time to each
competing event, at a vector of times supplied by the user.

## Usage

``` r
pdf_flexsurvmix(x, newdata = NULL, t = NULL)
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

- t:

  Vector of times at which to evaluate the probability density

## Value

A data frame with each row giving the fitted density `dens` for a
combination of covariate values, time and competing event.
