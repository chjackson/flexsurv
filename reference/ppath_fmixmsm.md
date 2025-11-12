# Probability of each pathway taken through a mixture multi-state model

Probability of each pathway taken through a mixture multi-state model

## Usage

``` r
ppath_fmixmsm(x, newdata = NULL, final = FALSE, B = NULL)
```

## Arguments

- x:

  Object returned by
  [`fmixmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmixmsm.md),
  representing a multi-state model built from piecing together mixture
  models fitted by
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- newdata:

  Data frame or list of covariate values. If omitted for a model with
  covariates, a default is used, defined by all combinations of factors
  if the only covariates in the model are factors, or all covariate
  values of zero if there are any non-factor covariates in the model.

- final:

  If `TRUE` then the probabilities of pathways with the same final state
  are added together, to produce the probability of each ultimate
  outcome or absorbing state from the multi-state model.

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

## Value

Data frame of pathway probabilities by covariate value and pathway.
