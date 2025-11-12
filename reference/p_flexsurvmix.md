# Transition probabilities from a flexsurvmix model

These quantities are variously known as transition probabilities, or
state occupancy probabilities, or values of the "cumulative incidence"
function, or values of the "subdistribution" function. They are the
probabilities that an individual has experienced an event of a
particular kind by time `t`.

## Usage

``` r
p_flexsurvmix(x, newdata = NULL, startname = "start", t = 1, B = NULL)
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

- startname:

  Name of the state where individuals start. This considers the model as
  a multi-state model where people start in this state, and may
  transition to one of the competing events.

- t:

  Vector of times `t` to calculate the probabilities of transition by.

- B:

  Number of simulations to use to compute 95% confidence intervals,
  based on the asymptotic multivariate normal distribution of the basic
  parameter estimates. If `B=NULL` then intervals are not computed.

## Value

A data frame with transition probabilities by time, covariate value and
destination state.

## Details

Note that "cumulative incidence" is a misnomer, as "incidence" typically
means a hazard, and the quantities computed here are not cumulative
hazards, but probabilities.
