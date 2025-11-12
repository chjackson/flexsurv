# Probabilities of final states in a flexible parametric competing risks model

This requires the model to be Markov, and is not valid for semi-Markov
models, as it works by wrapping
[`pmatrix.fs`](http://chjackson.github.io/flexsurv-dev/reference/pmatrix.fs.md)
to calculate the transition probability over a very large time. As it
also works on a `fmsm` object formed from transition-specific
time-to-event models, it therefore only works on competing risks models,
defined by just one starting state with multiple destination states
representing competing events. For these models, this function returns
the probability governing which competing event happens next. However
this function simply wraps `pmatrix.fs`, so for other models,
`pmatrix.fs` or `pmatrix.simfs` can be used with a large forecast time
`t`.

## Usage

``` r
pfinal_fmsm(x, newdata = NULL, fromstate, maxt = 1e+05, B = 0, cores = NULL)
```

## Arguments

- x:

  Object returned by
  [`fmsm`](http://chjackson.github.io/flexsurv-dev/reference/fmsm.md),
  representing a multi-state model formed from transition-specific
  time-to-event models fitted by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- newdata:

  Data frame of covariate values, with one column per covariate, and one
  row per alternative value.

- fromstate:

  State from which to calculate the transition probability state. This
  should refer to the name of a row of the transition matrix
  `attr(x,trans)`.

- maxt:

  Large time to use for forecasting final state probabilities. The
  transition probability from zero to this time is used. Note `Inf` will
  not work. The default is `100000`.

- B:

  Number of simulations to use to calculate 95% confidence intervals
  based on the asymptotic normal distribution of the basic parameter
  estimates. If `B=0` then no intervals are calculated.

- cores:

  Number of processor cores to use. If `NULL` (the default) then a
  single core is used.

## Value

A data frame with one row per covariate value and destination state,
giving the state in column `state`, and probability in column `val`.
Additional columns `lower` and `upper` for the confidence limits are
returned if `B=0`.
