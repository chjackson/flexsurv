# Simulate and summarise final outcomes from a flexible parametric multi-state model

Estimates the probability of each final outcome ("absorbing" state), and
the mean and quantiles of the time to that outcome for people who
experience it, by simulating a large sample of individuals from the
model. This can be used for both Markov and semi-Markov models.

## Usage

``` r
simfinal_fmsm(
  x,
  newdata = NULL,
  probs = c(0.025, 0.5, 0.975),
  t = 1000,
  M = 1e+05,
  B = 0,
  cores = NULL
)
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

- probs:

  Quantiles to calculate, by default, `c(0.025, 0.5, 0.975)` for a
  median and 95% interval.

- t:

  Maximum time to simulate to, passed to
  [`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md),
  so that the summaries are taken from the subset of individuals in the
  simulated data who are in the absorbing state at this time.

- M:

  Number of individuals to simulate.

- B:

  Number of simulations to use to calculate 95% confidence intervals
  based on the asymptotic normal distribution of the basic parameter
  estimates. If `B=0` then no intervals are calculated.

- cores:

  Number of processor cores to use. If `NULL` (the default) then a
  single core is used.

## Value

A tidy data frame with rows for each combination of covariate values and
quantity of interest. The quantity of interest is identified in the
column `quantity`, and the value of the quantity is in `val`, with
additional columns `lower` and `upper` giving 95% confidence intervals
for the quantity, if `B>0`.

## Details

For a competing risks model, i.e. a model defined by just one starting
state and multiple destination states representing competing events,
this returns the probability governing the next event that happens, and
the distribution of the time to each event conditionally on that event
happening.
