# Forms a tidy data frame for plotting the fit of parametric mixture multi-state models against nonparametric estimates

This computes Aalen-Johansen estimates of the probability of occupying
each state at a series of times, using
[`ajfit`](http://chjackson.github.io/flexsurv-dev/reference/ajfit.md).
The equivalent estimates from the parametric model are then produced
using
[`p_flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/p_flexsurvmix.md),
and concatenated with the nonparametric estimates to form a tidy data
frame. This data frame can then simply be plotted using
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Usage

``` r
ajfit_flexsurvmix(x, maxt = NULL, startname = "Start", B = NULL)
```

## Arguments

- x:

  Fitted model returned by
  [`flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvmix.md).

- maxt:

  Maximum time to produce parametric estimates for. By default this is
  the maximum event time in the data, the maximum time we have
  nonparametric estimates for.

- startname:

  Label to give the state corresponding to "no event happened yet". By
  default this is `"Start"`.

- B:

  Number of simulation replications to use to calculate a confidence
  interval for the parametric estimates in
  [`p_flexsurvmix`](http://chjackson.github.io/flexsurv-dev/reference/p_flexsurvmix.md).
  Comparable intervals for the Aalen-Johansen estimates are returned if
  this is set. Otherwise if `B=NULL` then no intervals are returned.
