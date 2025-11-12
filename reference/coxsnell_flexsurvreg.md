# Cox-Snell residuals from a parametric survival model

Cox-Snell residuals from a parametric survival model

## Usage

``` r
coxsnell_flexsurvreg(x)
```

## Arguments

- x:

  Object returned by
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  or
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  representing a fitted survival model

## Value

A data frame with a column called `est` giving the Cox-Snell residual,
defined as the fitted cumulative hazard at each data point. fitted
cumulative hazard at the given observed data point, and other columns
indicating the observation time, observed event status, and covariate
values defining the data at this point.

The cumulative hazards `est` should form a censored sample from an
Exponential(1). Therefore to check the fit of the model, plot a
nonparametric estimate of the cumulative hazard curve against a diagonal
line through the origin, which is the theoretical cumulative hazard
trajectory of the Exponential(1).

## Examples

``` r
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
  cs <- coxsnell_flexsurvreg(fitg)
  
  ## Model appears to fit well, with some small sample noise 
  surv <- survfit(Surv(cs$est, ovarian$fustat) ~ 1)
  plot(surv, fun="cumhaz")
  abline(0, 1, col="red")

  
```
