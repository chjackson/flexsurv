# Flexible parametric mixture models for times to competing events

In a mixture model for competing events, an individual can experience
one of a set of different events. We specify a model for the probability
that they will experience each event before the others, and a model for
the time to the event conditionally on that event occurring first.

## Usage

``` r
flexsurvmix(
  formula,
  data,
  event,
  dists,
  pformula = NULL,
  anc = NULL,
  partial_events = NULL,
  initp = NULL,
  inits = NULL,
  fixedpars = NULL,
  dfns = NULL,
  method = "direct",
  em.control = NULL,
  optim.control = NULL,
  aux = NULL,
  sr.control = survreg.control(),
  integ.opts,
  hess.control = NULL,
  ...
)
```

## Arguments

- formula:

  Survival model formula. The left hand side is a `Surv` object
  specified as in
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  This may define various kinds of censoring, as described in
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html). Any covariates
  on the right hand side of this formula will be placed on the location
  parameter for every component-specific distribution. Covariates on
  other parameters of the component-specific distributions may be
  supplied using the `anc` argument.

  Alternatively, `formula` may be a list of formulae, with one component
  for each alternative event. This may be used to specify different
  covariates on the location parameter for different components.

  A list of formulae may also be used to indicate that for particular
  individuals, different events may be observed in different ways, with
  different censoring mechanisms. Each list component specifies the data
  and censoring scheme for that mixture component.

  For example, suppose we are studying people admitted to hospital,and
  the competing states are death in hospital and discharge from
  hospital. At time t we know that a particular individual is still
  alive, but we do not know whether they are still in hospital, or have
  been discharged. In this case, if the individual were to die in
  hospital, their death time would be right censored at t. If the
  individual will be (or has been) discharged before death, their
  discharge time is completely unknown, thus interval-censored on
  (0,Inf). Therefore, we need to store different event time and status
  variables in the data for different alternative events. This is
  specified here as

  `formula = list("discharge" = Surv(t1di, t2di, type="interval2"), "death" = Surv(t1de, status_de))`

  where for this individual, `(t1di, t2di) = (0, Inf)` and
  `(t1de, status_de) = (t, 0)`.

  The "dot" notation commonly used to indicate "all remaining variables"
  in a formula is not supported in `flexsurvmix`.

- data:

  Data frame containing variables mentioned in `formula`, `event` and
  `anc`.

- event:

  Variable in the data that specifies which of the alternative events is
  observed for which individual. If the individual's follow-up is
  right-censored, or if the event is otherwise unknown, this variable
  must have the value `NA`.

  Ideally this should be a factor, since the mixture components can then
  be easily identified in the results with a name instead of a number.
  If this is not already a factor, it is coerced to one. Then the levels
  of the factor define the required order for the components of the list
  arguments `dists`, `anc`, `inits` and `dfns`. Alternatively, if the
  components of the list arguments are named according to the levels of
  `event`, then the components can be arranged in any order.

- dists:

  Vector specifying the parametric distribution to use for each
  component. The same distributions are supported as in
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- pformula:

  Formula describing covariates to include on the component membership
  proabilities by multinomial logistic regression. The first component
  is treated as the baseline.

  The "dot" notation commonly used to indicate "all remaining variables"
  in a formula is not supported.

- anc:

  List of component-specific lists, of length equal to the number of
  components. Each component-specific list is a list of formulae
  representing covariate effects on parameters of the distribution.

  If there are covariates for one component but not others, then a list
  containing one null formula on the location parameter should be
  supplied for the component with no covariates, e.g `list(rate=~1)` if
  the location parameter is called `rate`.

  Covariates on the location parameter may also be supplied here instead
  of in `formula`. Supplying them in `anc` allows some components but
  not others to have covariates on their location parameter. If a
  covariate on the location parameter was provided in `formula`, and
  there are covariates on other parameters, then a null formula should
  be included for the location parameter in `anc`, e.g `list(rate=~1)`

- partial_events:

  List specifying the factor levels of `event` which indicate knowledge
  that an individual will not experience particular events, but may
  experience others. The names of the list indicate codes that indicate
  partial knowledge for some individuals. The list component is a
  vector, which must be a subset of `levels(event)` defining the events
  that a person with the corresponding event code may experience.

  For example, suppose there are three alternative events called
  `"disease1"`,`"disease2"` and `"disease3"`, and for some individuals
  we know that they will not experience `"disease2"`, but they may
  experience the other two events. In that case we must create a new
  factor level, called, for example `"disease1or3"`, and set the value
  of `event` to be `"disease1or3"` for those individuals. Then we use
  the `"partial_events"` argument to tell `flexsurvmix` what the
  potential events are for individuals with this new factor level.

  `partial_events = list("disease1or3" = c("disease1","disease3"))`

- initp:

  Initial values for component membership probabilities. By default,
  these are assumed to be equal for each component.

- inits:

  List of component-specific vectors. Each component-specific vector
  contains the initial values for the parameters of the
  component-specific model, as would be supplied as the `inits` argument
  of
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
  By default, a heuristic is used to obtain initial values, which
  depends on the parametric distribution being used, but is usually
  based on the empirical mean and/or variance of the survival times.

- fixedpars:

  Indexes of parameters to fix at their initial values and not optimise.
  Arranged in the order: baseline mixing probabilities, covariates on
  mixing probabilities, time-to-event parameters by mixing component.
  Within mixing components, time-to-event parameters are ordered in the
  same way as in
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

  If `fixedpars=TRUE` then all parameters will be fixed and the function
  simply calculates the log-likelihood at the initial values.

  Not currently supported when using the EM algorithm.

- dfns:

  List of lists of user-defined distribution functions, one for each
  mixture component. Each list component is specified as the `dfns`
  argument of
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

- method:

  Method for maximising the likelihood. Either `"em"` for the EM
  algorithm, or `"direct"` for direct maximisation.

- em.control:

  List of settings to control EM algorithm fitting. The only options
  currently available are

  `trace` set to 1 to print the parameter estimates at each iteration of
  the EM algorithm

  `reltol` convergence criterion. The algorithm stops if the log
  likelihood changes by a relative amount less than `reltol`. The
  default is the same as in
  [`optim`](https://rdrr.io/r/stats/optim.html), that is,
  `sqrt(.Machine$double.eps)`.

  `var.method` method to compute the covariance matrix. `"louis"` for
  the method of Louis (1982), or `"direct"`for direct numerical
  calculation of the Hessian of the log likelihood.

  `optim.p.control` A list that is passed as the `control` argument to
  `optim` in the M step for the component membership probability
  parameters. The optimisation in the M step for the time-to-event
  parameters can be controlled by the `optim.control` argument to
  `flexsurvmix`.

  For example, `em.control = list(trace=1, reltol=1e-12)`.

- optim.control:

  List of options to pass as the `control` argument to
  [`optim`](https://rdrr.io/r/stats/optim.html), which is used by
  `method="direct"` or in the M step for the time-to-event parameters in
  `method="em"`. By default, this uses `fnscale=10000` and
  `ndeps=rep(1e-06,p)` where `p` is the number of parameters being
  estimated, unless the user specifies these options explicitly.

- aux:

  A named list of other arguments to pass to custom distribution
  functions. This is used, for example, by
  [`flexsurvspline`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md)
  to supply the knot locations and modelling scale (e.g. hazard or
  odds). This cannot be used to fix parameters of a distribution — use
  `fixedpars` for that.

- sr.control:

  For the models which use
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) to find the
  maximum likelihood estimates (Weibull, exponential, log-normal), this
  list is passed as the `control` argument to
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

- integ.opts:

  List of named arguments to pass to
  [`integrate`](https://rdrr.io/r/stats/integrate.html), if a custom
  density or hazard is provided without its cumulative version. For
  example,

  `integ.opts = list(rel.tol=1e-12)`

- hess.control:

  List of options to control covariance matrix computation. Available
  options are:

  `numeric`. If `TRUE` then numerical methods are used to compute the
  Hessian for models where an analytic Hessian is available. These
  models include the Weibull (both versions), exponential, Gompertz and
  spline models with hazard or odds scale. The default is to use the
  analytic Hessian for these models. For all other models, numerical
  methods are always used to compute the Hessian, whether or not this
  option is set.

  `tol.solve`. The tolerance used for
  [`solve`](https://rdrr.io/r/base/solve.html) when inverting the
  Hessian (default `.Machine$double.eps`)

  `tol.evalues` The accepted tolerance for negative eigenvalues in the
  covariance matrix (default `1e-05`).

  The Hessian is positive definite, thus invertible, at the maximum
  likelihood. If the Hessian computed after optimisation convergence
  can't be inverted, this is either because the converged result is not
  the maximum likelihood (e.g. it could be a "saddle point"), or because
  the numerical methods used to obtain the Hessian were inaccurate. If
  you suspect that the Hessian was computed wrongly enough that it is
  not invertible, but not wrongly enough that the nearest valid inverse
  would be an inaccurate estimate of the covariance matrix, then these
  tolerance values can be modified (reducing `tol.solve` or increasing
  `tol.evalues`) to allow the inverse to be computed.

- ...:

  Optional arguments to the general-purpose optimisation routine
  [`optim`](https://rdrr.io/r/stats/optim.html). For example, the BFGS
  optimisation algorithm is the default in
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
  but this can be changed, for example to `method="Nelder-Mead"` which
  can be more robust to poor initial values. If the optimisation fails
  to converge, consider normalising the problem using, for example,
  `control=list(fnscale = 2500)`, for example, replacing 2500 by a
  number of the order of magnitude of the likelihood. If 'false'
  convergence is reported with a non-positive-definite Hessian, then
  consider tightening the tolerance criteria for convergence. If the
  optimisation takes a long time, intermediate steps can be printed
  using the `trace` argument of the control list. See
  [`optim`](https://rdrr.io/r/stats/optim.html) for details.

## Value

List of objects containing information about the fitted model. The
important one is `res`, a data frame containing the parameter estimates
and associated information.

## Details

This differs from the more usual "competing risks" models, where we
specify "cause-specific hazards" describing the time to each competing
event. This time will not be observed for an individual if one of the
competing events happens first. The event that happens first is defined
by the minimum of the times to the alternative events.

The `flexsurvmix` function fits a mixture model to data consisting of a
single time to an event for each individual, and an indicator for what
type of event occurs for that individual. The time to event may be
observed or censored, just as in
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
and the type of event may be known or unknown. In a typical application,
where we follow up a set of individuals until they experience an event
or a maximum follow-up time is reached, the event type is known if the
time is observed, and the event type is unknown when follow-up ends and
the time is right-censored.

The model is fitted by maximum likelihood, either directly or by using
an expectation-maximisation (EM) algorithm, by wrapping
[`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
to compute the likelihood or to implement the E and M steps.

Some worked examples are given in the package vignette about multi-state
modelling, which can be viewed by running
`vignette("multistate", package="flexsurv")`.

## References

Jackson, C. H. and Tom, B. D. M. and Kirwan, P. D. and Mandal, S. and
Seaman, S. R. and Kunzmann, K. and Presanis, A. M. and De Angelis, D.
(2022) A comparison of two frameworks for multi-state modelling, applied
to outcomes after hospital admissions with COVID-19. Statistical Methods
in Medical Research 31(9) 1656-1674.

Larson, M. G., & Dinse, G. E. (1985). A mixture model for the regression
analysis of competing risks data. Journal of the Royal Statistical
Society: Series C (Applied Statistics), 34(3), 201-211.

Lau, B., Cole, S. R., & Gange, S. J. (2009). Competing risk regression
models for epidemiologic data. American Journal of Epidemiology, 170(2),
244-256.
