# Version ?? (??)

* `B` argument added to `predict.flexsurvreg` (#197).


# Version 2.3.2 (2024-08-16)

* Vignette source included to satisfy CRAN checks.

* Fix of prediction with `spline="splines2ns"`.


# Version 2.3.1 (2024-07-14)

* Fix for simulating from semi-Markov models with `tcovs` option,
  which was not completing for some model configurations.

* Fix for `simulate.flexsurvreg()` with vector `start` argument
  (#192).
  

# Version 2.3 (2024-04-21)

* Analytic Hessian calculation for models where this is possible, that
  is, Weibull, Gompertz, exponential, and spline models in hazard and
  odds scales.

* Analytic gradient calculation for Weibull proportional hazards models.

* Vignette now firmly warns against using flexsurv with time-dependent
  covariates (#176).

* New argument `spline="splines2ns"` can now be specified to use an
  orthogonal spline basis in `flexsurvspline()`.

* Weighted likelihood for relative survival models now implemented
  consistently with other models, as a weighted sum of individual
  log-likelihoods.

* `standsurv` now returns results in the same order of times `t`
  as given by the user, for consistency with `summary.flexsurvreg`.

* Quantiles of standardised survival now available in `standsurv`.

* Non-default factor contrasts now handled.

* `pmatrix.simfs` can now accept a vector of times `t` and has a
  `tidy` output option.

* `BIC` and `AICc` functions added.

* Column name of `predict()` output changed from `"time"` to `"eval_time"`,
  for consistency with tidymodels update.

* Default value for `t` now chosen in `hr_flexsurvreg`.

* `coxsnell.flexsurvreg` now handles delayed entry.

* Warning given if the name of a location parameter is included in the
  ancillary part of the model specification.

* Fix for computing quantiles for custom distributions (#187).

* Default spline knots for models with interval censoring are now based on quantiles of times that include the interval midpoints, as well as the event times.

* Estimation for models with left and interval censoring now uses numeric derivatives, rather than analytic derivatives that ignore left or interval censored observations.

Thank you to all who have contributed code for this version: @mikesweeting @stephematician @ndunnewind @mattwarkentin @hfrick @kkmann; or reported issues: @anddis @irtimmins @sbihorel @zou-ims @aghaynes @huftis @mafed @hezht3 @sebffischer (and anyone else who reported issues via email).


# Version 2.2.2 (2023-01-31)

* Allow unicode characters in vignette, to satisfy R CMD check on
  r-devel.


# Version 2.2.1 (2022-12-22)

* New `simulate.flexsurvreg` method to simulate data from a fitted
  `flexsurvreg()` or `flexsurvspline()` model.  Thanks to Mark Clements for
  help with this.

* Fix of bug for `summary()` method with type = `"quantile"` or `"median"`
  and left-truncation in the prediction (`start` > 0). 

* Correction to the examples and interpretation of Cox-Snell residuals.


# Version 2.2 (2022-06-16)

* New function `standsurv` for survival and hazards
  standardised over an observed distribution for covariates. 
  Contributed by Michael Sweeting <mikesweeting79@gmail.com>.
  See the new vignette about it.

* New function `hr_flexsurvreg` to compute the hazard ratio from a
  fitted `flexsurvreg()` or `flexsurvspline()` model as a function of time,
  with confidence intervals.

* `summary.flexsurvreg` has been rewritten to make it cleaner and
  faster.  User-visible changes: 

  - Custom summary functions in `summary.flexsurvreg` must now be
  vectorised.

  - A new argument `cross` specifies whether to compute summaries for
  all combinations of times t and covariate values, or for covariate
  values matched with corresponding times t, in custom summary
  functions.
  
  - Covariate names when `tidy=FALSE` now consistently don't include
  spaces.

* Better handling of NAs in summary and prediction functions. Thanks
  to Matthew Warkentin.

* New functions for Cox-Snell residuals: `coxsnell_flexsurvreg` or
  residuals(..., type="coxsnell").

* `rmst_generic`, `mean_survspline`, `rmst_survspline` and related functions
  (e.g. `mean_survspline1`) handle alternative parameter values in a
  vectorised way.

* Allow output functions to work on models that have been stripped of
  data with `x$data <- NULL`, if possible.


# Version 2.1 (2021-09-13)

* Fix of a bug that affected models with baseline hazard offsets and
  exponential, Weibull, Gompertz and hazard/odds-based spline models,
  where convergence may have been falsely reported due to incorrect
  likelihood derivatives. 
  
* New vignette section and better help page documentation about relative
  survival models.

* `start` parameter added to `predict.flexsurvreg`.

* New function `pdf_flexsurmix` for the fitted density function in a
  `flexsurvmix` model.

* Fix for tidy method with one-parameter exponential models.

* Fixes for `h/Hgamma` and `h/Hlnorm` with `log = TRUE`.

* Minor numerical improvements to `h/H` functions for some distributions.

* Bug fix for `simfs_bytrans` when some transitions don't happen.

* `NaNs produced` warnings should occur less often during fitting.


# Version 2.0 (2021-02-22)

* A new class of multi-state models based on mixtures (Larson and
  Dinse 1985).  A new vignette on multi-state modelling describes this
  model class and contrasts it with standard (cause-specific hazards)
  multi-state models.

* Different parametric families are now supported for different
  transitions in multi-state models.

* New function `fmsm` allows a list of flexsurvreg objects to be
  associated with a particular transition structure matrix, to create
  a multi-state model. 

* New function `simfinal_fmsm` to summarise times and probabilities of
  final absorbing events in multi-state models, using simulation.

* More features for right-truncated data.  Individual-level
  right-truncation times supported with new "rtrunc" argument to
  flexsurvreg and flexsurvspline.  A comparable non-parametric
  estimator for right-truncated data is provided in a new function,
  `survrtrunc`.  Alternative parametric estimators, which make use
  of the time of an initiating event, are provided in a new function,
  `flexsurvrtrunc`.

* A new vignette describes properties of the different built-in
   parametric distributions in more detail. 

* `qgeneric` is now vectorised, thanks to the `vuniroot` function imported
  from the `rstpm2` package by Mark Clements.  This massively boosts the
  speed of `rsurvspline`, hence speeding up simulations from
  spline-based multi-state models.

* `pmatrix.fs` can now calculate transition probabilities conditionally
  on the the transition being to one of a subset of the states.

* "tidy" argument to `pmatrix.fs` for tidy data frame output. 

* "tidy" argument to `sim.fmsm` for returning simulations in tidy data
  frame format with one row per transition, and associated function
  `simfs_bytrans`.

* Bootstrapping function `bootci.fmsm` made available for users to get
  confidence intervals / distributions for their own flexsurv model
  output functions.

* Parallel processing capability for bootstrap confidence intervals.

* Distribution and mean functions for the Royston-Parmar model named
  like `dsurvspline2`, `psurvspline4` and so on, with one argument per
  parameter, rather than all parameters collected in a single
  argument, going up to 7 knots / 9 parameters.

* Return value from `pars.fmsm` is now a list rather than a matrix, even
  if the model family is the same for each transition.

* `summary.flexsurvreg` given a new argument "na.action" to control
  whether missing values in "newdata" are dropped.  Defaults to
  producing summaries of `NA` when there are missing values, while
  previously missing values were dropped.

* S3 methods have been added for the generics defined in the broom
  package.  These functions create tidy data frames containing the
  results of fitted models.  The new functions are `tidy.flexsurvreg`,
  `glance.flexsurvreg`, and `augment.flexsurvreg`.

* S3 methods have been added for the predict and residuals
  generics. `predict.flexsurvreg` has full support for all model
  outcomes supported by `summary.flexsurvreg`.  `residuals.flexsurvreg`
  currently only supports a naive difference between observed survival and
  predicted mean, neglecting censoring.

* Case weights accounted for in nonparametric survival and cumulative
   hazard estimates in `plot.flexsurvreg`.  Thanks to https://github.com/andbe. 


# Version 1.1.1 (2019-03-18)

* New type="quantiles" and type="link" for summary.flexsurvreg.
  Thanks to Leonardo Marques for the contribution.

* Allow different covariates per transition in multi-state models
  supplied as a list of flexsurvreg objects.  Thanks to David
  McAllister for the contribution.

* Bug fix for qllogis with lower.tail=FALSE.

* Bug fix for likelihood when all events are observed.

* Bug fix for "ylim" argument in plot method for survival or
  cumulative hazard.

* Allow dynamic symbols in C code.

* Various other minor code and doc fixes, see github commit history.


# Version 1.1 (2017-03-27)

* Substantial speed improvements in fitting most of the built-in
  models, from implementing their PDFs and CDFs in C++.  Thanks to
  Paul Metcalfe for contributing this.  Results may therefore differ
  from previous versions in edge cases.

* As a result the package now depends on Rcpp.

* Mean, median and restricted mean included as built-in functions in
  summary.flexsurvreg. Thanks to Jordan Amdahl for the contribution.

* Documentation migrated to roxygen.  Thanks to Paul Metcalfe for the
  contribution.

* Various minor bug fixes, see github commits.


# Version 1.0.2 (2016-09-26)

* Bug fix: "start" was being ignored by plot.flexsurvreg. Thanks to
  Ruth Keogh.

* Built-in distribution names are now case-insensitive.  Thanks to
  Jordan Amdahl.

* Fix for Weibull hazard function to avoid numeric instability.
  Thanks to Jordan Amdahl.

* Fix to hsurvspline when t includes 0.  Thanks to Jordan Amdahl.

* Vectorised parameters supported in qgeneric.


# Version 1.0.1 (2016-05-31)

* Bug fix: covariates were labelled wrongly in summary.flexsurvreg
  tidy output.  Thanks to Owain Saunders.


# Version 1.0.0 (2016-05-10)

* Version number bumped to 1.0.0 to accompany the publication of the
  vignette in Journal of Statistical Software.


# Version 0.7.1 (2016-03-24)

* Slightly more efficient likelihood calculations, and removed
  spurious warning from likelihood with interval censoring.

* Tests modified to work with the latest (and current) testthat.


# Version 0.7 (2015-11-13)

* flexsurvspline now allows the log cumulative hazard (or its
  alternatives) to be modelled as a spline function of time instead of
  log time.

* The routine for generating initial values in flexsurvspline has been
  improved.  This now obeys the constraint that the log cumulative
  hazard is increasing, thus avoiding errors from optim() when this
  wasn't satisfied.  Cox regression is used as a fallback to
  initialise covariate effects if this fails.

* As a result, flexsurv now depends on the "quadprog" package.

* dweibullPH and related functions give the Weibull distribution in
  proportional hazards parameterisation, and "weibullPH" is supported
  as a built-in model for flexsurvreg.

* Option to summary.flexsurvreg to return a tidy data frame.

* New "logliki" component in model objects, containing vector of
  log-likelihoods for each observation at the estimated or fixed
  parameters.

* Fix of various bugs with supplying "newdata" to summary functions
  (github issue #7).  The behaviour here should now be like predict.lm,
  e.g. variables in newdata that were originally factors should be
  supplied as factor or character, not numeric.

* Fix of bug that prevented plots being drawn by categorical
  covariates by default.

* Fix of bugs with spline models and no data censored, or all data
  interval censored (github issue #3).

* Fix of bugs with subsetting in flexsurvspline (github issue #6).


# Version 0.6 (2015-04-13)

* CRAN release.  Also includes the changes from Version 0.5.1.

* Full support for multi-state models fitted as a list of
  independent transition-specific models.

* New function pars.fmsm to return transition-specific parameters in
  multi-state models.

* Bug fix for empirical hazard plots with categorical covariates.
  Thanks to Milan Bouchet-Valat.


# Version 0.5.1 (2015-02-24)

* github-only release.

* Log-logistic distribution built in to flexsurvreg, and distribution
  functions provided.

* Bug fix in tcovs option in semi-Markov model simulation.

* "digits" argument supported by default model print function.  This
  is passed to "format" to format the parameter estimates, and
  defaults to 3.

* Bug fix in "events" printed output for interval censored
  data. Thanks to Sabrina Russo.

* pgompertz returns Inf for q=Inf, even for parameters denoting
  "living forever", since the CDF is P(X <= q) not P(X < q).  This
  affected some fits of the Gompertz distribution.


# Version 0.5 (2014-09-22)

* Major new release, so version number bumped from 0.3 to 0.5.

* New package vignettes: a user guide and a vignette of examples.

* Development moved from r-forge to
  https://github.com/chjackson/flexsurv-dev.

Spline models and ancillary covariates:

* Major rewrite of flexsurvspline.  This now works by calling
  flexsurvreg with a custom distribution written dynamically.  Models
  can now include covariate effects that vary as spline functions of
  time, by including covariates on "gamma1" or on any further
  parameters.

* New argument "anc" to flexsurvreg, as an alternative and preferred
  way of modelling covariates on ancillary parameters.

* New general utility "unroll.function" which converts a function with
  matrix arguments to the equivalent function with vector arguments.
  The new flexsurvspline works by unrolling "dsurvspline".

* Quantile, random number, hazard and cumulative hazard functions for
  spline distribution.

* Autogeneration of initial values for flexsurvspline now accounts for
  left-truncation. Thanks to Ana Borges for the report.

Other new modelling features:

* Several utilities for parametric multi-state modelling, including
  transition probabilities and simulation ("pmatrix.fs", "totlos.fs",
  "sim.fmsm","pmatrix.simfs", "totlos.simfs").  An "msfit.flexsurvreg"
  method also gives cumulative transition-specific hazards in the
  format of the "mstate" package.

* Interval censoring supported in the Surv() response.

* Relative survival models, using the "bhazard" argument to
  flexsurvreg to specify the expected mortality rate.

* flexsurvreg now uses survreg internally to fit Weibull, exponential
  and log-normal models, unless there is left-truncation.

Custom distributions:

* Custom distributions can be defined through the hazard function.
  This can be optionally supplemented with the cumulative hazard
  function, otherwise this is obtained by numerical integration.

* In custom distributions specified by the density, the cumulative
  distribution can now be omitted, and it will be calculated by
  numerical integration.

* New arguments "dfns" and "aux" to flexsurvreg, which can be used to
  supply custom distribution functions and arguments to pass to them.

* Document that density functions for custom distributions need "log"
  argument.

* Documented how to supply derivatives of custom distributions for use
  in optimisation.

Output functions:

* New "newdata" argument to summary.flexsurvreg and plot.flexsurvreg
  for an easier way of supplying covariate values.

* User-defined summary functions can be used in summary.flexsurvreg
  and plot.flexsurvreg as an alternative to survival, hazard or
  cumulative hazard.

* New function "normboot.flexsurvreg" to simulate parameters from the
  asymptotic normal distribution of their estimates.  Used for
  representing uncertainty in any function of the parameters.

* summary.flexsurvreg can be called with ci=FALSE to omit confidence
  intervals.

* "start" argument defaults to 0 for all prediction times in
  summary.flexsurvreg.

* Bug fix in summary.flexsurvreg for left-truncated models, which was
  returning probabilities > 1 before the truncation time.

* New model.frame and model.matrix methods to extract the data from
  fitted flexsurvreg objects.

* Accept vector X in summary.flexsurvreg, and give informative error
  if X in wrong format.  Thanks to Mark Danese.

* Extra arguments can be passed to both muhaz and plot.muhaz in
  plot.flexsurvreg(...,type="hazard",...)

Printed output:

* Use format() not signif() for printing flexsurvreg objects, to avoid
  spurious zero significant figures. Thanks to Kenneth Chen.

* Standard errors included in printed output of flexsurvreg, but only
  for parameters optimised on natural or log scales (which includes
  all built-in distributions).

Distribution functions:

* Bug fix in rgengamma for Q=0 (log-normal) and sigma not equal to 1.

* Don't warn for shape parameters being exactly zero in generalized
  gamma and F, just give NaN.

* basis() and fss() functions for the natural cubic spline basis made
  available to users.



# Version 0.3.1 (2014-02-14)

* R-forge only release.

* Distribution functions tidied up, making special value handling and
  vectorisation consistent. Hazard and cumulative hazard functions for
  all supported distributions.

* Vectors of different col, lwd and lty can be passed to
  plot.flexsurvreg for multiple fitted lines.  Thanks to Julia
  Sandberg for the report.


# Version 0.3 (2014-01-19)

* CRAN release. Includes changes from 0.2.1 to 0.2.3.


# Version 0.2.3 (2013-10-09)

* R-forge only release.

* Parameters other than the location parameter can now have covariates
  on them in flexsurvreg.  Thanks to Milan Bouchet-Valat for help with
  this.

* subset and na.action arguments in flexsurvreg and flexsurvspline.

* coef, vcov and confint methods for all fitted model objects.

* Distribution functions for generalized gamma, generalized F, and
  Gompertz, now allow all parameters to be vectorised.

* Bug fix in analytic derivatives for Weibull.

* Restored print output introduced in 0.1.2 which had been
  accidentally removed in 0.1.5.


# Version 0.2.2 (2013-07-26)

* R-forge only release.

* Case weights supported in flexsurvreg and flexsurvspline.


# Version 0.2.1 (2013-07-03)

* R-forge only release.

* Default left truncation times were being set wrongly for
  user-supplied times in summary.flexsurvreg, giving wrong confidence
  intervals.  These now default to 0.

* Confidence intervals set to 1 for t=0 under spline models. Thanks to
  Paul Pynsent.

* dgompertz,dgengamma,dgengamma.orig,dgenf,dgenf.orig fixed to return
  -Inf instead of 0 when density is zero and log=TRUE.  Thanks to Gao
  Zheng.


# Version 0.2 (2013-05-13)

* New summary() method for fitted flexsurvreg and flexsurvspline model
  objects gives fitted survival, cumulative hazard or hazard curves,
  with confidence intervals mosly computed by a simulation method.

* This allows plot.flexsurvreg to plot confidence intervals for the
  fitted survival, hazard or cumulative hazard.

* Left-truncated survival observations are supported in flexsurvreg
  and flexsurvspline.

* New psurvspline and dsurvspline functions giving distribution and
  density function for the spline model.

* Analytic derivatives used in optimisation for spline (odds and
  hazard scale, not normal), exponential, Weibull and Gompertz models.

* Default to BFGS optimisation method, which uses derivatives where
  available and should be much faster, instead of Nelder-Mead.

* Work around NaN warnings from spline models presumably due to
  parameters violating implicit constraints.

* If "knots" specified, boundary knots set to min/max of uncensored
  times, not all times, to match results when "k" is specified.
  Thanks to Paul Pynsent.


# Version 0.1.5 (2012-08-29)

* Data are now stored in fitted flexsurvreg and flexsurvspline model
  objects, avoiding environment search errors and allowing package
  functions to be called within other functions.  Thanks to Hanna
  Daniel for the report.

* Gompertz documentation clarified for the case when there is a chance
  of living forever.  qgompertz and rgompertz now return Inf in these
  cases, with no warning, instead of NaN.  Thanks to Michael Sweeting.


# Version 0.1.4 (2012-03-22)

* maxt argument in plot.flexsurvreg.

* Plots no longer complain if data named "dat".

* Corrected wrong bug fix from Version 0.1.3 for transforming
  parameter estimates in output when fixedpars=TRUE.

* AIC penalty corrected for models with some fixed parameters.

* qgengamma corrected for parameter Q<0.  Thanks to Benn Ackley.


# Version 0.1.3 (2012-01-17)

* No longer complains about invalid initial values when there are zero
  survival times.

* Don't transform parameter estimates in output when fixedpars=TRUE.

* Checking functions for distribution utilities don't complain about
  vectorised parameter values.


# Version 0.1.2 (2011-11-08)

* Initial CRAN release.

* More features in print output for flexsurvreg and flexsurvspline
  models.


# Version 0.1.1 (2011-04-19)

* Fix of drop=FALSE bug in flexsurvspline.inits which caused
  flexsurvspline to fail with single covariates.


# Version 0.1 (2011-03-14)

Initial release
