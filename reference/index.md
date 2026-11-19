# Package index

## Fitting survival models

- [`flexsurvreg()`](https://chjackson.github.io/flexsurv/reference/flexsurvreg.md)
  : Flexible parametric regression for time-to-event data
- [`flexsurvspline()`](https://chjackson.github.io/flexsurv/reference/flexsurvspline.md)
  : Flexible survival regression using the Royston/Parmar spline model.
- [`flexsurvrtrunc()`](https://chjackson.github.io/flexsurv/reference/flexsurvrtrunc.md)
  : Flexible parametric models for right-truncated, uncensored data
  defined by times of initial and final events.
- [`survrtrunc()`](https://chjackson.github.io/flexsurv/reference/survrtrunc.md)
  : Nonparametric estimator of survival from right-truncated, uncensored
  data

## Getting results from survival models

- [`summary(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/summary.flexsurvreg.md)
  : Summaries of fitted flexible survival models
- [`standsurv()`](https://chjackson.github.io/flexsurv/reference/standsurv.md)
  : Marginal survival and hazards of fitted flexsurvreg models
- [`coef(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/coef.flexsurvreg.md)
  : Extract model coefficients from fitted flexible survival models
- [`normboot.flexsurvreg()`](https://chjackson.github.io/flexsurv/reference/normboot.flexsurvreg.md)
  : Simulate from the asymptotic normal distribution of parameter
  estimates.
- [`hr_flexsurvreg()`](https://chjackson.github.io/flexsurv/reference/hr_flexsurvreg.md)
  : Hazard ratio as a function of time from a parametric survival model
- [`vcov(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/vcov.flexsurvreg.md)
  : Variance-covariance matrix from a flexsurvreg model
- [`simulate(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/simulate.flexsurvreg.md)
  : Simulate censored time-to-event data from a fitted flexsurvreg model
- [`plot(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/plot.flexsurvreg.md)
  : Plots of fitted flexible survival models
- [`lines(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/lines.flexsurvreg.md)
  : Add fitted flexible survival curves to a plot
- [`plot(`*`<standsurv>`*`)`](https://chjackson.github.io/flexsurv/reference/plot.standsurv.md)
  : Plot standardized metrics from a fitted flexsurv model
- [`summary(`*`<flexsurvrtrunc>`*`)`](https://chjackson.github.io/flexsurv/reference/summary.flexsurvrtrunc.md)
  : Summarise quantities of interest from fitted flexsurvrtrunc models
- [`plot(`*`<survrtrunc>`*`)`](https://chjackson.github.io/flexsurv/reference/plot_survtrunc.md)
  [`lines(`*`<survrtrunc>`*`)`](https://chjackson.github.io/flexsurv/reference/plot_survtrunc.md)
  : Plot nonparametric estimates of survival from right-truncated data.

## Multi-state models based on transition-specific hazards

- [`fmsm()`](https://chjackson.github.io/flexsurv/reference/fmsm.md) :
  Construct a multi-state model from a set of parametric survival models
- [`pars.fmsm()`](https://chjackson.github.io/flexsurv/reference/pars.fmsm.md)
  : Transition-specific parameters in a flexible parametric multi-state
  model
- [`pmatrix.fs()`](https://chjackson.github.io/flexsurv/reference/pmatrix.fs.md)
  : Transition probability matrix from a fully-parametric,
  time-inhomogeneous Markov multi-state model
- [`pmatrix.simfs()`](https://chjackson.github.io/flexsurv/reference/pmatrix.simfs.md)
  : Transition probability matrix from a fully-parametric, semi-Markov
  multi-state model
- [`totlos.fs()`](https://chjackson.github.io/flexsurv/reference/totlos.fs.md)
  : Total length of stay in particular states for a fully-parametric,
  time-inhomogeneous Markov multi-state model
- [`totlos.simfs()`](https://chjackson.github.io/flexsurv/reference/totlos.simfs.md)
  : Expected total length of stay in specific states, from a
  fully-parametric, semi-Markov multi-state model
- [`sim.fmsm()`](https://chjackson.github.io/flexsurv/reference/sim.fmsm.md)
  : Simulate paths through a fully parametric semi-Markov multi-state
  model
- [`simfs_bytrans()`](https://chjackson.github.io/flexsurv/reference/simfs_bytrans.md)
  : Reformat simulated multi-state data with one row per simulated
  transition
- [`bootci.fmsm()`](https://chjackson.github.io/flexsurv/reference/bootci.fmsm.md)
  : Bootstrap confidence intervals for flexsurv output functions
- [`pfinal_fmsm()`](https://chjackson.github.io/flexsurv/reference/pfinal_fmsm.md)
  : Probabilities of final states in a flexible parametric competing
  risks model
- [`simfinal_fmsm()`](https://chjackson.github.io/flexsurv/reference/simfinal_fmsm.md)
  : Simulate and summarise final outcomes from a flexible parametric
  multi-state model
- [`ajfit_fmsm()`](https://chjackson.github.io/flexsurv/reference/ajfit_fmsm.md)
  : Check the fit of Markov flexible parametric multi-state models
  against nonparametric estimates.
- [`msfit.flexsurvreg()`](https://chjackson.github.io/flexsurv/reference/msfit.flexsurvreg.md)
  : Cumulative intensity function for parametric multi-state models

## Competing risks and multi-state models based on mixtures

- [`flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/flexsurvmix.md)
  : Flexible parametric mixture models for times to competing events
- [`get_basepars()`](https://chjackson.github.io/flexsurv/reference/get_basepars.md)
  : Evaluate baseline time-to-event distribution parameters given
  covariate values in a flexsurvmix model
- [`mean_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/mean_flexsurvmix.md)
  : Mean times to events from a flexsurvmix model
- [`pdf_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/pdf_flexsurvmix.md)
  : Fitted densities for times to events in a flexsurvmix model
- [`probs_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/probs_flexsurvmix.md)
  : Probabilities of competing events from a flexsurvmix model
- [`p_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/p_flexsurvmix.md)
  : Transition probabilities from a flexsurvmix model
- [`quantile_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/quantile_flexsurvmix.md)
  : Quantiles of time-to-event distributions in a flexsurvmix model
- [`rmst_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/rmst_flexsurvmix.md)
  : Restricted mean times to events from a flexsurvmix model
- [`simt_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/simt_flexsurvmix.md)
  : Simulate times to competing events from a mixture multi-state model
- [`fmixmsm()`](https://chjackson.github.io/flexsurv/reference/fmixmsm.md)
  : Constructor for a mixture multi-state model based on flexsurvmix
- [`meanfinal_fmixmsm()`](https://chjackson.github.io/flexsurv/reference/meanfinal_fmixmsm.md)
  : Mean time to final state in a mixture multi-state model
- [`ppath_fmixmsm()`](https://chjackson.github.io/flexsurv/reference/ppath_fmixmsm.md)
  : Probability of each pathway taken through a mixture multi-state
  model
- [`qfinal_fmixmsm()`](https://chjackson.github.io/flexsurv/reference/qfinal_fmixmsm.md)
  : Quantiles of the distribution of the time until reaching a final
  state in a mixture multi-state model
- [`ajfit()`](https://chjackson.github.io/flexsurv/reference/ajfit.md) :
  Aalen-Johansen nonparametric estimates comparable to a fitted
  flexsurvmix model
- [`ajfit_flexsurvmix()`](https://chjackson.github.io/flexsurv/reference/ajfit_flexsurvmix.md)
  : Forms a tidy data frame for plotting the fit of parametric mixture
  multi-state models against nonparametric estimates

## Model comparison and checking

- [`logLik(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/logLik.flexsurvreg.md)
  : Log likelihood from a flexsurvreg model
- [`AICc(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/AICc.flexsurvreg.md)
  [`AICC(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/AICc.flexsurvreg.md)
  : Second-order Akaike information criterion
- [`AICc()`](https://chjackson.github.io/flexsurv/reference/AICc.md)
  [`AICC()`](https://chjackson.github.io/flexsurv/reference/AICc.md) :
  Second-order Akaike information criterion
- [`BIC(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/BIC.flexsurvreg.md)
  : Bayesian Information Criterion (BIC) for comparison of flexsurvreg
  models
- [`AIC(`*`<fmsm>`*`)`](https://chjackson.github.io/flexsurv/reference/AIC.fmsm.md)
  : Akaike's information criterion from a flexible parametric multistate
  model
- [`nobs(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/nobs.flexsurvreg.md)
  : Number of observations contributing to a fitted flexible survival
  model
- [`residuals(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/residuals.flexsurvreg.md)
  : Calculate residuals for flexible survival models
- [`coxsnell_flexsurvreg()`](https://chjackson.github.io/flexsurv/reference/coxsnell_flexsurvreg.md)
  : Cox-Snell residuals from a parametric survival model

## Probability distributions

- [`dsurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`psurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`qsurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`rsurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`Hsurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`hsurvspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`rmst_survspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  [`mean_survspline()`](https://chjackson.github.io/flexsurv/reference/Survspline.md)
  : Royston/Parmar spline survival distribution
- [`mean_survspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`mean_survspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rmst_survspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`dsurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`psurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`qsurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`rsurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`hsurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline0()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline1()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline2()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline3()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline4()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline5()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline6()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  [`Hsurvspline7()`](https://chjackson.github.io/flexsurv/reference/Survsplinek.md)
  : Royston/Parmar spline survival distribution functions with one
  argument per parameter
- [`dweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  [`pweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  [`qweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  [`hweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  [`HweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  [`rweibullPH()`](https://chjackson.github.io/flexsurv/reference/WeibullPH.md)
  : Weibull distribution in proportional hazards parameterisation
- [`dgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  [`pgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  [`Hgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  [`hgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  [`qgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  [`rgengamma()`](https://chjackson.github.io/flexsurv/reference/GenGamma.md)
  : Generalized gamma distribution
- [`dgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  [`pgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  [`Hgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  [`hgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  [`qgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  [`rgengamma.orig()`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md)
  : Generalized gamma distribution (original parameterisation)
- [`dgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md)
  [`pgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md)
  [`Hgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md)
  [`hgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md)
  [`qgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md)
  [`rgenf()`](https://chjackson.github.io/flexsurv/reference/GenF.md) :
  Generalized F distribution
- [`dgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  [`pgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  [`Hgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  [`hgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  [`qgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  [`rgenf.orig()`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
  : Generalized F distribution (original parameterisation)
- [`dgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  [`pgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  [`qgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  [`rgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  [`hgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  [`Hgompertz()`](https://chjackson.github.io/flexsurv/reference/Gompertz.md)
  : The Gompertz distribution
- [`dllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  [`pllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  [`qllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  [`rllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  [`hllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  [`Hllogis()`](https://chjackson.github.io/flexsurv/reference/Llogis.md)
  : The log-logistic distribution
- [`hexp()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`Hexp()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`hgamma()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`Hgamma()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`hlnorm()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`Hlnorm()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`hweibull()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  [`Hweibull()`](https://chjackson.github.io/flexsurv/reference/hazard.md)
  : Hazard and cumulative hazard functions
- [`mean_exp()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_exp()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_gamma()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_gamma()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_genf()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_genf()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_genf.orig()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_genf.orig()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_gengamma()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_gengamma()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_gengamma.orig()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_gengamma.orig()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_gompertz()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_gompertz()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_llogis()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_llogis()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_lnorm()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_lnorm()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_weibull()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_weibull()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`rmst_weibullPH()`](https://chjackson.github.io/flexsurv/reference/means.md)
  [`mean_weibullPH()`](https://chjackson.github.io/flexsurv/reference/means.md)
  : Mean and restricted mean survival functions
- [`custom_pexp()`](https://chjackson.github.io/flexsurv/reference/custom_pexp.md)
  : Piecewise exponential distribution for use with flexsurvreg()
- [`qgeneric()`](https://chjackson.github.io/flexsurv/reference/qgeneric.md)
  : Generic function to find quantiles of a distribution

## Tools for using flexsurv with tidyverse

- [`predict(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/predict.flexsurvreg.md)
  : Predictions from flexible survival models
- [`tidy(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/tidy.flexsurvreg.md)
  : Tidy a flexsurv model object
- [`tidy(`*`<standsurv>`*`)`](https://chjackson.github.io/flexsurv/reference/tidy.standsurv.md)
  : Tidy a standsurv object.
- [`augment(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/augment.flexsurvreg.md)
  : Augment data with information from a flexsurv model object
- [`glance(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/glance.flexsurvreg.md)
  : Glance at a flexsurv model object

## Datasets

- [`bc`](https://chjackson.github.io/flexsurv/reference/bc.md) : Breast
  cancer survival data
- [`bos`](https://chjackson.github.io/flexsurv/reference/bos.md)
  [`bosms3`](https://chjackson.github.io/flexsurv/reference/bos.md)
  [`bosms4`](https://chjackson.github.io/flexsurv/reference/bos.md) :
  Bronchiolitis obliterans syndrome after lung transplants

## Miscellaneous / advanced tools

- [`model.frame(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/model.frame.flexsurvreg.md)
  [`model.matrix(`*`<flexsurvreg>`*`)`](https://chjackson.github.io/flexsurv/reference/model.frame.flexsurvreg.md)
  :

  Extract original data from `flexsurvreg` objects.

- [`model.frame(`*`<flexsurvmix>`*`)`](https://chjackson.github.io/flexsurv/reference/model.frame.flexsurvmix.md)
  : Model frame from a flexsurvmix object

- [`basis()`](https://chjackson.github.io/flexsurv/reference/basis.md) :
  Natural cubic spline basis

- [`rmst_generic()`](https://chjackson.github.io/flexsurv/reference/rmst_generic.md)
  : Generic function to find restricted mean survival time for some
  distribution

- [`unroll.function()`](https://chjackson.github.io/flexsurv/reference/unroll.function.md)
  : Convert a function with matrix arguments to a function with vector
  arguments.

## Package overview

- [`flexsurv-package`](https://chjackson.github.io/flexsurv/reference/flexsurv-package.md)
  [`flexsurv`](https://chjackson.github.io/flexsurv/reference/flexsurv-package.md)
  : flexsurv: Flexible parametric survival and multi-state models
