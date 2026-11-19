# flexsurv: Flexible parametric survival and multi-state models

flexsurv: Flexible parametric models for time-to-event data, including
the generalized gamma, the generalized F and the Royston-Parmar spline
model, and extensible to user-defined distributions.

## Details

[`flexsurvreg`](https://chjackson.github.io/flexsurv/reference/flexsurvreg.md)
fits parametric models for time-to-event (survival) data. Data may be
right-censored, and/or left-censored, and/or left-truncated. Several
built-in parametric distributions are available. Any user-defined
parametric model can also be employed by supplying a list with basic
information about the distribution, including the density or hazard and
ideally also the cumulative distribution or hazard.

Covariates can be included using a linear model on any parameter of the
distribution, log-transformed to the real line if necessary. This
typically defines an accelerated failure time or proportional hazards
model, depending on the distribution and parameter.

[`flexsurvspline`](https://chjackson.github.io/flexsurv/reference/flexsurvspline.md)
fits the flexible survival model of Royston and Parmar (2002) in which
the log cumulative hazard is modelled as a natural cubic spline function
of log time. Covariates can be included on any of the spline parameters,
giving either a proportional hazards model or an arbitrarily-flexible
time-dependent effect. Alternative proportional odds or probit
parameterisations are available.

Output from the models can be presented as survivor, cumulative hazard
and hazard functions
([`summary.flexsurvreg`](https://chjackson.github.io/flexsurv/reference/summary.flexsurvreg.md)).
These can be plotted against nonparametric estimates
([`plot.flexsurvreg`](https://chjackson.github.io/flexsurv/reference/plot.flexsurvreg.md))
to assess goodness-of-fit. Any other user-defined function of the
parameters may be summarised in the same way.

Multi-state models for time-to-event data can also be fitted with the
same functions. Predictions from those models can then be made using the
functions
[`pmatrix.fs`](https://chjackson.github.io/flexsurv/reference/pmatrix.fs.md),
[`pmatrix.simfs`](https://chjackson.github.io/flexsurv/reference/pmatrix.simfs.md),
[`totlos.fs`](https://chjackson.github.io/flexsurv/reference/totlos.fs.md),
[`totlos.simfs`](https://chjackson.github.io/flexsurv/reference/totlos.simfs.md),
or
[`sim.fmsm`](https://chjackson.github.io/flexsurv/reference/sim.fmsm.md),
or alternatively by
[`msfit.flexsurvreg`](https://chjackson.github.io/flexsurv/reference/msfit.flexsurvreg.md)
followed by `mssample` or `probtrans` from the package mstate.

Distribution (“dpqr”) functions for the generalized gamma and F
distributions are given in
[`GenGamma`](https://chjackson.github.io/flexsurv/reference/GenGamma.md),
[`GenF`](https://chjackson.github.io/flexsurv/reference/GenF.md)
(preferred parameterisations) and
[`GenGamma.orig`](https://chjackson.github.io/flexsurv/reference/GenGamma.orig.md),
[`GenF.orig`](https://chjackson.github.io/flexsurv/reference/GenF.orig.md)
(original parameterisations). `flexsurv` also includes the standard
Gompertz distribution with unrestricted shape parameter, see
[`Gompertz`](https://chjackson.github.io/flexsurv/reference/Gompertz.md).

## User guide

The **flexsurv user guide** vignette explains the methods in detail, and
gives several worked examples. A further vignette **flexsurv-examples**
gives a few more complicated examples, and users are encouraged to
submit their own.

## References

Jackson, C. (2016). flexsurv: A Platform for Parametric Survival
Modeling in R. Journal of Statistical Software, 70(8), 1-33.
doi:10.18637/jss.v070.i08

Royston, P. and Parmar, M. (2002). Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Statistics in Medicine 21(1):2175-2197.

Cox, C. (2008). The generalized \\F\\ distribution: An umbrella for
parametric survival analysis. Statistics in Medicine 27:4301-4312.

Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007). Parametric
survival analysis and taxonomy of hazard functions for the generalized
gamma distribution. Statistics in Medicine 26:4252-4374

## See also

Useful links:

- <https://github.com/chjackson/flexsurv>

- <http://chjackson.github.io/flexsurv/>

- Report bugs at <https://github.com/chjackson/flexsurv/issues>

## Author

Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
