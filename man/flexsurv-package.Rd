\encoding{latin1}
\name{flexsurv-package}
\alias{flexsurv-package}
\alias{flexsurv}
\docType{package}
\title{
flexsurv: Flexible parametric models for time-to-event data
}
\description{
  flexsurv: Flexible parametric models for time-to-event data,
  including the generalized gamma, the generalized F and the
  Royston-Parmar spline model, and extensible to user-defined
  distributions. 
}
\details{

  \code{\link{flexsurvreg}} fits parametric models for time-to-event
  (survival) data.  Data may be right-censored or left-truncated.
  Several built-in parametric distributions are available.  Any
  user-defined parametric model can also be employed by supplying a list
  with basic information about the distribution, including the density
  or hazard and ideally also the cumulative distribution or hazard.

  Covariates can be included using a linear model on any parameter of
  the distribution, log-transformed to the real line if necessary.  This
  typically defines an accelerated failure time or proportional hazards
  model, depending on the distribution and parameter.

  \code{\link{flexsurvspline}} fits the flexible survival model of
  Royston and Parmar (2002) in which the log cumulative hazard is
  modelled as a natural cubic spline function of log time.   Covariates
  can be included on any of the spline parameters, giving either a
  proportional hazards model or an arbitrarily-flexible time-dependent
  effect.  Alternative proportional odds or probit
  parameterisations are available.

  Output from the models can be presented as survivor, cumulative
  hazard and hazard functions (\code{\link{summary.flexsurvreg}}).
  These can be plotted against nonparametric estimates
  (\code{\link{plot.flexsurvreg}}) to assess goodness-of-fit.

  Distribution (``dpqr'') functions for the generalized gamma and F distributions are
  given in \code{\link{GenGamma}}, \code{\link{GenF}} (preferred
  parameterisations) and \code{\link{GenGamma.orig}},
  \code{\link{GenF.orig}} (original parameterisations).

  \code{\link{flexsurv}} also includes the standard Gompertz distribution
  with unrestricted shape parameter, see \code{\link{Gompertz}}.

}
\section{Related R packages}{
  \code{\link{flexsurv}} was written to encourage the use of flexible
  distributions to account for model uncertainty in survival analysis,
  initially the three-parameter generalized gamma, four-parameter
  generalized F and the Royston-Parmar spline
  model. However it was straightforward to modularise the design of the
  code to accept any generic parametric distribution.

  \code{\link[survival]{survreg}} from the \pkg{survival} package, the
  recommended R package for survival analysis, supports two-parameter
  location-scale parametric models.

  The \pkg{eha} package includes functions \code{\link[eha]{phreg}} and
  \code{\link[eha]{aftreg}} for parametric survival modelling under a
  variety of distributions and proportional hazards or accelerated
  failure time parameterisations.

  Other facilities for generic maximum likelihood model fitting exist,
  for example \code{\link[MASS]{fitdistr}} in the \pkg{MASS} package.
  \code{\link{flexsurvreg}} is intended to provide typical outputs and
  summaries of interest to survival analysts, particularly in medical
  applications.  Feature requests along these lines are welcome.

  Note that if an R package provides density and probability functions
  for a parametric distribution, it can then be used easily in
  \code{\link{flexsurvreg}}.  For instance, several ``reliability''
  distributions used in industrial statistics are available in the
  \pkg{VGAM} package.  Please report unexplained inconsistencies in
  results between \pkg{flexsurv} and other software.
}
\author{
  Christopher Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}
}
\references{
  Royston, P. and Parmar, M. (2002).  Flexible parametric
  proportional-hazards and proportional-odds models for censored survival
  data, with application to prognostic modelling and estimation of
  treatment effects. Statistics in Medicine 21(1):2175-2197.

  Cox, C. (2008). The generalized \eqn{F} distribution: An umbrella for
  parametric survival analysis.  Statistics in Medicine 27:4301-4312.

  Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007).
  Parametric survival analysis and taxonomy of hazard functions for the
  generalized gamma distribution.  Statistics in Medicine 26:4252-4374
}
\keyword{ package }
