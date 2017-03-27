##' flexsurv: Flexible parametric survival and multi-state models
##' 
##' flexsurv: Flexible parametric models for time-to-event data, including the
##' generalized gamma, the generalized F and the Royston-Parmar spline model,
##' and extensible to user-defined distributions.
##' 
##' \code{\link{flexsurvreg}} fits parametric models for time-to-event
##' (survival) data.  Data may be right-censored, and/or left-censored, and/or
##' left-truncated.  Several built-in parametric distributions are available.
##' Any user-defined parametric model can also be employed by supplying a list
##' with basic information about the distribution, including the density or
##' hazard and ideally also the cumulative distribution or hazard.
##' 
##' Covariates can be included using a linear model on any parameter of the
##' distribution, log-transformed to the real line if necessary.  This
##' typically defines an accelerated failure time or proportional hazards
##' model, depending on the distribution and parameter.
##' 
##' \code{\link{flexsurvspline}} fits the flexible survival model of Royston
##' and Parmar (2002) in which the log cumulative hazard is modelled as a
##' natural cubic spline function of log time.  Covariates can be included on
##' any of the spline parameters, giving either a proportional hazards model or
##' an arbitrarily-flexible time-dependent effect.  Alternative proportional
##' odds or probit parameterisations are available.
##' 
##' Output from the models can be presented as survivor, cumulative hazard and
##' hazard functions (\code{\link{summary.flexsurvreg}}).  These can be plotted
##' against nonparametric estimates (\code{\link{plot.flexsurvreg}}) to assess
##' goodness-of-fit.  Any other user-defined function of the parameters may be
##' summarised in the same way.
##' 
##' Multi-state models for time-to-event data can also be fitted with the same
##' functions.  Predictions from those models can then be made using the
##' functions \code{\link{pmatrix.fs}}, \code{\link{pmatrix.simfs}},
##' \code{\link{totlos.fs}}, \code{\link{totlos.simfs}}, or
##' \code{\link{sim.fmsm}}, or alternatively by \code{\link{msfit.flexsurvreg}}
##' followed by \code{mssample} or \code{probtrans} from the package
##' \pkg{mstate}.
##' 
##' Distribution (``dpqr'') functions for the generalized gamma and F
##' distributions are given in \code{\link{GenGamma}}, \code{\link{GenF}}
##' (preferred parameterisations) and \code{\link{GenGamma.orig}},
##' \code{\link{GenF.orig}} (original parameterisations).
##' \code{\link{flexsurv}} also includes the standard Gompertz distribution
##' with unrestricted shape parameter, see \code{\link{Gompertz}}.
##' 
##' @name flexsurv-package
##' @aliases flexsurv-package flexsurv
##' @docType package
##' @section User guide: The \bold{flexsurv user guide} vignette explains the
##' methods in detail, and gives several worked examples.  A further vignette
##' \bold{flexsurv-examples} gives a few more complicated examples, and users
##' are encouraged to submit their own.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @references Jackson, C. (2016). flexsurv: A Platform for Parametric
##' Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33.
##' doi:10.18637/jss.v070.i08
##' 
##' Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' 
##' Cox, C. (2008). The generalized \eqn{F} distribution: An umbrella for
##' parametric survival analysis.  Statistics in Medicine 27:4301-4312.
##' 
##' Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007).  Parametric
##' survival analysis and taxonomy of hazard functions for the generalized
##' gamma distribution.  Statistics in Medicine 26:4252-4374
##' @keywords package
##' @importFrom Rcpp sourceCpp
##' @useDynLib flexsurv, .registration = TRUE
##' @import stats
##' @importFrom mstate msfit probtrans
##' @importFrom graphics plot lines
##' @importFrom survival Surv survfit coxph survreg survreg.control
##' @importFrom muhaz muhaz
##' @importFrom mvtnorm rmvnorm
##' @importFrom deSolve ode
##' @importFrom quadprog solve.QP
"_PACKAGE"

.onUnload <- function(libpath) {
    library.dynam.unload("flexsurv", libpath)
}

##' Hazard and cumulative hazard functions
##' 
##' Hazard and cumulative hazard functions for distributions which are built
##' into flexsurv, and whose distribution functions are in base R.
##' 
##' For the exponential and the Weibull these are available analytically, and
##' so are programmed here in numerically stable and efficient forms.
##' 
##' For the gamma and log-normal, these are simply computed as minus the log of
##' the survivor function (cumulative hazard) or the ratio of the density and
##' survivor function (hazard), so are not expected to be robust to extreme
##' values or quick to compute.
##' 
##' @aliases hexp Hexp hweibull Hweibull hgamma Hgamma hlnorm Hlnorm
##' @param x Vector of quantiles
##' @param rate Rate parameter (exponential and gamma)
##' @param shape Shape parameter (Weibull and gamma)
##' @param scale Scale parameter (Weibull)
##' @param meanlog Mean on the log scale (log normal)
##' @param sdlog Standard deviation on the log scale (log normal)
##' @param log Compute log hazard or log cumulative hazard
##' @return Hazard (functions beginning 'h') or cumulative hazard (functions
##' beginning 'H').
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso
##' \code{\link{dexp}},\code{\link{dweibull}},\code{\link{dgamma}},\code{\link{dlnorm}},\code{\link{dgompertz}},\code{\link{dgengamma}},\code{\link{dgenf}}
##' @keywords distribution
##' @name hazard
NULL

##' Mean and restricted mean survival functions
##' 
##' Mean and restricted mean survival time functions for distributions which are built
##' into flexsurv.
##' 
##' For the exponential, Weibull, log-logistic, lognormal, and gamma, mean survival is
##' provided analytically.  Restricted mean survival for the exponential distribution
##' is also provided analytically.  Mean and restricted means for other distributions
##' are calculated via numeric integration.
##' 
##' @aliases mean_exp rmst_exp mean_weibull rmst_weibull mean_weibullPH rmst_weibullPH
##' mean_llogis rmst_llogis mean_lnorm rmst_lnorm mean_gamma rmst_gamma mean_gompertz
##' rmst_gompertz mean_gengamma rmst_gengamma mean_gengamma.orig rmst_gengamma.orig
##' mean_genf rmst_genf mean_genf.orig rmst_genf.orig
##' @param t Vector of times to which restricted mean survival time is evaluated
##' @param start Optional left-truncation time or times.  The returned
##' restricted mean survival will be conditioned on survival up to
##' this time.
##' @param rate Rate parameter (exponential and gamma)
##' @param shape Shape parameter (Weibull, gamma, log-logistic, generalized gamma [orig],
##' generalized F [orig])
##' @param scale Scale parameter (Weibull, log-logistic, generalized gamma [orig],
##' generalized F [orig])
##' @param meanlog Mean on the log scale (log normal)
##' @param sdlog Standard deviation on the log scale (log normal)
##' @param mu Mean on the log scale (generalized gamma, generalized F)
##' @param sigma Standard deviation on the log scale (generalized gamma, generalized F)
##' @param Q Vector of first shape parameters (generalized gamma, generalized F)
##' @param P Vector of second shape parameters (generalized F)
##' @param k vector of shape parameters (generalized gamma [orig]).
##' @param s1 Vector of first F shape parameters (generalized F [orig])
##' @param s2 vector of second F shape parameters (generalized F [orig])
##' @return mean survival (functions beginning 'mean') or restricted mean survival
##' (functions beginning 'rmst_').
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso
##' \code{\link{dexp}},\code{\link{dweibull}},\code{\link{dgamma}},\code{\link{dlnorm}},\code{\link{dgompertz}},\code{\link{dgengamma}},\code{\link{dgenf}}
##' @keywords distribution
##' @name means
NULL

##' Breast cancer survival data
##' 
##' Survival times of 686 patients with primary node positive breast cancer.
##' 
##' 
##' @format A data frame with 686 rows.  \tabular{rll}{ \code{censrec} \tab
##' (numeric) \tab 1=dead, 0=censored \cr \code{rectime} \tab (numeric) \tab
##' Time of death or censoring in days\cr \code{group} \tab (numeric) \tab
##' Prognostic group: \code{"Good"},\code{"Medium"} or \code{"Poor"}, \cr \tab
##' \tab from a regression model developed by Sauerbrei and Royston (1999).\cr
##' }
##' @seealso \code{\link{flexsurvspline}}
##' @references Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' 
##' Sauerbrei, W. and Royston, P. (1999). Building multivariable prognostic and
##' diagnostic models: transformation of the predictors using fractional
##' polynomials.  Journal of the Royal Statistical Society, Series A 162:71-94.
##' @source German Breast Cancer Study Group, 1984-1989.  Used as a reference
##' dataset for the spline-based survival model of Royston and Parmar (2002),
##' implemented here in \code{\link{flexsurvspline}}.  Originally provided with
##' the \code{stpm} (Royston 2001, 2004) and \code{stpm2} (Lambert 2009, 2010)
##' Stata modules.
##' @keywords datasets
"bc"


##' Bronchiolitis obliterans syndrome after lung transplants
##' 
##' A dataset containing histories of bronchiolitis obliterans syndrome (BOS)
##' from lung transplant recipients. BOS is a chronic decline in lung function,
##' often observed after lung transplantation.
##' 
##' The entry time of each patient into each stage of BOS was estimated by
##' clinicians, based on their history of lung function measurements and acute
##' rejection and infection episodes.  BOS is only assumed to occur beyond six
##' months after transplant.  In the first six months the function of each
##' patient's new lung stabilises.  Subsequently BOS is diagnosed by comparing
##' the lung function against the "baseline" value.
##' 
##' The same data are provided in the \pkg{msm} package, but in the
##' native format of \pkg{msm} to allow Markov models to be fitted.
##' In \pkg{flexsurv}, much more flexible models can be fitted.
##' @name bos
##' @aliases bosms3 bosms4
##' @docType data
##' @format A data frame containing a sequence of observed or censored
##' transitions to the next stage of severity or death.  It is grouped
##' by patient and includes histories of 204 patients.  All patients
##' start in state 1 (no BOS) at six months after transplant, and may
##' subsequently develop BOS or die.
##' 
##' \code{bosms3} contains the data for a three-state model: no BOS, BOS or
##' death. \code{bosms4} uses a four-state representation: no BOS, mild BOS,
##' moderate/severe BOS or death.  \tabular{rll}{ \code{id} \tab (numeric) \tab
##' Patient identification number \cr \code{from} \tab (numeric) \tab Observed
##' starting state of the transition \cr \code{to} \tab (numeric) \tab Observed
##' or potential ending state of the transition \cr \code{Tstart} \tab
##' (numeric) \tab Time at the start of the interval \cr \code{Tstop} \tab
##' (numeric) \tab Time at the end of the interval \cr \code{time} \tab
##' (numeric) \tab Time difference \code{Tstart}-\code{Tstop} \cr \code{status}
##' \tab (numeric) \tab 1 if the transition to state \code{to} was observed, or
##' 0 if the transition to state \code{to} was censored (for example, if the
##' patient was observed to move to a competing state) \cr \code{trans} \tab
##' (factor) \tab Number of the transition \code{from}-\code{to} in the set of
##' all \code{ntrans} allowed transitions, numbered from 1 to \code{ntrans}.  }
##' @references Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome:
##' Incidence, Natural History, Prognosis, and Risk Factors.  Journal of Heart
##' and Lung Transplantation 17(12)1255--1263.
##' @source Papworth Hospital, U.K.
##' @keywords datasets
NULL





