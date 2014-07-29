\name{bc}
\alias{bc}
\title{Breast cancer survival data}
\description{
  Survival times of 686 patients with primary node positive breast
  cancer. 
}
\usage{
data(bc)
}
\format{
  A data frame with 686 rows. 
    \tabular{rll}{
    \code{censrec} \tab (numeric) \tab 1=dead, 0=censored \cr
    \code{rectime} \tab (numeric) \tab Time of death or censoring in days\cr
    \code{group} \tab (numeric) \tab Prognostic group:
    \code{"Good"},\code{"Medium"} or \code{"Poor"}, \cr
    \tab \tab from a regression model developed by
              Sauerbrei and Royston (1999).\cr
  }
}
\source{
  German Breast Cancer Study Group, 1984-1989.  Used as a reference
  dataset for the spline-based survival model of Royston and Parmar
  (2002), implemented here in \code{\link{flexsurvspline}}.  Originally
  provided with the \code{stpm} (Royston 2001, 2004) and \code{stpm2}
  (Lambert 2009, 2010) Stata modules.
}
\references{
  Royston, P. and Parmar, M. (2002).  Flexible parametric
  proportional-hazards and proportional-odds models for censored survival
  data, with application to prognostic modelling and estimation of
  treatment effects. Statistics in Medicine 21(1):2175-2197. 

  Sauerbrei, W. and Royston, P. (1999). Building multivariable prognostic
  and diagnostic models: transformation of the predictors using
  fractional polynomials.  Journal of the Royal Statistical Society,
  Series A 162:71-94.   
}
\seealso{
  \code{\link{flexsurvspline}}
}
\keyword{datasets}
