\name{lines.flexsurvreg}
\title{Add fitted flexible survival curves to a plot}
\alias{lines.flexsurvreg}
\concept{Survival}
\description{
  Add fitted survival (or hazard or cumulative hazard) curves from a
  \code{\link{flexsurvreg}} model fit to an existing plot.

  Equivalent to \code{\link{plot.flexsurvreg}(...,add=TRUE)}.
}
\usage{
\method{lines}{flexsurvreg}(x, X=NULL, type="survival", t=NULL,
        est=TRUE, ci=NULL, B=1000, cl=0.95,
        col="red", lty=1, lwd=2, col.ci=NULL, lty.ci=2, lwd.ci=1, ...)
}
\arguments{

  \item{x}{Output from \code{\link{flexsurvreg}}, representing a fitted
    survival model object.}
  \item{X}{Matrix of covariate values to produce fitted curves
    for, as described in \code{\link{plot.flexsurvreg}}. }
  \item{type}{ \code{"survival"} for survival,
    \code{"cumhaz"} for cumulative hazard, or \code{"hazard"} for
    hazard, as in \code{\link{plot.flexsurvreg}}.
  }
  \item{t}{ Vector of times to plot fitted values for. }
  \item{est}{ Plot fitted curves (\code{TRUE} or \code{FALSE}.)}
  \item{ci}{ Plot confidence intervals for fitted curves.}
  \item{B}{Number of simulations controlling accuracy of confidence
    intervals, as used in \code{\link[=summary.flexsurvreg]{summary}}.}
  \item{cl}{ Width of confidence intervals, by default 0.95 for 95\% intervals.}

  \item{col}{ Colour of the fitted curve(s). }
  \item{lty}{ Line type of the fitted curve(s). }
  \item{lwd}{ Line width of the fitted curve(s). }

  \item{col.ci}{ Colour of the confidence limits, defaulting to the same
  as for the fitted curve.}
  \item{lty.ci}{ Line type of the confidence limits. }
  \item{lwd.ci}{ Line width of the confidence limits, defaulting to the same
  as for the fitted curve. }

  \item{...}{Other arguments to be passed to the generic
    \code{\link{plot}} and \code{\link{lines}} functions.}
}
\seealso{
  \code{\link{flexsurvreg}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
\keyword{aplot}
