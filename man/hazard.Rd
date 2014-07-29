\name{Hazard and cumulative hazard functions for standard distributions}
\alias{hexp}
\alias{Hexp}
\alias{hweibull}
\alias{Hweibull}
\alias{hgamma}
\alias{Hgamma}
\alias{hlnorm}
\alias{Hlnorm}
\title{Hazard and cumulative hazard functions}
\description{Hazard and cumulative hazard functions for distributions
  which are built into flexsurv, and whose distribution functions are in base R.
}
\usage{
hexp(x, rate=1, log=FALSE)
Hexp(x, rate=1, log=FALSE)
hweibull (x, shape, scale = 1, log = FALSE) 
Hweibull (x, shape, scale = 1, log = FALSE) 
hgamma(x, shape, rate=1, log=FALSE)
Hgamma(x, shape, rate=1, log=FALSE)
hlnorm(x, meanlog=0, sdlog=1, log=FALSE)
Hlnorm(x, meanlog=0, sdlog=1, log=FALSE)
}
\arguments{
  \item{x}{Vector of quantiles}
  \item{rate}{Rate parameter (exponential and gamma)}
  \item{shape}{Shape parameter (Weibull and gamma)}
  \item{scale}{Scale parameter (Weibull)}
  \item{meanlog}{Mean on the log scale (log normal)}
  \item{sdlog}{Standard deviation on the log scale (log normal)}
  \item{log}{Compute log hazard or log cumulative hazard}
}
\value{
  Hazard (functions beginning 'h') or cumulative hazard (functions beginning 'H').
}
\details{
  For the exponential and the Weibull these are available analytically,
  and so are programmed here in numerically stable and efficient forms.

  For the gamma and log-normal, these are simply computed as minus the
  log of the survivor function (cumulative hazard) or the ratio of the
  density and survivor function (hazard), so are not expected to be
  robust to extreme values or quick to compute.  
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
  }
\seealso{
  \code{\link{dexp}},\code{\link{dweibull}},\code{\link{dgamma}},\code{\link{dlnorm}},\code{\link{dgompertz}},\code{\link{dgengamma}},\code{\link{dgenf}}
}
\keyword{distribution}
