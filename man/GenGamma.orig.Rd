\name{GenGamma.orig}
\alias{GenGamma.orig}
\alias{dgengamma.orig}
\alias{pgengamma.orig}
\alias{qgengamma.orig}
\alias{rgengamma.orig}
\alias{Hgengamma.orig}
\alias{hgengamma.orig}
\title{Generalized gamma distribution (original parameterisation)}
\description{
     Density, distribution function, hazards, quantile function and random
     generation for the generalized gamma distribution, using the
     original parameterisation from Stacy (1962).
}
\usage{
     dgengamma.orig(x, shape, scale=1, k, log = FALSE)
     pgengamma.orig(q, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE)
     qgengamma.orig(p, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE)
     rgengamma.orig(n, shape, scale=1, k)
     Hgengamma.orig(x, shape, scale=1, k)
     hgengamma.orig(x, shape, scale=1, k)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape}{vector of ``Weibull'' shape parameters.}
  \item{scale}{vector of scale parameters.}
  \item{k}{vector of ``Gamma'' shape parameters.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as
    log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}
}
\details{
  If \eqn{w \sim Gamma(k,1)}{w ~ Gamma(k, 1)}, then \eqn{x = \exp(w/shape + \log(scale))}{x = exp(w/shape + log(scale))} follows the
  original generalised gamma distribution with the parameterisation given here
  (Stacy 1962).  Defining \code{shape}\eqn{=b>0}, \code{scale}\eqn{=a>0},
  \eqn{x} has probability density

  \deqn{f(x | a, b, k) = \frac{b}{\Gamma(k)} \frac{x^{bk - 1}}{a^{bk}}
    \exp(-(x/a)^b)}{
    f(x | a, b, k) = (b / \Gamma(k)) (x^{bk -1} / a^{bk}) exp(-(x/a)^b)}

  The original generalized gamma distribution simplifies to the gamma, exponential
  and Weibull distributions with the following parameterisations: 

  \tabular{lcl}{
    \code{dgengamma.orig(x, shape, scale, k=1)} \tab \code{=} \tab \code{\link{dweibull}(x, shape, scale)} \cr
    \code{dgengamma.orig(x, shape=1, scale, k)} \tab \code{=} \tab \code{\link{dgamma}(x, shape=k, scale)} \cr
    \code{dgengamma.orig(x, shape=1, scale, k=1)} \tab \code{=} \tab \code{\link{dexp}(x, rate=1/scale)} \cr
  }

  Also as k tends to infinity, it tends to the log normal (as in
  \code{\link{dlnorm}}) with the following parameters (Lawless, 1980):

  \code{dlnorm(x, meanlog=log(scale) + log(k)/shape, sdlog=1/(shape*sqrt(k)))}
  
  For more stable behaviour as the distribution tends to the
  log-normal, an alternative parameterisation was developed by Prentice
  (1974).  This is given in \code{\link{dgengamma}}, and is now
  preferred for statistical modelling.  It is also more flexible,
  including a further new class of distributions with negative shape
  \code{k}. 
  
  The generalized F distribution \code{\link{GenF.orig}}, and its
  similar alternative parameterisation \code{\link{GenF}}, extend the
  generalized gamma to four parameters.
}
\value{
  \code{dgengamma.orig} gives the density, \code{pgengamma.orig} gives the distribution
  function, \code{qgengamma.orig} gives the quantile function, \code{rgengamma.orig}
  generates random deviates,  \code{Hgengamma.orig} retuns the cumulative hazard
  and \code{hgengamma.orig} the hazard. 
}
\references{
  Stacy, E. W. (1962). A generalization of the gamma distribution.
  Annals of Mathematical Statistics 33:1187-92. 

  Prentice, R. L. (1974). A log gamma model and its maximum likelihood
  estimation. Biometrika 61(3):539-544.

  Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
  distributions.  Technometrics 22(3):409-419. 
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\seealso{
\code{\link{GenGamma}}, \code{\link{GenF.orig}}, \code{\link{GenF}},
\code{\link{Lognormal}}, \code{\link{GammaDist}}, \code{\link{Weibull}}.
}
\keyword{ distribution }
