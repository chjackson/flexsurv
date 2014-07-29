\encoding{latin1}
\name{GenGamma}
\alias{GenGamma}
\alias{dgengamma}
\alias{pgengamma}
\alias{qgengamma}
\alias{rgengamma}
\alias{Hgengamma}
\alias{hgengamma}
\title{Generalized gamma distribution}
\description{
     Density, distribution function, hazards, quantile function and random
     generation for the generalized gamma distribution, using the
     parameterisation originating from Prentice (1974). Also known as the
     (generalized) log-gamma distribution. 
}
\usage{
     dgengamma(x, mu=0, sigma=1, Q, log = FALSE)
     pgengamma(q, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE)
     qgengamma(p, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE)
     rgengamma(n, mu=0, sigma=1, Q)
     Hgengamma(x, mu=0, sigma=1, Q)
     hgengamma(x, mu=0, sigma=1, Q)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mu}{Vector of ``location'' parameters.}
  \item{sigma}{Vector of ``scale'' parameters.   Note the inconsistent
    meanings of the term ``scale'' - this parameter is analogous to the
    (log-scale) standard deviation of the log-normal distribution,
    ``sdlog'' in \code{\link{dlnorm}}, rather than the ``scale'' parameter
    of the gamma distribution \code{\link{dgamma}}. Constrained to be
    positive. 
  }
  \item{Q}{Vector of shape parameters.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as
    log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}
}
\details{ 
  If \eqn{\gamma \sim Gamma(Q^{-2}, 1)}{g ~ Gamma(Q^{-2}, 1)} , and 
  \eqn{w  = log(Q^2 \gamma) / Q}{w  = log(Q^2*g) / Q}, then \eqn{x = \exp(\mu + \sigma w)}{x = exp(mu + sigma w)}
  follows the generalized gamma distribution with probability density function

\deqn{f(x | \mu, \sigma, Q) = \frac{|Q|(Q^{-2})^{Q^{-2}}}{\sigma x \Gamma(Q^{-2})} \exp(Q^{-2}(Qw - \exp(Qw)))}{
  f(x | mu, sigma, Q) = |Q| (Q^{-2})^{Q^{-2}} / / (sigma * x * Gamma(Q^{-2})) exp(Q^{-2}*(Q*w - exp(Q*w)))}
  
  This parameterisation is preferred to the original parameterisation of
  the generalized gamma by Stacy (1962) since it is more numerically
  stable near to \eqn{Q=0} (the log-normal distribution), and allows \eqn{Q<=0}.
  The original is available in this package as
  \code{\link{dgengamma.orig}}, for the sake of completion and
  compatibility with other software - this is implicitly restricted to
  \code{Q}>0 (or \code{k}>0 in the original notation). The parameters of
  \code{\link{dgengamma}} and \code{\link{dgengamma.orig}} are related as
  follows.

  \code{dgengamma.orig(x, shape=shape, scale=scale, k=k) = }

  \code{dgengamma(x, mu=log(scale) + log(k)/shape, sigma=1/(shape*sqrt(k)), Q=1/sqrt(k))}

  The generalized gamma distribution simplifies to the gamma, log-normal
  and Weibull distributions with the following parameterisations: 

  \tabular{lcl}{
    \code{dgengamma(x, mu, sigma, Q=0)} \tab \code{=} \tab \code{dlnorm(x, mu, sigma)} \cr
    \code{dgengamma(x, mu, sigma, Q=1)} \tab \code{=} \tab  \code{dweibull(x, shape=1/sigma, scale=exp(mu))} \cr
    \code{dgengamma(x, mu, sigma, Q=sigma)} \tab \code{=} \tab \code{dgamma(x, shape=1/sigma^2, rate=exp(-mu) / sigma^2)} \cr
  }
  The properties of the generalized gamma and its applications to
  survival analysis are discussed in detail by Cox (2007). 
  
  The generalized F distribution \code{\link{GenF}} extends the
  generalized gamma to four parameters. 
}
\value{
  \code{dgengamma} gives the density, \code{pgengamma} gives the distribution
  function, \code{qgengamma} gives the quantile function, \code{rgengamma}
  generates random deviates,   \code{Hgengamma} retuns the cumulative hazard
  and \code{hgengamma} the hazard. 
}
\references{
  Prentice, R. L. (1974). A log gamma model and its maximum likelihood
  estimation. Biometrika 61(3):539-544.

  Farewell, V. T. and Prentice, R. L. (1977). A study of distributional
  shape in life testing. Technometrics 19(1):69-75. 
  
  Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
  distributions.  Technometrics 22(3):409-419. 

  Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007). 
  Parametric survival analysis and taxonomy of hazard functions for the
  generalized gamma distribution.  Statistics in Medicine 26:4252-4374

  Stacy, E. W. (1962). A generalization of the gamma distribution.
  Annals of Mathematical Statistics 33:1187-92
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\seealso{
\code{\link{GenGamma.orig}}, \code{\link{GenF}},
\code{\link{Lognormal}}, \code{\link{GammaDist}}, \code{\link{Weibull}}.
}
\keyword{ distribution }
