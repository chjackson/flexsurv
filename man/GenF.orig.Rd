\name{GenF.orig}
\alias{GenF.orig}
\alias{dgenf.orig}
\alias{pgenf.orig}
\alias{qgenf.orig}
\alias{rgenf.orig}
\alias{Hgenf.orig}
\alias{hgenf.orig}
\title{Generalized F distribution (original parameterisation)}
\description{
     Density, distribution function, quantile function and random
     generation for the generalized F distribution, using the
     less flexible original parameterisation described
     by Prentice (1975). 
}
\usage{
     dgenf.orig(x, mu=0, sigma=1, s1, s2, log = FALSE)
     pgenf.orig(q, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
     qgenf.orig(p, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
     rgenf.orig(n, mu=0, sigma=1, s1, s2)
     Hgenf.orig(x, mu=0, sigma=1, s1, s2)
     hgenf.orig(x, mu=0, sigma=1, s1, s2)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mu}{Vector of location parameters.}
  \item{sigma}{Vector of scale parameters.}
  \item{s1}{Vector of first F shape parameters.}
  \item{s2}{vector of second F shape parameters.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as
    log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}
}
\details{
  If \eqn{y \sim F(2s_1, 2s_2)}{y ~ F(2*s1, 2*s2)}, and \eqn{w =
    \log(y)}{w = log(y)} then \eqn{x = \exp(w\sigma + \mu)}{x =
    exp(w*sigma + mu)} has the original generalized F distribution with location
    parameter \eqn{\mu}{mu}, scale parameter \eqn{\sigma>0}{sigma>0} and
    shape parameters \eqn{s_1>0,s_2>0}{s1>0,s2>0}.  The probability density
    function of \eqn{x} is

    \deqn{f(x | \mu, \sigma, s_1, s_2) =
\frac{(s_1/s_2)^{s_1} e^{s_1 w}}{\sigma x (1 + s_1 e^w/s_2) ^ {(s_1 +
    s_2)} B(s_1, s_2)}
}{f(x | mu, sigma, s_1, s_2) = ((s1/s2)^{s1} e^{s1 w}) / (sigma x (1 + s1 e^w/s2) ^ (s1 + s2) B(s1, s2)) }
    where \eqn{w = (\log(x) - \mu)/\sigma}{w = (log(x) - mu)/sigma} ,
    \eqn{B(s_1,s_2) = \Gamma(s_1)\Gamma(s_2)/\Gamma(s_1+s_2)}{B(s1,s2) = \Gamma(s1)\Gamma(s2)/\Gamma(s1+s2)} is the beta function.
    
  As \eqn{s_2 \rightarrow \infty}{s2 -> infinity}, the distribution of
  \eqn{x} tends towards an original generalized gamma distribution with
  the following parameters:
  
  \code{\link{dgengamma.orig}(x, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1)}

  See \code{\link{GenGamma.orig}} for how this includes several other
  common distributions as special cases.
  
  The alternative parameterisation of the generalized F distribution,
  originating from Prentice (1975) and given in this package as
  \code{\link{GenF}},  is
  preferred for statistical modelling, since it is more stable as
  \eqn{s_1}{s1} tends to infinity, and includes a further new class
  of distributions with negative first shape parameter.   The original
  is provided here for the sake of completion and compatibility. 
}
\value{
  \code{dgenf.orig} gives the density, \code{pgenf.orig} gives the distribution
  function, \code{qgenf.orig} gives the quantile function, \code{rgenf.orig}
  generates random deviates, \code{Hgenf.orig} retuns the cumulative hazard
  and \code{hgenf.orig} the hazard. 
}
\references{
  R. L. Prentice (1975). Discrimination among some parametric
  models. Biometrika 62(3):607-614.
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\seealso{
  \code{\link{GenF}}, \code{\link{GenGamma.orig}}, \code{\link{GenGamma}}
}
\keyword{ distribution }
