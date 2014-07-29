\name{GenF}
\alias{GenF}
\alias{dgenf}
\alias{pgenf}
\alias{qgenf}
\alias{rgenf}
\alias{Hgenf}
\alias{hgenf}
\title{Generalized F distribution}
\description{
     Density, distribution function, hazards, quantile function and random
     generation for the generalized F distribution, using the
     reparameterisation by Prentice (1975).
}
\usage{
     dgenf(x, mu=0, sigma=1, Q, P, log = FALSE)
     pgenf(q, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
     qgenf(p, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
     rgenf(n, mu=0, sigma=1, Q, P)
     Hgenf(x, mu=0, sigma=1, Q, P)
     hgenf(x, mu=0, sigma=1, Q, P)
}
\arguments{
  \item{x,q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mu}{Vector of location parameters.}
  \item{sigma}{Vector of scale parameters.}
  \item{Q}{Vector of first shape parameters.}
  \item{P}{Vector of second shape parameters.}
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
    shape parameters \eqn{s_1,s_2}{s1,s2}.

    In this more stable version described by Prentice (1975),
    \eqn{s_1,s_2}{s1,s2} are replaced by shape parameters \eqn{Q,P},
    with \eqn{P>0}, and

    \deqn{s_1 = 2(Q^2 + 2P + Q\delta)^{-1}, \quad s_2 = 2(Q^2 + 2P - Q\delta)^{-1}}{s1 = 2 / (Q^2 + 2P + Q*delta), s2 = 2 / (Q^2 + 2P - Q*delta)}
    equivalently
    \deqn{Q = \left(\frac{1}{s_1} - \frac{1}{s_2}\right)\left(\frac{1}{s_1} + \frac{1}{s_2}\right)^{-1/2}, \quad P = \frac{2}{s_1 + s_2} }{Q = (1/s1 - 1/s2) / (1/s1 + 1/s2)^{-1/2}, P = 2 / (s1 + s2) }

    Define \eqn{\delta = (Q^2 + 2P)^{1/2}}{delta = (Q^2 + 2P)^{1/2}},
  and \eqn{w = (\log(x) - \mu)\delta /\sigma}{w = (log(x) - mu)delta /
  sigma}, then the probability density function of \eqn{x} is \deqn{ f(x)  =  \frac{\delta (s_1/s_2)^{s_1} e^{s_1 w}}{\sigma x (1 + s_1 e^w/s_2) ^ {(s_1 + s_2)} B(s_1, s_2)}
  }{ f(x)  =  (delta (s1/s2)^{s1} e^{s1 w}) / (sigma t (1 + s1 e^w/s2) ^ {(s1 + s2)} B(s1, s2)) }
  The original parameterisation is available in this package as
  \code{\link{dgenf.orig}}, for the sake of completion / compatibility.
  With the above definitions,

  \code{dgenf(x, mu=mu, sigma=sigma, Q=Q, P=P) = dgenf.orig(x, mu=mu, sigma=sigma/delta, s1=s1, s2=s2)}

  The generalized F distribution with \code{P=0} is equivalent to the
  generalized gamma distribution \code{\link{dgengamma}}, so that
  \code{dgenf(x, mu, sigma, Q, P=0)} equals \code{dgengamma(x, mu, sigma, Q)}.
  The generalized gamma reduces further to several
  common distributions, as described in the \code{\link{GenGamma}}
  help page.

  The generalized F distribution includes the log-logistic distribution
  (see \code{\link[eha]{Llogis}}) as a further special case:

  \code{dgenf(x, mu=mu, sigma=sigma, Q=0, P=1) = \link[eha]{dllogis}(x, shape=sqrt(2)/sigma, scale=exp(mu))}

  The range of hazard trajectories available under this distribution are
  discussed in detail by Cox (2008).  Jackson et al. (2010) give an
  application to modelling oral cancer survival for use in a health
  economic evaluation of screening.
}
\value{
  \code{dgenf} gives the density, \code{pgenf} gives the distribution
  function, \code{qgenf} gives the quantile function, \code{rgenf}
  generates random deviates, \code{Hgenf} retuns the cumulative hazard
  and \code{hgenf} the hazard.
}
\references{
  R. L. Prentice (1975). Discrimination among some parametric
  models. Biometrika 62(3):607-614.

  Cox, C. (2008). The generalized \eqn{F} distribution: An umbrella for
  parametric survival analysis.  Statistics in Medicine 27:4301-4312.

  Jackson, C. H. and Sharples, L. D. and Thompson, S. G. (2010).
  Survival models in health economic evaluations: balancing fit and
  parsimony to improve prediction.   International Journal of
  Biostatistics 6(1):Article 34.
}
\author{
  Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\note{
  The parameters \code{Q} and \code{P} are usually called \eqn{q} and
  \eqn{p} in the literature - they were made upper-case in these R
  functions to avoid clashing with the conventional arguments \code{q}
  in the probability function and \code{p} in the quantile function.
}
\seealso{
  \code{\link{GenF.orig}}, \code{\link{GenGamma}}
}
\keyword{ distribution }
