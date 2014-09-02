\name{sim.fmsm}
\alias{sim.fmsm}
\title{Simulate paths through a fully parametric multi-state model}
\description{
  Simulate changes of state and transition times from a fully parametric
  multi-state model fitted using \code{\link{flexsurvreg}}.
}
\usage{
sim.fmsm(x, trans, t, newdata=NULL, start=1, M=10, tvar="trans", debug=FALSE)
}
\arguments{
  \item{x}{A model fitted with \code{\link{flexsurvreg}}. See
\code{\link{msfit.flexsurvreg}} for the required form of the model and the data.}
  \item{trans}{Matrix indicating allowed transitions.  See \code{\link{msfit.flexsurvreg}}.}
  \item{t}{Time, or vector of times for each of the \code{M} individuals, to simulate trajectories until.}
  \item{newdata}{A data frame specifying the values of covariates in the
    fitted model, other than the transition number.  See \code{\link{msfit.flexsurvreg}}.}
  \item{start}{Starting state, or vector of starting states for each of the \code{M} individuals.}
  \item{M}{Number of individual trajectories to simulate.}
  \item{tvar}{Variable in the data representing the transition type.}
  \item{debug}{(for development use)}
}
\value{
A list of two matrices named \code{st} and \code{t}.  The rows of each
matrix represent simulated individuals.  The columns of \code{t} contain the
times when the individual changes state, to the corresponding states in
\code{st}.

The first columns will always contain the starting states and the
starting times. The last column of \code{t} represents either the
time when the individual moves to an absorbing state, or right-censoring
in a transient state at the time given in the \code{t} argument to \code{sim.fmsm}.
}
\details{
\code{sim.fmsm} relies on the presence of a function to sample random
numbers from the parametric survival distribution used in the fitted
model \code{x}.If \code{x} was fitted using a custom distribution,
called \code{dist} say, then there must be a function called (something like)
\code{rdist} either in the working environment, or supplied through the
\code{dfns} argument to \code{flexsurvreg}.  This must be in the same
format as standard R functions such as \code{rnorm}, with first argument
\code{n} and remaining arguments giving the parameters of the distribution.
}
\seealso{
\code{\link{pmatrix.simfs}}
}
\examples{
bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
sim.fmsm(bexp, M=10, t=5, trans=tmat)
}
\author{Christopher Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}.}
\keyword{models,survival}
