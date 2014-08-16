\name{msfit.flexsurvreg}
\title{Cumulative intensity function for parametric multi-state
  models}
\alias{msfit}
\alias{msfit.flexsurvreg}
\concept{Survival}
\description{
  Cumulative intensity function for fully-parametric multi-state or
  competing risks models, using a piecewise-constant approximation
  that will allow prediction using the functions in the \pkg{mstate}
  package. 
}
\usage{
msfit.flexsurvreg(object, t, newdata=NULL, variance=TRUE,
                 tvar="trans", trans, B=1000)
}
\arguments{
  \item{object}{
    Output from \code{\link{flexsurvreg}} or
    \code{\link{flexsurvspline}}, representing a fitted survival model
    object.

    The model should have been fitted to data consisting of one row for
    each observed transition and additional rows corresponding to
    censored times to competing transitions.  This is the "long" format,
    or counting process format, as explained in
    \code{\link[mstate]{msprep}}.   The \pkg{mstate} package
    documentation describes how to fit Cox models to this form of data. 

    The model should contain a categorical variable indicating the
    transition.  In \code{flexsurv} this variable can have any name,
    indicated here by the \code{tvar} argument.  In the Cox
    models demonstrated by \pkg{mstate} it is usually included in model
    formulae as \code{strata(trans)}, but note that the \code{strata}
    function does not do anything in \pkg{flexsurv}.  The formula
    supplied to \code{\link{flexsurvreg}} should be precise about which
    parameters are assumed to vary with the transition type.
  }
  
  \item{t}{ Vector of times.  These do not need to be the same as the
    observed event times, and since the model is parametric, they can be
    outside the range of the data.  A grid of more frequent times will
    provide a better approximation to the cumulative hazard trajectory
    for prediction with \code{\link[mstate]{probtrans}} or
    \code{\link[mstate]{mssample}}, at the cost of greater computational
    expense.}
  
  \item{newdata}{ A data frame specifying the values of covariates in
    the fitted model, other than the transition number.  This must be
    specified if there are other covariates. The variable names should
    be the same as those in the fitted model formula, but it should only
    have one row for each transition in the multi-state model.  }

  \item{variance}{ Calculate the variances and covariances of the
    transition cumulative hazards (\code{TRUE} or \code{FALSE}).  This
    is based on simulation from the normal asymptotic distribution
    of the estimates, which is computationally-expensive. }

  \item{tvar}{Name of the categorical variable in the model formula
    that represents the transition number. The values of this variable
    should correspond to elements of \code{trans}, conventionally a
    sequence of integers starting from 1.}

  \item{trans}{Matrix indicating allowed transitions in the multi-state
    model, in the format understood by \pkg{mstate}.}

  \item{B}{Number of simulations from the normal asymptotic distribution
    used to calculate variances.  Decrease for greater speed at the expense of accuracy.}
}
\value{
  An object of class \code{"msfit"}, in the same form as the objects
  used in the \pkg{mstate} package.  The \code{\link[mstate]{msfit}}
  method from \pkg{mstate} returns the equivalent cumulative intensities
  for Cox regression models fitted with \code{\link{coxph}}.

}
\examples{

## 3 state illness-death model for bronchiolitis obliterans
## Compare clock-reset / semi-Markov multi-state models

## Simple exponential model (reduces to Markov)
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans,
                    data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
mexp <- msfit.flexsurvreg(bexp, t=seq(0,12,by=0.1), trans=tmat, tvar="trans")

## Cox semi-parametric model within each transition
bcox <- coxph(Surv(years, status) ~ strata(trans), data=bosms3)
require(mstate)
mcox <- mstate::msfit(bcox, trans=tmat)

## Flexible parametric spline-based model 
bspl <- flexsurvspline(Surv(years, status) ~ trans + gamma1(trans),
                       data=bosms3, k=3)
mspl <- msfit.flexsurvreg(bspl, t=seq(0,12,by=0.1), trans=tmat, tvar="trans")

## Compare fit: exponential model is OK but the spline is better
plot(mcox, lwd=1, xlim=c(0, 12), ylim=c(0,4))
cols <- c("black","red","green")
for (i in 1:3){
    lines(mexp$Haz$time[mexp$Haz$trans==i], mexp$Haz$Haz[mexp$Haz$trans==i],
             col=cols[i], lwd=2, lty=2)
    lines(mspl$Haz$time[mspl$Haz$trans==i], mspl$Haz$Haz[mspl$Haz$trans==i],
             col=cols[i], lwd=3)
}
legend("topright", lwd=c(1,2,3), lty=c(1,2,1),
   c("Cox", "Exponential", "Flexible parametric"), bty="n")

}
\seealso{
  \code{\link{flexsurvreg}}, \code{\link{flexsurvspline}}.
}
\references{
   Liesbeth C. de Wreede, Marta Fiocco, Hein Putter (2011). \pkg{mstate}: An R
  Package for the Analysis of Competing Risks and Multi-State
  Models. \emph{Journal of Statistical Software}, 38(7), 1-30. 
  \url{http://www.jstatsoft.org/v38/i07}  

  Mandel, M. (2013). "Simulation based confidence intervals for
  functions with complicated derivatives." The American Statistician
  67(2):76-81
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
