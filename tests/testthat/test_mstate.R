context("Multi-state modelling and prediction")

require(mstate) # masks flexsurv's generic msfit, so call as msfit.flexsurvreg below

# library(testthat); library(devtools); 
# load_all("../..")
#unload("../..")
#library(flexsurv, lib.loc="~/work/flexsurv/lib/0.3")
#load(file="~/work/flexsurv/flexsurv/data/bosms3.rda")
#detach("package:flexsurv")

bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
bexp2 <- flexsurvreg(Surv(time, status) ~ trans, data=bosms3, dist="exp") 
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
mexp <- msfit.flexsurvreg(bexp, t=seq(0,150,1), trans=tmat, tvar="trans")

bcox <- coxph(Surv(time, status) ~ strata(trans), data=bosms3)
mcox <- mstate::msfit(bcox, trans=tmat)

bwei <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data=bosms3, dist="weibull", control=list(reltol=1e-16))
## NOTE NaN warnings: from zero shape or scale being visited by optimiser
## Is there a more sensible transform than log to optimise these on? 
mwei <- msfit.flexsurvreg(bwei, t=seq(0,150,1), trans=tmat, tvar="trans")

bgg <- flexsurvreg(Surv(time, status) ~ trans + sigma(trans) + Q(trans), data=bosms3, dist="gengamma")
mgg <- msfit.flexsurvreg(bgg, t=seq(0,150,1), trans=tmat, tvar="trans")

## how to constrain to exponential for 1-2 trans:  first shape par=1
bweifix <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + shape(trans), data=bosms3, dist="weibull", fixedpars=1)
mweifix <- msfit.flexsurvreg(bweifix, t=seq(0,150,1), trans=tmat, tvar="trans")

## Spline
bspl <- flexsurvspline(Surv(time, status) ~ trans + gamma1(trans), data=bosms3, k=3)
plot(bspl)
mspl <- msfit.flexsurvreg(bspl, t=seq(0,150,1), trans=tmat, tvar="trans")

sapply(list(bexp, bweifix, bwei, bgg, bspl),
       function(x)c("-2LL"=-2*logLik(x), p=attr(logLik(x),"df"), AIC=x$AIC))
## Weibull probably good enough

plot(mcox)
for (i in 1:3){
    lines(mexp$Haz$time[mexp$Haz$trans==i], mexp$Haz$Haz[mexp$Haz$trans==i], type="s")
    lines(mwei$Haz$time[mwei$Haz$trans==i], mwei$Haz$Haz[mwei$Haz$trans==i], type="s", col="red")
    lines(mweifix$Haz$time[mweifix$Haz$trans==i], mweifix$Haz$Haz[mweifix$Haz$trans==i], type="s", col="blue")
    lines(mgg$Haz$time[mgg$Haz$trans==i], mgg$Haz$Haz[mgg$Haz$trans==i], type="s", col="purple")
    lines(mspl$Haz$time[mspl$Haz$trans==i], mspl$Haz$Haz[mspl$Haz$trans==i], type="s", col="orange")
}

## Plot transition probabilities
## (though not valid for these semi-Markov models)
pt <- probtrans(mcox, predt=0, direction="forward")
plot(pt, from=1, xlim=c(0,150),ylim=c(0,1))
pt <- probtrans(mwei, predt=0, direction="forward")
par(new=TRUE)
plot(pt, xlim=c(0,150),ylim=c(0,1), col=rep("blue",3))
pt <- probtrans(mgg, predt=0, direction="forward")
par(new=TRUE)
plot(pt, xlim=c(0,150),ylim=c(0,1), col=rep("purple",3))


## With covariates

bosms3$x <- rnorm(nrow(bosms3))
bosms3$x2 <- factor(sample(0:2, size=nrow(bosms3), replace=TRUE))
bexp2 <- flexsurvreg(Surv(time, status) ~ trans + x + x2, data=bosms3, dist="exp")

test_that("Errors in msfit.flexsurvreg",{
    expect_error(mexp2 <- msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat), "\"tvar\" not supplied and variable .* not in model")
    expect_error(mexp2 <- msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat, tvar="foo"), "variable .* not in model")
    expect_error(mexp2 <- msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat, tvar="trans"), "Values of covariates .* not supplied")
    msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat, tvar="trans", newdata=data.frame(x=c(0,1,2), x2=c(0,1,2))) #works
})
