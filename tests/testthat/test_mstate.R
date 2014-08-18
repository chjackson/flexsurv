context("Multi-state modelling and prediction")

require(mstate) # masks flexsurv's generic msfit, so call as msfit.flexsurvreg below


### TODO change times in msfit


# library(testthat); library(devtools); 
# load_all("../..")
#unload("../..")
#library(flexsurv, lib.loc="~/work/flexsurv/lib/0.3")
#load(file="~/work/flexsurv/flexsurv/data/bosms3.rda")
#detach("package:flexsurv")

bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp") 
bexp2 <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

tgrid <- seq(0,14,by=0.1)

mexp <- msfit.flexsurvreg(bexp, t=tgrid, trans=tmat, tvar="trans")

bcox <- coxph(Surv(years, status) ~ strata(trans), data=bosms3)
mcox <- mstate::msfit(bcox, trans=tmat)

bwei <- flexsurvreg(Surv(years, status) ~ trans + shape(trans), data=bosms3, dist="weibull")
mwei <- msfit.flexsurvreg(bwei, t=tgrid, trans=tmat, tvar="trans")




bgg <- flexsurvreg(Surv(years, status) ~ trans + sigma(trans) + Q(trans), data=bosms3, dist="gengamma")
mgg <- msfit.flexsurvreg(bgg, t=tgrid, trans=tmat, tvar="trans")

## how to constrain to exponential for 1-2 trans:  first shape par=1
bweifix <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + shape(trans), data=bosms3, dist="weibull", fixedpars=1)
mweifix <- msfit.flexsurvreg(bweifix, t=tgrid, trans=tmat, tvar="trans")

## Spline
bspl <- flexsurvspline(Surv(years, status) ~ trans + gamma1(trans), data=bosms3, k=3)
plot(bspl)
mspl <- msfit.flexsurvreg(bspl, t=tgrid, trans=tmat, tvar="trans")

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
plot(pt, from=1, xlim=range(tgrid),ylim=c(0,1))
# pt <- probtrans(mwei, predt=0, direction="forward")
par(new=TRUE)
plot(pt, xlim=range(tgrid),ylim=c(0,1), col=rep("blue",3))
# pt <- probtrans(mgg, predt=0, direction="forward")
par(new=TRUE)
plot(pt, xlim=range(tgrid),ylim=c(0,1), col=rep("purple",3))

## 

## With covariates
bosms3$x <- rnorm(nrow(bosms3))
bosms3$x2 <- factor(sample(0:2, size=nrow(bosms3), replace=TRUE))
bexp2 <- flexsurvreg(Surv(years, status) ~ trans + x + x2, data=bosms3, dist="exp")

test_that("Errors in msfit.flexsurvreg",{
    expect_error(mexp2 <- msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat), "Values of covariates .+ not supplied")
    expect_error(mexp2 <- msfit.flexsurvreg(bexp2, t=seq(0,150,10), trans=tmat, tvar="foo"), "variable .* not in model")
})



tgrid <- seq(0,14,by=0.1)
mfwei <- msfit.flexsurvreg(bwei, t=tgrid, trans=tmat, tvar="trans")
ptw <- probtrans(mfwei, predt=0, direction="forward")[[1]]
ptw[ptw$time %in% c(5,10),]

tgrid <- seq(0,14,by=0.01)
mfwei <- msfit.flexsurvreg(bwei, t=tgrid, trans=tmat, tvar="trans")
ptw <- probtrans(mfwei, predt=0, direction="forward")[[1]]
ptw[ptw$time %in% c(5,10),]

tgrid <- seq(0,14,by=0.001)
mfwei <- msfit.flexsurvreg(bwei, t=tgrid, trans=tmat, tvar="trans")
ptw <- probtrans(mfwei, predt=0, direction="forward")[[1]]
ptw[ptw$time %in% c(5,10),]
