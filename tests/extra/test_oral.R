load("~/work/oral/tcrb.rda")
## Without covs
fit <- flexsurvreg(Surv(survtime, dead2) ~ 1, data=tcrb, dist="genf")
fit.gg <- flexsurvreg(Surv(survtime, dead2) ~ 1, data=tcrb, dist="gengamma")
fit.sp2 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=2, data=tcrb, control=list(maxit=10000))
fit.sp3 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=3, data=tcrb, control=list(maxit=10000))
fit.sp4 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=4, data=tcrb, control=list(maxit=10000))
fit.sp5 <- flexsurvspline(Surv(survtime, dead2) ~ 1, k=5, data=tcrb, method="BFGS", control=list(trace=1,REPORT=1,maxit=10000))
## min AIC with 4 knots, better than GF
save(fit, fit.gg, fit.sp2, fit.sp3, fit.sp4, fit.sp5, file="~/work/flexsurv/devel/tcr.rda")
plot(fit, ci=FALSE)
plot(fit, ci=FALSE, type="hazard")
plot(fit, ci=FALSE)
lines(fit.gg, col="blue", ci=FALSE)
lines(fit.sp4, col="green", ci=FALSE)

## With covs -- can't get arbitrarily good fit unless model 3 way interaction
fitc.f <- flexsurvreg(Surv(survtime, dead2) ~ age_10 + sex + stage, data=tcrb, dist="genf")
fitc.g <- flexsurvreg(Surv(survtime, dead2) ~ age_10 + sex + stage, data=tcrb, dist="gengamma")
fitc.sp2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, data=tcrb) # doesnt fit as well as gf/gg.  PH assumption vs AFT?
fitc.sp3 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=3, data=tcrb)
fitc.sp4 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=4, data=tcrb)
fitc.spo2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, scale="odds", data=tcrb) # PO doesn't fit better
fitc.spn2 <- flexsurvspline(Surv(survtime, dead2) ~ age_10 + sex + stage, k=2, scale="normal", data=tcrb) # normal model is worse
save(fitc.f, fitc.g, fitc.sp2, fitc.sp3, fitc.sp4, fitc.spo2, fitc.spn2, file="~/work/flexsurv/tests/tcrcov.rda")
## GF fits best, AFT assumption must be better than PH.

plot(survfit(Surv(survtime, dead2) ~ stage, data=tcrb, subset=(tcrb$age_10=="50-59" & tcrb$sex=="female")))
lines(fitc.f, X=rbind(c(1,0,0,0,1,0,0,0),
              c(1,0,0,0,1,1,0,0),
              c(1,0,0,0,1,0,1,0),
              c(1,0,0,0,1,0,0,1)))

fitc.f <- flexsurvreg(Surv(survtime, dead2) ~ stage, data=tcrb, dist="genf")
plot(fitc.f)
plot(fitc.f, type="hazard", min.time=0, max.time=9)
plot(survfit(Surv(survtime, dead2) ~ stage, data=tcrb))
lines(fitc.f, X=rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1)))

fitc.f <- flexsurvreg(Surv(survtime, dead2) ~ sex, data=tcrb, dist="genf")
plot(fitc.f, type="hazard", min.time=0, max.time=10)
