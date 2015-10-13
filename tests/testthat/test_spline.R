context("flexsurvspline model fits")

set.seed(11082012)
bc$foo <- factor(sample(1:3, nrow(bc), replace=TRUE))

spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=0)

test_that("Basic flexsurvspline, Weibull",{
    expect_equal(spl$loglik,  -810.926859076725, tol=1e-06)
})

test_that("flexsurvspline summary method",{
    summ <- summary(spl, B=3)$"group=Good,foo=1"
    expect_equal(summ$est[1],0.999838214156694, tol=1e-05)
    expect_true(all(summ$est > summ$lcl))
    expect_true(all(summ$est < summ$ucl))
    expect_true(all(summ$ucl < 1))
    expect_true(all(summ$lcl > 0))
    summ <- summary(spl, type="survival", B=3, t=1:5)$`group=Good,foo=1`
    expect_equal(summ$est[1], 0.9684014494284117, tol=1e-05)
    summ <- summary(spl, type="cumhaz", B=3, t=1:5)$`group=Good,foo=1`
    expect_equal(summ$est[1],  0.03210855719443461, tol=1e-05)
    summ <- summary(spl, type="hazard", B=3, t=1:5)$`group=Good,foo=1`
    expect_equal(summ$est[1],  0.04446356223043756, tol=1e-05)
})

if (interactive()){
    test_that("flexsurvspline plot method",{
        plot(spl, col=c("red","blue","green"))
        plot(spl, ci=TRUE, B=40)
        plot(spl, type="cumhaz")
        plot(spl, type="hazard")
        plot(spl, type="hazard", ci=TRUE, B=40)
        ## multicoloured plots
        plot(spl, col=c("red","purple","blue","green","brown","black","orange","yellow","pink"))
        plot(spl, lwd=1)
        lines(spl, X=rbind(c(0,0,0,0)), col="red")
        lines(spl, X=rbind(c(1,0,0,0)),  col="purple")
        lines(spl, X=rbind(c(0,1,0,0)), col="blue")
        legend("bottomleft", levels(bc$group), col=c("red","purple","blue"), lwd=2)
    })
}

test_that("Basic flexsurvspline, Weibull, no covs",{
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0)
    expect_equal(-873.207054864145, spl$loglik, tol=1e-06)
})

test_that("Basic flexsurvspline, one knot, best fitting in paper",{
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="odds")
    expect_equal(spl$loglik, -788.981901798638, tol=1e-06)
    expect_equal(spl$loglik  +   sum(log(bc$recyrs[bc$censrec==1])), -615.49431514184, tol=1e-06)
    expect_equal(spl$AIC  +   sum(log(bc$recyrs[bc$censrec==1])), 1761.45139087546, tol=1e-06)
    #results from the paper
    expect_equal(as.numeric(coef(spl))[1:3], c(-3.451, 2.915, 0.191), tol=1e-03)
    expect_equal(as.numeric(spl$res[1:3,"se"]), c(0.203, 0.298, 0.044), tol=1e-03)
    if (interactive()) {
        plot(spl)
        plot(spl, type="cumhaz")
        plot(spl, type="haz")
    }
})

test_that("Spline models with hazard and normal scales",{
    splh <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="hazard")
    expect_equal(splh$loglik, -792.863797823674, tol=1e-05)
    spln <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="normal")
    expect_equal(spln$loglik, -785.623256840396, tol=1e-05)
    if (interactive()){
        plot(spl, ci=TRUE, lwd.ci=1, B=30)
        lines(splh, col="blue", ci=TRUE, B=30)
        lines(spln, col="green", ci=TRUE, B=30)
    }
})

### log(H(t)) = g0 + g1 log(t) + Bz
### H(t) = exp(g0) t^g1 exp(Bz)
### S(t)  = exp(-H(t))   = exp( - exp(g0) t^g1 exp(Bz) )
### pweibull has par exp( - (x/b) ^ a)
### a = g1, 1/b^a = exp(g0 + Bz) = b ^ -a  = exp(-a log b)
### g0 + Bz = - a (log b + Bz)

test_that("Spline proportional hazards models reduce to Weibull",{
    wei <- survreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
    wei.base <- survreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
    a <- 1/wei$scale
    b1 <- exp(coef(wei)[1]); b2 <- exp(coef(wei)[1]+coef(wei)[2]); b3 <- exp(coef(wei)[1]+coef(wei)[3])
    a.base <- 1/wei.base$scale
    b.base <- exp(coef(wei.base[1]))
    ## Compare three implementations of the Weibull, with and without covs
    fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull", fixedpars=FALSE,
                       inits=c(a,b1,coef(wei)[2:3]))
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0,
                          inits=c(-a*log(b1), a, -a*coef(wei)[2:3]), fixedpars=FALSE)
    expect_equal(fit$loglik, spl$loglik)
    expect_equal(fit$loglik, wei$loglik[2])

    fit <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull", fixedpars=FALSE,
                       inits=c(a,b1)) # NaNs here dues to zeroes for dweibull being visited by optimiser

    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0,
                          inits=c(log(1 / b.base^a.base), a.base), fixedpars=FALSE)
    expect_equal(fit$loglik, spl$loglik)
    expect_equal(fit$loglik, wei.base$loglik[1])
})

test_that("Spline proportional odds models reduce to log-logistic",{
    if (is.element("eha", installed.packages()[,1])) {
        library(eha)
        custom.llogis <- list(name="llogis",
                              pars=c("shape","scale"),
                              location="scale",
                              transforms=c(log, log),
                              inv.transforms=c(exp, exp),
                              inits=function(t){ c(1, median(t)) })
        fitll <- flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc, dist=custom.llogis)
        fitsp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, scale="odds")
        expect_equal(fitsp$loglik, fitll$loglik)
        expect_equal(1/fitll$res["scale",1]^fitll$res["shape",1], exp(fitsp$res["gamma0",1]), tol=1e-02)
        expect_equal(fitsp$res["gamma1",1], fitll$res["shape",1], tol=1e-02)
        detach("package:eha")
    }
})

test_that("Spline normal models reduce to log-normal",{
    fitln <- flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc, dist="lnorm")
    fitsp <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, scale="normal")
    expect_equal(fitsp$res["gamma0",1], -fitln$res["meanlog",1]/fitln$res["sdlog",1], tol=1e-02)
    expect_equal(fitsp$res["gamma1",1], 1 /fitln$res["sdlog",1], tol=1e-02)
})

test_that("Spline models with time-varying covariate effects",{
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group), data=bc, knots=1)
    expect_equal(spl$loglik, -789.002098222827, tol=1e-04)
    expect_equal(spl$res["gamma1(groupMedium)","est"],  -0.410859123437426, tol=1e-04)
    summ.g <- summary(spl, ci=FALSE, t=c(1,5,10))$`group=Good`
    expect_equal(summ.g$est[1], 0.9853637097617495, tol=1e-04)
    summ.m <- summary(spl, ci=FALSE, t=c(1,5,10))$`group=Medium`
    expect_equal(summ.m$est[1], 0.9391492993851021, tol=1e-04)
    summ2 <- summary(spl, ci=FALSE, t=c(1,5,10), newdata=data.frame(group=c("Good","Medium")))
    expect_equal(summ2[[1]]$est[1], summ.g$est[[1]])
    expect_equal(summ2[[2]]$est[1], summ.m$est[[1]])
    spl2 <- flexsurvspline(Surv(recyrs, censrec) ~ group, anc=list(gamma1=~group), data=bc, knots=1)
    expect_equal(spl$loglik, spl2$loglik)
    if (interactive()) plot(spl)
    spl3 <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group), data=bc, k=2)
    if (interactive()) lines(spl3, col="blue")
})

test_that("Spline models with left-truncation",{
    bc2 <- bc[bc$recyrs>2,]
    (spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc2, k=0, fixedpars=TRUE))
    expect_equal(spl$loglik, -951.4567686517325, tol=1e-06)
    expect_equal(spl$loglik, sum(spl$logliki), tol=1e-06)
    (spl <- flexsurvspline(Surv(rep(0, nrow(bc2)), recyrs, censrec) ~ 1, data=bc2, k=0))
    expect_equal(spl$loglik, -432.9048860461514, tol=1e-06)
    spl <- flexsurvspline(Surv(rep(1.9, nrow(bc2)), recyrs, censrec) ~ 1, data=bc2, k=0)
    expect_equal(spl$loglik, -400.3556493114977, tol=1e-06)
#    if(interactive()) lines(spl, col="blue") # truncated model fits much better
})

test_that("Spline models with weighting",{
    wt <- rep(1, nrow(bc)); wt[1:30] <- 2
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, weights=wt)
    wei <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull", weights=wt)
    expect_equal(spl$loglik, wei$loglik, tol=1e-06)
})

test_that("Spline models with relative survival",{
    bc$bh <- rep(0.01, nrow(bc))
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0, bhazard=bh)
#    wei <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibullPH", bhazard=bh)
#    expect_equal(spl$loglik, wei$loglik, tol=1e-06)
})

test_that("flexsurvspline results match stpm in Stata",{
    ## Numbers copied from Stata output for equivalent stpm commands
    ## see ~/flexsurv/stpm/do1.do
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0)
    expect_equal(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])), -638.45432, tol=1e-05)
    expect_equal(as.numeric(spl$res[,"est"]), c(-3.360303, 1.379652, .8465394, 1.672433), tol=1e-04)
    
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=2)
    expect_equal(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])), -674.75128, tol=1e-04)
    expect_equal(as.numeric(spl$res[,"est"]), c(-1.728339,  3.476179,   .567432, -.3420749), tol=1e-02)
    
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=2, scale="odds")
    expect_equal(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])), -675.271, tol=1e-04)

    spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=2, scale="normal")
    expect_equal(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])), -675.73591 , tol=1e-04)

    ## stpm needs all or no ancillary pars to depend on covs
    ## TODO stpm2: our Stata 13.0 is too old for this, need 13.1
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group) + gamma2(group) + gamma3(group), data=bc, k=2, scale="hazard")
    expect_equal(spl$loglik  +  sum(log(bc$recyrs[bc$censrec==1])), -607.47942, tol=1e-04)
    ## coefficients are the same to about 2sf
})

test_that("Expected survival",{
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1)
    gamma <- coef(spl)[1:3]
    beta <- coef(spl)[4:5]
    surv <- function(x,...)psurvspline(q=x, gamma=gamma, beta=beta, knots=spl$knots, scale=spl$scale, lower.tail=FALSE, ...)
    expect_equal(integrate(surv, 0, 5, X=c(0,0))$value, 4.341222955052117, tol=1e-04)# For group="good"
    expect_equal(integrate(surv, 0, 5, X=c(1,0))$value, 3.664826479659649, tol=1e-04) # For group="medium"
    expect_equal(integrate(surv, 0, 5, X=c(0,1))$value, 2.713301623208948, tol=1e-04) # For group="poor"
})

test_that("gamma in d/psurvspline can be matrix or vector",{
    require(mvtnorm)
    boot <- rmvnorm(10, spl$opt$par, spl$cov)
    psurvspline(5, gamma=coef(spl)[1:2], knots=spl$knots, scale=spl$scale, lower.tail=FALSE)
    expect_equal(psurvspline(5, gamma=boot[3,1:2], knots=spl$knots, scale=spl$scale, lower.tail=FALSE),
                 psurvspline(5, gamma=boot[,1:2], knots=spl$knots, scale=spl$scale, lower.tail=FALSE)[3])
    dsurvspline(5, gamma=coef(spl)[1:2], knots=spl$knots, scale=spl$scale)
    expect_equal(dsurvspline(5, gamma=boot[3,1:2], knots=spl$knots, scale=spl$scale),
                 dsurvspline(5, gamma=boot[,1:2], knots=spl$knots, scale=spl$scale)[3])
})

test_that("Errors in spline",{
expect_error(flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, knots="foo"), "\"knots\" must be a numeric vector")
expect_error(flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, knots=c(-5, 2)), "knot -5 less than or equal to minimum log time")
expect_error(flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, knots=c(2)), "knot 2 greater than or equal to maximum log time")
expect_error(flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=0, bknots=c("foo")), "boundary knots should be")
expect_error(flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=0, bknots=c(0,1,2)), "boundary knots should be")
expect_warning(
    flexsurvspline(Surv(recyrs, censrec) ~ foo, data=bc, subset=1:10, k=0, inits=c(1,1,0,0,0,0), fixedpars=TRUE)
  , "minimum and maximum log death times are the same")
})

test_that("supplying knots",{
    ldtimes <- log(bc$recyrs[bc$censrec==1])
    expect_equal(flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1)$loglik,
                 flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, knots=median(ldtimes), bknots=range(ldtimes))$loglik, tol=1e-08)   
    expect_true(flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=1)$loglik !=
                flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=1, bknots=c(-5, 2))$loglik)
})

test_that("subset",{
    subflex <- flexsurvspline(Surv(time = time, event = status) ~ sex + cut(age, 4), data = survival::lung,
                              subset = !is.na(wt.loss))
    subflex2 <- flexsurvspline(Surv(time = time, event = status) ~ sex + cut(age, 4), data = survival::lung[!is.na(survival::lung$wt.loss),])
    expect_equal(subflex$loglik, subflex2$loglik, tol=1e-08)
    subflex <- flexsurvspline(Surv(time = time, event = status) ~ sex + cut(age, 4), data = survival::lung, 
                              subset = survival::lung$age > 60) # empty factor level in subset, should be dropped since 0.7
})
