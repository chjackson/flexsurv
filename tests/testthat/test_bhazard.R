
test_that("Background hazards / relative survival", { 
    bc$bh <- rep(0.01, nrow(bc))

    ## Compare with stata stgenreg, using Weibull PH model
    fs6b <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibullPH", bhazard=bh)
      
    expect_equal(log(fs6b$res[1,"est"]), 0.3268327417773233, tolerance=1e-05)
    expect_equal(log(fs6b$res[2,"est"]), -3.5308925743338038, tolerance=1e-05)
    expect_equal(fs6b$res["groupMedium","est"], 0.9343799681269026, tolerance=1e-04)
    expect_equal(fs6b$res["groupPoor","est"], 1.799204192587765, tolerance=1e-04)

    ## Check fit from 3 par model reduces to 1 par 
    ## Deriv calculation bug causing false convergence fixed in 2.1
    fshgg <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="gengamma", 
                         inits=c(1,1,1), fixedpars=2:3, bhazard=bh)
    fshe <-  flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="exponential", bhazard=bh)
    fshw <-  flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull", 
                         inits = c(1,10), fixedpars=1, bhazard=bh)
    expect_equal(fshgg$loglik, fshe$loglik, tolerance=1e-06)
    expect_equal(fshgg$loglik, fshw$loglik, tolerance=1e-06)
    
    ## same results as 
    ## cd /home/chris/flexsurv/stata
    ## use stpm/bc
    ## gen rec = 1 - censrec
    ## gen recyrs = rectime / 365
    ## gen bh = 0.01
    ## stset recyrs, failure(censrec)
    ## stgenreg, loghazard([ln_lambda] :+ [ln_gamma] :+ (exp([ln_gamma]) :- 1) :* log(#t)) nodes(100) ln_lambda(group2 group3) bhazard(bh)

    ## Check we can convert from partial to full likelihood by adding the
    ## sum of the cumulative hazards 
    mdl_0 <- flexsurvreg(Surv(time/365, status == 2) ~ 1, data = lung, dist = "exp")
    bhaz <- 0.1
    mdl_1 <- flexsurvreg(Surv(time/365, status == 2) ~ 1, dist = "exp", data = lung, 
                         inits = mdl_0$res[, "est"] - bhaz, 
                         fixedpars = TRUE, bhazard = rep(bhaz, nrow(lung)))
  
    expect_equal(mdl_0$loglik, mdl_1$loglik - sum(bhaz * lung$time/365))
})

test_that("background hazards with delayed entry", {
    old <- options(flexsurv.test.analytic.derivatives=TRUE)
    on.exit(options(old))
    start <- rep(5,10)
    stop <- rep(7,10)
    fit <- flexsurvreg(Surv(start,stop,rep(1,10))~1, dist="exp", bhazard=rep(0.1, 10))
    expect_equal(unname(fit$res["rate","est"]), 0.4, tolerance=1e-4)
    expect_equal(fit$loglik, 10*(log(.5) - 0.4*2), tolerance=1e-7)
    expect_lt(fit$deriv.test$error, 1e-5)
    expect_lt(fit$hess.test$error, 1e-5)
})
