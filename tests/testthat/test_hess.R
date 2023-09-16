if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()
if (!require("numDeriv")) return()

library(numDeriv)
library(testthat)
pars <- log(c(shape=1.2,scale=1.1))

test_that("Weibull AFT Hessian",{
  lds <- function(pars){dweibull(2,exp(pars[1]),exp(pars[2]),log=TRUE)}
  expect_equivalent(grad(lds, pars),
                    flexsurv:::DLdweibull(2, exp(pars[1]), exp(pars[2])))
  expect_equivalent(hessian(lds, pars),
                    D2Ldweibull(2, exp(pars[1]), exp(pars[2]))[1,,])

  lss <- function(pars){pweibull(2,exp(pars[1]),exp(pars[2]),log=TRUE, lower.tail = FALSE)}
  expect_equivalent(grad(lss, pars),
                    flexsurv:::DLSweibull(2, exp(pars[1]), exp(pars[2])))
  expect_equivalent(hessian(lss, pars),
                    D2LSweibull(2, exp(pars[1]), exp(pars[2]))[1,,])
})

test_that("Weibull PH Hessian",{
  lds <- function(pars){dweibullPH(2,exp(pars[1]),exp(pars[2]),log=TRUE)}
  expect_equivalent(grad(lds, pars),
                    flexsurv:::DLdweibullPH(2, exp(pars[1]), exp(pars[2])))
  lss <- function(pars){pweibullPH(2,exp(pars[1]),exp(pars[2]),log=TRUE, lower.tail = FALSE)}
  expect_equivalent(grad(lss, pars),
                    flexsurv:::DLSweibullPH(2, exp(pars[1]), exp(pars[2])))

  expect_equivalent(hessian(lds, pars),
                    D2LdweibullPH(2, exp(pars[1]), exp(pars[2]))[1,,])
  hess <- hessian(lss, pars)
  expect_equivalent(hessian(lss, pars),
                    D2LSweibullPH(2, exp(pars[1]), exp(pars[2]))[1,,])
})

test_that("Gompertz Hessian",{
  pars <- c(shape=1.2,scale=log(1.1))

  ## nonzero shape
  lds <- function(pars){dgompertz(2,pars[1],exp(pars[2]),log=TRUE)}
  lss <- function(pars){pgompertz(2,pars[1],exp(pars[2]),log=TRUE, lower.tail = FALSE)}
  expect_equivalent(grad(lds, pars),
                    flexsurv:::DLdgompertz(2, pars[1], exp(pars[2])))
  expect_equivalent(grad(lss, pars),
                    flexsurv:::DLSgompertz(2, pars[1], exp(pars[2])))
  expect_equivalent(hessian(lss, pars),
                    D2LSgompertz(2, pars[1], exp(pars[2]))[1,,])
  expect_equivalent(hessian(lss, pars),
                    D2Ldgompertz(2, pars[1], exp(pars[2]))[1,,])

  ## zero shape
  pars <- c(shape=0,scale=log(1.1))
  lds <- function(pars){dgompertz(2,0,exp(pars[2]),log=TRUE)}
  lss <- function(pars){pgompertz(2,0,exp(pars[2]),log=TRUE, lower.tail = FALSE)}
  expect_equivalent(grad(lds, pars),
                    flexsurv:::DLdgompertz(2, pars[1], exp(pars[2])))
  expect_equivalent(grad(lss, pars),
                    flexsurv:::DLSgompertz(2, pars[1], exp(pars[2])))
  expect_equivalent(hessian(lss, pars),
                    D2LSgompertz(2, pars[1], exp(pars[2]))[1,,])
  expect_equivalent(hessian(lds, pars),
                    D2Ldgompertz(2, pars[1], exp(pars[2]))[1,,])
})


hess_error <- function(object){
    if (!isTRUE(getOption("flexsurv.test.analytic.derivatives")))
        stop("flexsurv.test.analytic.derivatives option not set")
    object$hess.test$error
}
err <- 1e-03

options(flexsurv.test.analytic.derivatives=TRUE)

test_that("flexsurvreg fit hessian",{
  fl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="exp")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data=bc, dist="exp")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="weibull")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, dist="weibull")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="weibullPH")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data=bc, dist="weibullPH")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, dist="weibullPH")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gompertz")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data=bc, dist="gompertz")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, dist="gompertz")
  expect_lt(hess_error(fl), err)

  set.seed(1)
  simt <- rweibull(1000, 2, 0.5)
  status <- ifelse(simt>0.6, 0, 1)
  simt[status==0] <- 0.6
  tmin <- simt
  tmax <- ifelse(status==1, simt, Inf)
  tmax.sr <- ifelse(status==1, simt, NA)
  fl <- flexsurvreg(Surv(tmin, tmax.sr, type="interval2") ~ 1, dist="weibull")
  expect_lt(hess_error(fl), err)
  
  set.seed(1)
  sim <- rgenf(3000, 1.5, 1, -0.4, 0.6)
  dead <- as.numeric(sim<=30)
  simt <- ifelse(sim<=30, sim, 30)
  obs <- simt>3; simt <- simt[obs]; dead <- dead[obs]
  fl <- flexsurvreg(Surv(simt, dead) ~ 1, dist="weibull")
  expect_lt(hess_error(fl), err)

  wts <- rep(1, length(simt)); wts[1:200] <- 1.2
  fl <- flexsurvreg(Surv(simt, dead) ~ 1, weights=wts, dist="weibull")
  expect_lt(hess_error(fl), err)
  
  bc$bhaz <- rep(0.1, nrow(bc))
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, bhazard=bhaz, 
                    dist="weibull")
  expect_lt(hess_error(fl), err)
  
  bc$wts <- 1; bc$wts[1:300] <- 1.1
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, 
                    weights=wts, dist="weibull")
  expect_lt(deriv_error(fl), err)
  expect_lt(hess_error(fl), err)
  
  bc$bhaz <- rep(0.1, nrow(bc))
  bc$wts <- 1; bc$wts[1:300] <- 1.1
  fl <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    anc=list(shape=~group), data=bc, bhazard=bhaz, 
                    weights=wts,
                    dist="weibull")
  expect_lt(deriv_error(fl), err)
  expect_lt(hess_error(fl), err)
})

## TODO documentation 


test_that("flexsurvspline fit hessian",{
  fl <- flexsurvspline(formula = Surv(recyrs, censrec) ~ group,
                     k = 1, data=bc, scale = "hazard")
  expect_lt(hess_error(fl), err)
  fl <- flexsurvspline(formula = Surv(recyrs, censrec) ~ group,
                     k = 1, data=bc, scale="odds")
  expect_lt(hess_error(fl), err)
})
