context("Local multi-state modelling tests")

library(msm)
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
Q3 <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))
bos3$years <- round((bos3$time - 6)*30) / 365.25
bmsm <- msm(state ~ years, subject=ptnum, data=bos3, exacttimes=TRUE, qmatrix=Q3)
bwei <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + shape(trans), data=bosms3, dist="weibull") 
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

test_that("flexsurv exponential model matches msm",{
    fs.q <- as.numeric(exp(coef(bexp)[1] + c(0, coef(bexp)[2:3])))
    msm.q <- qmatrix.msm(bmsm, ci="none")[rbind(c(1,2),c(1,3),c(2,3))]
    expect_equal(fs.q, msm.q, tol=1e-03)
    ## markov + "semi-markov" models same in this case
    expect_equal(bexp$loglik, bmsm$minus2loglik*-0.5, tol=1e-06)
    bexp2 <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp") 
    expect_equal(bexp2$loglik, bmsm$minus2loglik*-0.5, tol=1e-06)
})

test_that("P matrix using numerical ODE integrator",{
    ## one or more times, CIs or no CIs
    (P <- pmatrix.fs(bwei, tmat))
    (P <- pmatrix.fs(bwei, tmat, ci=TRUE, B=5))
    (P <- pmatrix.fs(bwei, t=c(5,10), tmat))
    (P <- pmatrix.fs(bwei, t=c(5,10), tmat, ci=TRUE, B=5))

    ## with covariates: newdata argument
    bosms3$x <- rbinom(nrow(bosms3), 1, 0.4)
    bxwei <- flexsurvreg(Surv(Tstart, Tstop, status) ~ (trans + shape(trans))*x, data=bosms3, dist="weibull") 
    (P <- pmatrix.fs(bxwei, t=c(5,10), newdata=list(x=rep(0,3)), tmat, ci=TRUE, B=5))
    (P <- pmatrix.fs(bxwei, t=c(5,10), newdata=list(x=0), tmat, ci=TRUE, B=5))
    (P <- pmatrix.fs(bxwei, t=c(5,10), newdata=list(x=1), tmat, ci=TRUE, B=5))
  
})

test_that("P matrix for exponential model matches msm",{
    expect_equal(pmatrix.msm(bmsm, t=1)[1,1], pmatrix.fs(bexp, tmat)[1,1], tol=1e-04)
    expect_equal(pmatrix.msm(bmsm, t=1)[1,2], pmatrix.fs(bexp, tmat)[1,2], tol=1e-03)
    expect_equal(pmatrix.msm(bmsm, t=1)[2,3], pmatrix.fs(bexp, tmat)[2,3], tol=1e-03)
    pmatrix.fs(bexp, tmat, ci=TRUE, B=5)
    pmatrix.fs(bexp, tmat, t=c(5,10))
    pmatrix.fs(bexp, tmat, t=c(5,10), ci=TRUE, B=5)
})
    
test_that("ODE method matches Aalen-Johansen",{
    tgrid <- seq(0,14,by=0.01)
    mwei <- msfit.flexsurvreg(bwei, t=tgrid, trans=tmat, tvar="trans")
    pt <- probtrans(mwei, predt=0)
    pa <- as.numeric(pt[[1]][pt[[1]]$time==1,-1][1:3])
    po <- pmatrix.fs(bwei, tmat)[1,]
    expect_equal(pa, po, tol=1e-03)
    ## With a finer time grid, AJ approaches the ODE solution.  faster
    ## than accurate 0.01, though slower than less accurate 0.1.  Though
    ## probtrans returns standard errors efficiently, ODE would need
    ## normal resampling.
})

test_that("bootstrap CIs in multi state models",{
    bosms3$x <- rnorm(nrow(bosms3))
    bexp.markov.cov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + x, data=bosms3, dist="exp")
    pmatrix.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1), ci=TRUE, B=3)
    totlos.fs(bexp.markov, t=c(5), trans=tmat, ci=TRUE, B=5)
    tl <- totlos.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1), ci=TRUE, B=5)
    pmatrix.simfs(bexp.markov.cov, t=5, trans=tmat, newdata=list(x=1), M=10, ci=TRUE, B=3)
    totlos.simfs(bexp.markov.cov, t=5, trans=tmat, newdata=list(x=1), M=10, ci=TRUE, B=3)

    bwei.list <- bweic.list <- bweim.list <- bexpc.list <- vector(3, mode="list")
    for (i in 1:3) {
        bwei.list[[i]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==i),
                                      data = bosms3, dist = "weibull")
        bweic.list[[i]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==i),
                                       data = bosms3, dist = "weibull")
        bweim.list[[i]] <- flexsurvreg(Surv(Tstart, Tstop, status) ~ 1, subset=(trans==i),
                                       data=bosms3, dist="weibull")
        bexpc.list[[i]] <- flexsurvreg(Surv(years, status) ~ x, subset=(trans==i), data=bosms3, dist="exp")
    }

    totlos.simfs(bwei.list, t=5, trans=tmat, M=10, ci=TRUE, B=10)
    totlos.simfs(bweic.list, t=5, trans=tmat, M=100, newdata=list(x=0), ci=TRUE, B=10)
    pmatrix.simfs(bwei.list, t=5, trans=tmat, M=100, ci=TRUE, B=10)
    pmatrix.simfs(bweic.list, t=5, trans=tmat, M=100, newdata=list(x=0), ci=TRUE, B=10)
    pmatrix.fs(bweim.list, t=5, trans=tmat, ci=TRUE, B=10)
    pmatrix.fs(bweim.list, t=c(5,10), trans=tmat, ci=TRUE, B=10)
    totlos.fs(bweim.list, t=5, trans=tmat, ci=TRUE, B=10)
})

## Qualitative comparisons for msfit variance between list, non-list
if (0) { 
    ms1 <- msfit.flexsurvreg(bexpci, newdata=list(x=1), trans=tmat, t=1:10, variance=TRUE, B=1000)
    ms2 <- msfit.flexsurvreg(bexpc.list, newdata=list(x=1), trans=tmat, t=1:10, variance=TRUE, B=1000)
    ms1$varHaz[1:10,]
    ms2$varHaz[1:10,]
    ms1$varHaz[21:30,] # these cross-correlations take B=1000 to visibly converge
    ms2$varHaz[21:30,]
}
