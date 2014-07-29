
bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
bexp2 <- flexsurvreg(Surv(time, status) ~ trans, data=bosms3, dist="exp") 
test_that("flexsurv exponential model matches msm",{
    Q3 <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))
    bmsm <- msm(state ~ time, subject=ptnum, data=bos3, exacttimes=TRUE, qmatrix=Q3)
    fs.q <- as.numeric(exp(coef(bexp)[1] + c(0, coef(bexp)[2:3])))
    msm.q <- qmatrix.msm(bmsm, ci="none")[rbind(c(1,2),c(1,3),c(2,3))]
    expect_equal(fs.q, msm.q, tol=1e-05)
    ## markov + "semi-markov" models same in this case
    expect_equal(bexp$loglik, bmsm$minus2loglik*-0.5, tol=1e-06)
    expect_equal(bexp2$loglik, bmsm$minus2loglik*-0.5, tol=1e-06)
})

