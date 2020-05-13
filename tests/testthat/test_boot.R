bosms3$x <- rnorm(nrow(bosms3))
bexp.markov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")
bexp.markov.cov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + x, data=bosms3, dist="exp")
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))

test_that("bootstrap CIs in multi state models",{
    expect_error(
        pmatrix.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1), ci=TRUE, B=3)
      , NA)
    expect_error(
        totlos.fs(bexp.markov, t=c(5), trans=tmat, ci=TRUE, B=5)
      , NA)
})

test_that("bootstrap CIs with multicore",{
    expect_error(
        pmatrix.simfs(bexp.markov.cov, t=5, trans=tmat, newdata=list(x=1), M=1000, ci=TRUE, B=3, cores=2)
      , NA)
    expect_error(
        totlos.simfs(bexp.markov, t=5, trans=tmat, ci=TRUE, M=1000, B=5, cores=2)
      , NA)
})

test_that("bootstrapping a function within a function",{

    fun <- function(x, t, trans, newdata, M, ci, B=10, cores=1){
        pmatrix.simfs(x, t=t, trans=trans, newdata=newdata, M=M, ci=ci, B=B, cores=cores)
    }

    expect_error(
        fun(bexp.markov.cov, t=5, trans=tmat, newdata=list(x=1), M=1000, ci=TRUE, B=3),
        NA)

    expect_error(
        fun(bexp.markov.cov, t=5, trans=tmat, newdata=list(x=1), M=1000, ci=TRUE, B=3, cores=2)
        ,NA)
})
