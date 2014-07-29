test_that("normboot.flexsurvreg",{
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
    set.seed(1); b2 <- normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
    expect_equal(b1, b2)
})
