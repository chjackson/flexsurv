context("Testing generalized F (original)")

tol <- 1e-06

test_that("Generalized F (original)",{
    expect_equal(dgenf.orig(c(-1,1,2,3,4), mu=0, sigma=1, s1=1, s2=1),
         c(0, 0.25, 0.111111111111111, 0.0625, 0.04))
    x <- c(-1,0,1,2,3,4); mu <- 0.1; sigma <- 1.2; s1 <- 1.7; s2 <- 10000000
    dgenf.orig(x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    dgengamma.orig(x, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1) # equal for large s2

    expect_equal(pgenf.orig(c(-1,0,1,2,3,4), mu=0, sigma=1, s1=1, s2=1),
         c(0, 0, 0.5, 0.666666666666667, 0.75, 0.8))
    x <- c(-1,0,1,2,3,4); mu <- 0.1; sigma <- 1.2; s1 <- 1.7; s2 <- 10000000
    pgenf.orig(x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    pgengamma.orig(x, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1) # equal for large s2

    expect_equal(qgenf.orig(p=0.25, mu=0, sigma=1, s1=1, s2=1), 0.333333333333333)
    expect_equal(qgenf.orig(p=0.25, mu=0, sigma=1, s1=1, s2=1), qgeneric(pgenf.orig, p=0.25, mu=0, sigma=1, s1=1, s2=1))
    expect_equal(qgenf.orig(p=0, mu=0, sigma=1, s1=1, s2=1), 0)
    expect_equal(qgenf.orig(pgenf.orig(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, s1=1, s2=1), mu=0, sigma=1, s1=1, s2=1),
         c(0,0,0,1,2,3,4))
    x <- c(0.1, 0.4, 0.7); mu <- 0.1; sigma <- 1.2; s1 <- 1.7; s2 <- 10000000
    qgenf.orig(x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    hgenf.orig(x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    qgengamma.orig(x, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1) # equal for large s2
    expect_error(Hgenf.orig(x, mu=mu, sigma=sigma, s1=s1), "argument \"s2\" is missing")

    rgenf.orig(n=10, mu=0, sigma=1, s1=1, s2=1)
    if (interactive())  {
        mu <- 0.1; sigma <- 1.2; s1 <- 1.7; s2 <- 100000000
        plot(density(rgenf.orig(10000, mu=mu, sigma=sigma, s1=s1, s2=s2)))
        lines(density(rgengamma.orig(10000, shape=1/sigma, scale=exp(mu) / s1^sigma, k=s1)), lty=2)
    }
})
