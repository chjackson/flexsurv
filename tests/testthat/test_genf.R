context("Generalized F distribution")

tol <- 1e-06

test_that("Generalized F",
          {
    expect_equal(dgenf(c(-1,1,2,3,4), mu=0, sigma=1, Q=0, P=1),
                                        # FIXME add limiting value for x=0
                 c(0, 0.353553390593274, 0.140288989252053, 0.067923038519582, 0.038247711235678), tol=tol)
    expect_error(Hgenf(c(-1,1,2,3,4), mu=0, sigma=1, P=1), "argument \"Q\" is missing")
    expect_error(Hgenf(c(-1,1,2,3,4), mu=0, sigma=1, Q=1), "argument \"P\" is missing")
})

test_that("Generalized F reduces to generalized gamma: d",{
    expect_equal(dgenf(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0, P=0),
                 dgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0))
    x <- c(-1,0,1,2,3,4); mu <- 2.2; sigma <- 1.6; Q <- 0.2; P <- 1.2
    delta <- (Q^2 + 2*P)^{1/2}
    s1 <- 2 / (Q^2 + 2*P + Q*delta); s2 <- 2 / (Q^2 + 2*P - Q*delta)
    expect_equal(dgenf(x, mu=mu, sigma=sigma, Q=Q, P=P),
                 dgenf.orig(x, mu=mu, sigma=sigma/delta, s1=s1, s2=s2))
    expect_equal(dgenf(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=c(0,1,2), P=0),
                 dgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=c(0,1,2)))
})

x <- c(-1,1,2,3,4); mu <- 2.2; sigma <- 1.6; Q <- 0; P <- 1
delta <- (Q^2 + 2*P)^{1/2}
s1 <- 2 / (Q^2 + 2*P + Q*delta); s2 <- 2 / (Q^2 + 2*P - Q*delta)

test_that("Generalized F reduces to log logistic",
          {
    expect_equal(dgenf(x, mu=mu, sigma=sigma, Q=Q, P=P),
                 dllogis(x, shape=sqrt(2)/sigma, scale=exp(mu)), tol=tol)
})

test_that("Generalized F reduces to gamma",{
    expect_equal(dgengamma(x, mu=mu, sigma=sigma, Q=sigma),
         dgamma(x, shape=1/sigma^2, scale=exp(mu)*sigma^2))
})

test_that("Generalized F reduces to generalized gamma: p",{
    expect_equal(pgenf(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0, P=1),
         c(0, 0, 0.5, 0.727159434644773, 0.825443507527843, 0.876588815661789))
    expect_equal(pgenf(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0, P=0),
         pgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0))
    expect_equal(pgenf(x, mu=mu, sigma=sigma, Q=Q, P=P),
         pgenf.orig(x, mu=mu, sigma=sigma/delta, s1=s1, s2=s2))
})

test_that("Generalized F reduces to generalized gamma: q",{
    expect_equal(qgenf(p=0.25, mu=0, sigma=1, Q=0, P=1), 0.459858613264917)
    expect_equal(qgenf(p=0.25, mu=0, sigma=1, Q=0, P=1), qgeneric(pgenf, p=0.25, mu=0, sigma=1, Q=0, P=1))
    expect_equal(qgenf(p=0, mu=0, sigma=1, Q=0, P=1), 0)
    expect_equal(qgenf(p=0.25, mu=0, sigma=1, Q=0, P=0),  qgengamma(p=0.25, mu=0, sigma=1, Q=0))
    expect_equal(qgenf(pgenf(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=0, P=1), mu=0, sigma=1, Q=0, P=1),
         c(0,0,0,1,2,3,4))
    expect_equal(qgenf(pgenf(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=-1, P=1), mu=0, sigma=1, Q=-1, P=1),
         c(0,0,0,1,2,3,4))
    expect_equal(qgengamma(pgenf(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=0, P=0), mu=0, sigma=1, Q=0),
         c(0,0,0,1,2,3,4))
    x <- c(0.1,0.4,0.6)
    expect_equal(qgenf(x, mu=mu, sigma=sigma, Q=Q, P=P),
         qgenf.orig(x, mu=mu, sigma=sigma/delta, s1=s1, s2=s2))
})

test_that("Generalized F reduces to generalized gamma: r",{
    rgenf(n=10, mu=0, sigma=1, Q=0, P=1)
    set.seed(22061976)
    x <- rgenf(n=10, mu=0, sigma=1, Q=0, P=0)
    set.seed(22061976)
    y <- rgengamma(n=10, mu=0, sigma=1, Q=0)
    expect_equal(x, y)
    if (interactive())  {
        x <- c(-1,0,1,2,3,4); mu <- 2.2; sigma <- 0.6; Q <- 0.2; P=0.1
        delta <- (Q^2 + 2*P)^{1/2}
        s1 <- 2 / (Q^2 + 2*P + Q*delta); s2 <- 2 / (Q^2 + 2*P - Q*delta)
        plot(density(rgenf(10000, mu=mu, sigma=sigma, Q=Q, P=P)))
        lines(density(rgenf.orig(10000, mu=mu, sigma=sigma/delta, s1=s1, s2=s2)), lty=2)
    }
})

test_that("Gives errors",{
    expect_warning(dgenf(c(1,1), 1, -2, 1, -1), "Negative scale")
    expect_warning(dgenf(c(1,1), 1, -2, 1, -1), "Negative shape")
    expect_warning(pgenf(1, 1, -2, 1, 1), "Negative")
    expect_warning(qgenf(0.1, 1, -2, 1, 0), "Negative")
    expect_warning(rgenf(4, 1, -2, 1, 1), "Negative")
})
