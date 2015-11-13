context("Distribution functions and utilities")

## note - standard q fns in R return zero for p=0 for positive dists, but -Inf for real dists.

tol <- 1e-06

test_that("Generalized F",{
    expect_equal(dgenf(c(-1,1,2,3,4), mu=0, sigma=1, Q=0, P=1), # FIXME add limiting value for x=0
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

test_that("Generalized F reduces to log logistic",{
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

x <- c(-1,1,2,3,4); shape <- 2.2; scale <- 1.6; k <- 1.9

test_that("Generalized gamma reductions: d",{
    expect_equal(dgengamma(c(-1,1,2,3,4), mu=0, sigma=1, Q=1),  # FIXME value for x=0 and add here
         c(0, 0.367879441171442, 0.135335283236613, 0.0497870683678639, 0.0183156388887342))
    expect_equal(dgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0),
         dlnorm(c(-1,0,1,2,3,4), meanlog=0, sdlog=1))
    expect_equal(dgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1:3, Q=0),
         dlnorm(c(-1,0,1,2,3,4), meanlog=0, sdlog=1:3))
    expect_equal(dgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1:3, Q=0, log=TRUE),
         dlnorm(c(-1,0,1,2,3,4), meanlog=0, sdlog=1:3, log=TRUE))
    expect_equal(dgengamma(c(1,2,3,4), mu=0.1, sigma=1.2, Q=1),
         dweibull(c(1,2,3,4), shape=1/1.2, scale=exp(0.1)))  # only defined for x>0 anyway
    x <- c(1,2,3,4); mu <- 0.4; sigma <- 1.2
    expect_equal(dgengamma(x, mu=mu, sigma=sigma, Q=sigma),
         dgamma(x, shape=1/sigma^2, scale=exp(mu)*sigma^2))
# FIXME add limiting value for x=0
    expect_equal(dgengamma.orig(x, shape=shape, scale=scale, k=k),
         dgengamma(x, mu=log(scale) + log(k)/shape, sigma=1/(shape*sqrt(k)), Q=1/sqrt(k)))
})

test_that("Generalized gamma reductions: p",{
    expect_equal(pgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=1),
         c(0, 0, 0.632120558828558, 0.864664716763387, 0.950212931632136, 0.981684361111266))
    expect_equal(pgengamma(c(-1,0,1,2,3,4), mu=0, sigma=1, Q=0),
         plnorm(c(-1,0,1,2,3,4), meanlog=0, sdlog=1))
    expect_equal(pgengamma(c(-1,0,1,2,3,4), mu=0.1, sigma=1.2, Q=1),
         pweibull(c(-1,0,1,2,3,4), shape=1/1.2, scale=exp(0.1)))
    x <- c(1,2,3,4); mu <- 0.4; sigma <- 1.2
    expect_equal(pgengamma(x, mu=mu, sigma=sigma, Q=sigma),
         pgamma(x, shape=1/sigma^2, scale=exp(mu)*sigma^2))
})

test_that("Generalized gamma reductions: q",{
    expect_equal(qgengamma(p=0.25, mu=0, sigma=1, Q=1), 0.287682072451781)
    expect_equal(qgengamma(p=0.25, mu=0, sigma=1, Q=1), qgeneric(pgengamma, p=0.25, mu=0, sigma=1, Q=1))
    expect_equal(qgengamma(p=0, mu=0, sigma=1, Q=1), 0)
    expect_equal(qgengamma(p=0.25, mu=0, sigma=1, Q=0), qlnorm(p=0.25, meanlog=0, sdlog=1))
    expect_equal(qgengamma(p=0.25, mu=0.1, sigma=1.2, Q=1), qweibull(p=0.25, scale=exp(0.1), shape=1/1.2))
    expect_equal(qgengamma(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=1), mu=0, sigma=1, Q=1),
         c(0,0,0,1,2,3,4))
    expect_equal(qgengamma(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=-1), mu=0, sigma=1, Q=-1),
         c(0,0,0,1,2,3,4))
    expect_equal(qlnorm(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=0), meanlog=0, sdlog=1),
         c(0,0,0,1,2,3,4))
    expect_equal(qlnorm(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=0, log=TRUE), meanlog=0, sdlog=1, log=TRUE),
         c(0,0,0,1,2,3,4))
    expect_equal(qgengamma(p=0.25, mu=0, sigma=1, Q=1, lower.tail=TRUE), qgeneric(pgengamma, p=0.25, mu=0, sigma=1, Q=1, lower.tail=TRUE))
})

test_that("Generalized gamma reductions: r",{
    rgengamma(n=10, mu=0, sigma=1, Q=0)
    set.seed(22061976)
    x <- rgengamma(n=10, mu=0, sigma=1.1, Q=0)
    set.seed(22061976)
    y <- rlnorm(n=10, meanlog=0, sdlog=1.1)
    expect_equal(x, y)
})

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

test_that("Generalized gamma (original)",{
    expect_equal(dgengamma.orig(c(-1,0,1,2,3,4), shape=1.2, scale=1.3, k=1.4),
         c(0, 0, 0.419477559803262, 0.260699967439176, 0.120081193404263, 0.0474236822588797))
    expect_equal(dgengamma.orig(c(1,2,3,4), shape=1.2, scale=1.3, k=1),
         dweibull(c(1,2,3,4), shape=1.2, scale=1.3))
    expect_equal(dgengamma.orig(c(1,2,3,4), shape=1, scale=1.3, k=1),
         dexp(c(1,2,3,4), rate=1/1.3))
    expect_equal(dgengamma.orig(c(1,2,3,4), shape=1, scale=1.3, k=1.4),
         dgamma(c(1,2,3,4), shape=1.4, scale=1.3))

    shape <- 1.2; scale <- 1.3; k <- 10000
    pgengamma.orig(2800 + 1:4, shape=shape, scale=scale, k=k)
    plnorm(2800 + 1:4, log(scale) + log(k)/shape, 1/(shape*sqrt(k)))

    expect_equal(pgengamma.orig(c(1,2,3,4), shape=1.2, scale=1.3, k=1),
         pweibull(c(1,2,3,4), shape=1.2, scale=1.3))
    expect_equal(pgengamma.orig(c(1,2,3,4), shape=1, scale=1.3, k=1),
         pexp(c(1,2,3,4), rate=1/1.3))
    expect_equal(pgengamma.orig(c(1,2,3,4), shape=1, scale=1.3, k=1.4),
         pgamma(c(1,2,3,4), shape=1.4, scale=1.3))

    expect_equal(qgengamma.orig(p=0.25, shape=1.2, scale=1.3, k=1), qgeneric(pgengamma.orig, p=0.25, shape=1.2, scale=1.3, k=1))
    expect_equal(qgengamma.orig(c(0.1, 0.4, 0.7), shape=1.2, scale=1.3, k=1),
         qweibull(c(0.1, 0.4, 0.7), shape=1.2, scale=1.3))
    expect_equal(qgengamma.orig(c(0.1, 0.4, 0.7), shape=1, scale=1.3, k=1),
         qexp(c(0.1, 0.4, 0.7), rate=1/1.3))
    expect_equal(qgengamma.orig(c(0.1, 0.4, 0.7), shape=1, scale=1.3, k=1.4),
         qgamma(c(0.1, 0.4, 0.7), shape=1.4, scale=1.3))

    if (interactive()){
        plot(density(rgengamma.orig(100000, shape=1.2, scale=1.3, k=1)))
        lines(density(rweibull(100000, shape=1.2, scale=1.3)), lty=2)
        plot(density(rgengamma.orig(100000, shape=1, scale=1.5, k=1)))
        lines(density(rexp(100000, rate=1/1.5)), lty=2)
        plot(density(rgengamma.orig(100000, shape=1, scale=3.3, k=1.2)))
        lines(density(rgamma(100000, shape=1.2, scale=3.3)), lty=2)
    }
})

test_that("Errors in generalized gamma and F",{
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative shape")
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative scale")
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=c(-1.2, 1), scale=1.3, k=1), "Negative shape")
    expect_warning(pgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative")
    expect_warning(qgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative")
    expect_warning(rgengamma.orig(3, shape=-1.2, scale=-1.3, k=-1), "Negative")
    expect_warning(dgengamma(1, 1, -2, 1), "Negative")
    expect_warning(pgengamma(1, 1, -2, 1), "Negative")
    expect_warning(qgengamma(0.1, 1, -2, 1), "Negative")
    expect_warning(rgengamma(1, 1, -2, 1), "Negative")
    expect_warning(dgenf(c(1,1), 1, -2, 1, -1), "Negative scale")
    expect_warning(dgenf(c(1,1), 1, -2, 1, -1), "Negative shape")
    expect_warning(pgenf(1, 1, -2, 1, 1), "Negative")
    expect_warning(qgenf(0.1, 1, -2, 1, 0), "Negative")
    expect_warning(rgenf(4, 1, -2, 1, 1), "Negative")
})

## test haz and cum haz functions
x <- seq(0.1, 100, by=0.1)
if (interactive())
    plot(x, hgengamma.orig(x, 1, 1, 1), type="l")
## Blows up at 35, since num and denom both converge to 0.

test_that("dgompertz",{
    x <- c(-1,0,1,2,3,4)
    expect_equal(dgompertz(x, shape=0.1, rate=0.2), c(0, 0.2, 0.179105591827508, 0.156884811322895, 0.134101872197705, 0.111571759992743))
    dgompertz(x, shape=0.0001, rate=0.2)
    dgompertz(x, shape=-0.0001, rate=0.2)
    dexp(x, rate=0.2)
    expect_equal(dgompertz(x, shape=0, rate=0.2), dexp(x, rate=0.2))
    d <- numeric(6); for (i in 1:6) d[i] <- dgompertz(x[i], shape=-0.0001, rate=i/5)
    expect_equal(d, dgompertz(x, shape=-0.0001, rate=1:6/5))
})

test_that("pgompertz",{
    x <- c(-1,0,1,2,3,4)
    pgompertz(x, shape=0, rate=0.2)
    pgompertz(x, shape=0.001, rate=0.2)
    pgompertz(x, shape=-0.001, rate=0.2)
    expect_equal(pgompertz(x, shape=0, rate=0.2), pexp(x, rate=0.2))
    p <- numeric(6); for (i in 1:6) p[i] <- pgompertz(x[i], shape=-0.0001, rate=i/5)
    expect_equal(p, pgompertz(x, shape=-0.0001, rate=1:6/5))
    expect_equal(p, 1 - exp(-Hgompertz(x, shape=-0.0001, rate=1:6/5)))
})

test_that("qgompertz",{
    x <- c(0.1, 0.2, 0.7)
    expect_equal(qgompertz(x, shape=0.1, rate=0.2), qgeneric(pgompertz, p=x, shape=0.1, rate=0.2))
    expect_equal(qgompertz(x, shape=0, rate=0.2), qexp(x, rate=0.2))
    expect_equal(x, pgompertz(qgompertz(x, shape=0.1, rate=0.2), shape=0.1, rate=0.2))
    q <- numeric(3); for (i in 1:3) q[i] <- qgompertz(x[i], shape=-0.0001, rate=i/5)
    expect_equal(q, qgompertz(x, shape=-0.0001, rate=1:3/5))
    x <- c(0.5, 1.06, 4.7)
    expect_equal(x, qgompertz(pgompertz(x, shape=0.1, rate=0.2), shape=0.1, rate=0.2))
    q <- numeric(3); for (i in 1:3) q[i] <- qgompertz(x[i], shape=-0.0001, rate=i/5)
    expect_equal(q, qgompertz(x, shape=-0.0001, rate=1:3/5))
    qgompertz(p=c(-1, 0, 1, 2), 1, 1)
})

test_that("rgompertz",{
    set.seed(1); rgompertz(4, shape=shape, rate=1:4)
    set.seed(1);
    rgompertz(1, shape=shape, rate=1)
    rgompertz(1, shape=shape, rate=2)
    rgompertz(1, shape=shape, rate=3)
    rgompertz(1, shape=shape, rate=4)
    rgompertz(1:2, shape=1, rate=1)
    rgompertz(3, shape=c(1,NaN, NA), rate=1)
    rgompertz(3, shape=c(NaN), rate=1)
    rgompertz(4, shape=c(1,NA), rate=1)
})

if (interactive()) {
    plot(density(rgompertz(10000, shape=0.1, rate=0.2)))
    x <- seq(0, 20, by=0.001)
    lines(x, dgompertz(x, shape=0.1, rate=0.2), lty=2)
}

test_that("Gompertz with chance of living forever",{
    shape <- -0.6; rate <- 1.8
    x <- c(0.8, 0.9, 0.97, 0.99)
    expect_equal(qgompertz(x, shape=shape, rate=rate), c(1.28150707286845, 2.4316450975351, Inf, Inf))
                                        # qgeneric(pgompertz, p=x, shape=shape, rate=rate) # won't work - needs smoothness
    expect_equal(pgompertz(Inf, shape=shape, rate=rate), 1)
})

test_that("Spline distribution functions",{
    regscale <- 0.786; cf <- 1.82
    a <- 1/regscale; b <- exp(cf)
    d1 <- dweibull(1, shape=a, scale=b)
    d2 <- dsurvspline(1, gamma=c(log(1 / b^a), a))
    expect_equal(d1, d2)
    p1 <- pweibull(1, shape=a, scale=b)
    p2 <- psurvspline(1, gamma=c(log(1 / b^a), a))
    expect_equal(p1, p2)
    meanlog <- 1.52; sdlog <- 1.11
    d1 <- dlnorm(1, meanlog, sdlog)
    d2 <- dsurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
    expect_equal(d1, d2)
    p1 <- plnorm(1, meanlog, sdlog)
    p2 <- psurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
    expect_equal(p1, p2)

    ## other way round 
    gamma <- c(0.1, 0.2)
    d1 <- dweibull(1, shape=gamma[2], scale= exp(-gamma[1]/gamma[2]))
    d2 <- dsurvspline(1, gamma=gamma)
    expect_equal(d1, d2)

    d1 <- dllogis(1, shape=gamma[2], scale= exp(-gamma[1]/gamma[2]))
    d2 <- dsurvspline(1, gamma=gamma, scale="odds")
    expect_equal(d1, d2)

    d1 <- dlnorm(1, meanlog=-gamma[1]/gamma[2], sdlog=1/gamma[2])
    d2 <- dsurvspline(1, gamma=gamma, scale="normal")
    expect_equal(d1, d2)

                                        # TODO document 
                                        #H1 <- Hlnorm(1, meanlog, sdlog)
                                        #H2 <- Hsurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
                                        #expect_equal(H1, H2)
    g <- c(0.1, 0.2, 0.3); k <- c(2,3,4)
    expect_equal(dsurvspline(1,g,knots=k)/(1 - psurvspline(1,g,knots=k)),  hsurvspline(1,g,knots=k))
    expect_equal(dsurvspline(1,g,knots=k,scale="odds")/(1 - psurvspline(1,g,knots=k,scale="odds")),  hsurvspline(1,g,knots=k,scale="odds"))
    expect_equal(dsurvspline(1,g,knots=k,scale="normal")/(1 - psurvspline(1,g,knots=k,scale="normal")),  hsurvspline(1,g,knots=k,scale="normal"))
    expect_equal(-log(1 - psurvspline(0.2,g,knots=k)), Hsurvspline(0.2,g,knots=k))
    expect_equal(-log(1 - psurvspline(0.2,g,knots=k,scale="odds")), Hsurvspline(0.2,g,knots=k,scale="odds"))
    expect_equal(-log(1 - psurvspline(0.2,g,knots=k,scale="normal")), Hsurvspline(0.2,g,knots=k,scale="normal"))

    expect_equal(dsurvspline(c(-1,0,NA,NaN,Inf), g, knots=k), c(0,0,NA,NA,NaN))
    expect_equal(psurvspline(qsurvspline(0.2,g,knots=k), g, knots=k), 0.2)
    expect_equal(qsurvspline(psurvspline(0.2,g,knots=k), g, knots=k), 0.2)
    expect_equal(psurvspline(c(NA,NaN,-1,0), gamma=c(1,1), knots=c(-10, 10)), c(NA,NA,0,0))
    expect_equal(qsurvspline(c(0,1), gamma=c(1,1), knots=c(-10, 10)), c(-Inf, Inf))

    expect_equal(1 - psurvspline(1, g, knots=k),  psurvspline(1, g, knots=k, lower.tail=FALSE))
    expect_equal(log(psurvspline(c(-1,NA,1), g, knots=k)),  psurvspline(c(-1,NA,1), g, knots=k, log.p=TRUE))

    ## single x=0: fully defined in in dbase.survspline
    expect_equal(dsurvspline(0, g, knots=k), 0)
    ## value for x=0?  currently zero, should it be limit as x reduces to 0? 
    expect_equal(hsurvspline(0, g, knots=k), 0)
    
    ## TODO special value handling and vectorisation for d function
    if(0){
        bc$foo <- factor(sample(1:3, nrow(bc), replace=TRUE))
        spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=0)
        system.time(hist(rsurvspline(10000, gamma=c(log(1 / b^a), a)), prob=TRUE))
        x <- 1:50
        lines(x, dweibull(x, a, b))
    }
})

test_that("Exponential hazards",{
    expect_equal(hexp(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,2,3)), c(0,NaN,3,NA,0,3,1))
    expect_equal(hexp(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,2,3), log=TRUE), log(c(0,NaN,3,NA,0,3,1)))
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3)), c(0,NaN,Inf, NA,0,0, 1, 4, 12))
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3), log=TRUE), log(c(0,NaN,Inf, NA,0,0, 1, 4, 12)))
    expect_equal(dexp(c(1,2,3), c(2,3,4)) / (1 - pexp(c(1,2,3), c(2,3,4))), hexp(c(1,2,3), c(2,3,4)))
    expect_equal(-log(1 - pexp(c(1,2,3), c(2,3,4))), Hexp(c(1,2,3), c(2,3,4)))
    expect_equal(log(-log(1 - pexp(c(1,2,3), c(2,3,4)))), Hexp(c(1,2,3), c(2,3,4), log=TRUE))
    expect_warning(hexp(c(1,1,1), c(-1, 0, 2)), "Negative rate")
    expect_warning(dexp(c(1,1,1), c(-1, 0, 2)), "NaNs produced")
})

test_that("Weibull hazards",{
    expect_equal(hweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)), 
         dweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)) / (1 - pweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3))))
    expect_equal(Hweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)),
         -pweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3), lower.tail=FALSE, log.p=TRUE))
    ## exponential reduction
    expect_equal(hweibull(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,1,1), c(1,2,3)), c(0,NaN,1/3,NA,0,1/3,1))
### positive shape, increasing 
    hweibull(c(-Inf, -1, 0, 1, 2, Inf), 2, 2.5)
### neg shape, decreasing 
    hweibull(c(-Inf, -1, 0, 1, 2, Inf), 0.5, 1)
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3)), c(0,NaN,Inf, NA,0,0, 1, 4, 12))
})

test_that("Gompertz hazards",{
    expect_equal(hgompertz(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)),
         dgompertz(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)) / (1 - pgompertz(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3))))
    expect_equal(Hgompertz(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)), 
         -pgompertz(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3), lower.tail=FALSE, log.p=TRUE))
    ## reduction to exponential
    expect_equal(hgompertz(c(2,4), c(0,0), c(2,2)), hexp(c(2,4), c(2,2)))
    expect_equal(Hgompertz(c(2,4), c(0,0), c(2,2)), Hexp(c(2,4), c(2,2)))
})

test_that("Generalized gamma definition in Stata manual",{
    Sgg <- function(t, mu, sigma, kappa){
        IGF <- function(a, x){ pgamma(x, a) } # incomplete gamma function   
        gamma <- 1 / kappa^2
        z <- (log(t) - mu)/sigma
        z[kappa<0] <- -z[kappa<0] # Stata manual uses sign(0) = 1
        u <- gamma*exp(abs(kappa)*z)
        ifelse(kappa > 0,
               1 - IGF(gamma, u),
               ifelse(kappa==0,
                      1 - pnorm(z),
                      IGF(gamma, u)))
    }
    expect_equal(
        Sgg(c(1,2,3),c(-1,2,0.2),c(1,2,1),c(-1, 0, 1)),
        pgengamma(c(1,2,3),c(-1,2,0.2),c(1,2,1),c(-1, 0, 1), lower.tail=FALSE)
        )
})

test_that("llogis",{
    x <- c(0.1, 0.2, 0.7)
    if (require("eha"))
        expect_equal(dllogis(x, shape=0.1, scale=0.2), eha::dllogis(x, shape=0.1, scale=0.2))
    expect_equal(qllogis(x, shape=0.1, scale=0.2), qgeneric(pllogis, p=x, shape=0.1, scale=0.2))
    expect_equal(x, pllogis(qllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2))
    expect_equal(x, 1 - exp(-Hllogis(qllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2)))
    expect_equal(hllogis(x, shape=0.1, scale=0.2),  dllogis(x, shape=0.1, scale=0.2) / (1 - pllogis(x, shape=0.1, scale=0.2)))
    q <- numeric(3); for (i in 1:3) q[i] <- qllogis(x[i], shape=0.0001, scale=i/5)
    expect_equal(q, qllogis(x, shape=0.0001, scale=1:3/5))
    x <- c(0.5, 1.06, 4.7)
    expect_equal(x, qllogis(pllogis(x, shape=0.1, scale=0.2), shape=0.1, scale=0.2))
    if (interactive()) {
        rl <- rllogis(10000, shape=1.5, scale=1.2)
        plot(density(rl[rl<100]), xlim=c(0,10))
        x <- seq(0, 10, by=0.001)
        lines(x, dllogis(x, shape=1.5, scale=1.2), lty=2)
    }
    expect_equal(mean.llogis(shape=0.1, scale=0.2), NaN)
    expect_equal(var.llogis(shape=1.1, scale=0.2), NaN)
    mean.llogis(shape=1.1, scale=0.2)
    var.llogis(shape=2.1, scale=0.2)
})

test_that("WeibullPH",{
    a <- 0.1; m <- 2
    b <- m^(-1/a)
    x <- c(-Inf, NaN, NA, -1, 0, 1, 2, Inf)
    expect_equal(dweibullPH(x, a, m), dweibull(x, a, b))
    expect_equal(dweibullPH(x, a, m), dweibull(x, a, b), log=TRUE)
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b))
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b), log.p=TRUE)
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b), lower.tail=FALSE)
    qq <- c(0, 0.5, 0.7, 1)
    expect_equal(qweibullPH(qq, a, m), qweibull(qq, a, b))
    expect_equal(hweibullPH(x, a, m), hweibull(x, a, b))
    expect_equal(hweibullPH(x, a, m), hweibull(x, a, b), log=TRUE)
    expect_equal(HweibullPH(x, a, m), Hweibull(x, a, b))
    expect_equal(HweibullPH(x, a, m), Hweibull(x, a, b), log=TRUE)
    set.seed(1); x1 <- rweibull(10, a, b)
    set.seed(1); x2 <- rweibullPH(10, a, m)
    expect_equal(x1, x2)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ rx, data = ovarian, dist="weibull")
    fitwp <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ rx, data = ovarian, dist="weibullPH")
    expect_equal(fitw$res["shape","est"], fitwp$res["shape","est"], tol=1e-06)
    expect_equal(fitw$res["scale","est"], fitwp$res["scale","est"]^(-1/fitwp$res["shape","est"]), tol=1e-05)
    expect_equal(coef(fitw)["rx"], -coef(fitwp)["rx"] / fitwp$res["shape","est"])
    
})
