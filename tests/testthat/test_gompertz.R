context("Testing Gompertz")

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

test_that("Gompertz with chance of living forever",{
    shape <- -0.6; rate <- 1.8
    x <- c(0.8, 0.9, 0.97, 0.99)
    expect_equal(qgompertz(x, shape=shape, rate=rate), c(1.28150707286845, 2.4316450975351, Inf, Inf))
                                        # qgeneric(pgompertz, p=x, shape=shape, rate=rate) # won't work - needs smoothness
    expect_equal(pgompertz(Inf, shape=shape, rate=rate, lower.tail = F), exp(rate/shape))
    expect_equal(
      pgompertz(Inf, shape=shape, rate=rate, lower.tail = F),
      pgompertz(9999999, shape=shape, rate=rate, lower.tail = F)
    )
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
