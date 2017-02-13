context("Generalized gamma (original)")

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


test_that("Gives errors",
          {
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative shape")
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative scale")
    expect_warning(dgengamma.orig(c(1,2,3,4), shape=c(-1.2, 1), scale=1.3, k=1), "Negative shape")
    expect_warning(pgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative")
    expect_warning(qgengamma.orig(c(1,2,3,4), shape=-1.2, scale=-1.3, k=-1), "Negative")
    expect_warning(rgengamma.orig(3, shape=-1.2, scale=-1.3, k=-1), "Negative")
    
})
