context("Generalized gamma distribution")

tol <- 1e-06

x <- c(-1,1,2,3,4); mu <- 2.2; sigma <- 1.6; Q <- 0; P <- 1
delta <- (Q^2 + 2*P)^{1/2}
s1 <- 2 / (Q^2 + 2*P + Q*delta); s2 <- 2 / (Q^2 + 2*P - Q*delta)

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
    expect_equal(qgengamma(p=0.25, mu=0, sigma=1.2, Q=0), qlnorm(p=0.25, meanlog=0, sdlog=1.2))
    expect_equal(qgengamma(p=0.25, mu=0.1, sigma=1.2, Q=1), qweibull(p=0.25, scale=exp(0.1), shape=1/1.2))
    expect_equal(qgengamma(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=1), mu=0, sigma=1, Q=1),
         c(0,0,0,1,2,3,4))
    expect_equal(qgengamma(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1, Q=-1), mu=0, sigma=1, Q=-1),
         c(0,0,0,1,2,3,4))
    expect_equal(qlnorm(pgengamma(q=c(-2,-1,0,1,2,3,4), mu=0, sigma=1.2, Q=0), meanlog=0, sdlog=1.2),
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

test_that("Gives errors",
          {
    expect_warning(dgengamma(1, 1, -2, 1), "Negative")
    expect_warning(pgengamma(1, 1, -2, 1), "Negative")
    expect_warning(qgengamma(0.1, 1, -2, 1), "Negative")
    expect_warning(rgengamma(1, 1, -2, 1), "Negative")
    
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
