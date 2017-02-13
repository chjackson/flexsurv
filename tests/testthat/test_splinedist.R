context("Spline distribution functions")

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

    qsurvspline(0.2,g,knots=k)

    expect_equal(psurvspline(qsurvspline(0.2,g,knots=k), g, knots=k), 0.2)
    expect_equal(qsurvspline(psurvspline(0.2,g,knots=k), g, knots=k), 0.2)

    ## Vectorised q and p functions, through qgeneric
    kvec <- rbind(k, k+1)
    gvec <- rbind(g, g*1.1)
    expect_equal(psurvspline(qsurvspline(0.2,gvec,knots=kvec), gvec, knots=kvec), c(0.2,0.2))
    expect_equal(qsurvspline(psurvspline(0.2,gvec,knots=kvec), gvec, knots=kvec), c(0.2,0.2))

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
