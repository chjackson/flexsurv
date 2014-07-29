
## SIMULATION TESTS. Simulate data from r and refit model to test
set.seed(12082012)
if (0) {
    sim <- rgenf(3000, 1.5, 1, -0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="genf", control=list(trace=1,REPORT=1))
    fit$res # OK

    sim <- rgengamma(3000, 1.5, 1, -0.4)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma", control=list(trace=1,REPORT=1))
    fit$res # OK

    sim <- rgenf.orig(3000, 1.5, 1, 0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="genf.orig", control=list(trace=1,REPORT=1,maxit=10000))
    fit$res # OK

    sim <- rgengamma.orig(3000, 1.5, 1, 0.4)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma.orig", control=list(trace=1,REPORT=1))
    fit$res # OK

    xg <- rgompertz(1000, 0.12, 4); hist(xg)
    flexsurvreg(Surv(xg, rep(1,1000)) ~ 1, dist="gompertz") ## OK. robust to starting values

    if (is.element("eha", installed.packages()[,1])) {
        library(eha)
        foo <- phreg(Surv(xg, rep(1,1000)) ~ 1, dist="gompertz") ## OK - names of parameters other way round, see dgompertz help
        xl <- rllogis(1000, 5.4, 0.1); hist(xl)
        flexsurvreg(Surv(xl, rep(1,1000)) ~ 1, dist=custom.llogis) ## OK. robust to starting values
    }
}
