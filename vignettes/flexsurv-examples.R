### R code from vignette source 'flexsurv-examples.Rnw'

###################################################
### code chunk number 1: flexsurv-examples.Rnw:39-59
###################################################
library(flexsurv)

hgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

HgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * Hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

custom.gengammaPH <- list(name="gengammaPH", 
                          pars=c("dummy","mu","sigma","Q"), location="dummy",
                          transforms=c(log, identity, log, identity),
                          inv.transforms=c(exp, identity, exp, identity),
                          inits=function(t){
                              lt <- log(t[t>0])
                              c(1, mean(lt), sd(lt), 0)
                          })
fs7 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, 
                   dist=custom.gengammaPH, fixedpars=1)


###################################################
### code chunk number 2: flexsurv-examples.Rnw:84-127
###################################################

fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ group + sigma(group), 
                   data=bc, dist="gengamma")
B <- 5000
t <- seq(0.1, 8, by=0.1)

hrAFT.est <-
    summary(fs2, t=t, type="hazard",
            newdata=data.frame(group="Medium"),ci=FALSE)[[1]][,"est"] /
    summary(fs2, t=t, type="hazard",
            newdata=data.frame(group="Good"),ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(fs2, B=B, newdata=data.frame(group=c("Good","Medium")))
hrAFT <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){ 
    haz.medium.rep <- do.call(hgengamma, c(list(t[i]), as.data.frame(pars[[2]])))
    haz.good.rep <- do.call(hgengamma, c(list(t[i]), as.data.frame(pars[[1]])))
    hrAFT[,i] <- haz.medium.rep / haz.good.rep
}
hrAFT <- apply(hrAFT, 2, quantile, c(0.025, 0.975))

hrPH.est <-
    summary(fs7, t=t, type="hazard",
            newdata=data.frame(group="Medium"),ci=FALSE)[[1]][,"est"] /
    summary(fs7, t=t, type="hazard",
            newdata=data.frame(group="Good"),ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(fs7, B=B, newdata=data.frame(group=c("Good","Medium")))
hrPH <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){ 
    haz.medium.rep <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[2]])))
    haz.good.rep <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[1]])))
    hrPH[,i] <- haz.medium.rep / haz.good.rep
}
hrPH <- apply(hrPH, 2, quantile, c(0.025, 0.975))

plot(t, hrAFT[1,], type="l", ylim=c(0, 10), col="red", xlab="Years",
     ylab="Hazard ratio (Medium / Good)", lwd=1, lty=2)
lines(t, hrAFT[2,], col="red", lwd=1, lty=2)
lines(t, hrPH[1,], col="darkgray", lwd=1, lty=2)
lines(t, hrPH[2,], col="darkgray", lwd=1, lty=2)
lines(t, hrAFT.est, col="red", lwd=2)
lines(t, hrPH.est, col="darkgray", lwd=2)
legend("topright", lwd=c(2,2), col=c("red","darkgray"), bty="n",
       c("Generalized gamma: standard AFT", "Generalized gamma: proportional hazards"))


###################################################
### code chunk number 3: flexsurv-examples.Rnw:150-157
###################################################
mean.gengamma <- function(mu, sigma, Q, horizon=100, ...){
    surv <- function(t, ...) {  1 - pgengamma(q=t, mu=mu, sigma=sigma, Q=Q, ...) }
    integrate(surv, 0, horizon, ...)$value
}
summary(fs2, newdata=list(group="Good"), t=1, fn=mean.gengamma)
summary(fs2, newdata=list(group="Medium"), t=1, fn=mean.gengamma)
summary(fs2, newdata=list(group="Poor"), t=1, fn=mean.gengamma)


###################################################
### code chunk number 4: flexsurv-examples.Rnw:160-164
###################################################
median.gengamma <- function(mu, sigma, Q) { 
    qgengamma(0.5, mu=mu, sigma=sigma, Q=Q) 
}
summary(fs2, newdata=list(group="Good"), t=1, fn=median.gengamma)


###################################################
### code chunk number 5: flexsurv-examples.Rnw:189-209
###################################################
if (require("TH.data")){

GBSG2 <- transform(GBSG2,
                   X1a=(age/50)^-2,
                   X1b=(age/50)^-0.5,
                   X4=tgrade %in% c("II","III"),
                   X5=exp(-0.12*pnodes),
                   X6=(progrec+1)^0.5
                   )
(progc <- coxph(Surv(time, cens) ~ horTh + X1a + X1b + X4 + 
                  X5 + X6, data=GBSG2))
(prog3 <- flexsurvspline(Surv(time, cens) ~ horTh + X1a + X1b + X4 + 
                           X5 + X6, k=3, data=GBSG2))
predc <- predict(progc, type="lp")
progc <- cut(predc, quantile(predc, 0:3/3))
predf <- model.matrix(prog3) %*% prog3$res[-(1:5),"est"]
progf <- cut(predf, quantile(predf, 0:3/3))
table(progc, progf)

}


###################################################
### code chunk number 6: flexsurv-examples.Rnw:229-233
###################################################
set.seed(1)
nsim <- 10000
onsetday <- runif(nsim, 0, 30) 
deathday <- onsetday + rgamma(nsim, shape=1.5, rate=1/10)


###################################################
### code chunk number 7: flexsurv-examples.Rnw:243-247
###################################################
datt <- data.frame(delay = deathday - onsetday,
                   event = rep(1, nsim),
                   rtrunc = 40 - onsetday)
datt <- datt[datt$delay < datt$rtrunc, ] 


###################################################
### code chunk number 8: flexsurv-examples.Rnw:257-260
###################################################
fitt <- flexsurvreg(Surv(delay, event) ~ 1, data=datt, rtrunc = rtrunc, dist="gamma")
fitt
summary(fitt, t=1, fn = mean_gamma)


