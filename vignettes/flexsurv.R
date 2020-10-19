## ----echo=FALSE-------------------------------------------
options(width=60)
options(prompt="R> ")
library(knitr)
opts_chunk$set(fig.path="flexsurv-")
render_sweave()

## ----results='hide'---------------------------------------
library("flexsurv")

## ---------------------------------------------------------
head(bc, 2)

## ---------------------------------------------------------
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                   dist = "weibull")

## ---------------------------------------------------------
fs1

## ---------------------------------------------------------
survreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "weibull")

## ----results='hide',eval=FALSE----------------------------
#  library(eha)
#  aftreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "weibull")

## ---------------------------------------------------------
fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, 
                   dist = "gengamma")
fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ group + sigma(group), data = bc, 
                   dist = "gengamma")

## ----eval=FALSE-------------------------------------------
#  fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
#                     anc = list(sigma = ~ group), dist = "gengamma")

## ----surv,include=FALSE-----------------------------------
cols <- c("#E495A5", "#86B875", "#7DB0DD") # from colorspace::rainbow_hcl(3)
plot(fs1, col = cols[2], lwd.obs = 2, xlab = "Years", ylab = "Recurrence-free survival")
lines(fs2, col = cols[3], lty = 2)
lines(fs3, col = cols[3])
text(x=c(2,3.5,4), y=c(0.4, 0.55, 0.7), c("Poor","Medium","Good"))
legend("bottomleft", col = c("black", cols[2], cols[3], cols[3]), 
       lty = c(1, 1, 2, 1), bty = "n", lwd = rep(2, 4),
       c("Kaplan-Meier", "Weibull", "Generalized gamma (AFT)",
         "Generalized gamma (time-varying)"))

## ----haz,include=FALSE------------------------------------
plot(fs1, type = "hazard", col = cols[2], lwd.obs = 2, max.time=6, 
     xlab = "Years", ylab = "Hazard")
lines(fs2, type = "hazard", col = cols[3], lty = 2)
lines(fs3, type = "hazard", col = cols[3])
text(x=c(2,2,2), y=c(0.3, 0.13, 0.05), c("Poor","Medium","Good"))
legend("topright", col = c("black", cols[2], cols[3], cols[3]),
       lty = c(1, 1, 2, 1),  bty = "n", lwd = rep(2, 4),
       c("Kernel density estimate", "Weibull", "Gen. gamma (AFT)",
         "Gen. gamma (time-varying)"))

## ---------------------------------------------------------
median.weibull <- function(shape, scale) { 
    qweibull(0.5, shape = shape, scale = scale) 
}
summary(fs1, fn = median.weibull, t = 1, B = 10000)

## ----eval=FALSE-------------------------------------------
#  custom.llogis <- list(name = "llogis",  pars = c("shape", "scale"),
#                        location = "scale",
#                        transforms = c(log, log),
#                        inv.transforms = c(exp, exp),
#                        inits = function(t){ c(1, median(t)) })
#  fs4 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
#                     dist = custom.llogis)

## ----eval=FALSE-------------------------------------------
#  dmakeham3 <- function(x, shape1, shape2, scale, ...)  {
#      dmakeham(x, shape = c(shape1, shape2), scale = scale, ...)
#  }
#  pmakeham3 <- function(q, shape1, shape2, scale, ...)  {
#      pmakeham(q, shape = c(shape1, shape2), scale = scale, ...)
#  }

## ----eval=FALSE-------------------------------------------
#  dmakeham3 <- Vectorize(dmakeham3)
#  pmakeham3 <- Vectorize(pmakeham3)

## ----eval=FALSE-------------------------------------------
#  pmakeham3(c(0, 1, 1, Inf), 1, c(1, 1, 2, 1), 1)

## ----echo=FALSE-------------------------------------------
options(warn=-1)

## ---------------------------------------------------------
hweibullPH <- function(x, shape, scale = 1, log = FALSE){
    hweibull(x, shape = shape, scale = scale ^ {-1 / shape}, log = log)
}
HweibullPH <- function(x, shape, scale = 1, log = FALSE){
    Hweibull(x, shape = shape, scale = scale ^ {-1 / shape}, log = log)
}
custom.weibullPH <- list(name = "weibullPH", 
                         pars = c("shape", "scale"), location = "scale",
                         transforms = c(log, log),
                         inv.transforms = c(exp, exp),
                         inits = function(t){
                             c(1, median(t[t > 0]) / log(2))
                         })
fs6 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                   dist = custom.weibullPH)
fs6$res["scale", "est"] ^ {-1 / fs6$res["shape", "est"]}
- fs6$res["groupMedium", "est"] / fs6$res["shape", "est"]
- fs6$res["groupPoor", "est"] / fs6$res["shape", "est"]

## ----echo=FALSE-------------------------------------------
options(warn=0)

## ---------------------------------------------------------
sp1 <- flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k = 1, 
                      scale = "odds")

## ---------------------------------------------------------
sp2 <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group),
                      data = bc, k = 1, scale = "odds")

## ----splinehaz,include=FALSE------------------------------
plot(sp1, type = "hazard", col=cols[3], ylim = c(0, 0.5), xlab = "Years", 
     ylab = "Hazard")
lines(sp2, type = "hazard", col = cols[3], lty = 2)
lines(fs2, type = "hazard", col = cols[2])
text(x=c(2,2,2), y=c(0.3, 0.15, 0.05), c("Poor","Medium","Good"))
legend("topright", col = c("black",cols[c(3,3,2)]), 
       lty = c(1,1,2,1), lwd = rep(2,4),
       c("Kernel density estimate","Spline (proportional odds)",
         "Spline (time-varying)","Generalized gamma (time-varying)"))

## ---------------------------------------------------------
sp3 <- flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k = 1, 
                      scale = "hazard")
sp3$res[c("groupMedium", "groupPoor"), c("est", "se")]
cox3 <- coxph(Surv(recyrs, censrec) ~ group, data = bc)
coef(summary(cox3))[ , c("coef", "se(coef)")]

## ---------------------------------------------------------
sp4 <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group) + 
                        gamma2(group), data = bc, k = 1, scale = "hazard")

## ----eval=FALSE-------------------------------------------
#  flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(group) +
#                 gamma2(group) + treat + gamma1(treat),
#                 data = bc, k = 1, scale = "hazard")

## ---------------------------------------------------------
res <- t(sapply(list(fs1, fs2, fs3, sp1, sp2, sp3, sp4), 
                function(x)rbind(-2 * round(x$loglik,1), x$npars, 
                                 round(x$AIC,1))))
rownames(res) <- c("Weibull (fs1)", "Generalized gamma (fs2)",
                   "Generalized gamma (fs3)", 
                   "Spline (sp1)", "Spline (sp2)", "Spline (sp3)", 
                   "Spline (sp4)")
colnames(res) <- c("-2 log likelihood", "Parameters", "AIC")

## ----size="scriptsize"------------------------------------
res

## ---------------------------------------------------------
gamma <- sp1$res[c("gamma0", "gamma1", "gamma2"), "est"]
1 - psurvspline(5, gamma = gamma, knots = sp1$knots)

## ---------------------------------------------------------
pfn <- unroll.function(psurvspline, gamma = 0:2)
1 - pfn(5, gamma0 = gamma[1], gamma1 = gamma[2], gamma2 = gamma[3], 
        knots = sp1$knots)

## ---------------------------------------------------------
hsurvspline.lh <- function(x, gamma, knots){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow = 1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x, length = nret))
    x <- rep(x, length = nret)
    loghaz <- rowSums(basis(knots, log(x)) * gamma)
    exp(loghaz)
}

## ---------------------------------------------------------
hsurvspline.lh3 <- unroll.function(hsurvspline.lh, gamma = 0:2)

## ---------------------------------------------------------
custom.hsurvspline.lh3 <- list(
    name = "survspline.lh3",
    pars = c("gamma0", "gamma1", "gamma2"),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms = rep(c(identity), 3)
    )
dtime <- log(bc$recyrs)[bc$censrec == 1]
ak <- list(knots = quantile(dtime, c(0, 0.5, 1)))

## ----eval=FALSE-------------------------------------------
#  sp5 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, aux = ak,
#                     inits = c(0, 0, 0, 0, 0),
#                     dist = custom.hsurvspline.lh3,
#                     method = "L-BFGS-B", lower = c(-Inf, -Inf, -0.5),
#                     upper = c(Inf, Inf, 0.5),
#                     control = list(trace = 1, REPORT = 1))

