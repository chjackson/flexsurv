## ----echo=FALSE-------------------------------------------
options(width=60)
options(prompt="R> ")
library(knitr)
opts_chunk$set(fig.path="flexsurv-")
render_sweave()

## ---------------------------------------------------------
library(flexsurv)
bosms3[18:22, ]

## ---------------------------------------------------------
crexp <- flexsurvreg(Surv(years, status) ~ trans, data = bosms3, 
                     dist = "exp")
crwei <- flexsurvreg(Surv(years, status) ~ trans + shape(trans), 
                     data = bosms3, dist = "weibull")
cfwei <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + shape(trans), 
                     data = bosms3, dist = "weibull")

## ---------------------------------------------------------
crcox <- coxph(Surv(years, status) ~ strata(trans), data = bosms3)
cfcox <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = bosms3)

## ---------------------------------------------------------
mod_nobos_bos <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==1),  
                             data = bosms3, dist = "weibull")
mod_nobos_death <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==2),  
                               data = bosms3, dist = "weibull")
mod_bos_death <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==3), 
                             data = bosms3, dist = "weibull")

## ---------------------------------------------------------
tmat <- rbind(c(NA, 1, 2), c(NA, NA, 3), c(NA, NA, NA))
crfs <- fmsm(mod_nobos_bos, mod_nobos_death, mod_bos_death, trans = tmat)

## ---------------------------------------------------------
require("mstate")
mrcox <- msfit(crcox, trans = tmat)
mfcox <- msfit(cfcox, trans = tmat)

## ---------------------------------------------------------
tgrid <- seq(0, 14, by = 0.1)
mrwei <- msfit.flexsurvreg(crwei, t = tgrid, trans = tmat)
mrexp <- msfit.flexsurvreg(crexp, t = tgrid, trans = tmat)
mfwei <- msfit.flexsurvreg(cfwei, t = tgrid, trans = tmat)

## ----cumhaz,include=FALSE---------------------------------
cols <- c("black", "#E495A5", "#86B875", "#7DB0DD") # colorspace::rainbow_hcl(3)
plot(mrcox, xlab = "Years after baseline", lwd = 3, xlim = c(0, 14), cols = cols[1:3])
for (i in 1:3){
    lines(tgrid, mrexp$Haz$Haz[mrexp$Haz$trans == i], col = cols[i], lty = 2, lwd = 2)
    lines(tgrid, mrwei$Haz$Haz[mrwei$Haz$trans == i], col = cols[i], lty = 3, lwd = 2)
}
lines(mfcox$Haz$time[mfcox$Haz$trans == 3], mfcox$Haz$Haz[mfcox$Haz$trans == 3],
      type = "s", col = cols[4], lty = 1, lwd = 2)
lines(tgrid, mfwei$Haz$Haz[mfwei$Haz$trans == 3], col = cols[4], lty = 3, lwd = 2)
legend("topleft", inset = c(0, 0.2), lwd = 2, col = cols[4], 
       c("2 -> 3 (clock-forward)"), bty = "n")
legend("topleft", inset = c(0, 0.3), c("Non-parametric", "Exponential", "Weibull"),
       lty = c(1, 2, 3), lwd = c(3, 2, 2), bty = "n")

## ---------------------------------------------------------
pmatrix.fs(cfwei, t = c(5, 10), trans = tmat)

## ----eval=FALSE-------------------------------------------
#  pmatrix.simfs(crwei, trans = tmat, t = 5)
#  pmatrix.simfs(crwei, trans = tmat, t = 10)

## ---------------------------------------------------------
ptc <- probtrans(mfcox, predt = 0, direction = "forward")[[1]]
round(ptc[c(165, 193),], 3)

## ---------------------------------------------------------
ptw <- probtrans(mfwei, predt = 0, direction = "forward")[[1]]
round(ptw[ptw$time %in% c(5, 10),], 3)

## ----eval=FALSE-------------------------------------------
#  mssample(mrcox$Haz, trans = tmat, clock = "reset", M = 1000,
#           tvec = c(5, 10))
#  mssample(mrwei$Haz, trans = tmat, clock = "reset", M = 1000,
#           tvec = c(5, 10))

## ---------------------------------------------------------
tmat_nobos <- rbind("No BOS"=c(NA,1,2),
"BOS"=c(NA,NA,NA),"Death"=c(NA,NA,NA))
crfs_nobos <- fmsm(mod_nobos_bos, mod_nobos_death, trans = tmat_nobos)
pmatrix.fs(crfs_nobos, from=1, t=100000)["No BOS",c("BOS","Death")]

## ---------------------------------------------------------
modexp_nobos_bos <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==1),  
                                data = bosms3, dist = "exponential")
modexp_nobos_death <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==2),
                                  data = bosms3, dist = "exponential")
crfs_nobos <- fmsm(modexp_nobos_bos, modexp_nobos_death, trans=tmat_nobos)
pmatrix.fs(crfs_nobos, from=1, t=100000)["No BOS",c("BOS","Death")]
rate12 <- modexp_nobos_bos$res["rate","est"]
rate13 <- modexp_nobos_death$res["rate","est"]
rate12 / (rate12  + rate13)

## ---------------------------------------------------------
crfs_nobos <- fmsm(mod_nobos_bos, mod_nobos_death, trans = tmat_nobos)
simfinal_fmsm(crfs_nobos, probs = c(0.25, 0.5, 0.75), M=10000)

## ---------------------------------------------------------
simfinal_fmsm(crfs, probs = c(0.25, 0.5, 0.75), M=10000)

## ---------------------------------------------------------
bosms3[bosms3$id  %in% c(4,5),]
bosmx3 <- bosms3
bosmx3$Tstart <- bosmx3$Tstop <- bosmx3$trans <- NULL
bosmx3 <- bosmx3[!(bosmx3$status==0 & duplicated(paste(bosmx3$id, bosmx3$from))),]
bosmx3$event <- ifelse(bosmx3$status==0, NA, bosmx3$to)
bosmx3$event <- factor(bosmx3$event, labels=c("BOS","Death"))
bosmx3$to <- NULL
bosmx3[bosmx3$id %in% c(4,5),]

## ---------------------------------------------------------
bosfs <- flexsurvmix(Surv(years, status) ~ 1, event=event, 
                     data=bosmx3[bosmx3$from==1,],
            dists = c(BOS="weibull", Death="exponential"))
bosfs

## ---------------------------------------------------------
set.seed(1)
bosmx3$x <- rnorm(nrow(bosmx3),0,1)
bosfsp <- flexsurvmix(Surv(years, status) ~ 1, event=event,
                      data=bosmx3[bosmx3$from==1,],
            dists = c(BOS="weibull", Death="exponential"),
            pformula = ~ x)
bosfsp

## ---------------------------------------------------------
bosfst <- flexsurvmix(Surv(years, status) ~ x, event=event, 
                      data=bosmx3[bosmx3$from==1,],
            dists = c(BOS="weibull", Death="exponential"))

## ---------------------------------------------------------
flexsurvmix(list(Surv(years, status) ~ x, 
                 Surv(years, status) ~ 1),
                 event=event, data=bosmx3[bosmx3$from==1,],
            dists = c(BOS="weibull", Death="exponential"))

## ---------------------------------------------------------
flexsurvmix(Surv(years, status) ~ 1,
            event=event, data=bosmx3[bosmx3$from==1,],
            dists = c(BOS="weibull", Death="exponential"),
            anc = list(BOS=list(shape=~x), 
                       Death=list(rate=~1)))

## ---------------------------------------------------------
nd <-  data.frame(x=c(0,1))
probs_flexsurvmix(bosfsp, newdata=nd, B=50)

## ---------------------------------------------------------
mean_flexsurvmix(bosfst, newdata=nd)
quantile_flexsurvmix(bosfst, B=50, newdata=nd, probs=c(0.25, 0.5, 0.75))

## ---------------------------------------------------------
p_flexsurvmix(bosfst, t=c(5, 10), B=10, startname="No BOS")

## ---------------------------------------------------------
bdat <- bosmx3[bosmx3$from==2,]
bdat$event <- factor(bdat$event)
bosfst_bos <- flexsurvmix(Surv(years, status) ~ x, event=event, data=bdat,
                          dists = c(Death="exponential"))
bosfmix <- fmixmsm("No BOS"=bosfst, "BOS"=bosfst_bos)

## ---------------------------------------------------------
ppath_fmixmsm(bosfmix, B=20)
ppath_fmixmsm(bosfmix, B=20, final=TRUE)

## ---------------------------------------------------------
meanfinal_fmixmsm(bosfmix, B=10)
qfinal_fmixmsm(bosfmix, B=10)
meanfinal_fmixmsm(bosfmix, B=10, final=TRUE)

## ----fig.height=3-----------------------------------------
aj <- ajfit_flexsurvmix(bosfs, start="No BOS")
bosfs_bos <- flexsurvmix(Surv(years, status) ~ 1, event=event, data=bdat, 
                          dists = c(Death="exponential"))
aj2 <- ajfit_flexsurvmix(bosfs_bos, start="BOS")
library(ggplot2)
ggplot(aj, aes(x=time, y=val, lty=model, col=state)) + 
  geom_line() +  
  xlab("Years after transplant") + ylab("Probability of having moved to state")
ggplot(aj2, aes(x=time, y=val, lty=model, col=state)) + 
  xlab("Years after transplant") + ylab("Probability of having moved to state") +
  geom_line()

## ----message=FALSE----------------------------------------
library(dplyr)

aj3 <- ajfit_fmsm(crfs_nobos) %>%
  filter(model=="Parametric") %>%
  mutate(model = "Cause specific hazards") %>%
  mutate(state = factor(state, labels=c("No BOS","BOS","Death"))) %>%
  full_join(aj) 

ggplot(aj3, aes(x=time, y=val, lty=model, col=state)) + 
  xlab("Years after transplant") + ylab("Probability of having moved to state") +
  geom_line()

