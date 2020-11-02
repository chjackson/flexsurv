## No covariates 

n <- 1000
p <- 0.5 
set.seed(1)
death <- rbinom(n, 1, p)
t <- numeric(n)
t[death==0] <- rgamma(sum(death==0), 1, 3.2)
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 3)
status <- 1 - cens
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
t[cens] <- 3
dat <- data.frame(t, status, event)
dat$evname <- c("cure", "death")[dat$event]
dat$evnamef <- factor(dat$evname)

flexsurvmix(Surv(t, status) ~ 1, data=dat, 
            event=event, dists=c("gamma","gamma"), method="em", 
            em.control = list(var.method="direct"))
flexsurvmix(Surv(t, status) ~ 1, data=dat, 
            event=event, dists=c("gamma","gamma"), method="em", 
            em.control = list(var.method="louis"))

# complete data
flexsurvmix(Surv(t, status) ~ 1, data=dat[!is.na(dat$event),], 
            event=event, dists=c("gamma","gamma"), method="em", 
            em.control = list(var.method="direct"))
flexsurvmix(Surv(t, status) ~ 1, data=dat[!is.na(dat$event),], 
            event=event, dists=c("gamma","gamma"), method="em", 
            em.control = list(var.method="louis"))


# With covariates

n <- 1000
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
p <- plogis(qlogis(0.5) + 2*x - 3*y)
death <- rbinom(n, 1, p)
t <- numeric(n)
t[death==0] <- rgamma(sum(death==0), 1, 3.2)
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
cens <- as.numeric(t > 3)
status <- 1 - cens
t[cens] <- 3
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x, y)
dat$evname <- c("cure", "death")[dat$event]

flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                 dists=c("gamma","gamma"), pformula = ~ x + y, 
                 method="em")
flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
                 dists=c("gamma","gamma"), pformula = ~ x + y, 
                 method="em",
                 em.control = list(var.method="louis"))  

datc <- dat[!is.na(dat$evname),]
datc$evname <- factor(datc$evname)
flexsurvmix(Surv(t, status) ~ 1, data=datc, event=evname, 
            dists=c("gamma","gamma"), pformula = ~ y, 
            method="em")
flexsurvmix(Surv(t, status) ~ 1, data=datc, event=evname, 
            dists=c("gamma","gamma"), pformula = ~ y, 
            method="em",
            em.control = list(var.method="louis"))  


## 3 categories, and covariates
n <- 1000
set.seed(1)
x <- rnorm(n)
y <- rbinom(n, 1,  0.5)
alpha1 <- qlogis(0.5) + 2*x - 3*y
alpha2 <- qlogis(0.5) + 1*x - 2*y
alpha <- cbind(0, alpha1, alpha2)
p <- exp(alpha)/(rowSums(exp(alpha)))
death <- numeric(n)
for (i in 1:n)
  death[i] <- sample(x=0:2, size=1, prob=p[i,])
t <- numeric(n)
t[death==0] <- rgamma(sum(death==0), 1, 3.2)
t[death==1] <- rgamma(sum(death==1), 2.5, 1.2)
t[death==2] <- rgamma(sum(death==2), 3.5, 2.2)
cens <- as.numeric(t > 3)
status <- 1 - cens
t[cens] <- 3
event <- ifelse(cens, NA, death+1) # 1 is cure, 2 is death
dat <- data.frame(t, status, event, x, y)
dat$evname <- c("cure", "icu","death")[dat$event]

## complete data 
datc <- dat[!is.na(dat$evname),]
datc$evname <- factor(datc$evname)

flexsurvmix(Surv(t, status) ~ 1, data=datc, event=evname, 
            dists=c("gamma","gamma","gamma"), pformula = ~ y, 
            method="em")

flexsurvmix(Surv(t, status) ~ 1, data=datc, event=evname, 
            dists=c("gamma","gamma","gamma"), pformula = ~ y, 
            method="em",
            em.control = list(var.method="louis")) 


## incomplete data 
flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
            dists=c("gamma","gamma","gamma"), pformula = ~ y, 
            method="em")

flexsurvmix(Surv(t, status) ~ 1, data=dat, event=evname, 
            dists=c("gamma","gamma","gamma"), pformula = ~ y, 
            method="em",
            em.control = list(var.method="louis")) 
