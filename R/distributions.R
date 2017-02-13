sr.weib.inits <- function(t,aux){
    if (aux$counting){
        lt <- log(t[t>0])
###       c(1, exp(median(lt)) / log(2))
        c(1.64/var(lt), exp(mean(lt)+0.572)) # from survreg
    } else {
        aux$formula <- aux$forms[[1]]
        aux$forms <- NULL
        aux$dist <- "weibull"
        sr <- do.call(survreg, aux)
        sr2fswei(sr)
    }
}

sr.weibPH.inits <- function(aux){
    if (aux$counting){
        lt <- log(t[t>0])
        shape <- 1.64/var(lt)
        scale <- exp(mean(lt)+0.572)
        c(shape, scale^{-shape})
    } else {
        aux$formula <- aux$forms[[1]]
        aux$forms <- NULL
        aux$dist <- "weibull"
        sr <- do.call(survreg, aux)
        sr2fswei(sr, ph=TRUE)
    }
}

sr.exp.inits <- function(t,aux){
    if (aux$counting){
        1 / mean(t)
    } else {
        aux$formula <- aux$forms[[1]]
        aux$forms <- NULL
        aux$dist <- "exponential"
        sr <- do.call(survreg, aux)
        sr2fsexp(sr)
    }
}

sr.ln.inits <- function(t,aux){
    if (aux$counting){
        lt <- log(t[t>0])
        c(mean(lt), sd(lt))
    } else {
        aux$formula <- aux$forms[[1]]
        aux$forms <- NULL
        aux$dist <- "lognormal"
        sr <- do.call(survreg, aux)
        sr2fsln(sr)
    }
}

sr.llog.inits <- function(t,aux){
    if (aux$counting){
        scale <- median(t)
        shape <- 1 / log(quantile(t, 0.25)/scale, base=3)
        if (shape < 0) shape <- 1
        c(shape, scale)
    } else {
        aux$formula <- aux$forms[[1]]
        aux$forms <- NULL
        aux$dist <- "lognormal"
        sr <- do.call(survreg, aux)
        sr2fsllog(sr)
    }
}

## Convert parameters of survreg models to flexsurvreg
## parameterisation, for use as initial values

sr2fswei <- function(sr, ph=FALSE){
    scale <- exp(coef(sr)[1])
    beta.scale <- coef(sr)[-1]
    shape <- mean(1/sr$scale)
    beta.shape <- if (length(sr$scale)>1) log(sr$scale[1]/sr$scale[-1]) else numeric()
    if (ph) c(shape, scale^{-shape}, -beta.scale*shape, beta.shape)
    else c(shape, scale, beta.scale, beta.shape)
}

sr2fsexp <- function(sr){
    rate <- exp(-coef(sr)[1])
    beta <- -coef(sr)[-1]
    c(rate, beta)
}

sr2fsln <- function(sr){
    meanlog <- coef(sr)[1]
    sdlog <- sr$scale
    beta <- coef(sr)[-1]
    c(meanlog, sdlog, beta)
}

sr2fsllog <- function(sr){
    shape <- 1/sr$scale
    scale <- exp(coef(sr)[1])
    beta <- coef(sr)[-1]
    c(shape, scale, beta)
}

##' @export
flexsurv.dists <-
    list(genf = list(
             name="genf",
             pars=c("mu","sigma","Q","P"),
             location="mu",
             transforms=c(identity, log, identity, log),
             inv.transforms=c(identity, exp, identity, exp),
             inits=function(t){
                 lt <- log(t[t>0])
                 c(mean(lt), sd(lt), 0, 1)
             }
        ),
         genf.orig = list(
             name="genf.orig",
             pars=c("mu","sigma","s1","s2"),
             location="mu",
             transforms=c(identity, log, log, log),
             inv.transforms=c(identity, exp, exp, exp),
             inits=function(t){
                 lt <- log(t[t>0])
                 c(mean(lt), sd(lt), 1, 1)
             }
             ),
         gengamma = list(
             name="gengamma",
             pars=c("mu","sigma","Q"),
             location="mu",
             transforms=c(identity, log, identity),
             inv.transforms=c(identity, exp, identity),
             inits=function(t){
                 lt <- log(t[t>0])
                 c(mean(lt), sd(lt), 0)
             }
             ),
         gengamma.orig = list(
             name="gengamma.orig",
             pars=c("shape","scale","k"),
             location="scale",
             transforms=c(log, log, log),
             inv.transforms=c(exp, exp, exp),
             inits=function(t){c(1, mean(t), 1)}
             ),
         exp = list(
             name="exp",
             pars=c("rate"),
             location="rate",
             transforms=c(log),
             inv.transforms=c(exp),
             inits=sr.exp.inits
             ),
         weibull = list(
             name = "weibull.quiet",
             pars=c("shape","scale"),
             location="scale",
             transforms=c(log, log),
             inv.transforms=c(exp, exp),
             inits=sr.weib.inits
             ),
         weibullPH = list(
             name="weibullPH",
             pars=c("shape","scale"),
             location="scale",
             transforms=c(log, log),
             inv.transforms=c(exp, exp),
             inits = sr.weibPH.inits
             ),
         lnorm = list(
             name="lnorm",
             pars=c("meanlog","sdlog"),
             location="meanlog",
             transforms=c(identity, log),
             inv.transforms=c(identity, exp),
             inits=sr.ln.inits
             ),
         gamma = list(
             name="gamma",
             pars=c("shape","rate"),
             location="rate",
             transforms=c(log, log),
             inv.transforms=c(exp, exp),
             inits=function(t){
                 m=mean(t); v=var(t);
                 c(m^2/v, m/v)
             }
             ),
         gompertz = list(
             name="gompertz",
             pars=c("shape","rate"),
             location="rate",
             transforms=c(identity, log),
             inv.transforms=c(identity, exp),
             inits=function(t){c(0.001,1 / mean(t))}
             ),
         llogis = list(
             name="llogis",
             pars=c("shape","scale"),
             location="scale",
             transforms=c(log, log),
             inv.transforms=c(exp, exp),
             inits=sr.llog.inits
             )
         )
flexsurv.dists$exponential <- flexsurv.dists$exp
flexsurv.dists$lognormal <- flexsurv.dists$lnorm
