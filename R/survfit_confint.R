### From the R "survival" package
## Therneau T (2020). _A Package for Survival Analysis in R_. R package
### version 3.2-3., <URL: https://CRAN.R-project.org/package=survival>.

survfit_confint <- function (p, se, logse = TRUE, conf.type, conf.int, selow, ulimit = TRUE) 
{
    zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
    if (missing(selow)) 
        scale <- 1
    else scale <- ifelse(selow == 0, 1, selow/se)
    if (!logse) 
        se <- ifelse(se == 0, 0, se/p)
    if (conf.type == "plain") {
        se2 <- se * p * zval
        if (ulimit) 
            list(lower = pmax(p - se2 * scale, 0), upper = pmin(p + 
                se2, 1))
        else list(lower = pmax(p - se2 * scale, 0), upper = p + 
            se2)
    }
    else if (conf.type == "log") {
        xx <- ifelse(p == 0, NA, p)
        se2 <- zval * se
        temp1 <- exp(log(xx) - se2 * scale)
        temp2 <- exp(log(xx) + se2)
        if (ulimit) 
            list(lower = temp1, upper = pmin(temp2, 1))
        else list(lower = temp1, upper = temp2)
    }
    else if (conf.type == "log-log") {
        xx <- ifelse(p == 0 | p == 1, NA, p)
        se2 <- zval * se/log(xx)
        temp1 <- exp(-exp(log(-log(xx)) - se2 * scale))
        temp2 <- exp(-exp(log(-log(xx)) + se2))
        list(lower = temp1, upper = temp2)
    }
    else if (conf.type == "logit") {
        xx <- ifelse(p == 0, NA, p)
        se2 <- zval * se * (1 + xx/(1 - xx))
        temp1 <- 1 - 1/(1 + exp(log(p/(1 - p)) - se2 * scale))
        temp2 <- 1 - 1/(1 + exp(log(p/(1 - p)) + se2))
        list(lower = temp1, upper = temp2)
    }
    else if (conf.type == "arcsin") {
        xx <- ifelse(p == 0, NA, p)
        se2 <- 0.5 * zval * se * sqrt(xx/(1 - xx))
        list(lower = (sin(pmax(0, asin(sqrt(xx)) - se2 * scale)))^2, 
            upper = (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
    }
    else stop("invalid conf.int type")
}
