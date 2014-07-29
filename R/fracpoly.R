### simplified version of gamlss::bfp with no shift or scale

bfp <- function (x, powers = c(1, 2)) 
{
    nobs <- length(x)
    npoly <- length(powers)
    X <- matrix(0, nrow = nobs, ncol = npoly)
    x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
    X[, 1] <- x1
    if (npoly >= 2) {
        for (i in 2:npoly) {
            if (powers[i] == powers[(i - 1)]) 
                x2 <- log(x) * x1
            else x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], 
                log(x))
            X[, i] <- x2
            x1 <- x2
        }
    }
    X
}

dbfp <- function (x, powers = c(1, 2)) 
{
    nobs <- length(x)
    npoly <- length(powers)
    X <- matrix(0, nrow = nobs, ncol = npoly)
    x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
    dx1 <- ifelse(powers[1] != rep(0, nobs), powers[1]*x^{powers[1]-1}, 1/x)
    X[, 1] <- dx1
    if (npoly >= 2) {
        for (i in 2:npoly) {
            if (powers[i] == powers[(i - 1)]) {
                x2 <- log(x) * x1
                dx2 <- log(x) * dx1 + 1/x * x1
            }
            else {
                x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], log(x))
                dx2 <- ifelse(powers[i] != rep(0, nobs), powers[i]*x^{powers[i]-1}, 1/x)
            }
            X[, i] <- dx2
            x1 <- x2
            dx1 <- dx2
        }
    }
    X
}
