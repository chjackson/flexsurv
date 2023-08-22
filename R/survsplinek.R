##' Royston/Parmar spline survival distribution functions with one argument per parameter
##' 
##' Probability density, distribution, quantile, random generation, hazard, 
##' cumulative hazard, mean and restricted mean functions for the Royston/Parmar
##' spline model, with one argument per parameter.   For the equivalent functions with all parameters collected together in a single argument, see \code{\link{Survspline}}.
##'
##' These functions go up to 7 spline knots, or 9 parameters.  If you'd like higher-dimension versions, just submit an issue at \url{https://github.com/chjackson/flexsurv-dev/issues}.
##'
##' @aliases dsurvspline0 dsurvspline1 dsurvspline2 dsurvspline3 dsurvspline4 dsurvspline5 dsurvspline6 dsurvspline7     psurvspline0 psurvspline1 psurvspline2 psurvspline3 psurvspline4 psurvspline5 psurvspline6 psurvspline7     qsurvspline0 qsurvspline1 qsurvspline2 qsurvspline3 qsurvspline4 qsurvspline5 qsurvspline6 qsurvspline7   rsurvspline0 rsurvspline1 rsurvspline2 rsurvspline3 rsurvspline4 rsurvspline5 rsurvspline6 rsurvspline7  hsurvspline0 hsurvspline1 hsurvspline2 hsurvspline3 hsurvspline4 hsurvspline5 hsurvspline6 hsurvspline7  Hsurvspline0 Hsurvspline1 Hsurvspline2 Hsurvspline3 Hsurvspline4 Hsurvspline5 Hsurvspline6 Hsurvspline7  mean_survspline0 mean_survspline1 mean_survspline2 mean_survspline3 mean_survspline4 mean_survspline5 mean_survspline6 mean_survspline7    rmst_survspline0 rmst_survspline1 rmst_survspline2 rmst_survspline3 rmst_survspline4 rmst_survspline5 rmst_survspline6 rmst_survspline7      
##'
##'
##' @param gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8 Parameters describing the baseline spline function, as
##' described in \code{\link{flexsurvspline}}. 
##'
##' @inheritParams Survspline
##'
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##'
##' @name Survsplinek
NULL

##' @rdname Survsplinek
##' @export
mean_survspline0 <- function(gamma0, gamma1, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline1 <- function(gamma0, gamma1, gamma2, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline2 <- function(gamma0, gamma1, gamma2, gamma3, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline3 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline4 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline5 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline6 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
mean_survspline7 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots=c(-10, 10), scale="hazard", timescale="log"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale)
}




##' @rdname Survsplinek
##' @export
rmst_survspline0 <- function(t, gamma0, gamma1, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline1 <- function(t, gamma0, gamma1, gamma2, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline2 <- function(t, gamma0, gamma1, gamma2, gamma3, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline3 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline4 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline5 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline6 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, start=start)
}

##' @rdname Survsplinek
##' @export
rmst_survspline7 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots=c(-10, 10), scale="hazard", timescale="log", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, start=start)
}



##' @rdname Survsplinek
##' @export
dsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, log=log)
}

##' @rdname Survsplinek
##' @export
dsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, log=log)
}


##' @rdname Survsplinek
##' @export
psurvspline0 <- function(q, gamma0, gamma1, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline1 <- function(q, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline2 <- function(q, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline3 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline4 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline5 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline6 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
psurvspline7 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}




##' @rdname Survsplinek
##' @export
qsurvspline0 <- function(p, gamma0, gamma1, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline1 <- function(p, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline2 <- function(p, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline3 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline4 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline5 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline6 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}

##' @rdname Survsplinek
##' @export
qsurvspline7 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, lower.tail=lower.tail, log.p=log.p)
}




##' @rdname Survsplinek
##' @export
rsurvspline0 <- function(n, gamma0, gamma1, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline1 <- function(n, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline2 <- function(n, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline3 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline4 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline5 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline6 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
rsurvspline7 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale)
}


##' @rdname Survsplinek
##' @export
hsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
hsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale)
}


##' @rdname Survsplinek
##' @export
Hsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale)
}

##' @rdname Survsplinek
##' @export
Hsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale)
}
