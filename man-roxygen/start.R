##' @param start Optional left-truncation time or times.  The returned
##' survival, hazard or cumulative hazard will be conditioned on survival up to
##' this time.   Predicted times returned with \code{"rmst"}, \code{"mean"}, \code{"median"} or \code{"quantile"}
##' will be times since time zero, not times since the \code{start} time.
##'
##' A vector of the same length as \code{t} can be supplied to allow different
##' truncation times for each prediction time, though this doesn't make sense
##' in the usual case where this function is used to calculate a predicted
##' trajectory for a single individual.  This is why the default \code{start}
##' time was changed for version 0.4 of \pkg{flexsurv} - this was previously a
##' vector of the start times observed in the data.
