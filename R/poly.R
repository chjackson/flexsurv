dpolyhaz <- function(base, ncause, pars) {
    ## which pars vary between causes?
    ## Different for different distributions.
    ## Exponential: sum of hazards gives another exponential
    ## Weibull: shape has to change. same shape reduces to Weibull
    ## What others? See Louzada-Neto. GG, loglogistic, lognormal.

    ## What about covariates?  For our app, has prop haz for one cause and identical for other.
    ## see Louzada-Neto. various ways of doing it.

    ## can people just use a custom dist for this? e.g. for biweibull.
    ## yes if just one covariate effect parameter.
}
