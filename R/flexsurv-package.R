## unload the shared object
.onUnload <- function(libpath) {
    library.dynam.unload("flexsurv", libpath)
}
