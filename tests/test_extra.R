#library(devtools)
#load_all("~/work/flexsurv/flexsurv")
library(flexsurv) # package needs to be installed to test these, unless put devtools in Suggests
if (require("testthat") && !covr::in_covr()){
    test_dir("extra")
}
