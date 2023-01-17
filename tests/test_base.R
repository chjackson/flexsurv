test_partial <- TRUE

if (test_partial)
  options(
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE,
    warnPartialMatchDollar = TRUE
  )

if (require("testthat")){
    test_check("flexsurv")
}
