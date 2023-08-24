gamma <- seq(1, 1.2, length.out=9)
names(gamma) <- 0:8

test_that("mean_survsplinek",{
  expect_equal(mean_survspline0(gamma["0"], gamma["1"]), mean_survspline(gamma[1:2]))
  expect_equal(mean_survspline1(gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,0)),
               mean_survspline(gamma[1:3], knots = c(-10,0,0)))
  expect_equal(mean_survspline2(gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,0,0)), 
               mean_survspline(gamma[1:4], knots=c(-10,0,0,0)))
  expect_equal(mean_survspline3(gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], knots=c(-10,0,0,0,0)),
               mean_survspline(gamma[1:5], knots=c(-10,0,0,0,0)))
  expect_equal(mean_survspline4(gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], knots=c(-10,0,0,0,0,0)),
               mean_survspline(gamma[1:6], knots=c(-10,0,0,0,0,0)))
  expect_equal(mean_survspline5(gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], knots=c(-10,0,0,0,0,0,0)),
               mean_survspline(gamma[1:7], knots=c(-10,0,0,0,0,0,0)))
  expect_equal(mean_survspline6(gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], knots=c(-10,0,0,0,0,0,0,0)),
               mean_survspline(gamma[1:8], knots=c(-10,0,0,0,0,0,0,0)))
  expect_equal(mean_survspline7(gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], knots=c(-10,0,0,0,0,0,0,0,0)),
               mean_survspline(gamma[1:9], knots=c(-10,0,0,0,0,0,0,0,0)))
})

test_that("rmst_survsplinek",{
  t <- 1
  expect_equal(rmst_survspline0(t=t, gamma["0"], gamma["1"]), rmst_survspline(t=t, gamma[1:2]))
  expect_equal(rmst_survspline1(t=t, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               rmst_survspline(t=t, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(rmst_survspline2(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               rmst_survspline(t=t, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(rmst_survspline3(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               rmst_survspline(t=t, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(rmst_survspline4(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               rmst_survspline(t=t, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(rmst_survspline5(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               rmst_survspline(t=t, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(rmst_survspline6(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               rmst_survspline(t=t, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(rmst_survspline7(t=t, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               rmst_survspline(t=t, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})

test_that("dsurvsplinek",{
  x <- 1
  expect_equal(dsurvspline0(x=x, gamma["0"], gamma["1"], knots=c(-10,10)), dsurvspline(x=x, gamma[1:2], knots=c(-10,10)))
  expect_equal(dsurvspline1(x=x, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               dsurvspline(x=x, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(dsurvspline2(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               dsurvspline(x=x, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(dsurvspline3(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               dsurvspline(x=x, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(dsurvspline4(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               dsurvspline(x=x, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(dsurvspline5(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               dsurvspline(x=x, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(dsurvspline6(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               dsurvspline(x=x, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(dsurvspline7(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               dsurvspline(x=x, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})


test_that("psurvsplinek",{
  x <- 1
  expect_equal(psurvspline0(q=x, gamma["0"], gamma["1"], knots=c(-10,10)), psurvspline(q=x, gamma[1:2], knots=c(-10,10)))
  expect_equal(psurvspline1(q=x, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               psurvspline(q=x, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(psurvspline2(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               psurvspline(q=x, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(psurvspline3(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               psurvspline(q=x, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(psurvspline4(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               psurvspline(q=x, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(psurvspline5(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               psurvspline(q=x, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(psurvspline6(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               psurvspline(q=x, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(psurvspline7(q=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               psurvspline(q=x, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})

if (covr::in_covr()){ # these are slow
  
test_that("qsurvsplinek",{
  p <- 0.4
  expect_equal(qsurvspline0(p=p, gamma["0"], gamma["1"], knots=c(-10,10)), qsurvspline(p=p, gamma[1:2], knots=c(-10,10)))
  expect_equal(qsurvspline1(p=p, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               qsurvspline(p=p, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(qsurvspline2(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               qsurvspline(p=p, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(qsurvspline3(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               qsurvspline(p=p, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(qsurvspline4(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               qsurvspline(p=p, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(qsurvspline5(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               qsurvspline(p=p, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(qsurvspline6(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               qsurvspline(p=p, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(qsurvspline7(p=p, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               qsurvspline(p=p, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})

test_that("rsurvsplinek",{
  expect_equal({set.seed(1); rsurvspline0(n=1, gamma["0"], gamma["1"], knots=c(-10,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:2], knots=c(-10,10))})
  expect_equal({set.seed(1); rsurvspline1(n=1, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:3], knots = c(-10,0,10))})
  expect_equal({set.seed(1); rsurvspline2(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10))}, 
               {set.seed(1); rsurvspline(n=1, gamma[1:4], knots=c(-10,0,1,10))})
  expect_equal({set.seed(1); rsurvspline3(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:5], knots=c(-10,0,1,2,10))})
  expect_equal({set.seed(1); rsurvspline4(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:6], knots=c(-10,0,1,2,3,10))})
  expect_equal({set.seed(1); rsurvspline5(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:7], knots=c(-10,0,1,2,3,4,10))})
  expect_equal({set.seed(1); rsurvspline6(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10))})
  expect_equal({set.seed(1); rsurvspline7(n=1, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10))},
               {set.seed(1); rsurvspline(n=1, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10))})
})

}

test_that("hsurvsplinek",{
  x <- 1
  expect_equal(hsurvspline0(x=x, gamma["0"], gamma["1"], knots=c(-10,10)), hsurvspline(x=x, gamma[1:2], knots=c(-10,10)))
  expect_equal(hsurvspline1(x=x, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               hsurvspline(x=x, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(hsurvspline2(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               hsurvspline(x=x, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(hsurvspline3(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               hsurvspline(x=x, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(hsurvspline4(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               hsurvspline(x=x, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(hsurvspline5(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               hsurvspline(x=x, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(hsurvspline6(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               hsurvspline(x=x, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(hsurvspline7(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               hsurvspline(x=x, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})

test_that("Hsurvsplinek",{
  x <- 1
  expect_equal(Hsurvspline0(x=x, gamma["0"], gamma["1"], knots=c(-10,10)), Hsurvspline(x=x, gamma[1:2], knots=c(-10,10)))
  expect_equal(Hsurvspline1(x=x, gamma["0"], gamma["1"], gamma["2"], knots = c(-10,0,10)),
               Hsurvspline(x=x, gamma[1:3], knots = c(-10,0,10)))
  expect_equal(Hsurvspline2(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], knots=c(-10,0,1,10)), 
               Hsurvspline(x=x, gamma[1:4], knots=c(-10,0,1,10)))
  expect_equal(Hsurvspline3(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], 
                                knots=c(-10,0,1,2,10)),
               Hsurvspline(x=x, gamma[1:5], knots=c(-10,0,1,2,10)))
  expect_equal(Hsurvspline4(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], 
                                knots=c(-10,0,1,2,3,10)),
               Hsurvspline(x=x, gamma[1:6], knots=c(-10,0,1,2,3,10)))
  expect_equal(Hsurvspline5(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], 
                                knots=c(-10,0,1,2,3,4,10)),
               Hsurvspline(x=x, gamma[1:7], knots=c(-10,0,1,2,3,4,10)))
  expect_equal(Hsurvspline6(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], gamma["4"], gamma["5"], gamma["6"], gamma["7"], 
                                knots=c(-10,0,1,2,3,4,5,10)),
               Hsurvspline(x=x, gamma[1:8], knots=c(-10,0,1,2,3,4,5,10)))
  expect_equal(Hsurvspline7(x=x, gamma["0"], gamma["1"], gamma["2"], gamma["3"], 
                                gamma["4"], gamma["5"], gamma["6"], gamma["7"], gamma["8"], 
                                knots=c(-10,0,1,2,3,4,5,6,10)),
               Hsurvspline(x=x, gamma[1:9], knots=c(-10,0,1,2,3,4,5,6,10)))
})
