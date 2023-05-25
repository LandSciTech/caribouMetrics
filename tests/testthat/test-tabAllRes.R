test_that("default works", {
  mod <- caribouBayesianIPM(Nchains = 1, Niter = 100, Nburn = 10,
                     Nthin = 2)
  
  expect_is(tabAllRes(mod$result, 2009, 2043), "data.frame")
  
  expect_is(tabAllRes(mod$result, 2009, 2043, doSummary = F), "data.frame")
  # TODO: why are there 45 rows per year and parameter?
})
