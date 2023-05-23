test_that("default works", {
  mod <- caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                     Nthin = 2)
  
  expect_is(tabAllRes(mod$result, 2009, 2023), "data.frame")
  
  expect_is(tabAllRes(mod$result, 2009, 2023, doSummary = F), "data.frame")
  # TODO: why are there 45 rows per year and parameter?
})
