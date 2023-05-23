test_that("default works", {
  mod <- caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                     Nthin = 2)
  
  tabAllRes(mod$result, 2009, 2023)
})
