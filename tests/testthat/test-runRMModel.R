
test_that("Runs with defaults", {
  # reduce some defaults to make fast
  # note that the default csv does not match the default startYear
  expect_warning(runRMModel(Nchains = 1, Niter = 100, Nburn = 10, Nthin = 2))
  expect_is(runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                       Nthin = 2),
            "list")
})
