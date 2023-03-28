# To make a new test use
# shinytest2::record_test("tests/testthat/app")
# make sure to copy the code created since it will not be saved

# To run tests interactively use
# app <- shinytest2::AppDriver$new("tests/testthat/app", name = "app-test")

# To compare snapshots after testing:
# testthat::snapshot_review('demographicProjectionApp/', path = "tests/testthat/app/")

# Tests use INSTALLED version of the package so make sure it is updated before
# doing tests

test_that("app starts and runs with defaults", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  app <- shinytest2::AppDriver$new(".", name = "app-test", seed = 7)

  app$set_inputs(Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2)

  app$set_inputs(Run.model = "click")

  app$set_inputs(graphPanel = "Female population size")
  app$set_inputs(graphPanel = "Female recruitment")
  app$set_inputs(graphPanel = "Adult female survival")
  app$set_inputs(graphPanel = "Population growth rate")
  app$set_inputs(graphPanel = "Recruitment")

  app$expect_values()
})

test_that("app starts and runs KS calculated", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  app <- shinytest2::AppDriver$new(".", name = "app-test2", seed = 7)

  app$set_inputs(Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2)

  app$set_inputs(adjustR = TRUE)
  app$set_inputs(getKSDists = TRUE)

  app$set_inputs(Run.model = "click")

  app$set_inputs(graphPanel = "Recruitment")
  app$set_inputs(graphPanel = "Adult female survival")
  app$set_inputs(graphPanel = "Recruitment Kolmogorov-Smirnov Distance")
  app$set_inputs(graphPanel = "Adult female survival Kolmogorov-Smirnov Distance")
  app$set_inputs(graphPanel = "Population growth rate Kolmogorov-Smirnov Distance")

  app$expect_values()
})

