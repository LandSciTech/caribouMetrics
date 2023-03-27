# To make a new test use
# shinytest2::record_test("tests/testthat/app")
# make sure to copy the code created since it will not be saved

# To run tests interactively use
# app <- shinytest2::AppDriver$new("tests/testthat/app", name = "hello")

# Tests use INSTALLED version of the package so make sure it is updated before
# doing tests

test_that("app starts and runs with defaults", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  app <- shinytest2::AppDriver$new(".", name = "hello")

  app$set_inputs(Run.model = "click")

  app$set_inputs(graphPanel = "Female population size")
  app$set_inputs(graphPanel = "Female recruitment")
  app$set_inputs(graphPanel = "Adult female survival")
  app$set_inputs(graphPanel = "Population growth rate")
  app$set_inputs(graphPanel = "Recruitment")

  app$expect_values()
})
