test_that("app starts ok", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  shiny_app <- demographicProjectionApp(n = 10)
  app <- shinytest2::AppDriver$new(shiny_app, name = "hello")

  app$expect_values()
})
