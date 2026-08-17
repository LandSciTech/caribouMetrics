# Annual series tested in simulateObservations
test_that("monthly data works", {
  plt <- bboudata::bbousurv_a %>% plotSurvivalSeries()
  expect_is(plt, "ggplot2::ggplot")
})
