test_that("defaults work", {
 expect_is(getScenarioDefaults(), "data.frame")
})

test_that("can change with individual args", {
  expect_equal(getScenarioDefaults(obsAnthroSlope = c(1:3))$obsAnthroSlope, 1:3)
})

test_that("can change with data.frame", {
  out <- getScenarioDefaults(paramTable = data.frame(obsYears = 1:5, curYear = 2010))
  expect_equal(out$obsYears, 1:5)
  expect_equal(out$startYear[5], 2006)
})
