test_that("defaults work", {
 expect_is(getScenarioDefaults(), "data.frame")
})

test_that("can change with individual args", {
  expect_equal(getScenarioDefaults(aS = c(1:3))$aS, 1:3)
})

test_that("can change with data.frame", {
  out <- getScenarioDefaults(scns = data.frame(P = 1:5, curYear = 2010))
  expect_equal(out$P, 1:5)
  expect_equal(out$iYr[5], 2006)
})
