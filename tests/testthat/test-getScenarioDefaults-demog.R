test_that("defaults work", {
 expect_is(getScenarioDefaults(), "data.frame")
})

test_that("can change with individual args", {
  defs <- getScenarioDefaults(obsAnthroSlope = c(1:3))
  expect_equal(defs$obsAnthroSlope, 1:3)
  # ensure labels are unique
  expect_equal(n_distinct(defs$label), 3)
})

test_that("can change with data.frame", {
  out <- getScenarioDefaults(paramTable = data.frame(obsYears = 1:5, curYear = 2010))
  expect_equal(out$obsYears, 1:5)
  expect_equal(out$startYear[5], 2006)
})

test_that("N0 variation works", {
  scns <- getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = 1, 
                              N0 = data.frame(N0 = 5000, N.sd = 100))
  expect_is(scns, "data.frame")
})
