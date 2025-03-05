test_that("default works", {
  mod <- caribouBayesianPM(niters=100)
  
  expect_is(tabAllRes(mod$result, 2009, 2043), "data.frame")
  
  expect_is(tabAllRes(mod$result, 2009, 2043, doSummary = F), "data.frame")
  # TODO: why are there 45 rows per year and parameter?
})
