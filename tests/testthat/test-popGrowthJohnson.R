test_that("popGrowthJohnson works", {
  expect_is(popGrowthJohnson(200, 20, 0.5, 0.8), "data.frame")
})
