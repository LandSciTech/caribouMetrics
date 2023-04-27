test_that("caribouPopGrowth works", {
  expect_is(caribouPopGrowth(200, 20, 0.5, 0.8, progress = FALSE), "data.frame")
})

test_that("caribouPopGrowth options work", {
  # vector for N
  res1 <- caribouPopGrowth(c(200, 400, 600) , 20, 0.5, 0.8, progress = FALSE)
  
  # vector for R_bar
  expect_error(caribouPopGrowth(1000, 20, c(rep(0.2, 5), rep(0.8, 5)),
                                0.8, progress = FALSE), 
               "must have length")
  
  res2 <- caribouPopGrowth(1:10*100, 20, c(rep(0.2, 5), rep(0.8, 5)),
                           0.8, progress = FALSE) 
  
  expect_true(all(which(res2$lambda < 1) == 1:5))
  
  # vector for S_bar
  expect_error(caribouPopGrowth(1000, 20, 0.8, 
                                c(rep(0.7, 5), rep(0.8, 5)), progress = FALSE), 
               "must have length")
  
  res3 <- caribouPopGrowth(1:10*100, 20, 0.8, 
                           c(rep(0.7, 5), rep(0.8, 5)), progress = FALSE) 
  
  expect_true(all(which(res3$lambda < 1) == 1:5))
  
  # there is an error for extreme values of R_bar or S_bar
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.6, S_bar = 0.2),
               "expected survival S_bar")
  
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.9, S_bar = 0.7),
               "expected recruitment R_bar")
  
  res4 <- caribouPopGrowth(1000, 200, R_bar = 0.7, S_bar = 0.9, progress = FALSE)
  # lower carrying capacity
  res5 <- caribouPopGrowth(1000, 200, R_bar = 0.7, S_bar = 0.9, P_K = 0.4,
                           progress = FALSE)
  
  expect_lt(res5$N, res4$N)
})
