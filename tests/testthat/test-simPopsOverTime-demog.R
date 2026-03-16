test_that("Works for multiple populations with constant rates", {
  multPop <- simPopsOverTime(N0 = 100, R_samp = 1:10/10, S_samp = 1:10/10, numSteps = 10, 
                             interannualVar = FALSE, 
                             # set to avoid warnings
                             l_R = 0, h_R = 1, l_S = 0, h_S = 1, 
                             progress = FALSE)
  
  expect_true(all(unique(multPop$R_t) == 1:10/10))
  expect_equal(nrow(multPop), 100)
  
})

test_that("Works for single population with dynamic rates", {
  onePop <- simPopsOverTime(N0 = 100, R_samp = 1:10/10, S_samp = 1:10/10, numSteps = 10, 
                            dynamicRates = TRUE,
                            interannualVar = FALSE, 
                            # set to avoid warnings
                            l_R = 0, h_R = 1, l_S = 0, h_S = 1, 
                            progress = FALSE)
  expect_true(all(unique(onePop$R_t) == 1:10/10))
  expect_equal(nrow(onePop), 10)
})

test_that("Works for multiple populations with dynamic rates", {
  R_samp = matrix(rep(1:10/10, 10)*sort(rep(1:10/20, 10)*rep(c(1, 1.5), 50)), nrow = 10, byrow = TRUE)
  S_samp = matrix(rep(1:10/10, 10)*sort(rep(1:10/20, 10)) + 0.5, nrow = 10, byrow = TRUE)
  multDynamPop <- simPopsOverTime(N0 = 100, R_samp = R_samp, S_samp = S_samp, numSteps = 10, 
                            dynamicRates = TRUE,
                            interannualVar = FALSE, 
                            # set to avoid warnings
                            l_R = 0, h_R = 1, l_S = 0, h_S = 1, 
                            progress = FALSE)
  expect_equal(multDynamPop %>% filter(id == 10) %>% pull(R_t), R_samp[10,])
  expect_equal(multDynamPop %>% filter(id == 10) %>% pull(S_t), S_samp[10,])
  expect_equal(nrow(multDynamPop), 100)
})

