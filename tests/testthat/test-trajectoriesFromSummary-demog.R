test_that("summary gives expected trajectory", {
  trajs <- trajectoriesFromSummary(
    numSteps = 10, replicates = 5000, N0 = 100, R_bar = 0.3,
    S_bar = 0.8, R_sd = 0.05, S_sd = 0.1, R_iv_shape = 0.01, 
    R_iv_mean = 0.01, S_iv_mean = 0.05, S_iv_shape = 0.05, 
    scn_nm = "test")
  
  expect_equal(
    trajs %>% filter(type == "mean") %>% pull(R_t) %>% unique(),
    0.3
  )

  trajs_beta <- trajectoriesFromSummary(
    numSteps = 10, replicates = 5000, N0 = 100, R_bar = 0.3,
    S_bar = 0.8, R_sd = 0.05, S_sd = 0.1, R_iv_shape = 0.01, 
    R_iv_mean = 0.01, S_iv_mean = 0.05, S_iv_shape = 0.05, 
    scn_nm = "test", type = "beta")
  
  # beta and logistic give similar results
  expect_equal(
    trajs_beta %>% filter(type == "samp") %>% 
      summarise(mlambda = mean(lambda)),
    trajs %>% filter(type == "samp") %>% 
      summarise(mlambda = mean(lambda)),
    tolerance = 0.01
  )  
  
  trajs_w_sum <- trajectoriesFromSummary(
    numSteps = 10, replicates = 5000, N0 = 100, R_bar = 0.3,
    S_bar = 0.8, R_sd = 0.05, S_sd = 0.1, R_iv_shape = 0.01, 
    R_iv_mean = 0.01, S_iv_mean = 0.05, S_iv_shape = 0.05, 
    scn_nm = "test", doSummary = TRUE)
  
  # setting no Summary doesn't change mean lambda
  expect_equal(
    trajs_w_sum$summary %>% filter(MetricTypeID == "lambda") %>% 
      summarise(mlambda = mean(Mean)),
    trajs %>% filter(type == "samp") %>% 
      summarise(mlambda = mean(lambda)),
    tolerance = 0.01
  ) 
  
  # works with N0 length 2
  
  # can't have multiple R_bar because then sample from multiple distributions
  # which is confusing
  expect_error(trajectoriesFromSummary(
    numSteps = 10, replicates = 5000, N0 = c(100, 200), R_bar = c(1:3/9.01),
    S_bar = c(7:9/9.01), R_sd = 0.05, S_sd = 0.1, R_iv_shape = 0.01, 
    R_iv_mean = 0.01, S_iv_mean = 0.05, S_iv_shape = 0.05, 
    scn_nm = "test"), "length one")
  
  trajs_rng <- trajectoriesFromSummary(
    numSteps = 10, replicates = 5000, N0 = c(100, 200), R_bar = 0.3,
    S_bar = 0.8, R_sd = 0.05, S_sd = 0.1, R_iv_shape = 0.01, 
    R_iv_mean = 0.01, S_iv_mean = 0.05, S_iv_shape = 0.05, 
    scn_nm = "test")
  
  expect_equal(trajs_rng %>% filter(time == 1, type == "samp") %>% pull(N0) %>% n_distinct(), 
               101)
  
})