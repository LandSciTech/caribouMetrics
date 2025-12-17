test_that("bias correction works", {
  # number or reps
  nr <- 10
  set.seed(777)
  c_vals <- compositionBiasCorrection(w = 6,
                                      q = runif(nr, 0, 0.6),
                                      u = runif(nr, 0, 0.2),
                                      z = runif(nr, 0, 0.2),
                                      approx = FALSE)
  expect_length(c_vals, nr)
  
  set.seed(777)
  c_dist_table <- compositionBiasCorrection(w = 6,
                                            q = runif(nr, 0, 0.6),
                                            u = runif(nr, 0, 0.2),
                                            z = runif(nr, 0, 0.2),
                                            approx = TRUE)
  expect_is(c_dist_table, "data.frame")
  
  expect_equal(c_dist_table$m, mean(c_vals))
})

