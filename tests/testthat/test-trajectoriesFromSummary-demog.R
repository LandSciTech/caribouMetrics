test_that("trajectoriesFromSummary works with variation in N0", {
  N0df <- data.frame(N0 = 1000, N.sd =1000^0.5,N.lower=900,N.upper=1100)

  
  traj <- trajectoriesFromSummary(replicates = 35, N0 = N0df,
                          Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                                            upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                                            PopulationName = "A"),
                          Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                                            upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                                            PopulationName = "A"),
                          Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
                          Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
                          type = "bbou")
  
  traj$popInfo$N0 %>% range() %>% diff() %>% expect_gt(100)
  
  #Expect truncated poisson distribution
  Ndist <- traj$popInfo$N0
  
  expect_true((var(Ndist) - 1000) < 250)
  expect_true(abs(mean(Ndist)-1000)< 10)  
  expect_true(min(Ndist)>=900)
  expect_true(max(Ndist)<=1100)
})

test_that("summary gives expected trajectory", {
  trajs <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
    Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
    type = "bbou")
  
  testthat::expect_true(
    all(trajs$summary %>% filter(MetricTypeID == "recruitment") %>% pull(Mean) >
          trajs$summary %>% filter(MetricTypeID == "Rbar") %>% pull(lower))
  )
  
  testthat::expect_true(
    all(trajs$summary %>% filter(MetricTypeID == "recruitment") %>% pull(Mean) <
          trajs$summary %>% filter(MetricTypeID == "Rbar") %>% pull(upper))
  )
  
  trajs_beta <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_cv_min = 0.01, R_cv_max = 0.7),
    Siv = data.frame(S_cv_min = 0.01, S_cv_max = 0.13),
    type = "beta")
  
  trajs_beta %>% plotTrajectories(metrics = "Population growth rate")
  trajs %>% plotTrajectories(metrics = "Population growth rate")
  
  trajs_beta %>% plotTrajectories(metrics = "Expected growth rate")
  trajs %>% plotTrajectories(metrics = "Expected growth rate")
  
  # WIP not working
  # # beta and logistic give similar results
  # expect_equal(
  #   trajs_beta %>% filter(type == "samp") %>%
  #     summarise(mlambda = mean(lambda)),
  #   trajs %>% filter(type == "samp") %>%
  #     summarise(mlambda = mean(lambda)),
  #   tolerance = 0.01
  # )
  
  # varPersists works as expected
  trajs_beta_vPF <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_cv_min = 0.01, R_cv_max = 0.7),
    Siv = data.frame(S_cv_min = 0.01, S_cv_max = 0.13),
    type = "beta", varPersists = FALSE)

  trajs_beta_vPF %>% plotTrajectories(metrics = "Population growth rate")
  
  trajs_no_sum <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
    Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
    type = "bbou", doSummary = FALSE)
  
  # setting no Summary doesn't change mean lambda
  expect_equal(
    trajs$summary %>% filter(MetricTypeID == "lambda") %>%
      summarise(mlambda = mean(Mean)),
    trajs_no_sum %>%
      summarise(mlambda = mean(lambda)),
    tolerance = 0.01
  )

  # can have just one year in S/Rbar
  expect_is(trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Year = 2001, Annual = 2001,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Year = 2001, Annual = 2001,
                      PopulationName = "A"),
    Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
    Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
    type = "bbou"), "list")

})


test_that("adjust works as expected", {
  trajs <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
    Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
    type = "bbou")
  
  
  trajs_adjR <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A", adjust.mu = -0.05, adjust.sd = 0.01),
    Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
                      PopulationName = "A"),
    Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
    Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
    type = "bbou")
  
  plotTrajectories(trajs_adjR, metrics = c("Expected recruitment", "Recruitment"))+
    ggplot2::ylim(0, 0.5)+ggplot2::ggtitle("adjustedR trajs")
  plotTrajectories(trajs, metrics = c("Expected recruitment", "Recruitment"))+
    ggplot2::ylim(0, 0.5)+ggplot2::ggtitle("trajs")
  
  
  # Currently errors see # 165
  # trajs_adjR <- trajectoriesFromSummary(
  #   replicates = 5000, N0 = 100,
  #   Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
  #                     upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
  #                     PopulationName = "A", adjust.mu = -0.05, adjust.sd = 0.2),
  #   Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
  #                     upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
  #                     PopulationName = "A"),
  #   Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
  #   Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
  #   type = "bbou")
  

  # trajs_adjS <- trajectoriesFromSummary(
  #   replicates = 5000, N0 = 100,
  #   Rbar = data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
  #                     upper = 0.27, Annual = 2010:2015, Year = 2010:2015,
  #                     PopulationName = "A"),
  #   Sbar = data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
  #                     upper = 0.98, Annual = 2010:2015, Year = 2010:2015,
  #                     PopulationName = "A", adjust.mu = 0.1, adjust.sd = 0.01),
  #   Riv = data.frame(R_iv_mean = 0.36, R_iv_shape = 2),
  #   Siv = data.frame(S_iv_mean = 0.63, S_iv_shape = 1.4),
  #   type = "bbou")
})

test_that("works with multipop",{
  trajs_beta <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = bind_rows(data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                      upper = 0.27, Annual = 1:5, Year = 1:5,
                      PopulationName = "A"),
                     data.frame(mean = 0.3, sd = 0.23, lower = 0.2,
                                upper = 0.4, Annual = 1:5, Year = 1:5,
                                PopulationName = "B")), 
    Sbar = bind_rows(data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                      upper = 0.98, Annual = 1:5, Year = 1:5,
                      PopulationName = "A"),
                     data.frame(mean = 0.9, sd = 0.61, lower = 0.8,
                                upper = 0.95, Annual = 1:5, Year = 1:5,
                                PopulationName = "B")),
    Riv = data.frame(R_cv_min = 0.01, R_cv_max = 0.7),
    Siv = data.frame(S_cv_min = 0.01, S_cv_max = 0.13),
    type = "beta")
  
  # When Year is < 100 still can track timesteps 
  expect_contains(colnames(trajs_beta$summary), "Timestep")
})

