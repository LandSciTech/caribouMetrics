test_that("works for mult pops and fo app", {
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
  
  trajs_beta1 <- trajectoriesFromSummary(
    replicates = 5000, N0 = 100,
    Rbar = bind_rows(data.frame(mean = 0.19, sd = 0.23, lower = 0.13,
                                upper = 0.27, Annual = 1:5, Year = 1:5,
                                PopulationName = "A")), 
    Sbar = bind_rows(data.frame(mean = 0.94, sd = 0.61, lower = 0.86,
                                upper = 0.98, Annual = 1:5, Year = 1:5,
                                PopulationName = "A")),
    Riv = data.frame(R_cv_min = 0.01, R_cv_max = 0.7),
    Siv = data.frame(S_cv_min = 0.01, S_cv_max = 0.13),
    type = "beta")
  
  # When Year is < 100 still can track timesteps 
  expect_contains(colnames(trajs_beta$summary), "Timestep")
  
  # default works mult pop
  expect_is(plotTrajectories(trajs_beta), 
            "ggplot")

  
  # default works 1 pop
  plotTrajectories(trajs_beta1, metrics = "Female population size")
  
  # Works for app 1 pop
  plotTrajectories(trajs_beta1, metrics = "Female population size", 
                   multicolour_samples = FALSE, ribbon = FALSE)+
    ggplot2::labs(col = "Scenario", x = "Years in the future",
                  linewidth = "Trajectory")+
    ggplot2::guides(linewidth = ggplot2::guide_legend(override.aes = list(linewidth = c(1,2))),
                    color = ggplot2::guide_legend(override.aes = list(linewidth = c(2))))
  
  # Works for app 2 pop
  plotTrajectories(trajs_beta, metrics = "Female population size", 
                   multicolour_samples = FALSE, ribbon = FALSE, 
                   add_N0 = TRUE)+
    ggplot2::labs(col = "Scenario", x = "Years in the future",
                  linewidth = "Trajectory")+
    ggplot2::guides(linewidth = ggplot2::guide_legend(override.aes = list(linewidth = c(1,2))),
                    color = ggplot2::guide_legend(override.aes = list(linewidth = c(2))))
})
