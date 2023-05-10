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
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.6, S_bar = 0.2, progress = FALSE),
               "expected survival S_bar")
  
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.9, S_bar = 0.7, progress = FALSE),
               "expected recruitment R_bar")
  
  res4 <- caribouPopGrowth(1000, 200, R_bar = 0.7, S_bar = 0.9, progress = FALSE)
  # lower carrying capacity
  res5 <- caribouPopGrowth(1000, 200, R_bar = 0.7, S_bar = 0.9, P_K = 0.4,
                           progress = FALSE)
  
  expect_lt(res5$N, res4$N)
  
  # population gets larger over longer time with high R and S
  res6 <- caribouPopGrowth(1000, 2, R_bar = 0.7, S_bar = 0.99, progress = FALSE)
  res7 <- caribouPopGrowth(1000, 200, R_bar = 0.8, S_bar = 0.99, progress = FALSE)
  
  expect_gt(res7$N - res6$N, 200)
})

test_that("pop Growth matches Johnson figures", {
  johnsonCompare <- read.csv(file.path("data", "Johnson et al. figures5_6.csv"))
  
  covTableSim <- data.frame(Anthro = c(0, 10, 20, 30, 40, 50, 
                                       60, 70, 80, 90, 100), 
                            fire_excl_anthro = 5)
  covTableSim$polygon <- paste0("Missisa_", covTableSim$Anthro + 1)
  covTableSim$area <- "FarNorth"
  covTableSim$Total_dist <- covTableSim$Anthro + covTableSim$fire_excl_anthro
  
  popGrowthPars <- demographicCoefficients(500,
                                           modelVersion = "Johnson",
                                           survivalModelNumber = "M1",
                                           recruitmentModelNumber = "M4")
  
  rateSummaries <- demographicRates(
    covTable = covTableSim, popGrowthPars = popGrowthPars,
    ignorePrecision = F, returnSample = F, useQuantiles = F
  )

  johnsonCompare$Anthro=johnsonCompare$anthro
  johnsonCompare$S_bar=johnsonCompare$Sresp
  johnsonCompare$S_PIlow = johnsonCompare$Slow_95_pred
  johnsonCompare$S_PIhigh = johnsonCompare$Shigh_95_pred
  johnsonCompare$R_bar=johnsonCompare$Rresp/100
  johnsonCompare$R_PIlow = johnsonCompare$Rlow_95_pre/100
  johnsonCompare$R_PIhigh = johnsonCompare$Rhigh_95_pred/100
  johnsonCompare$method <- "johnson"
  
  rateSummaries$method <- "caribouMetrics"
  
  comp <- bind_rows(rateSummaries%>% select(any_of(colnames(johnsonCompare))),
                    johnsonCompare %>% 
                      select(any_of(colnames(rateSummaries))) %>% 
                      semi_join(rateSummaries, by = "Anthro")) %>% 
    arrange(Anthro) %>% 
    pivot_longer(cols = matches("R_|S_") , names_to = "var", values_to = "value") %>% 
    pivot_wider(names_from = "method", values_from = "value") %>% 
    mutate(dif = johnson - caribouMetrics)
  
  ggplot(comp, aes(Anthro, dif))+
    geom_point()+
    facet_wrap(~var)
  
  expected_difs <- dplyr::tribble(
    ~var, ~exp_dif,
    "S_bar", 0.005,
    "S_PIlow", 0.02,
    "S_PIhigh", 0.01,
    "R_bar", 0.005,
    "R_PIlow", 0.01,
    "R_PIhigh", 0.06
  )
  
  dif_comp <- comp %>% mutate(abs_dif = abs(dif)) %>% 
    group_by(var) %>% 
    summarise(mean_dif = mean(abs_dif)) %>% 
    left_join(expected_difs, by = "var")
  
  expect_true(all(dif_comp$mean_dif < dif_comp$exp_dif))
  
  dia_shp <- 23
  
  err_col <- "grey50"
  
  theme_set(theme_bw())
  
  # TODO: finish adapting this to work as a test
  
  # demography
  pars <- data.frame(N0 = c(round(745 / 2)))
  # increase to get a better sample size, or set interannualVar to NA
  pars <- merge(pars, data.frame(rrp = 1))
  pars <- merge(pars, rateSamplesLarge)
  numSteps <- 20
  pars1 <- cbind(pars, popGrowthJohnson(pars$N0,
                                        numSteps = numSteps, R_bar = pars$R_bar,
                                        S_bar = pars$S_bar, probOption = "binomial"
  ))
  
  pars <- data.frame(N0 = round(745 / 2))
  # increase to get a better sample size, or set interannualVar to NA
  pars <- merge(pars, data.frame(rrp = 1)) 
  pars <- merge(pars, rateSamples)
  numSteps <- 20
  pars2 <- cbind(pars, popGrowthJohnson(pars$N0,
                                        numSteps = numSteps, R_bar = pars$R_bar,
                                        S_bar = pars$S_bar, probOption = "binomial"
  ))
  
  
  check <- subset(pars1, Anthro == 5)
  testPars <- list(R_bar = mean(check$R_bar), S_bar = mean(check$S_bar))
  popGrowthJohnson(round(745 / 2), numSteps = 20, R_bar = testPars$R_bar, 
                   S_bar = testPars$S_bar, probOption = "continuous", 
                   interannualVar = F)
  # this is theoretical lambda - to confirm
  (testPars$S_bar) * (1 + testPars$R_bar * 0.5) 
  mean(check$lambda)
  testPars
  # The math seems ok here.
  
  oo <- pars2 %>%
    select(Anthro, lambda, fullGrp, rrp) %>%
    group_by(fullGrp, Anthro) %>%
    summarise(lambda = median(lambda))
  
  subset(oo, Anthro == 5)
  ooT <- pars1 %>%
    select(Anthro, lambda, fullGrp, rrp) %>%
    group_by(fullGrp, Anthro) %>%
    summarise(lambda = median(lambda))
  
  ooS <- ooT %>%
    select(Anthro, lambda) %>%
    group_by(Anthro) %>%
    summarise(lambdaH = max(lambda), lambdaL = min(lambda), lambda = median(lambda))
  
  oo$lambdaH <- oo$lambda
  oo$lambdaL <- oo$lambda
  str(ooS)
  plot_lambda <- ggplot(oo, 
                        aes(x = Anthro, y = lambda, ymin = lambdaH, ymax = lambdaL)) +
    geom_line(size = 0.5,
              aes(x = Anthro, y = lambda, group = fullGrp, color = fullGrp), 
              alpha = 0.5) +
    geom_errorbar(data = subset(ooS, Anthro < 10), width = 2, 
                  size = 0.7, col = err_col) +
    geom_point(data = data.frame(Anthro = 0.27, lambda = 0.86, fullGrp = 0,
                                 lambdaH = 0.86, lambdaL = 0.86), 
               size = 2, shape = dia_shp, fill = "black") +
    scale_x_continuous(limits = c(-1, 90), breaks = c(0, 20, 40, 60, 80)) +
    xlab("Anthropogenic disturbance (%)") +
    ylab(expression("Average Growth Rate " * lambda)) +
    theme(legend.position = "none", plot.margin = margin(l = 0.6, unit = "cm"))
  
})
