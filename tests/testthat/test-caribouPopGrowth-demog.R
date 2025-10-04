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
                           0.85, progress = FALSE) 
  
  expect_true(all(which(res2$lambdaE < 1) == 1:5))
  
  # vector for S_bar
  expect_error(caribouPopGrowth(1000, 20, 0.8, 
                                c(rep(0.7, 5), rep(0.8, 5)), progress = FALSE), 
               "must have length")
  
  res3 <- expect_warning(caribouPopGrowth(1:10*100, 20, 0.8, 
                           c(rep(0.6, 5), rep(0.95, 5)), progress = FALSE))
  
  expect_true(all(which(res3$lambda < 1) == 1:5))
  
  # there is an error for extreme values of R_bar or S_bar
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.6, S_bar = 0.2, progress = FALSE),
               "expected survival S_bar")
  
  expect_warning(caribouPopGrowth(1000, 20, R_bar = 0.9, S_bar = 0.7, progress = FALSE),
               "expected recruitment R_bar")
  
  res4 <- caribouPopGrowth(9000, 200, R_bar = 0.7, S_bar = 0.9, progress = FALSE)
  # lower recruitment at carrying capacity
  res5 <- caribouPopGrowth(9000, 200, R_bar = 0.7, S_bar = 0.9, P_K = 0.4,
                           progress = FALSE)
  
  expect_lt(res5$N, res4$N)
  
  res4_2 <- caribouPopGrowth(9000, 200, R_bar = 0.7, S_bar = 0.9, progress = FALSE, K = 50000)
  # lower actual carrying capacity
  res5_2 <- caribouPopGrowth(9000, 200, R_bar = 0.7, S_bar = 0.9, progress = FALSE, K = 10000)
  expect_lt(res5_2$N, res4_2$N)
  
  # population gets larger over longer time with high R and S
  res6 <- caribouPopGrowth(1000, 2, R_bar = 0.7, S_bar = 0.99, progress = FALSE)
  res7 <- caribouPopGrowth(1000, 200, R_bar = 0.8, S_bar = 0.99, progress = FALSE)
  
  expect_gt(res7$N - res6$N, 200)
})

test_that("pop Growth matches Johnson figures", {
  
  johnsonCompare <- read.csv(file.path("data", "Johnson_figures5_6.csv"))
  # johnsonCompare <- read.csv(file.path("tests/testthat/data", "Johnson_figures5_6.csv"))
  
  covTableSim <- data.frame(Anthro = c(0, 10, 20, 30, 40, 50, 
                                       60, 70, 80, 90, 100), 
                            fire_excl_anthro = 4.27)
  covTableSim$polygon <- paste0("Missisa_", covTableSim$Anthro + 1)
  covTableSim$area <- "FarNorth"
  covTableSim$Total_dist <- covTableSim$Anthro + covTableSim$fire_excl_anthro
  
  popGrowthPars <- demographicCoefficients(500,
                                           modelVersion = "Johnson",
                                           survivalModelNumber = "M1",
                                           recruitmentModelNumber = "M4"
  )
  
  rateSamplesLarge <- demographicRates(
    covTable = covTableSim,
    popGrowthPars = popGrowthPars,
    ignorePrecision = F, returnSample = T, useQuantiles = T
  )
  
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
  
  comp_plot <- ggplot2::ggplot(comp, ggplot2::aes(Anthro, dif))+
    ggplot2::geom_point()+
    ggplot2::facet_wrap(~var, scales = "free")
  
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
  
  # demographic rates match Johnson
  expect_true(all(dif_comp$mean_dif < dif_comp$exp_dif))
  
  # TODO: investigate why intercepts of R_bar and S_bar seem to be right but
  # slopes are slightly off. JH says I suspect it has something to do with
  # compounding slight deviations over time. I bet if we looked at 1 yr ahead
  # projections (each year independent of the next) the errors would remain more
  # constant as anthropogenic disturbance changes.
  
  if(F){
    # reproduce comparison plots used in the MS
    library(ggplot2)
    popGrowthParsSmall <- demographicCoefficients(35,
                                                  modelVersion = "Johnson",
                                                  survivalModelNumber = "M1",
                                                  recruitmentModelNumber = "M4"
    )
    
    rateSamples <- demographicRates(
      covTable = covTableSim,
      popGrowthPars = popGrowthParsSmall,
      ignorePrecision = F, returnSample = T, useQuantiles = T
    )
    
    rateSamples$S_PIlow <- 1
    rateSamples$S_PIhigh <- 1
    rateSamples$rep <- as.factor(rateSamples$replicate)
    levels(rateSamples$rep) <- sample(unique(rateSamples$replicate), replace = F)
    rateSamples$rep <- as.character(rateSamples$rep)
    rateSamplesLarge$S_PIlow <- 1
    rateSamplesLarge$S_PIhigh <- 1

    rateSamples$R_PIlow <- 1
    rateSamples$R_PIhigh <- 1
    rateSamples$fullGrp <- paste(rateSamples$rep) # ,rateSamples$Fire)
    rateSamplesLarge$R_PIlow <- 1
    rateSamplesLarge$R_PIhigh <- 1
    rateSamplesLarge$fullGrp <- paste(rateSamplesLarge$rep) # ,rateSamplesLarge$Fire)

    theme_set(theme_bw())
    base1 <- ggplot(data = rateSummaries,
                    aes(x = Anthro, y = S_bar, ymin = S_PIlow, ymax = S_PIhigh)) +
      geom_ribbon(fill="grey95",colour="grey95",data=johnsonCompare)+
      geom_line(data = subset(rateSamples), size = 0.5, alpha = 0.5,
                aes(x = Anthro, y = S_bar, group = rep, colour = rep)) +
      geom_line(colour = "grey50", linewidth = 2, linetype = "dotted") +
      geom_line(colour = "black", linewidth = 2, linetype = "dotted",data=johnsonCompare) +
      xlab("Anthropogenic Disturbance (%)") +
      ylab("Adult Female Survival") +
      scale_x_continuous(limits = c(-1, 90), breaks = c(0, 20, 40, 60, 80)) +
      scale_y_continuous(limits = c(0.65, 1)) +
      theme(legend.position = "none", plot.margin = margin(l = 0.6, unit = "cm"))

    plot_recruitment3 <- ggplot(data = rateSummaries,
                                aes(x = Anthro, y = R_bar * 100,
                                    ymin = R_PIlow * 100, ymax = R_PIhigh * 100)) +
      geom_ribbon(fill="grey95",colour="grey95",data=johnsonCompare)+
      geom_line(data = rateSamples, size = 0.5,
                aes(x = Anthro, y = R_bar * 100, group = fullGrp, color = fullGrp),
                alpha = 0.5) +
      geom_line(colour = "grey50", linewidth = 2, linetype = "dotted") +
      geom_line(colour = "black", linewidth = 2, linetype = "dotted",data=johnsonCompare) +
      scale_x_continuous(limits = c(-1, 90), breaks = c(0, 20, 40, 60, 80)) +
      scale_y_continuous(limits = c(0, 60), breaks = c(0, 10, 20, 30, 40, 50, 60)) +
      xlab("Anthropogenic disturbance (%)") +
      ylab("Recruitment (calves/100 cows)") +
      theme(legend.position = "none", plot.margin = margin(l = 0.6, unit = "cm"))
  }
  
  # demography
  rates <- data.frame(N0 = c(round(745 / 2)), rrp = 1)
  numSteps <- 20
  ratesLarge <- merge(rates, rateSamplesLarge)
  set.seed(1234)
  ratesLarge <- cbind(ratesLarge, caribouPopGrowth(ratesLarge$N0,
                                        numSteps = numSteps, R_bar = ratesLarge$R_bar,
                                        S_bar = ratesLarge$S_bar, probOption = "binomial", 
                                        progress = FALSE
  ))
  
  check <- subset(ratesLarge, Anthro == 10)
  testRates <- list(R_bar = mean(check$R_bar), S_bar = mean(check$S_bar))
  set.seed(1234)
  testCheck <- caribouPopGrowth(round(745 / 2), numSteps = 20, R_bar = testRates$R_bar, 
                   S_bar = testRates$S_bar, probOption = "continuous", 
                   interannualVar = F, progress = FALSE)
  # this is theoretical lambda - to confirm
  theor_lam <- (testRates$S_bar) * (1 + testRates$R_bar * 0.5) 
  
  # these should all be similar
  expect_equal(mean(check$lambdaE), theor_lam, tolerance = 0.003)
  expect_equal(testCheck$lambdaE, theor_lam, tolerance = 0.003)
  expect_equal(testCheck$lambdaE, mean(check$lambdaE), tolerance = 0.003)
})

test_that("interannualVar works as expected", {
  def_IV <- caribouPopGrowth(rep(1000, 10), 200, R_bar = 0.7, S_bar = 0.99, progress = FALSE)
  expect_s3_class(def_IV, "data.frame")
  
  phi_IV <- caribouPopGrowth(rep(1000, 10), 200, R_bar = 0.7, S_bar = 0.99, progress = FALSE,
                             interannualVar = list(R_phi = 0.5, S_phi = 0.5))
  expect_s3_class(phi_IV, "data.frame")
  
  annual_IV <- caribouPopGrowth(rep(1000, 10), 200, R_bar = 0.7, S_bar = 0.99, progress = FALSE,
                                interannualVar = list(R_annual = 0.5, S_annual = 0.5))
  expect_s3_class(annual_IV, "data.frame")
})
