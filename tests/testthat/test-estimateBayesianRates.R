s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)

s_data <- s_data %>% getCaribouYear() %>% filter(CaribouYear >= 2005)
r_data <- r_data %>% getCaribouYear() %>% filter(CaribouYear >= 2005)

test_that("multipop works", {

  multPop <- estimateBayesianRates(s_data, r_data, N0 = 500, niters = 20)
  
  # N0 is not really used until passed into trajectories
  wN0Var <- estimateBayesianRates(s_data, r_data, 
                                  N0 = data.frame(PopulationName = c("A", "B"), 
                                                  N0 = c(500, 5000), 
                                                  N.sd = c(20, 1000)),
                                  niters = 20, return_mcmc = TRUE)
  trajwN0 <- trajectoriesFromBayesian(wN0Var)
  
  # check N0 is different in two pops
  N0Pops <- trajwN0$samples %>% 
    filter(MetricTypeID == "N", Timestep == min(Timestep)) %>% 
    group_by(PopulationName) %>% 
    summarise(meanN = mean(Amount), minN = min(Amount), maxN = max(Amount))
  
  expect_lt(abs(N0Pops$meanN[1] - 500), 50)
  expect_lt(abs(N0Pops$meanN[2] - 5000), 500)
  
})

test_that("No survival works", {
  
  s_data <- s_data %>% 
    mutate(MortalitiesCertain = ifelse(Year > 2013, StartTotal, MortalitiesCertain))
  
  r_data <- r_data %>% 
    mutate(Cows = ifelse(Year > 2013, 0, Cows),
           Calves = ifelse(Year > 2013, 0, Calves)) 
  
  lowRates <- estimateBayesianRates(s_data, r_data, N0 = 500, niters = 20, return_mcmc = TRUE)
  
  lowSims <- simulateObservations(getScenarioDefaults(obsAnthroSlope = 0, projAnthroSlope = 0, 
                                                      curYear = 2016, collarCount = 20),
                                  trajectories = trajectoriesFromBayesian(lowRates)$samples %>%
                                    filter(Replicate == "x1"))
  # Gives error see issue #163
  # lowRates2 <- estimateBayesianRates(lowSims$simSurvObs %>% filter(!is.na(StartTotal)),
  #                                    lowSims$simRecruitObs%>% filter(!is.na(Cows)),
  #                                    N0 = 500, niters = 20, return_mcmc = TRUE)
  
})
