s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)

s_data <- s_data %>% getCaribouYear() %>% filter(CaribouYear >= 2010)
r_data <- r_data %>% getCaribouYear() %>% filter(CaribouYear >= 2010)

test_that("multipop works", {

  estimateBayesianRates(s_data, r_data, N0 = 500, niters = 20)
})

test_that("No survival works", {
  
  s_data <- s_data %>% 
    mutate(MortalitiesCertain = ifelse(Year > 2013, StartTotal, MortalitiesCertain))%>% 
    filter(PopulationName == "A")
  
  r_data <- r_data %>% 
    mutate(Cows = ifelse(Year > 2013, 0, Cows),
           Calves = ifelse(Year > 2013, 0, Calves)) %>% 
    filter(PopulationName == "A")
  
  lowRates <- estimateBayesianRates(s_data, r_data, N0 = 500, niters = 20, return_mcmc = TRUE)
  
  lowSims <- simulateObservations(getScenarioDefaults(obsAnthroSlope = 0, projAnthroSlope = 0, 
                                                      curYear = 2016, collarCount = 20),
                                  trajectories = trajectoriesFromBayesian(lowRates)$samples %>%
                                    filter(Replicate == "x1"))
  
  lowRates2 <- estimateBayesianRates(lowSims$simSurvObs %>% filter(!is.na(StartTotal)), 
                                     lowSims$simRecruitObs%>% filter(!is.na(Cows)),
                                     N0 = 500, niters = 20, return_mcmc = TRUE)
})
