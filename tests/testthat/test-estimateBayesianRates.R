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
  
  
  # when one pop is missing all years of data
  expect_warning({
    multPop2 <- estimateBayesianRates(bboudata::bbousurv_multi %>% getCaribouYear() %>%
                                        filter(CaribouYear %>% between(2010, 2015)), 
                                      bboudata::bbourecruit_multi %>% getCaribouYear() %>%
                                        filter(CaribouYear %>% between(2010, 2015)),
                                      N0 = 500, niters = 20)
  })
  
  expect_equal(multPop2$PopulationName, c("B", "C"))

})

test_that("No survival works", {
  
  s_data <- s_data %>% 
    mutate(MortalitiesCertain = ifelse(Year > 2013, StartTotal, MortalitiesCertain))
  
  r_data <- r_data %>% 
    mutate(Cows = ifelse(Year > 2010, 0, Cows),
           Calves = ifelse(Year > 2010, 0, Calves)) 
  
  lowRates <- estimateBayesianRates(s_data, r_data, N0 = 500, niters = 20, return_mcmc = TRUE)
  
  lowSims <- simulateObservations(getScenarioDefaults(obsAnthroSlope = 0, projAnthroSlope = 0, 
                                                      curYear = 2016, collarCount = 20),
                                  trajectories = trajectoriesFromBayesian(lowRates)$samples %>%
                                    filter(Replicate == "x1"))
  # Works now
  lowRates2 <- estimateBayesianRates(lowSims$simSurvObs %>%
                                       mutate(StartTotal = ifelse(is.na(StartTotal),
                                                                  10, StartTotal),
                                              MortalitiesCertain = ifelse(is.na(MortalitiesCertain),
                                                                          10, MortalitiesCertain),
                                              MortalitiesUncertain = 0, 
                                              Malfunctions = 0),
                                     lowSims$simRecruitObs %>% 
                                       mutate(Cows = ifelse(is.na(Cows), 10, Cows),
                                              Calves = ifelse(is.na(Calves), 0, Calves),
                                              Bulls = ifelse(is.na(Bulls), 0, Bulls),
                                              UnknownAdults = ifelse(is.na(UnknownAdults), 0, UnknownAdults),
                                              Yearlings = ifelse(is.na(Yearlings), 0, Yearlings),
                                              CowsBulls = ifelse(is.na(CowsBulls), 0, CowsBulls)),
                                     N0 = 500, niters = 20, return_mcmc = TRUE)
  
  expect_is(lowRates2, "list")
})
