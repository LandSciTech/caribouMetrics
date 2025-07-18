test_that("default works", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, cowMult = 3,
                              collarCount = 50)
  trajs <- getSimsInitial(replicates = 2)$samples
  expect_is(simulateObservations(trajs, paramTable = scns),
            "list")
})

test_that("multiple scenarios not allowed",{
  scns <- getScenarioDefaults(data.frame(iFire = 1:2), projYears = 10, obsYears = 10)
  trajs <- getSimsInitial(replicates = 2)$samples
  expect_error(simulateObservations(trajs, scns,
                                 freqStartsByYear = data.frame(Year = 2014:2023,
                                                               numStarts = 10),
                                 cowCounts = data.frame(Year = 2014:2023,
                                                        Cows = 10)),
            "cannot have multiple rows")
})

# TODO: add test for non-default popGrow table should use testPopGrowTable
# internally, do same in getPriors

test_that("collarCount and cowCount behave", {
  scns <- getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = NA)
  trajs <- getSimsInitial(replicates = 2)$samples
  
  simObs <- simulateObservations(trajs, scns)
  
  # if cowCount is 100 we observe 100
  expect_true(all(simObs$simRecruitObs$Cows == 100))
  
  expect_error(getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = 2),
               "not both")
  
  # if cowMult is 2 we observe max 2*collarCount but fewer when some deaths were observed
  scns2 <- getScenarioDefaults(collarCount = 30, cowMult = 2)
  
  simObs2 <- simulateObservations(trajs, scns2)
  
  simObs2$simRecruitObs %>% 
    left_join(simObs2$simSurvObs, 
              by = join_by(PopulationName, Replicate, Year)) %>% 
    mutate(pass = Cows == (30 - MortalitiesCertain) * 2) %>% 
    pull(pass) %>% all %>% 
    expect_true()
  
  # Test with months
  scns10 <- getScenarioDefaults(collarCount = 15, cowMult = 2)
  simObs_mon <- simulateObservations(trajs, scns10, 
                       surv_data = bboudata::bbousurv_a,
                       recruit_data = bboudata::bbourecruit_a)

  scns10$curYear
  simObs_mon$simSurvObs$Year %>% max()
  
  # should be different cow counts in simulated years
  expect_gt(simObs_mon$simRecruitObs %>% filter(Year < 2016) %>% pull(Cows) %>% mean, 
            simObs_mon$simRecruitObs %>% filter(Year > 2016) %>% pull(Cows) %>% mean)
  
  # Visual:
  # simObs_mon$simRecruitObs %>% ggplot(aes(Year, Cows, colour = Replicate))+
  #   geom_point()
  
  
  
  # Got a warning of NAs produced by: rbinom(n = nrow(simRecruitObs), size =
  # round(apparentCows), prob = simRecruitObs$recruitment)
  # Only happens sometimes with low collarCount. Need to figure it out and make a test
  # Haven't been able to reproduce...
  
  #TODO: This seems off... What is the point of exData? If we are simulating
  # observations why does this output go beyond the simulation time period?
  
  trajs_low <- getSimsInitial(replicates = 2,
                              cPars = getScenarioDefaults(lQuantile = 0.99))$samples

  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, projYears = 100, 
                                projAnthroSlope = 0)

  simObs8 <- simulateObservations(trajs_low, scns10,
                                  surv_data = bboudata::bbousurv_a,
                                  recruit_data = bboudata::bbourecruit_a)

  ggplot(simObs8$exData %>% filter(MetricTypeID == "N"),
         aes(Year, Amount, colour = Replicate))+
    geom_point()
  
  # if tables are supplied they should not be modified by cowCount or collarCount
  simObs3 <- simulateObservations(trajs, scns,
                       freqStartsByYear = data.frame(Year = 2009:2023,
                                                     numStarts = 10),
                       cowCounts = data.frame(Year = 2009:2023,
                                              Cows = 10))
  
  simObs3$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # Total n collars goes up at first and then balances with mortality
  expect_lt(
    simObs3$simSurvObs %>% filter(Year == 2010) %>% pull(StartTotal) %>% mean(),
    simObs3$simSurvObs %>% filter(Year == 2020) %>% pull(StartTotal) %>% mean()
  )
  # Visual test
  # simObs3$simSurvObs %>% 
  #   ggplot(aes(Year, StartTotal, colour = Replicate)) +
  #   geom_point()
  
  
  # cowMult doesn't affect cowCounts table
  simObs4 <- simulateObservations(trajs, scns2,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10))
  
  simObs4$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # can supply just freqStartsByYear and cowMult
  simObs4b <- simulateObservations(trajs, scns2,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10))
  scns2$cowMult <- 10
  simObs4c <- simulateObservations(trajs, scns2,
                                   freqStartsByYear = data.frame(Year = 2009:2023,
                                                                 numStarts = 10))
  
  # cowCounts is created from freqStartsByYear and cowMult 
  expect_true(all(
    simObs4b$simRecruitObs %>% pull(Cows) <
      simObs4c$simRecruitObs %>% pull(Cows))) 
  
  # for collarInterval 
  scns3 <- getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = NA,
                               collarInterval = 3)
  
  simObs5 <- simulateObservations(trajs, scns3)
  
  # simObs5$simSurvObs %>%
  #   ggplot(aes(Year, StartTotal, colour = Replicate)) +
  #   geom_point()
  
  # Should go back up to 30 collars every 3 years
  simObs5$simSurvObs %>% pull(StartTotal) %>% .[seq(1,15, by = 3)] %>% 
    {. == 30} %>% all() %>% 
    expect_true()
  
  simObs5$simRecruitObs %>% filter(Cows != 100) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # collarInterval doesn't affect tables
  simObs6 <- simulateObservations(trajs, scns3,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10))
  
  simObs6$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # confirm that if freqStartByYear table does skips year it still works
  simObs7 <- simulateObservations(trajs, scns,
                                  freqStartsByYear = data.frame(Year = seq(2009, 2023, by = 3),
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10))
  
  simObs7$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
})