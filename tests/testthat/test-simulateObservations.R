test_that("default works", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, cowCount = 100, collarCount = 50)
  expect_is(simulateObservations(scns),
            "list")
})

test_that("error messages are as expected", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10)
  expect_error(
    simulateObservations(scns,
                         freqStartsByYear = data.frame(Year = 2014:2023,
                                                       numStarts = 10),
                         cowCounts = data.frame(Year = 2014:2016,
                                                Count = 10,
                                                Class = "cow")),
    "Year is missing expected values")

  expect_error(
    simulateObservations(scns,
                         freqStartsByYear = data.frame(Year = 2014:2023,
                                                       numStarts = 10),
                         cowCounts = data.frame(Year = 2014:2023,
                                                Count = 10,
                                                Class = "cat")),
    "Class contains unexpected values")
})

test_that("multiple scenarios not allowed",{
  scns <- getScenarioDefaults(data.frame(iFire = 1:2), projYears = 10, obsYears = 10)
  expect_error(simulateObservations(scns,
                                 freqStartsByYear = data.frame(Year = 2014:2023,
                                                               numStarts = 10),
                                 cowCounts = data.frame(Year = 2014:2023,
                                                        Count = 10,
                                                        Class = "cow")),
            "must have length 1")
})

# TODO: add test for non-default popGrow table should use testPopGrowTable
# internally, do same in getPriors

test_that("collarCount and cowCount behave", {
  scns <- getScenarioDefaults(collarCount = 30, cowCount = 100)
  
  simObs <- simulateObservations(scns)
  
  # number of rows in each year should be 30 next years should add n
  # died plus n dropped in prev year
  simObs$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% 
    {pull(., ncollar) == 30} %>% all() %>% 
    expect_true()
  
  simObs$ageRatioOut %>% filter(Class == "cow", Count != 100) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # cowMult
  expect_error(getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = 2),
               "not both")
  
  scns2 <- getScenarioDefaults(collarCount = 30, cowMult = 2)
  
  simObs2 <- simulateObservations(scns2)
  
  collarSum <- simObs2$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7), 
              survsCalving = sum(exit >= 6))
    
  expect_true(all(
  simObs2$ageRatioOut %>% filter(Class == "cow") %>% pull(Count) ==
    collarSum$survsCalving*2))
  
  # if tables are supplied they should not be modified by cowCount or collarCount
  simObs3 <- simulateObservations(scns,
                       freqStartsByYear = data.frame(Year = 2009:2023,
                                                     numStarts = 10),
                       cowCounts = data.frame(Year = 2009:2023,
                                              Count = 10,
                                              Class = "cow"))
  
  simObs3$ageRatioOut %>% filter(Class == "cow", Count != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # number added remains constant rather than number collars
  simObs3$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% 
    {pull(., nadded) == 10} %>% all() %>% 
    expect_true()
  
  # cowMult doesn't affect cowCounts table
  simObs4 <- simulateObservations(scns2,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Count = 10,
                                                         Class = "cow"))
  
  simObs4$ageRatioOut %>% filter(Class == "cow", Count != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # number added remains constant rather than number collars
  simObs4$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% 
    {pull(., nadded) == 10} %>% all() %>% 
    expect_true()
  
  # can supply just freqStartsByYear and cowMult
  simObs4b <- simulateObservations(scns2,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10))
  
  collarSum2 <- simObs4b$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7), 
              survsCalving = sum(exit >= 6))
  
  # cowCounts is created from freqStartsByYear and cowMult 
  expect_true(all(
    simObs4b$ageRatioOut %>% filter(Class == "cow") %>% pull(Count) ==
      collarSum2$survsCalving*2)) 
  
  # number added remains constant rather than number collars
  collarSum2 %>% 
    {pull(., nadded) == 10} %>% all() %>% 
    expect_true()
  
  # for collarInterval 
  scns3 <- getScenarioDefaults(collarCount = 30, cowCount = 100, collarInterval = 3)
  
  simObs5 <- simulateObservations(scns3)
  
  # number of rows in each year should be 30 next years should add n
  # died plus n dropped in prev year
  simObs5$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% pull(ncollar) %>%
    .[seq(1,15, by = 3)] %>% 
    {. == 30} %>% all() %>% 
    expect_true()
  
  simObs5$ageRatioOut %>% filter(Class == "cow", Count != 100) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # collarInterval doesn't affect tables
  simObs6 <- simulateObservations(scns3,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Count = 10,
                                                         Class = "cow"))
  
  simObs6$ageRatioOut %>% filter(Class == "cow", Count != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # number added remains constant rather than number collars
  simObs6$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% 
    {pull(., nadded) == 10} %>% all() %>% 
    expect_true()
  
  # confirm that if freqStartByYear table does skips year it still works
  simObs7 <- simulateObservations(scns,
                                  freqStartsByYear = data.frame(Year = seq(2009, 2023, by = 3),
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Count = 10,
                                                         Class = "cow"))
  
  simObs7$ageRatioOut %>% filter(Class == "cow", Count != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # number added remains constant rather than number collars
  simObs7$simSurvObs %>% group_by(Year) %>%
    summarise(ncollar = n(), ndeaths = sum(event), 
              ndropped = sum(exit == 5 & event == 0),
              nadded = sum(enter == 7)) %>% pull(nadded) %>%
    .[seq(1,15, by = 3)] %>% 
    {. == 10} %>%
    all() %>% 
    expect_true()
})