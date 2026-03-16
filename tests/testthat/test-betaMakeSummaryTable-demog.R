# set values that work for test and are fast
niters <- 100
nthin <- 5
nc <- 3      # number of chains
ni <- niters * nthin * 2   # number of samples for each chain
nb <- ni / 2 

test_that("Basic inputs works", {
  out <- betaMakeSummaryTable(
    surv_data = bboudata::bbousurv_a %>% filter(Year > 2010),
    recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
    disturbance = data.frame(Year = 2010:2020, Anthro = 10:20, Fire_excl_anthro = 10:20),
    priors = betaNationalPriors(), 
    nc, nthin, ni, nb
  )
  expect_type(out, "list")
})

test_that("multiple populations works", {
 # disturbance doesn't have PopName so assumed to apply to both
    out <- betaMakeSummaryTable(
      surv_data = bboudata::bbousurv_a %>% bind_rows(bboudata::bbousurv_b) %>% filter(Year > 2010),
      recruit_data = bboudata::bbourecruit_a %>% bind_rows(bboudata::bbourecruit_b) %>% filter(Year > 2010),
      disturbance = data.frame(Year = 2010:2020, Anthro = 10:20, Fire_excl_anthro = 10:20),
      priors = betaNationalPriors(), 
      nc, nthin, ni, nb
    )
  
  expect_warning(
    out2 <- betaMakeSummaryTable(
      surv_data = bboudata::bbousurv_a %>% bind_rows(bboudata::bbousurv_b) %>% filter(Year > 2010),
      recruit_data = bboudata::bbourecruit_a %>% bind_rows(bboudata::bbourecruit_b) %>% filter(Year > 2010),
      disturbance = bind_rows(
        data.frame(PopulationName = "A", Year = 2010:2020, Anthro = 10:20, Fire_excl_anthro = 10:20),
        data.frame(PopulationName = "B", Year = 2010:2020, Anthro = 10:20, Fire_excl_anthro = 10:20)
      ),
      priors = betaNationalPriors(), 
      nc, nthin, ni, nb
    ),
    "no data for population"
  )
})

test_that("works with simulated data",{
  simO <- simulateObservations(getScenarioDefaults(collarCount = 10))
  
  out <- betaMakeSummaryTable(surv_data = simO$simSurvObs, 
                              recruit_data = simO$simRecruitObs, 
                              disturbance = simO$simDisturbance, 
                              priors = betaNationalPriors(), 
                              nc, nthin, ni, nb)
  expect_type(out, "list")
})
