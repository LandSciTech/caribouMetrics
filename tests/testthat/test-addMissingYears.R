test_that("missing years works", {
  zeroAdd <- addMissingYears(bboudata::bbousurv_a, integer(0))
 
  # when there are incomplete Caribou Years in the data this adds rows for those
  # months
  expect_equal(zeroAdd %>% filter(!is.na(MortalitiesCertain)),
          as.data.frame(bboudata::bbousurv_a))

  completeSurvData <- bboudata::bbousurv_a %>% getCaribouYear() %>% 
    group_by(CaribouYear) %>% filter(sum(!is.na(MortalitiesCertain)) == 12) %>% 
    ungroup() %>% 
    select(-CaribouYear)
  
  zeroAdd2 <- addMissingYears(completeSurvData, integer(0))
  
  # Nothing added when all years are complete
  expect_equal(zeroAdd2, as.data.frame(completeSurvData))
  
  add10S <- addMissingYears(completeSurvData, 2017:2026)
  
  expect_true(nrow(add10S) == 12*10+nrow(completeSurvData))
  
  # includes Jan-Mar 2027 because based on CaribouYear
  expect_equal(add10S %>% filter(Year == 2027) %>% pull(Month) %>% max(), 3)
  
  
  zeroAddR <- addMissingYears(bboudata::bbourecruit_a, integer(0))
  
  expect_equal(zeroAddR, as.data.frame(bboudata::bbourecruit_a))
  
  add10R <- addMissingYears(bboudata::bbourecruit_a, 2017:2026)

  # should have rec survey in 2027 since 2026 caribou year was added. 
  expect_equal(add10R %>% filter(Year == 2027) %>% nrow(), 1)
  
  # error from some simulated and some observed in same 
  add16S <- addMissingYears(bboudata::bbousurv_a %>% filter(Year > 2010), 2016) %>% getCaribouYear()
  add16R <- addMissingYears(bboudata::bbourecruit_a %>% filter(Year > 2010), 2016)
  disturbance <-  data.frame(Year = unique(add16S$CaribouYear), Anthro = 3, Fire_excl_anthro = 5)
  rates <- estimateBayesianRates(add16S, add16R,disturbance = disturbance, niters = 20)
  expect_is(rates, "list")
})
