test_that("works without disturbance", {
  mod_flc <- here::here("results/test_bbou_nodist_nomonitor.rds")
  if (file.exists(mod_flc)) {
    bbouNoDist <- readRDS(mod_flc)
  } else {
    bbouNoDist <- estimateBayesianRates(bboudata::bbousurv_a %>% filter(Year > 2010),
                                       bboudata::bbourecruit_a %>% filter(Year > 2010),
                                       N0 = NA, return_mcmc = T, niters = 3000
    )
    if(dir.exists(dirname(mod_flc))){
      saveRDS(bbouNoDist, mod_flc)
    }
  }
  
  trajNoDist <- trajectoriesFromBayesian(bbouNoDist)
  
  pltNoDist <- plotTrajectories(trajNoDist)
  
  expect_is(pltNoDist, "ggplot2::ggplot")
})

test_that("works with disturbance and N0 variation", {
  disturbance <- unique(subset(bboudata::bbourecruit_a %>% filter(Year > 2010), 
                               select = Year))
  disturbance$Anthro <- 0
  disturbance$Fire_excl_anthro <- 0
  
  mod_fld <- here::here("results/test_bbou_dist.rds")
  if (file.exists(mod_fld)) {
    bbouDist <- readRDS(mod_fld)
  } else {
    bbouDist <- estimateBayesianRates(bboudata::bbousurv_a %>% filter(Year > 2010),
                                      bboudata::bbourecruit_a %>% filter(Year > 2010),
                                      N0 = NA, return_mcmc = T, niters = 3000,
                                      disturbance = disturbance)
    if(dir.exists(dirname(mod_fld))){
      saveRDS(bbouDist, mod_fld)
    }
  }
  
  trajDist <- trajectoriesFromBayesian(bbouDist)
  
  expect_type(trajDist, "list")
  
  pltDist <- plotTrajectories(trajDist)
  
  expect_is(pltDist, "ggplot2::ggplot")
  
  # with variation in N0
  trajDistwN0 <- trajectoriesFromBayesian(bbouDist, 
                                          N0 = data.frame(N0 = 125,
                                                          N.lower = 100,
                                                          N.upper = 150))
  
  expect_type(trajDistwN0, "list")
  
  pltDistwN0 <- plotTrajectories(trajDistwN0)
  
  expect_is(pltDistwN0, "ggplot2::ggplot")
})

test_that("can project into future with missing years in data", {
  bbouInformativeFile <- here::here("results/vignetteBbbouExample.rds")
  
  surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
  surv_data_add <- expand.grid(Year = seq(2017, 2022), Month = seq(1:12),
                               PopulationName = unique(surv_data$PopulationName))
  surv_data <- merge(surv_data, surv_data_add, all.x = TRUE, all.y = TRUE)
  
  recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
  recruit_data_add <- expand.grid(Year = seq(2017, 2022), 
                                  PopulationName = unique(recruit_data$PopulationName))
  recruit_data <- merge(recruit_data, recruit_data_add, all.x = TRUE, all.y = TRUE)
  
  if (file.exists(bbouInformativeFile)) {
    bbouInformative <- readRDS(bbouInformativeFile)
  } else {
    bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
                                             return_mcmc = TRUE)
    if (dir.exists(dirname(bbouInformativeFile))) {
      saveRDS(bbouInformative, bbouInformativeFile)
    }
  }
  
  expect_true(all(seq(2017, 2022) %in% bbouInformative$parList$Rbar$Year))
  
  bbouTraj <- trajectoriesFromBayesian(bbouInformative)
  
  expect_true(all(seq(2017, 2022) %in% unique(bbouTraj$summary$Year)))
})
