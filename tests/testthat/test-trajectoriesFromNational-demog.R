test_that("cacheing happens", {
  # clear the cache
  rm(list = ls(envir = cacheEnv), envir = cacheEnv)
  expect_message((trajectoriesFromNational()), "Updating cached")
  expect_message(trajectoriesFromNational(), "saved object")
})

test_that("sample trajectories are not returned when the national model is used", {
  noDist <- trajectoriesFromNational(forceUpdate = TRUE)
  
  expect_null(noDist$samples)
})

test_that("Can set constant Anthro", {
  scnsConst <- getScenarioDefaults()
  distScn <- data.frame(Anthro = 20, Fire_excl_anthro = 0, Year = 1:10+2000)
  scnsConst <- merge(scnsConst, distScn)
  scnsConst$iAnthro <- NULL
  distConst <- trajectoriesFromNational(cPars = scnsConst, forceUpdate = TRUE)$summary
  
  expect_equal(unique(distConst$AnthroID), 20)
})


test_that("can specify multiple disturbance scenarios", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                 projYears = 100,iAnthro=c(0,5))
  
  summary2 <- trajectoriesFromNational(replicates = 2, cPars = scns10m)$summary
  
  # The first year will have both values for Anthro
  summary2 %>% filter(Year == min(Year), MetricTypeID == "lambda") %>% pull(AnthroID) %>% 
    range() %>% 
    expect_equal(c(0,5))
})

test_that("get samples by default and projects over time when Year supplied",{
  disturbance <- data.frame(Anthro = 40, Fire_excl_anthro = 2)
  disturbance <-  data.frame(step = 0:4) %>% bind_cols(disturbance) %>% 
    mutate(Anthro = Anthro + 10 * step, 
           Year = step * 10+2000)
  
  wDist <- trajectoriesFromNational(disturbance = disturbance, replicates = 35, 
                                    interannualVar = FALSE, useQuantiles = TRUE,
                                    N0 = 100, numSteps = 10)
  expect_named(wDist, c("summary", "samples"))
  
  if(F&interactive()){
    proj <- ggplot(data = wDist$samples, aes(x = Timestep, y = Amount, colour = Replicate,
                                             group = Replicate)) +
      geom_line() +
      facet_wrap(~MetricTypeID, scales = "free") +
      xlab("Time") +
      theme(legend.position = "none")
    proj
  }
})
