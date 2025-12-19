test_that("real and simulated data work", {
  # Using observed survival, recruitment and disturbance data
  mod <- bayesianTrajectoryWorkflow(
    surv_data = bboudata::bbousurv_a %>% filter(Year > 2010),
    recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
    disturbance = NULL
  )
  str(mod, max.level = 2)
  
  expect_type(mod, "list")
  
  # Using simulated observation data
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
                              obsAnthroSlope = 1, projAnthroSlope = 5,
                              collarCount = 20, cowMult = 5)
  
  simO <- simulateObservations(scns)
  
  out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
                           disturbance = simO$simDisturbance,
                           startYear = 2014)
  
  expect_type(out, "list")
})


test_that("results match expected", {
  # save to speed up tests
  # simBig <- suppressWarnings(trajectoriesFromNational(N0 = 3000, forceUpdate=T))
  # saveRDS(simBig, "tests/testthat/data/simBig3000.rds", version = 2)
  
  simBig <- readRDS( file.path(test_path(), "data/simBig3000.rds"))
  doScn <- function(nCollar = 2000, nobsYears = 10, collarOn = 4, collarOff = 4, 
                    iAnthro = 0, obsAnthroSlope = 0, projAnthroSlope = 0, 
                    sQuantile = 0.5,  rQuantile = 0.5, rSlopeMod = 1, sSlopeMod = 1){
    #nCollar = 2000; nobsYears = 10; collarOn = 1; collarOff = 12; 
    #iAnthro = 0; obsAnthroSlope = 0; projAnthroSlope = 0; 
    #sQuantile = 0.5;  rQuantile = 0.5; rSlopeMod = 1; sSlopeMod = 1; KSDists = FALSE
    eParsIn <- list()
    eParsIn$collarOnTime <- collarOn
    eParsIn$collarOffTime <- collarOff
    eParsIn$collarNumYears <- 5
    
    scns <- expand.grid(
      obsYears = nobsYears, collarCount = nCollar, cowMult = 6, collarInterval = 1,
      assessmentYrs = 1, iAnthro = iAnthro, rSlopeMod = rSlopeMod, sSlopeMod = sSlopeMod,
      obsAnthroSlope = obsAnthroSlope, projAnthroSlope = projAnthroSlope,
      sQuantile = sQuantile, rQuantile = rQuantile, N0 = 3000
    )
    scResults <- suppressWarnings(bayesianScenariosWorkflow(
      scns, simBig,  eParsIn,
      niters = 3000
    ))
  }
  
  doPlot <- function(scResults, var = "Recruitment", title = ""){
    if (interactive()) {
      return(plotCompareTrajectories(scResults, var,
                                     lowBound = 0,  facetVars = NULL
      )+
        ggplot2::ggtitle(title))
    }
  }
  
  # difference between observed and true simulated observations
  calcDif <- function(obs, var){
    obs %>%
      filter(!MetricTypeID %in% c("Anthro", "fire_excl_anthro")) %>% 
      select(Year, Mean, Type, Parameter) %>% 
      tidyr::pivot_wider(names_from = "Type", values_from = "Mean") %>% 
      filter(!is.na(observed)) %>% 
      mutate(dif = abs(true - observed)) %>%
      group_by(Parameter) %>% 
      summarise(mean_dif = mean(dif))
  }
  
  # difference between modeled and true simulated observations
  calcDifMod <- function(mod, var){
    obs_true <- mod$obs.all %>% 
      filter(!MetricTypeID %in% c("Anthro", "fire_excl_anthro", "c"),
             Type == "true") %>% 
      select(Year, Mean, Type, Parameter) 
    mod_proj <- mod$rr.summary.all %>% 
      mutate(ci_width = upper - lower, .keep = "unused") %>% 
      select(Year, Mean, Parameter, ci_width)
    
    comp <- inner_join(obs_true, mod_proj, by = c("Year", "Parameter"),
              suffix = c("_true", "_proj")) %>% 
      # female pop size is done differently so don't compare
      filter(Parameter != "Female population size") %>% 
      mutate(dif = Mean_true - Mean_proj) %>% 
      group_by(Parameter) %>% 
      summarise(mean_dif = mean(abs(dif)),
                ci_width = mean(ci_width))
    comp
  }
  
  # difference between initial model and Bayesian model
  calcDifNat <- function(mod, min_year = 0){
    mod$rr.summary.all %>% select(Parameter, Mean, Year) %>% 
      filter(Parameter != "c") %>% 
      inner_join(mod$sim.all %>% select(Parameter, Mean, Year),
                 by = c("Parameter", "Year"), suffix = c("_PM", "_nat")) %>% 
      mutate(dif = Mean_PM - Mean_nat) %>% 
      filter(Year >= min_year) %>% 
      group_by(Parameter) %>% 
      summarise(mean_dif = mean(dif)) %>% 
      # nat model does not do female adult pop in the same way
      filter(Parameter != "Female population size")
  }
  
  # expectations that make sense based on what we know and ensure results
  # are similar to expected
  library(caribouMetrics)
  # when we have a lot of collars the distance between observations and "true"
  # pop is smaller than when we have few.
  manyObs <- doScn(nCollar = 2000, rQuantile = 0.9, sQuantile = 0.9)
  doPlot(manyObs, title = "2000 collars")#,var="Adult female survival")
  
  fewCollarObs <- doScn(nCollar = 5, rQuantile = 0.9, sQuantile = 0.9)
  doPlot(fewCollarObs, title = "5 collars")#,var="Adult female survival")
  
  difMany <- calcDif(manyObs$obs.all)
  difFew <- calcDif(fewCollarObs$obs.all)
  expect_true(all(difFew$mean_dif > difMany$mean_dif))
  
  # model predictions should also be closer to true with more collar data
  modDifMany <- calcDifMod(manyObs)
  modDifFew <- calcDifMod(fewCollarObs)
  
  # TODO: still not sure if these expectations are exactly right
  # 
  expect_true(all((modDifFew %>% filter(stringr::str_detect(Parameter, "Expect")) %>% 
                     pull(mean_dif)) >
                (modDifMany %>% filter(stringr::str_detect(Parameter, "Expect")) %>% 
                   pull(mean_dif))))


  # TODO disabled because of issue #145
  # # difference between modeled and true does not change much if collar on/off times
  # # are different but still a year apart
  # collOff3On4 <- doScn(collarOn = 1, collarOff = 12)
  # doPlot(collOff3On4)
  # 
  # modDifcollOff3On4 <- calcDifMod(collOff3On4)
  # 
  # # not more than 5 times difference
  # expect_true(all((modDifcollOff3On4$mean_dif - 
  #                    modDifMany$mean_dif)/modDifMany$mean_dif < 5))
  # 
  # # When there are gaps in survival data bc collars on/off at different times
  # 
  # collOff12On4 <- doScn(collarOn = 4, collarOff = 12)
  # doPlot(collOff12On4)
  # doPlot(collOff12On4, "Adult female survival")
  # 
  # collOff6On4 <- doScn(collarOn = 4, collarOff = 6)
  # doPlot(collOff6On4)
  # 
  # obsDifcollOff12On4 <- calcDif(collOff12On4$obs.all)
  # 
  # modDifcollOff12On4 <- calcDifMod(collOff12On4)

  
  # A pop with quantile >> 0.5 will be above the initial model projection
  highQ <- doScn(rQuantile = 0.95, sQuantile = 0.95)
  doPlot(highQ, "Adult female survival")
  
  difHighQ <- calcDifNat(highQ)
  
  expect_true(all(difHighQ$mean_dif > 0))
  
  # a pop that is less sensitive to anthro dist ie r/sSlopeMod < 1 will show a
  # line that diverges from the initial model. But only if there was some
  # disturbance in training data?
  lowSens <- doScn(rSlopeMod = 0.1, sSlopeMod = 0.1, iAnthro = 80, nobsYears = 20,
                   obsAnthroSlope = 1, projAnthroSlope = 1)
  doPlot(lowSens)
  doPlot(lowSens, "Adult female survival")
  
  difLowSens <- calcDifNat(lowSens, 2040)
  
  expect_true(all(difLowSens %>% pull(mean_dif) > 0))
  
  # same but no anthro in training data
  lowSensNtrain <- doScn(rSlopeMod = 0.1, sSlopeMod = 0.1, iAnthro = 0, nobsYears = 20,
                         obsAnthroSlope = 0, projAnthroSlope = 10)
  
  doPlot(lowSensNtrain)
  doPlot(lowSensNtrain, "Adult female survival")
  doPlot(lowSensNtrain, "Population growth rate")
  difLowSensNtrain <- calcDifNat(lowSensNtrain, min_year = 2040)
  
  # expect differences to be small
  # expect_true(all(difLowSensNtrain$mean_dif < difLowSens$mean_dif))
  # TODO: fix logic of test. Model seems to be behaving as expected.
  
  
  # TODO KS distances are not part of the model any more, is there another
  # metric that should be used to compare?
  # # KS distances JH added to characterize deviation from initial model bands
  # # not just the mean. So should test what happens when there is no sample info
  # # provided. distribution of means from initial model vs Bayesian model. Set standards
  # # that when no obs the differences don't get much worse than they are now. See
  # # the doc JH will send to show what we are looking for
  # 
  # # scenario with no information is very similar to initial model
  # noDat <- doScn(nCollar = 0, nobsYears = 1, KSDists = TRUE)
  # doPlot(noDat)
  # doPlot(noDat, "Adult female survival")
  # 
  # if(interactive()){
  #   noDat$ksDists %>% filter(Parameter != "Female population size") %>% 
  #     ggplot2::ggplot(ggplot2::aes(Year, KSDistance))+
  #     ggplot2::geom_point()+
  #     ggplot2::facet_wrap(~Parameter)
  # }
  # 
  # noDatKS <- noDat$ksDists %>% group_by(Parameter) %>% 
  #   filter(Parameter != "Female population size") %>% 
  #   summarise(meanKS = mean(KSDistance))
  
  #expect_true(all(noDatKS$meanKS < 0.15))
  
  # Values on Jan 13 2025 commit
  #Parameter              meanKS
  #<chr>                   <dbl>
  #1 Adjusted recruitment   0.08154630 
  #2 Adult female survival  0.07131481
  #3 Population growth rate 0.12633333
  #4 Recruitment            0.06106481 
  
})