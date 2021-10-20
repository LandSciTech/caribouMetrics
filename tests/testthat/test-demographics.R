
# example disturbance table
covTableSim <- expand.grid(Anthro = seq(1, 90, by = 20), 
                           Fire = seq(1, 70, by = 20)) 
covTableSim$polygon = paste0("Missisa_", covTableSim$Anthro + 1)
covTableSim$area = "FarNorth"
covTableSim$Total_dist = covTableSim$Anthro + covTableSim$Fire
covTableSim$fire_excl_anthro = covTableSim$Fire - 1
covTableSim$fire_prop_dist = covTableSim$Fire/covTableSim$Total_dist
# these are just dummies for testing
covTableSim$ln_nn <- 1
covTableSim$hqh <- 1

test_that("basic example works", {
  demCoefs <- demographicCoefficients(replicates = 10)

  demRates <- demographicRates(covTable = covTableSim,
                               popGrowthPars = demCoefs)
  
  expect_equal(nrow(demRates), nrow(covTableSim))
  
  demRates <- demographicRates(covTable = covTableSim,
                               popGrowthPars = demCoefs,
                               returnSample = TRUE)
  
  expect_equal(nrow(demRates), nrow(covTableSim)*10)
  
})

test_that("useQuantiles works as expected", {
  demCoefswQ <- demographicCoefficients(replicates = 10)
  demCoefsnQ <- demographicCoefficients(replicates = 10, useQuantiles = FALSE)
  
  expect_gt(length(demCoefswQ$coefSamples_Survival), 
            length(demCoefsnQ$coefSamples_Survival))
  
  # custom quantiles
  demCoefsCustomQ <- demographicCoefficients(replicates = 10, useQuantiles = c(0.2, 0.8))
  
  expect_gt(demCoefswQ$coefSamples_Survival$quantiles %>% max(),
            demCoefsCustomQ$coefSamples_Survival$quantiles %>% max())
  
  # in demRates
  demRates1 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefswQ,
                                returnSample = TRUE,
                                ignorePrecision = FALSE, 
                                useQuantiles = FALSE)
  
  demRates2 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefsnQ,
                                returnSample = TRUE, 
                                ignorePrecision = FALSE,
                                useQuantiles = FALSE)
  
  demRates3 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefswQ,
                                returnSample = TRUE,
                                ignorePrecision = FALSE,
                                useQuantiles = TRUE)
  
  demRates4 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefsnQ,
                                returnSample = TRUE,
                                ignorePrecision = FALSE,
                                useQuantiles = TRUE)
  
  # if quantiles are used then the replicate is a predictor of S_bar
  lmQF <- lm(S_bar~Anthro+replicate, data = demRates1) %>% summary()

  lmQT <- lm(S_bar~Anthro+replicate, data = demRates3) %>% summary()
  
  lmQF2 <- lm(S_bar~Anthro+replicate, data = demRates2) %>% summary()
  
  lmQT2 <- lm(S_bar~Anthro+replicate, data = demRates4) %>% summary()
  
  # expect replicates signif if quantiles used
  expect_false(sum(lmQF[["coefficients"]][, 4] < 0.01) > 5)
  expect_true(sum(lmQT[["coefficients"]][, 4] < 0.01) > 5)
  expect_false(sum(lmQF2[["coefficients"]][, 4] < 0.01) > 5)
  expect_true(sum(lmQT2[["coefficients"]][, 4] < 0.01) > 5)
  
  # Try using different quantiles
  # in demRates
  demRates5 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefsnQ,
                                returnSample = TRUE,
                                useQuantiles = c(0.001, 0.999))
  
  expect_gt(max(demRates5$S_bar), max(demRates4$S_bar))
  
  # from demCoefs
  # should override quantiles added in demRates
  expect_warning(demRates6 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefsCustomQ,
                                returnSample = TRUE,
                                useQuantiles = c(0.001, 0.999)))
  
  expect_gt(max(demRates4$S_bar), max(demRates6$S_bar))
  
  demRates7 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefsCustomQ,
                                returnSample = TRUE)
  
  expect_gt(max(demRates4$S_bar), max(demRates6$S_bar))

  # currently causes error see issue #74
  
})

test_that("ignorePrecision works as expected", {
  demCoefs <- demographicCoefficients(replicates = 10)
  
  demRates1 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefs,
                                ignorePrecision = TRUE)
  demRates2 <- demographicRates(covTable = covTableSim,
                                popGrowthPars = demCoefs,
                                ignorePrecision = FALSE)
  
  expect_true(all(demRates2$S_stdErr > demRates1$S_stdErr))
})

test_that("all model versions in table work", {
  versJ <- distinct(popGrowthTableJohnsonECCC, modelVersion, ModelNumber) %>%
    filter(modelVersion == "Johnson")
  
  allJmods <- purrr::map2(versJ$modelVersion, versJ$ModelNumber,
                          ~demographicRates(covTableSim, 
                                            demographicCoefficients(10, modelVersion = .x, 
                                                                    survivalModelNumber = .y, 
                                                                    recruitmentModelNumber = .y),
                                            ignorePrecision = TRUE))
  expect_equal(nrow(versJ), length(allJmods))
  
  # only 1 femaleSurvival model
  versEC <- popGrowthTableJohnsonECCC %>% 
    filter(modelVersion == "ECCC", responseVariable == "recruitment") %>% 
    distinct(ModelNumber)
  
  allECmods <- purrr::map(versEC$ModelNumber,
                          ~demographicRates(covTableSim, 
                                            demographicCoefficients(10, modelVersion = "ECCC", 
                                                                    survivalModelNumber = "M1", 
                                                                    recruitmentModelNumber = .x),
                                            ignorePrecision = TRUE))
})

test_that("demoCoefs has reasonable errors", {
  # wrong model #
  expect_error(demographicCoefficients(10, modelVersion = "ECCC", 
                                         survivalModelNumber = "M2"), 
                "Model not available")

  
  # multiple models getCoefs can take multiple models but demCoefs can't
  expect_error(demographicCoefficients(10, modelVersion = c("ECCC", "Johnson"), 
                          survivalModelNumber = c("M1", "M2"), 
                          recruitmentModelNumber = c("M3", "M2")),
               "Multiple models")
})

test_that("demoRates has reasonable errors",{
  expect_error(demographicRates(rename(covTableSim, ant = Anthro),
                   demographicCoefficients(10)),
               "Covariates missing")
  
  expect_error(demographicRates(covTableSim,
                                demographicCoefficients(10, modelVersion = "ECCC", 
                                                        recruitmentModelNumber = "M7")),
               "Missing precision")
})

# Compare output to previous run. This will raise a flag if the result has
# changed. Update the stored result if the change was expected.
resultCompare <- readRDS(file.path("data", "demog_resultCompare.rds"))

demCoefs <- demographicCoefficients(replicates = 10)

demRates <- demographicRates(covTable = covTableSim,
                             popGrowthPars = demCoefs)

# To update
# saveRDS(demRates, file.path("tests/testthat/data", "demog_resultCompare.rds"))

testthat::test_that("results match previous results",{
  testthat::expect_equal(demRates$S_bar, resultCompare$S_bar)
  testthat::expect_equal(demRates$R_bar, resultCompare$R_bar)
})
