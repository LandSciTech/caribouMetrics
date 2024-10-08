test_that("testScript still works", {
  eParsIn <- list()
  eParsIn$collarOnTime <- 1
  eParsIn$collarOffTime <- 12
  eParsIn$collarNumYears <- 3

  ##########
  # Get full set of sims for comparison
  simBig <- suppressWarnings(getSimsNational()) # If called with default parameters, use saved object to speed things up.

  ###############
  # Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
  #################
  # source("CaribouDemoFns.R")
  # eParsIn$collarNumYears=1

  scns <- expand.grid(
    obsYears = 8, collarCount = 30, cowMult = 2, collarInterval = 2,
    assessmentYrs = 1, iAnthro = 0,
    tA = 0, obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.960908218594268, 
    rQuantile = 0.744425233039074, N0 = 1000
  )
  scResults <- suppressWarnings(runScnSet(scns, eParsIn, simBig, getKSDists = F,
                                          Niter = 100, Nburn = 10))

  expect_s3_class(scResults$rr.summary.all, "data.frame")

  if (interactive()) {
    print(plotRes(scResults$rr.summary.all, "Population growth rate",
      obs = scResults$obs.all,
      lowBound = 0, simRange = scResults$sim.all, facetVars = c("obsYears", "sQuantile")
    ))

    print(plotRes(scResults$rr.summary.all, "Recruitment",
      obs = scResults$obs.all,
      lowBound = 0, simRange = scResults$sim.all, facetVars = c("obsYears", "sQuantile")
    ))

    print(plotRes(scResults$rr.summary.all, "Adult female survival",
      obs = scResults$obs.all,
      lowBound = 0.65, simRange = scResults$sim.all, facetVars = c("obsYears", "sQuantile")
    ))

    print(plotRes(scResults$rr.summary.all, "Female population size",
      obs = scResults$obs.all,
      lowBound = 0, highBound = 2000, facetVars = c("obsYears", "sQuantile")
    ))
  }
})
