test_that("testScript still works", {
  # Use Eacker example data for collaring parameters
  eParsIn <- list()
  eParsIn$cowCounts <- data.frame(
    Year = 1981:2023,
    Count = 100,
    Class = "cow"
  )
  eParsIn$freqStartsByYear <- data.frame(
    Year = 1981:2023,
    numStarts = 30
  )
  eParsIn$collarOnTime <- 1
  eParsIn$collarOffTime <- 12
  eParsIn$collarNumYears <- 3

  adjustR <- T # adjust recruitment for delayed age of first reproduction or no.

  ##########
  # Get full set of sims for comparison
  simBig <- suppressWarnings(getSimsNational(adjustR = adjustR)) # If called with default parameters, use saved object to speed things up.

  ###############
  # Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
  #################
  # source("CaribouDemoFns.R")
  # eParsIn$collarNumYears=1

  scns <- expand.grid(
    P = 8, st = 30, cmult = 2, ri = 2, assessmentYrs = 1, iA = 0,
    tA = 0, aS = 0, aSf = 0, sQ = 0.960908218594268, rQ = 0.744425233039074, N0 = 1000
  )
  scResults <- suppressWarnings(runScnSet(scns, eParsIn, simBig, getKSDists = F))

  expect_s3_class(scResults$rr.summary.all, "data.frame")

  if (interactive()) {
    print(plotRes(scResults$rr.summary.all, "Population growth rate",
      obs = scResults$obs.all,
      lowBound = 0, simRange = scResults$sim.all, facetVars = c("P", "sQ")
    ))

    print(plotRes(scResults$rr.summary.all, "Recruitment",
      obs = scResults$obs.all,
      lowBound = 0, simRange = scResults$sim.all, facetVars = c("P", "sQ")
    ))

    print(plotRes(scResults$rr.summary.all, "Adult female survival",
      obs = scResults$obs.all,
      lowBound = 0.65, simRange = scResults$sim.all, facetVars = c("P", "sQ")
    ))

    print(plotRes(scResults$rr.summary.all, "Female population size",
      obs = scResults$obs.all,
      lowBound = 0, highBound = 2000, facetVars = c("P", "sQ")
    ))
  }
})
