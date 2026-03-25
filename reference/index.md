# Package index

## Overview

Overview of the package.

- [`caribouMetrics`](https://landscitech.github.io/caribouMetrics/reference/caribouMetrics-package.md)
  [`caribouMetrics-package`](https://landscitech.github.io/caribouMetrics/reference/caribouMetrics-package.md)
  : caribouMetrics: Models and Metrics of Boreal Caribou Demography and
  Habitat Selection

## Caribou demography

Functions for predicting boreal caribou demographic rates and projecting
population growth

- [`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianScenariosWorkflow.md)
  : Run the Bayesian population model for multiple parameter sets
- [`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianTrajectoryWorkflow.md)
  : Bayesian population model for boreal caribou
- [`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/reference/betaNationalPriors.md)
  : Get prior parameters for Bayesian Beta demographic rate models
- [`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/reference/caribouPopGrowth.md)
  : Caribou demographic model
- [`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/compareTrajectories.md)
  : Summarize results of Bayesian demographic model in tables
- [`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/reference/compositionBiasCorrection.md)
  : Calculate bias correction term for calf:cow composition survey.
- [`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/reference/dataFromSheets.md)
  : Demographic data from Google sheet
- [`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/reference/demographicProjectionApp.md)
  : Run the Bayesian caribou demographic projection app
- [`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateBayesianRates.md)
  : Create summary table of demographic rates from survival and
  recruitment surveys
- [`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/reference/estimateNationalRates.md)
  [`estimateNationalRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateNationalRates.md)
  : Sample demographic rates
- [`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md)
  [`sampleNationalCoefs()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md)
  [`subsetNationalCoefs()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md)
  : Sample demographic regression model coefficients
- [`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/reference/getScenarioDefaults.md)
  : Default parameters for simulation of example demographic
  trajectories.
- [`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotCompareTrajectories.md)
  : Plot Bayesian population model results
- [`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/reference/plotSurvivalSeries.md)
  : Plot survival time series
- [`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotTrajectories.md)
  : TO DO: document plotTrajectories
- [`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md)
  : Population growth model table for Johnson models
- [`simulateObservations()`](https://landscitech.github.io/caribouMetrics/reference/simulateObservations.md)
  : Simulate observations
- [`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md)
  [`summarizeTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md)
  [`simulateTrajectoriesFromPosterior()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md)
  : Format trajectory tables
- [`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromBayesian.md)
  : Get trajectories from a Bayesian model result
- [`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromNational.md)
  : Get a set of simulation results from the national demographic model
- [`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummary.md)
  : Projections of population growth from demographic model summaries.
- [`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummaryForApp.md)
  : Demographic projections for cases with no change in demographic
  rates over time. This is the method used (so far) in the demography
  app. TO DO: Consider removing and replacing with call to
  trajectoriesFromSummary.

## Summarizing disturbance

Functions for determining the percentage of a landscape that has been
disturbed

- [`DisturbanceMetrics-class`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md)
  [`DisturbanceMetrics`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md)
  : Disturbance Metrics class
- [`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/reference/disturbanceMetrics.md)
  : Calculate metrics of natural and anthropogenic disturbance
- [`loadSpatialInputs()`](https://landscitech.github.io/caribouMetrics/reference/loadSpatialInputs.md)
  : Load Spatial Input Data
- [`rasterizeLineDensity()`](https://landscitech.github.io/caribouMetrics/reference/rasterizeLineDensity.md)
  : Rasterize line density
- [`reclassDist()`](https://landscitech.github.io/caribouMetrics/reference/reclassDist.md)
  : Reclassify natural disturbance and harvest layers
- [`results()`](https://landscitech.github.io/caribouMetrics/reference/results.md)
  : Extract results
- [`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md)
  : Update an existing DisturbanceMetrics Object
