# Dynamically simulate a trajectory over time based on the national model

Dynamically simulate a trajectory over time based on the national model

## Usage

``` r
simTrajectory(
  numYears,
  covariates,
  survivalModelNumber = "M1",
  recruitmentModelNumber = "M4",
  popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
  recSlopeMultiplier = 1,
  sefSlopeMultiplier = 1,
  rQuantile = NULL,
  sQuantile = NULL,
  stepLength = 1,
  N0 = 1000,
  cowMult = 1,
  qMin = 0,
  qMax = 0,
  uMin = 0,
  uMax = 0,
  zMin = 0,
  zMax = 0,
  interannualVar = eval(formals(caribouPopGrowth)$interannualVar)
)
```

## Arguments

- numYears:
- covariates:
- survivalModelNumber:
- recruitmentModelNumber:
- popGrowthTable:
- recSlopeMultiplier:
- sefSlopeMultiplier:
- rQuantile:
- sQuantile:
- stepLength:
- N0:
- cowMult:
- qMin:
- qMax:
- uMin:
- uMax:
- zMin:
- zMax:
- interannualVar:

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`bbouNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/bbouNationalPriors.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)
