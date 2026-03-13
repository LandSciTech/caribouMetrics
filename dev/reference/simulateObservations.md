# Simulate observations

Parameters specify a monitoring program that is applied to simulate
observations from the example trajectories. Parameters for the caribou
monitoring program, disturbance scenario and the true population
trajectory can be specified with
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md).

## Usage

``` r
simulateObservations(
  paramTable,
  trajectories = NULL,
  cowCounts = NULL,
  freqStartsByYear = NULL,
  collarNumYears = 4,
  collarOffTime = 3,
  collarOnTime = 4,
  caribouYearStart = 4,
  topUp = F,
  recSurveyMonth = 3,
  recSurveyDay = 15,
  distScen = NULL,
  writeFilesDir = NULL,
  surv_data = NULL,
  recruit_data = NULL
)
```

## Arguments

- paramTable:

  data.frame. Parameters for the simulations. See
  [`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md)
  for details.

- trajectories:

  data.frame. Optional example demographic trajectory. If NULL the
  trajectory will be simulated from the national model.

- cowCounts:

  data.frame. Optional. Number of cows counted in aerial surveys each
  year. If NULL, and `paramTable` contains `cowMult` the number of cows
  that survive calving based on the collar data is multiplied by
  `cowMult` to determine the number of cows counted in aerial surveys.
  If `paramTable` does not contain `cowMult` `paramTable$cowCount` is
  used to set the number of cows counted in aerial surveys each year. If
  a data.frame is provided it must have columns "Year" and "Cows".

- freqStartsByYear:

  data.frame. Optional. Number of collars deployed in each year. If NULL
  `paramTable$collarCount` is used as the target number of collars and
  each year that collars are deployed they will be topped up to this
  number. If a data.frame is provided it must have 2 columns "Year" and
  "numStarts" or "numTarget" (but not both). "numStarts" is the absolute
  number of collars deployed in that year, and "numTarget" is the target
  number of collars.

- collarNumYears:

  integer. Number of years until collar falls off

- collarOffTime:

  integer. Month that collars fall off. A number from 1 (January) to 12
  (December)

- collarOnTime:

  integer. Month that collars are deployed. A number from 1 (January) to
  12 (December)

- caribouYearStart:

  integer. The first month of the year for caribou.

- recSurveyMonth:

  integer. The month of simulated recruitment surveys.

- recSurveyDay:

  integer. The day for simulated recruitment surveys.

- distScen:

  data.frame. Disturbance scenario. Must have columns "Year", "Anthro",
  and "Fire_excl_anthro" containing the year, percentage of the
  landscape covered by anthropogenic disturbance buffered by 500 m, and
  the percentage covered by fire that does not overlap anthropogenic
  disturbance. See
  [`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceMetrics.md).
  If NULL the disturbance scenario is simulated based on `paramTable`

- writeFilesDir:

  character. If not NULL `simSurvObs` and `simRecruitObs` results will
  be saved to csv files in the directory provided

- surv_data:

  data.frame. Optional existing survival data in bboudata format. Will
  be combined with simulated data if ... Otherwise ignored.

- recruit_data:

  data.frame. Optional existing recruitment data in bboudata format.
  Will be combined with simulated data if ... Otherwise ignored.

## Value

a list with elements:

- minYr: first year in the simulations,

- maxYr: last year in the simulations,

- simDisturbance: a data frame with columns Anthro, Fire_excl_anthro,
  Total_dist, and Year,

- simSurvObs: a data frame of survival data in bboutools format,

- simRecruitObs: a data frame of recruitment data in bboutools format,

- exData: a tibble of expected population metrics based on the initial
  model,

- paramTable: a data frame recording the input parameters for the
  simulation.

## Details

For a detailed description of the process for simulating data see the
[vignette](https://landscitech.github.io/caribouMetrics/articles/bayesian-model-outputs.html#simulation-of-local-population-dynamics-and-monitoring)
([`vignette("bayesian-model-outputs", package = "caribouMetrics")`](https://landscitech.github.io/caribouMetrics/articles/bayesian-model-outputs.html))
and [Hughes et al. 2025](https://doi.org/10.1016/j.ecoinf.2025.103095).

## References

Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
Integration of national demographic-disturbance relationships and local
data can improve caribou population viability projections and inform
monitoring decisions. Ecological Informatics, 87, p.103095.
<https://doi.org/10.1016/j.ecoinf.2025.103095>

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
                            collarCount = 20, cowMult = 5)

simO <- simulateObservations(scns)
```
