# Default parameters for simulating demographic trajectories.

Returns default parameter values for simulating any type of demographic
trajectories. See
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)
or
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

## Usage

``` r
demographyDefaults(
  paramTable = NULL,
  N0 = 1000,
  qMin = 0,
  qMax = 0,
  uMin = 0,
  uMax = 0,
  zMin = 0,
  zMax = 0,
  cowMult = 6,
  lQuantile = NA,
  correlateRates = F,
  ...
)
```

## Arguments

- paramTable:

  a data.frame with column names matching the arguments below. Any
  columns that are missing will be filled with the default values.

- N0:

  number. Initial number of adult females.

- qMin:

  number in 0, 1. Minimum ratio of bulls to cows in composition survey
  groups.

- qMax:

  number in 0, 1. Maximum ratio of bulls to cows in composition survey
  groups.

- uMin:

  number in 0, 1. Minimum probability of misidentifying young bulls as
  adult females and vice versa in composition survey.

- uMax:

  number in 0, 1. Maximum probability of misidentifying young bulls as
  adult females and vice versa in composition survey.

- zMin:

  number in 0, 1. Minimum probability of missing calves in composition
  survey.

- zMax:

  number in 0, \<1. Maximum probability of missing calves in composition
  survey.

- cowMult:

  number \>= 1. The apparent number of adult females per collared animal
  in composition survey. Set to NA to use `cowCount`.

- lQuantile:

  number in 0, 1. Lambda quantile

- correlateRates:

  logical. Set TRUE to force correlation between recruitment and
  survival.

## Value

a data.frame of parameter values.

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
[`disturbanceDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceDefaults.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`monitoringDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/monitoringDefaults.md),
[`nationalTrajectoryDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/nationalTrajectoryDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`timeDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/timeDefaults.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
demographyDefaults()
#> # A tibble: 1 × 10
#>      N0  qMin  qMax  uMin  uMax  zMin  zMax cowMult lQuantile correlateRates
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl> <lgl>     <lgl>         
#> 1  1000     0     0     0     0     0     0       6 NA        FALSE         
```
