# Default scenario parameters.

`getScenarioDefaults()`: Returns default parameters for scenarios. Use
this function to get a combination of `disturbanceDefaults()`,
`timeDefaults()`, `demographyDefaults()`,
`nationalTrajectoryDefaults()`, and `monitoringDefaults()`. If only one
of these sets of parameters is needed consider using the relevant
component function instead.

`timeDefaults():` Returns default parameter values for scenario
durations. See
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

`disturbanceDefaults()`: Returns default parameter values for
disturbance scenarios. See
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

`demographyDefaults()`: Returns default parameter values for simulating
any type of demographic trajectories. See
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)
or
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

`nationalTrajectoryDefaults()`: Returns default parameter values for
national demographic trajectories. See
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md)
and
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

`monitoringDefaults`: Returns default parameter values for monitoring.
See
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
and
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md)
for additional details.

## Usage

``` r
getScenarioDefaults(
  paramTable = NULL,
  includeDist = T,
  includeTime = T,
  includeDemography = T,
  includeNational = T,
  includeMonitoring = T,
  ...
)

timeDefaults(
  paramTable = NULL,
  projYears = 35,
  obsYears = 15,
  preYears = 0,
  curYear = 2023,
  startYear = NA,
  ...
)

disturbanceDefaults(
  paramTable = NULL,
  iFire = 0,
  iAnthro = 0,
  obsAnthroSlope = 2,
  projAnthroSlope = 2,
  ...
)

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

nationalTrajectoryDefaults(
  paramTable = NULL,
  sQuantile = NA,
  rQuantile = NA,
  rSlopeMod = 1,
  sSlopeMod = 1,
  interannualVar = list(eval(formals(caribouPopGrowth)$interannualVar)),
  ...
)

monitoringDefaults(
  paramTable = NULL,
  collarInterval = NA,
  cowCount = NA,
  collarCount = NA,
  ...
)
```

## Arguments

- paramTable:

  a data.frame with column names matching the arguments below. Any
  columns that are missing will be filled with the default values.

- includeDist:

  logical. Include `disturbanceDefaults()`?

- includeTime:

  logical. Include `timeDefaults()`?

- includeDemography:

  logical. Include `demographyDefaults()`?

- includeNational:

  logical. Include `nationalTrajectoryDefaults()`?

- includeMonitoring:

  logical. Include `monitoringDefaults()`?

- projYears:

  Number of years of projections

- obsYears:

  Number of years of observations

- preYears:

  Number of years before monitoring begins

- curYear:

  year. The current year. All years before are part of the observation
  period and years after are part of the projection period.

- startYear:

  year. First year in observation period. Optional, if not provided it
  will be calculated from `curYear` and `obsYears`

- iFire:

  number. Initial fire disturbance percentage.

- iAnthro:

  number. Initial anthropogenic disturbance percentage

- obsAnthroSlope:

  number. Percent change in anthropogenic disturbance per year in the
  observation period

- projAnthroSlope:

  number. Percent change in anthropogenic disturbance per year in the
  projection period

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

- sQuantile:

  number in 0,1. Survival quantile.

- rQuantile:

  number in 0,1. Recruitment quantile.

- rSlopeMod:

  number. Disturbance-recruitment slope multiplier

- sSlopeMod:

  number. Disturbance-survival slope multiplier

- interannualVar:

  list or logical. List containing interannual variability parameters.
  These can be either coefficients of variation (R_CV, S_CV), beta
  precision parameters (R_phi, S_phi), or random effects parameters from
  a logistic glmm (R_annual, S_annual). Set to `FALSE` to ignore
  interannual variability.

- collarInterval:

  number. Optional. Number of years between collar deployments. If
  missing assumed to be every year

- cowCount:

  Optional. Only used in
  [`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md)
  to set the number of cows per year in recruitment survey

- collarCount:

  number \>= 1. The target number of collars active each year. Set to NA
  to use `freqStartsPerYear` in
  [`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)

## Value

a data.frame of parameter values and for `getScenarioDefaults()`, a
label column that combines all the parameter names and values into a
string

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
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
getScenarioDefaults()
#> # A tibble: 1 × 24
#>      N0  qMin  qMax  uMin  uMax  zMin  zMax cowMult correlateRates rSlopeMod
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl> <lgl>              <dbl>
#> 1  1000     0     0     0     0     0     0       6 FALSE                  1
#> # ℹ 14 more variables: sSlopeMod <dbl>, interannualVar <list>, iFire <dbl>,
#> #   iAnthro <dbl>, obsAnthroSlope <dbl>, projAnthroSlope <dbl>, hasYear <lgl>,
#> #   projYears <dbl>, obsYears <dbl>, preYears <dbl>, curYear <dbl>,
#> #   startYear <dbl>, ID <int>, label <chr>

# paramTable list takes precedence over argument values
getScenarioDefaults(paramTable = data.frame(iFire = 10, iAnthro = 20, obsYears = 1), obsYears = 5)
#>     N0 qMin qMax uMin uMax zMin zMax cowMult correlateRates rSlopeMod sSlopeMod
#> 1 1000    0    0    0    0    0    0       6          FALSE         1         1
#>     interannualVar iFire iAnthro obsAnthroSlope projAnthroSlope hasYear
#> 1 0.46000, 0.08696    10      20              2               2   FALSE
#>   projYears obsYears preYears curYear startYear ID
#> 1        35        1        0    2023      2023  1
#>                                                                                                                                                                                                                                                                          label
#> 1 ID1_startYear2023_curYear2023_preYears0_obsYears1_projYears35_hasYearFALSE_projAnthroSlope2_obsAnthroSlope2_iAnthro20_iFire10_interannualVarlist(R_CV = 0.46, S_CV = 0.08696)_sSlopeMod1_rSlopeMod1_correlateRatesFALSE_cowMult6_zMax0_zMin0_uMax0_uMin0_qMax0_qMin0_N01000_

timeDefaults()
#> # A tibble: 1 × 6
#>   hasYear projYears obsYears preYears curYear startYear
#>   <lgl>       <dbl>    <dbl>    <dbl>   <dbl>     <dbl>
#> 1 FALSE          35       15        0    2023      2009

disturbanceDefaults()
#> # A tibble: 1 × 10
#>   iFire iAnthro obsAnthroSlope projAnthroSlope hasYear projYears obsYears
#>   <dbl>   <dbl>          <dbl>           <dbl> <lgl>       <dbl>    <dbl>
#> 1     0       0              2               2 FALSE          35       15
#> # ℹ 3 more variables: preYears <dbl>, curYear <dbl>, startYear <dbl>

demographyDefaults()
#> # A tibble: 1 × 10
#>      N0  qMin  qMax  uMin  uMax  zMin  zMax cowMult lQuantile correlateRates
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl> <lgl>     <lgl>         
#> 1  1000     0     0     0     0     0     0       6 NA        FALSE         

nationalTrajectoryDefaults()
#> # A tibble: 1 × 13
#>   rSlopeMod sSlopeMod interannualVar      N0  qMin  qMax  uMin  uMax  zMin  zMax
#>       <dbl>     <dbl> <list>           <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1         1         1 <named list [2]>  1000     0     0     0     0     0     0
#> # ℹ 3 more variables: cowMult <dbl>, lQuantile <lgl>, correlateRates <lgl>

monitoringDefaults()
#> # A tibble: 1 × 10
#>      N0  qMin  qMax  uMin  uMax  zMin  zMax cowMult lQuantile correlateRates
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl> <lgl>     <lgl>         
#> 1  1000     0     0     0     0     0     0       6 NA        FALSE         
```
