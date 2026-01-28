# Default parameters for simulation of example demographic trajectories.

Returns default parameter values for simulations of example demographic
trajectories. See
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)
for additional details.

## Usage

``` r
getScenarioDefaults(
  paramTable = NULL,
  iFire = 0,
  iAnthro = 0,
  obsAnthroSlope = 2,
  projAnthroSlope = 2,
  rSlopeMod = 1,
  sSlopeMod = 1,
  lQuantile = NA,
  sQuantile = NA,
  rQuantile = NA,
  correlateRates = F,
  projYears = 35,
  obsYears = 15,
  preYears = 0,
  N0 = 1000,
  qMin = 0,
  qMax = 0,
  uMin = 0,
  uMax = 0,
  zMin = 0,
  zMax = 0,
  cowMult = 6,
  collarInterval = NA,
  cowCount = NA,
  collarCount = NA,
  startYear = NA,
  interannualVar = list(eval(formals(caribouPopGrowth)$interannualVar)),
  curYear = 2023
)
```

## Arguments

- paramTable:

  a data.frame with column names matching the arguments below. Any
  columns that are missing will be filled with the default values.

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

- rSlopeMod:

  number. Disturbance-recruitment slope multiplier

- sSlopeMod:

  number. Disturbance-survival slope multiplier

- lQuantile:

  number in 0, 1. Lambda quantile

- sQuantile:

  number in 0,1. Survival quantile.

- rQuantile:

  number in 0,1. Recruitment quantile.

- correlateRates:

  logical. Set TRUE to force correlation between recruitment and
  survival.

- projYears:

  Number of years of projections

- obsYears:

  Number of years of observations

- preYears:

  Number of years before monitoring begins

- N0:

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \$\_t=S_t\*(1+cR_t)/s\$.

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

- startYear:

  year. First year in observation period. Optional, if not provided it
  will be calculated from `curYear` and `obsYears`

- interannualVar:

  list or logical. List containing interannual variability parameters.
  These can be either coefficients of variation (R_CV, S_CV), beta
  precision parameters (R_phi, S_phi), or random effects parameters from
  a logistic glmm (R_annual, S_annual). Set to `FALSE` to ignore
  interannual variability.

- curYear:

  year. The current year. All years before are part of the observation
  period and years after are part of the projection period.

## Value

a data.frame of parameter values including a label that combines all the
parameter names and values into a string

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
#> # A tibble: 1 × 23
#>   iFire iAnthro obsAnthroSlope projAnthroSlope rSlopeMod sSlopeMod
#>   <dbl>   <dbl>          <dbl>           <dbl>     <dbl>     <dbl>
#> 1     0       0              2               2         1         1
#> # ℹ 17 more variables: correlateRates <lgl>, projYears <dbl>, obsYears <dbl>,
#> #   preYears <dbl>, N0 <dbl>, qMin <dbl>, qMax <dbl>, uMin <dbl>, uMax <dbl>,
#> #   zMin <dbl>, zMax <dbl>, cowMult <dbl>, interannualVar <list>,
#> #   curYear <dbl>, ID <int>, label <chr>, startYear <dbl>

# paramTable list takes precedence over argument values
getScenarioDefaults(paramTable = data.frame(iFire = 10, iAnthro = 20, obsYears = 1), obsYears = 5)
#>   iFire iAnthro obsAnthroSlope projAnthroSlope rSlopeMod sSlopeMod
#> 1    10      20              2               2         1         1
#>   correlateRates projYears obsYears preYears   N0 qMin qMax uMin uMax zMin zMax
#> 1          FALSE        35        1        0 1000    0    0    0    0    0    0
#>   cowMult   interannualVar curYear ID
#> 1       6 0.46000, 0.08696    2023  1
#>                                                                                                                                                                                                                                               label
#> 1 ID1_curYear2023_interannualVarlist(R_CV = 0.46, S_CV = 0.08696)_cowMult6_zMax0_zMin0_uMax0_uMin0_qMax0_qMin0_N01000_preYears0_obsYears1_projYears35_correlateRatesFALSE_sSlopeMod1_rSlopeMod1_projAnthroSlope2_obsAnthroSlope2_iAnthro20_iFire10_
#>   startYear
#> 1      2023
```
