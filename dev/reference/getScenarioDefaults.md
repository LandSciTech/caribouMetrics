# Default scenario parameters.

Returns default parameters for scenarios. Use this function to get a
combination of
[`disturbanceDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceDefaults.md),
[`timeDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/timeDefaults.md),
[`demographyDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographyDefaults.md),
[`nationalTrajectoryDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/nationalTrajectoryDefaults.md),
and
[`monitoringDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/monitoringDefaults.md).
If only one of these sets of parameters is needed consider using the
relevant component function instead.

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
```

## Arguments

- paramTable:

  a data.frame with column names matching the arguments below. Any
  columns that are missing will be filled with the default values.

- includeDist:

  logical. Include
  [`disturbanceDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceDefaults.md)?

- includeTime:

  logical. Include
  [`timeDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/timeDefaults.md)?

- includeDemography:

  logical. Include
  [`demographyDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographyDefaults.md)?

- includeNational:

  logical. Include
  [`nationalTrajectoryDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/nationalTrajectoryDefaults.md)?

- includeMonitoring:

  logical. Include
  [`monitoringDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/monitoringDefaults.md)?

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
[`demographyDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographyDefaults.md),
[`disturbanceDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceDefaults.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
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
```
