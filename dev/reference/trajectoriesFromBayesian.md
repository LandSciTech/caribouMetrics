# Get trajectories from a Bayesian model result

Get trajectories from a Bayesian model result

## Usage

``` r
trajectoriesFromBayesian(
  bayesianResults,
  N0 = NULL,
  cPars = demographyDefaults(),
  returnSamples = TRUE,
  doSummary = TRUE,
  ...
)
```

## Arguments

- bayesianResults:

  A result from `estimateBayesianRates`

- N0:

- cPars:

- returnSamples:

- doSummary:

- ...:

## See also

Caribou demography functions:
[`addN0Variation()`](https://landscitech.github.io/caribouMetrics/dev/reference/addN0Variation.md),
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
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
                                         return_mcmc = TRUE)

trajB <- trajectoriesFromBayesian(bbouInformative)
str(trajB, max.level = 1)
#> List of 5
#>  $ summary     :'data.frame':    60 obs. of  8 variables:
#>  $ samples     :'data.frame':    198000 obs. of  8 variables:
#>  $ surv_data   :'data.frame':    84 obs. of  9 variables:
#>  $ recruit_data:'data.frame':    6 obs. of  9 variables:
#>  $ popInfo     :'data.frame':    3000 obs. of  4 variables:
plotTrajectories(trajB)
#> Warning: Removed 420 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 12 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).
#> Warning: Removed 12 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
