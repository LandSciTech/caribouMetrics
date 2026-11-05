# Demographic projections for cases with no change in demographic rates over time. This is the method used (so far) in the demography app. TO DO: Consider removing and replacing with call to trajectoriesFromSummary.

Demographic projections for cases with no change in demographic rates
over time. This is the method used (so far) in the demography app. TO
DO: Consider removing and replacing with call to
trajectoriesFromSummary.

## Usage

``` r
trajectoriesFromSummaryForApp(
  numSteps,
  replicates,
  N0,
  R_bar,
  S_bar,
  R_sd,
  S_sd,
  R_iv_mean,
  R_iv_shape,
  S_iv_mean,
  S_iv_shape,
  scn_nm,
  type = "logistic",
  addl_params = list(),
  doSummary = F,
  returnSamples = T
)
```

## Arguments

- numSteps:

  Number. Number of years to project.

- replicates:

  integer. Number of replicate populations.

- N0:

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \\\lambda_t=S_t\*(1+cR_t/s)\\.

- R_bar:

  Number or vector of numbers. Expected recruitment rate (calf:cow
  ratio) for one or more sample populations.

- S_bar:

  Number or vector of numbers. Expected adult female survival for one or
  more sample populations.

- R_sd, S_sd:

  standard deviation of R_bar and S_bar

- R_iv_mean, R_iv_shape, S_iv_mean, S_iv_shape:

  define the mean and shape of the interannual variation

- scn_nm:

  Scenario name

- type:

  "logistic" or "beta" defines how demographic rates are sampled from
  the given mean and standard deviation.

- addl_params:

  a list of additional parameters for `caribouPopGrowth`

- doSummary:

  logical. Default TRUE. If FALSE returns unprocessed outcomes from
  caribouPopGrowth. If TRUE returns summaries and (if returnSamples = T)
  sample trajectories from prepareTrajectories.

- returnSamples:

  logical. If FALSE returns only summaries. If TRUE returns example
  trajectories as well.

## Value

a data.frame

## See also

Caribou demography functions:
[`addMissingYears()`](https://landscitech.github.io/caribouMetrics/dev/reference/addMissingYears.md),
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
[`getCaribouYear()`](https://landscitech.github.io/caribouMetrics/dev/reference/getCaribouYear.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
 outParTab <- trajectoriesFromSummaryForApp(
   numSteps = 5, replicates = 2, N0 = NA, R_bar = 0.18, S_bar = 0.87,
   R_sd = 0.085, S_sd = 0.16,
   R_iv_mean = 0.34, S_iv_mean = 0.31,
   R_iv_shape = 18, S_iv_shape = 3.3,
   scn_nm = "base", addl_params = NULL, type = "logistic"
 )
 outParTab
#>    N0    lambda   lambdaE  N       R_t        X_t       S_t n_recruits
#> 1  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 2  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 3  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 4  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 5  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 6  NA 0.9556370 0.9027069 NA 0.1819505 0.09097526 0.8759474         NA
#> 7  NA 0.9514840 0.9627004 NA 0.1982733 0.09913664 0.8656649         NA
#> 8  NA 0.9057918 0.9027069 NA 0.1518784 0.07593922 0.8418615         NA
#> 9  NA 0.9603580 0.9627004 NA 0.2179213 0.10896064 0.8659983         NA
#> 10 NA 0.8786781 0.9027069 NA 0.1279517 0.06397583 0.8258440         NA
#> 11 NA 0.9740977 0.9627004 NA 0.1936647 0.09683234 0.8881008         NA
#> 12 NA 0.8917126 0.9027069 NA 0.1379252 0.06896261 0.8341850         NA
#> 13 NA 0.9705263 0.9627004 NA 0.2255463 0.11277314 0.8721691         NA
#> 14 NA 0.9408761 0.9027069 NA 0.1409463 0.07047313 0.8789348         NA
#> 15 NA 0.9398549 0.9627004 NA 0.1663683 0.08318415 0.8676778         NA
#>    surviving_adFemales id time type  scn
#> 1                   NA  1    1 mean base
#> 2                   NA  1    2 mean base
#> 3                   NA  1    3 mean base
#> 4                   NA  1    4 mean base
#> 5                   NA  1    5 mean base
#> 6                   NA  1    1 samp base
#> 7                   NA  2    1 samp base
#> 8                   NA  1    2 samp base
#> 9                   NA  2    2 samp base
#> 10                  NA  1    3 samp base
#> 11                  NA  2    3 samp base
#> 12                  NA  1    4 samp base
#> 13                  NA  2    4 samp base
#> 14                  NA  1    5 samp base
#> 15                  NA  2    5 samp base
```
