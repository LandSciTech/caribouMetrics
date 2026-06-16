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

- N0:

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \$\_t=S_t\*(1+cR_t)/s\$.

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
#> 6  NA 0.8234332 0.9677935 NA 0.2324943 0.11624715 0.7376800         NA
#> 7  NA 0.9668711 0.9351453 NA 0.1222831 0.06114157 0.9111613         NA
#> 8  NA 0.9626445 0.9677935 NA 0.2004618 0.10023090 0.8749477         NA
#> 9  NA 0.9908783 0.9351453 NA 0.1472652 0.07363260 0.9229212         NA
#> 10 NA 0.9856054 0.9677935 NA 0.1673664 0.08368318 0.9094959         NA
#> 11 NA 0.9602650 0.9351453 NA 0.1581016 0.07905079 0.8899164         NA
#> 12 NA 0.9883486 0.9677935 NA 0.1634996 0.08174981 0.9136573         NA
#> 13 NA 0.7770321 0.9351453 NA 0.1374835 0.06874174 0.7270532         NA
#> 14 NA 0.9134489 0.9677935 NA 0.1956025 0.09780123 0.8320713         NA
#> 15 NA 0.9092682 0.9351453 NA 0.2288713 0.11443565 0.8159002         NA
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
