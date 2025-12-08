# Do multiple caribouPopGrowth runs

Do multiple caribouPopGrowth runs

## Usage

``` r
trajectoriesFromSummary(
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

  Sceanrio name

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
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md)

## Examples

``` r
 outParTab <- trajectoriesFromSummary(
   numSteps = 5, replicates = 2, N0 = NA, R_bar = 0.18, S_bar = 0.87,
   R_sd = 0.085, S_sd = 0.16,
   R_iv_mean = 0.34, S_iv_mean = 0.31,
   R_iv_shape = 18, S_iv_shape = 3.3,
   scn_nm = "base", addl_params = NULL, type = "logistic"
 )
 outParTab
#>    N0    lambda  lambdaE  N        R_t        X_t       S_t n_recruits
#> 1  NA 0.9483000 0.948300 NA 0.18000000 0.09000000 0.8700000         NA
#> 2  NA 0.9483000 0.948300 NA 0.18000000 0.09000000 0.8700000         NA
#> 3  NA 0.9483000 0.948300 NA 0.18000000 0.09000000 0.8700000         NA
#> 4  NA 0.9483000 0.948300 NA 0.18000000 0.09000000 0.8700000         NA
#> 5  NA 0.9483000 0.948300 NA 0.18000000 0.09000000 0.8700000         NA
#> 6  NA 0.9906811 0.945523 NA 0.33758306 0.16879153 0.8476115         NA
#> 7  NA 0.9445214 0.944343 NA 0.14824558 0.07412279 0.8793421         NA
#> 8  NA 0.9697466 0.945523 NA 0.22407055 0.11203527 0.8720466         NA
#> 9  NA 0.9490204 0.944343 NA 0.16240332 0.08120166 0.8777460         NA
#> 10 NA 0.9428867 0.945523 NA 0.21493632 0.10746816 0.8513894         NA
#> 11 NA 0.9010762 0.944343 NA 0.19358218 0.09679109 0.8215568         NA
#> 12 NA 0.9325672 0.945523 NA 0.15200955 0.07600478 0.8666943         NA
#> 13 NA 0.8571384 0.944343 NA 0.22738776 0.11369388 0.7696355         NA
#> 14 NA 0.8628831 0.945523 NA 0.06635537 0.03317768 0.8351740         NA
#> 15 NA 0.9064271 0.944343 NA 0.14212023 0.07106011 0.8462897         NA
#>    surviving_adFemales id  time type  scn
#> 1                   NA  1 value mean base
#> 2                   NA  1 value mean base
#> 3                   NA  1 value mean base
#> 4                   NA  1 value mean base
#> 5                   NA  1 value mean base
#> 6                   NA  1 value samp base
#> 7                   NA  2 value samp base
#> 8                   NA  1 value samp base
#> 9                   NA  2 value samp base
#> 10                  NA  1 value samp base
#> 11                  NA  2 value samp base
#> 12                  NA  1 value samp base
#> 13                  NA  2 value samp base
#> 14                  NA  1 value samp base
#> 15                  NA  2 value samp base
```
