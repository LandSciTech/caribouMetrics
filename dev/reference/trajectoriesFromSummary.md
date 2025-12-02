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
  addl_params
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
#>    N0    lambda   lambdaE  N       R_t        X_t       S_t n_recruits
#> 1  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 2  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 3  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 4  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 5  NA 0.9483000 0.9483000 NA 0.1800000 0.09000000 0.8700000         NA
#> 6  NA 0.8732208 0.9291704 NA 0.2464673 0.12323363 0.7774169         NA
#> 7  NA 0.9592624 0.9732188 NA 0.1398352 0.06991761 0.8965760         NA
#> 8  NA 0.9104215 0.9291704 NA 0.1585313 0.07926566 0.8435565         NA
#> 9  NA 0.9779617 0.9732188 NA 0.1887095 0.09435474 0.8936423         NA
#> 10 NA 0.9700922 0.9291704 NA 0.1595875 0.07979376 0.8984051         NA
#> 11 NA 0.9755772 0.9732188 NA 0.1736364 0.08681822 0.8976452         NA
#> 12 NA 0.9618322 0.9291704 NA 0.1618140 0.08090702 0.8898381         NA
#> 13 NA 0.9816932 0.9732188 NA 0.2053741 0.10268704 0.8902737         NA
#> 14 NA 0.7646397 0.9291704 NA 0.2121796 0.10608981 0.6912998         NA
#> 15 NA 0.9572864 0.9732188 NA 0.1446053 0.07230267 0.8927390         NA
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
