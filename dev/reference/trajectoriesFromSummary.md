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
#>    N0    lambda  lambdaE  N       R_t        X_t       S_t n_recruits
#> 1  NA 0.9483000 0.948300 NA 0.1800000 0.09000000 0.8700000         NA
#> 2  NA 0.9483000 0.948300 NA 0.1800000 0.09000000 0.8700000         NA
#> 3  NA 0.9483000 0.948300 NA 0.1800000 0.09000000 0.8700000         NA
#> 4  NA 0.9483000 0.948300 NA 0.1800000 0.09000000 0.8700000         NA
#> 5  NA 0.9483000 0.948300 NA 0.1800000 0.09000000 0.8700000         NA
#> 6  NA 0.9305566 0.948335 NA 0.1436576 0.07182880 0.8681952         NA
#> 7  NA 1.0643953 0.986013 NA 0.3578052 0.17890260 0.9028696         NA
#> 8  NA 0.9339037 0.948335 NA 0.1737680 0.08688402 0.8592487         NA
#> 9  NA 0.9910432 0.986013 NA 0.2006741 0.10033703 0.9006724         NA
#> 10 NA 0.9345237 0.948335 NA 0.1642731 0.08213654 0.8635913         NA
#> 11 NA 0.9995650 0.986013 NA 0.2155775 0.10778874 0.9023065         NA
#> 12 NA 0.9525828 0.948335 NA 0.1703149 0.08515746 0.8778291         NA
#> 13 NA 1.0016319 0.986013 NA 0.2359655 0.11798276 0.8959278         NA
#> 14 NA 0.9377782 0.948335 NA 0.1753349 0.08766744 0.8621921         NA
#> 15 NA 0.9483883 0.986013 NA 0.1046146 0.05230728 0.9012465         NA
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
