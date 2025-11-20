# Get informative priors for bboutools logistic models

Returns intercept prior parameter values for bboutools logistic models
that are informed by national demographic-disturbance relationships.

## Usage

``` r
bbouNationalPriors(Anthro, fire_excl_anthro, month = T)
```

## Arguments

- Anthro:

  numeric. Percent non-overlapping buffered anthropogenic disturbance.

- fire_excl_anthro:

  numeric. Percent fire not overlapping with anthropogenic disturbance.

- month:

  logical. Does the survival data to be analyzed include month, or is it
  aggregated by year?

## Value

a list with values:

- priors_survival = c(b0_mu,b0_sd). Passed to
  [`bboutools::bb_fit_survival`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_survival.html).
  See `bb_priors_survival` for details.

- priors_recruitment=c(b0_mu,b0_sd). Passed to
  [`bboutools::bb_fit_recruitment`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_recruitment.html).
  See `bb_priors_recruitment` for details.

## Details

This is done by ... TO DO.

See `betaNationalPriors` for additional details.

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
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
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
bbouNationalPriors(Anthro=50,fire=5)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 0
#>    Unobserved stochastic nodes: 78
#>    Total graph size: 349
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 5
#>    Unobserved stochastic nodes: 24
#>    Total graph size: 120
#> 
#> Initializing model
#> 
#> $priors_recruitment
#>      b0_mu      b0_sd 
#> -1.7302119  0.5594394 
#> 
#> $priors_survival
#>     b0_mu     b0_sd 
#> 4.3290487 0.5171383 
#> 
```
