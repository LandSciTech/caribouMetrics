# Add variation to initial population size

Applies stochastic variation to initial population size estimates
(`N0`). If uncertainty columns are present, a new value of `N0` is
sampled for each row. When `N.sd` is provided, variation is sampled from
either a Poisson distribution (if `N.sd = sqrt(N0)`) or a truncated
Normal distribution. When `N.lower` and `N.upper` are provided,
variation is sampled from a Uniform distribution bounded by those
values. Simulated values are rounded to integers and constrained to any
specified lower and upper bounds.

## Usage

``` r
addN0Variation(popInfo, forceDataFrame = F)
```

## Arguments

- popInfo:

  numeric, list, or data.frame. Initial population size information.
  Must contain an `N0` column or value. Optional columns `N.sd`,
  `N.lower`, and `N.upper` specify uncertainty in abundance estimates.

- forceDataFrame:

  logical. If `TRUE` and the result is numeric, return a data frame with
  column `N0`. Default `FALSE`.

## Value

A modified version of `popInfo` with updated `N0` values. If no
uncertainty columns are present, the input is returned unchanged.

## Details

This function is intended to be called at the same stage of simulation
and projection workflows as
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md)
so that uncertainty in abundance estimates is propagated through
subsequent analyses.

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
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
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
addN0Variation(500)
#> [1] 500

addN0Variation(data.frame(PopulationName = rep("A", 10),
                          N0 = 500, N.sd = 50))
#>    PopulationName  N0
#> 1               A 430
#> 2               A 513
#> 3               A 378
#> 4               A 500
#> 5               A 531
#> 6               A 557
#> 7               A 409
#> 8               A 488
#> 9               A 488
#> 10              A 486

addN0Variation(data.frame(PopulationName = rep("A", 10),
                          N0 = 500,
                          N.lower = 400,
                          N.upper = 600))
#>    PopulationName  N0
#> 1               A 458
#> 2               A 536
#> 3               A 547
#> 4               A 439
#> 5               A 596
#> 6               A 548
#> 7               A 410
#> 8               A 506
#> 9               A 539
#> 10              A 538
```
