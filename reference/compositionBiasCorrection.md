# Calculate bias correction term for calf:cow composition survey.

When composition surveys are conducted there is a possibility of bias in
calf cow ratios due to misidentifying young bulls as adult females and
vice versa or missing calves. Here we address this gap with a bias term
derived from a simple model of the recruitment survey observation
process. See [Hughes et al. (2025) Section
2.2](https://doi.org/10.1016/j.ecoinf.2025.103095) for a detailed
description of the model.

## Usage

``` r
compositionBiasCorrection(w, q, u, z, approx = F)
```

## Arguments

- w:

  number. The apparent number of adult females per collared animal in
  composition survey.

- q:

  number in 0, 1. Ratio of bulls to cows in composition survey groups.

- u:

  number in 0, 1. Probability of misidentifying young bulls as adult
  females and vice versa in composition survey.

- z:

  number in 0, \<1. Probability of missing calves in composition survey.

- approx:

  logical. If TRUE approximate the uncertainty about the value of the
  composition bias correction value (c) with the log-normal distribution
  of c given all the supplied values of `q`, `u`, and `z`. If FALSE the
  composition bias correction value (c) is returned for each value of
  `q`, `u`, and `z`

## Value

number or tibble. If `approx = FALSE` a vector of composition bias
correction values (c) of the same length as `q`, `u`, and `z`. If
`approx = TRUE` a tibble with one row per unique value of `w` and
columns `w`, `m`, `v`, `sig2`, `mu` representing `w`, mean `c`, variance
of `c`, and parameters for a log-normal approximation of the
distribution of `c`.

## References

Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
Integration of national demographic-disturbance relationships and local
data can improve caribou population viability projections and inform
monitoring decisions. Ecological Informatics, 87, p.103095.
<https://doi.org/10.1016/j.ecoinf.2025.103095>

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/compareTrajectories.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummaryForApp.md)

## Examples

``` r
# number of reps
nr <- 10

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = FALSE)
#>  [1] 1.1781994 0.9761174 1.1199332 1.0712302 1.1224818 0.9939395 0.8997471
#>  [8] 0.9974560 1.0607055 1.0023524

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = TRUE)
#> # A tibble: 1 × 5
#>       w     m       v    sig2     mu
#>   <dbl> <dbl>   <dbl>   <dbl>  <dbl>
#> 1     6  1.05 0.00761 0.00687 0.0458

```
