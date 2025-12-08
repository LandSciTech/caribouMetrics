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
`approx = TRUE` a tibble with on row per unique value of `w` and columns
`w`, `m`, `v`, `sig2`, `mu` representing `w`, mean `c`, variance of `c`,
and parameters for a log-normal approximation of the distribution of
`c`.

## References

Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
Integration of national demographic-disturbance relationships and local
data can improve caribou population viability projections and inform
monitoring decisions. Ecological Informatics, 87, p.103095.
<https://doi.org/10.1016/j.ecoinf.2025.103095>

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`bbouNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/bbouNationalPriors.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
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
# number or reps
nr <- 10

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = FALSE)
#>  [1] 1.116720 1.143107 0.997899 1.046624 1.006478 1.059040 1.054061 1.123068
#>  [9] 1.040555 1.113304

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = TRUE)
#> # A tibble: 1 × 5
#>       w     m       v    sig2     mu
#>   <dbl> <dbl>   <dbl>   <dbl>  <dbl>
#> 1     6  1.07 0.00461 0.00405 0.0616

```
