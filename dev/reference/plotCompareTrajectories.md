# Plot Bayesian population model results

Plot Bayesian population model results, with (optionally) the
distribution of outcomes from the initial model, local observations, and
true local state for comparison.

## Usage

``` r
plotCompareTrajectories(
  modTables,
  parameter,
  lowBound = 0,
  highBound = 1,
  facetVars = NULL,
  labFontSize = 14,
  legendPosition = "right",
  breakInterval = 1,
  typeLabels = c("Bayesian", "initial")
)
```

## Arguments

- modTables:

  list. A list of model results tables created using
  `[compareTrajectories()]`.

- parameter:

  character. Which parameter to plot, if more than one, a list of plots
  is returned.

- lowBound, highBound:

  numeric. Lower and upper y axis limits

- facetVars:

  character. Optional. Vector of column names to facet by

- labFontSize:

  numeric. Optional. Label font size if there are not facets. Font size
  is 10 pt if facets are used.

- legendPosition:

  "bottom", "right", "left","top", or "none". Legend position.

- breakInterval:

  number. How many years between x tick marks?

- typeLabels:

  vector of two labels. Default c("Bayesian","initial"). Names of models
  to be compared.

## Value

a ggplot object or list of ggplot objects if a vector of parameters was
given.

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
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
                            obsAnthroSlope = 1, projAnthroSlope = 5,
                            collarCount = 20, cowMult = 5)

simO <- simulateObservations(scns)

out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
                          disturbance = simO$simDisturbance,
                          startYear = 2014, niters=10)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 10
#>    Unobserved stochastic nodes: 33
#>    Total graph size: 665
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 29
#>    Unobserved stochastic nodes: 70
#>    Total graph size: 694
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf

out_tbl <- compareTrajectories(out, exData = simO$exData, paramTable = simO$paramTable,
                           simInitial = trajectoriesFromNational())
#> Using saved object

plotCompareTrajectories(out_tbl, parameter = "Recruitment")
```
