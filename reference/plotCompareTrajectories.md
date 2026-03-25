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
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/reference/estimateNationalRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/reference/getScenarioDefaults.md),
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
#> Warning: Adaptation incomplete
#> NOTE: Stopping adaptation
#> 
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 20
#>    Unobserved stochastic nodes: 79
#>    Total graph size: 694
#> 
#> Initializing model
#> 
#> Warning: Adaptation incomplete
#> NOTE: Stopping adaptation
#> 
#> 
#> Warning: no non-missing arguments to max; returning -Inf

out_tbl <- compareTrajectories(out, exData = simO$exData, paramTable = simO$paramTable,
                           simInitial = trajectoriesFromNational())
#> Using saved object

plotCompareTrajectories(out_tbl, parameter = "Recruitment")
```
