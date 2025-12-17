# Run the Bayesian population model for multiple parameter sets

Define scenarios in a table and
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
run the
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md)
model and
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md)
for each scenario.

## Usage

``` r
bayesianScenariosWorkflow(
  scns,
  simInitial,
  ePars = list(collarOnTime = 4, collarOffTime = 4, collarNumYears = 4),
  Rep = NULL,
  printProgress = FALSE,
  priors = "default",
  niters = formals(bboutools::bb_fit_survival)$niters,
  nthin = formals(bboutools::bb_fit_survival)$nthin,
  returnSamples = F,
  ...
)
```

## Arguments

- scns:

  data.frame. Parameters for the simulations. See
  [`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md)
  for details.

- simInitial:

  Initial simulation results, produced by calling
  `trajectoriesFromAny()`

- ePars:

  list. Additional parameters passed on to
  [`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md)

- Rep:

  integer. Optional. If specified, select specified replicate
  trajectory.

- printProgress:

  logical. Should the scenario number and parameters be printed at each
  step?

- priors:

  a list of model priors. If disturbance is NA, this should be
  list(priors_survival=c(...),priors_recruitment=c(...)); see
  [`bboutools::bb_priors_survival`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_survival.html)
  and
  [`bboutools::bb_priors_recruitment`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_recruitment.html)
  for details. If disturbance is not NA, see
  [`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md)
  for details.

- niters:

  A whole number of the number of iterations per chain after thinning
  and burn-in.

- nthin:

  integer. The number of the thinning rate.

- returnSamples:

  logical. Optional. If true, return full results from
  [`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md).

## Value

A list similar to
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md)
where tables for each scenario have been appended together. Plus an
error log for any scenarios that failed to run.

## See also

Caribou demography functions:
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
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
scns <- expand.grid(
  obsYears =c(10, 20), collarCount = c(30, 300), cowMult = 2, collarInterval = 2,
  iAnthro = 0,
  obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.9,
  rQuantile = 0.7, N0 = 1000
)

eParsIn <- list(collarOnTime = 4, collarOffTime = 4, collarNumYears = 3)
simsIn <- trajectoriesFromNational()
#> Warning: Setting expected survival S_bar to be between l_S and h_S.
#> Updating cached initial simulations.
scResults <- bayesianScenariosWorkflow(scns, simsIn, eParsIn,
                       niters = 10)# only set to speed up example. Normally keep defaults.
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 10
#>    Unobserved stochastic nodes: 83
#>    Total graph size: 409
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 54
#>    Unobserved stochastic nodes: 170
#>    Total graph size: 666
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 20
#>    Unobserved stochastic nodes: 93
#>    Total graph size: 489
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 74
#>    Unobserved stochastic nodes: 200
#>    Total graph size: 806
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 10
#>    Unobserved stochastic nodes: 83
#>    Total graph size: 409
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 54
#>    Unobserved stochastic nodes: 170
#>    Total graph size: 666
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 20
#>    Unobserved stochastic nodes: 93
#>    Total graph size: 489
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 74
#>    Unobserved stochastic nodes: 200
#>    Total graph size: 806
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf
```
