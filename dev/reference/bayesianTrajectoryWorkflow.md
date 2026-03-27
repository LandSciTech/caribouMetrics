# Bayesian population model for boreal caribou

Bayesian population model for boreal caribou

## Usage

``` r
bayesianTrajectoryWorkflow(
  surv_data = bboudata::bbousurv_a,
  recruit_data = bboudata::bbourecruit_a,
  disturbance = NULL,
  priors = "default",
  startYear = NULL,
  endYear = NULL,
  N0 = NA,
  returnSamples = F,
  inputList = list(),
  niters = formals(bboutools::bb_fit_survival)$niters,
  nthin = formals(bboutools::bb_fit_survival)$nthin,
  ...
)
```

## Arguments

- surv_data:

  either a path to a csv file or a survival data table in bboutools
  format.

- recruit_data:

  either a path to a csv file or a recruitment data table in bboutools
  format.

- disturbance:

  either a path to a csv file or a dataframe containing the columns
  "Anthro","Fire_excl_anthro", and "Year".

- priors:

  a list of model priors. If disturbance is NA, this should be
  list(priors_survival=c(...),priors_recruitment=c(...)); see
  [`bboutools::bb_priors_survival`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_survival.html)
  and
  [`bboutools::bb_priors_recruitment`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_recruitment.html)
  for details. If disturbance is not NA, see
  [`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md)
  for details.

- startYear, endYear:

  year defining the beginning of the observation period and the end of
  the projection period.

- N0:

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \$\_t=S_t\*(1+cR_t)/s\$.

- returnSamples:

  logical. If F returns only summaries. If T returns example
  trajectories.

- inputList:

  an optional list of inputs with names matching the above. If an
  argument is included in this list it will override the named argument.

- niters:

  A whole number of the number of iterations per chain after thinning
  and burn-in.

- nthin:

  integer. The number of the thinning rate.

## Value

a list with elements:

- result: a list of model results:

  - summary: a data.frame

  - samples: a tibble providing the full range of MCMC trajectories from
    the model. It is in a long format where "Amount" gives the value for
    each metric in Anthro, Fire_excl_anthro, c, survival, recruitment,
    X, N, lambda, Sbar, Rbar, Xbar, Nbar, and lambda_bar, with a row for
    each combination of "MetricTypeID", "Replicate", "Year" and
    "LambdaPercentile"

  - surv_data: a data.frame

  - recruit_data: a tibble

  - popInfo: a data.frame

- inData: a list of data that is used as input to the jags model:

  - survDataIn: survival data

  - disturbanceIn: disturbance data

  - recruitDataIn: composition data

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
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
# \donttest{
  # Note these examples take a long time to run!
  
  # Using observed survival, recruitment and disturbance data
  mod <- bayesianTrajectoryWorkflow(
    surv_data = bboudata::bbousurv_a,
    recruit_data = bboudata::bbourecruit_a,
    disturbance = NULL
  )
#> Warning: requested year range: 1986 - 2016 does not match recruitment data year range:  1990 - 2016
#> Warning: missing years of recruitment data: 1986, 1987, 1988, 1989
#> Warning: no non-missing arguments to max; returning -Inf
  str(mod, max.level = 2)
#> List of 4
#>  $ result :List of 4
#>   ..$ summary     :'data.frame': 300 obs. of  8 variables:
#>   ..$ surv_data   :'data.frame': 363 obs. of  8 variables:
#>   ..$ recruit_data:'data.frame': 30 obs. of  8 variables:
#>   ..$ popInfo     :'data.frame': 3000 obs. of  4 variables:
#>  $ inData :List of 1
#>   ..$ disturbanceIn: NULL
#>  $ parTab :'data.frame': 1 obs. of  18 variables:
#>   ..$ PopulationName: chr "A"
#>   ..$ R_bar         : num 0.198
#>   ..$ R_sd          : num 0.0831
#>   ..$ R_iv_mean     : num 0.325
#>   ..$ R_iv_shape    : num 13.4
#>   ..$ R_bar_lower   : num 0.174
#>   ..$ R_bar_upper   : num 0.226
#>   ..$ S_bar         : num 0.873
#>   ..$ S_sd          : num 0.173
#>   ..$ S_iv_mean     : num 0.346
#>   ..$ S_iv_shape    : num 4.29
#>   ..$ S_bar_lower   : num 0.834
#>   ..$ S_bar_upper   : num 0.909
#>   ..$ N0            : logi NA
#>   ..$ nCollarYears  : int NA
#>   ..$ nSurvYears    : int 31
#>   ..$ nCowsAllYears : int NA
#>   ..$ nRecruitYears : int 31
#>  $ parList:List of 5
#>   ..$ Rbar:'data.frame': 27 obs. of  7 variables:
#>   ..$ Sbar:'data.frame': 27 obs. of  7 variables:
#>   ..$ Siv :'data.frame': 1 obs. of  2 variables:
#>   ..$ Riv :'data.frame': 1 obs. of  2 variables:
#>   ..$ type: chr "bbou"
  
  # Using simulated observation data
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
                              obsAnthroSlope = 1, projAnthroSlope = 5,
                              collarCount = 20, cowMult = 5)
  
  simO <- simulateObservations(scns)
  
  out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
                           disturbance = simO$simDisturbance,
                           startYear = 2014)
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
#>    Observed stochastic nodes: 20
#>    Unobserved stochastic nodes: 79
#>    Total graph size: 694
#> 
#> Initializing model
#> 
#> Warning: no non-missing arguments to max; returning -Inf
# }
```
