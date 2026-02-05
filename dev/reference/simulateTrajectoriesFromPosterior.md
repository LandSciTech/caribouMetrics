# Format trajectory tables

Format trajectory tables

Get 95% prediction intervals from trajectories

Get a set of simulation results from fitted demographic models in raw
form

Assumes that rec_pred and surv_pred each include the same years and
populations.TO DO: check this.

## Usage

``` r
convertTrajectories(pars)

summarizeTrajectories(pars, returnSamples = T)

simulateTrajectoriesFromPosterior(
  popInfo = NA,
  rec_pred,
  surv_pred,
  initYear = NULL,
  correlateRates = FALSE,
  returnExpected = FALSE,
  c = formals(caribouPopGrowth)$c,
  K = FALSE,
  ...
)
```

## Arguments

- pars:

- returnSamples:

- popInfo:

  If NA (default) predictions are made without populations size, density
  dependence, or demographic stochasticity. See
  [`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md)
  for details.

- rec_pred:

  results returned by
  [`bboutools::bb_fit_recruitment()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_recruitment.html)
  or
  [`bboutools::bb_predict_calf_cow_ratio()`](https://poissonconsulting.github.io/bboutools/reference/bb_predict_calf_cow_ratio.html)
  functions of bboutools R package, or recruit_fit returned by
  [`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md).

- surv_pred:

  bboufit object return by
  [`bboutools::bb_fit_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_survival.html)
  or
  [`bboutools::bb_predict_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_predict_survival.html)
  functions of bboutools R package, or surv_fit returned by
  [`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md).

- initYear:

  numeric. Initial year.

- correlateRates:

  logical. Set TRUE to force correlation between recruitment and
  survival. Ignored

- returnExpected:

  logical. Default FALSE. Set TRUE to return expected values of R, S,
  and lambda (without interannual variation). Ignored if
  rec_pred/surv_pred are
  [`bboutools::bb_predict_calf_cow_ratio()`](https://poissonconsulting.github.io/bboutools/reference/bb_predict_calf_cow_ratio.html)/[`bboutools::bb_predict_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_predict_survival.html)
  results.

- c:

  Number. Bias correction term.

- K:

  Number. Carrying capacity.

## Value

convertTrajectories: formatted tables

summarizeTrajectories:

simulateTrajectoriesFromPosterior: a data frame with results from
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md)
for each set of survival/recruitment predictions.

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
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

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
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
  mod <- estimateBayesianRates(
    surv_data = bboudata::bbousurv_a,
    recruit_data = bboudata::bbourecruit_a, 
    N0 = NA, return_mcmc = TRUE
  )
                           
  outmcmc = simulateTrajectoriesFromPosterior(popInfo=NA,mod$recruit_fit,mod$surv_fit)
  names(outmcmc)
#>  [1] "N0"                  "lambda"              "lambdaE"            
#>  [4] "N"                   "R_t"                 "X_t"                
#>  [7] "S_t"                 "n_recruits"          "surviving_adFemales"
#> [10] "lab"                 "Year"                "PopulationName"     
#> [13] "id"                 

  #get 95% prediction intervals from demographic trajectories
  PImcmc <- summarizeTrajectories(convertTrajectories(outmcmc))
  str(PImcmc, max.level = 3, give.attr = FALSE)
#> List of 2
#>  $ summary:'data.frame': 189 obs. of  8 variables:
#>   ..$ MetricTypeID  : chr [1:189] "N" "N" "N" "N" ...
#>   ..$ Year          : num [1:189] 1989 2012 1995 1997 2001 ...
#>   ..$ PopulationName: chr [1:189] "A" "A" "A" "A" ...
#>   ..$ Mean          : num [1:189] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
#>   ..$ lower         : Named num [1:189] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ upper         : Named num [1:189] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ probViable    : num [1:189] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
#>   ..$ Parameter     : chr [1:189] "Female population size" "Female population size" "Female population size" "Female population size" ...
#>  $ samples: tibble [567,000 × 7] (S3: tbl_df/tbl/data.frame)
#>   ..$ Replicate       : chr [1:567000] "x1" "x1" "x1" "x1" ...
#>   ..$ LambdaPercentile: logi [1:567000] NA NA NA NA NA NA ...
#>   ..$ Year            : num [1:567000] 1989 1989 1989 1989 1989 ...
#>   ..$ PopulationName  : chr [1:567000] "A" "A" "A" "A" ...
#>   ..$ Timestep        : num [1:567000] 1989 1989 1989 1989 1989 ...
#>   ..$ MetricTypeID    : chr [1:567000] "c" "survival" "recruitment" "X" ...
#>   ..$ Amount          : num [1:567000] NA 0.811 0.24 0.12 NA ...
# }

```
