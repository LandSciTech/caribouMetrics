# Create summary table of demographic rates from survival and recruitment surveys

Create summary table of demographic rates from survival and recruitment
surveys

## Usage

``` r
estimateBayesianRates(
  surv_data,
  recruit_data,
  N0 = NA,
  disturbance = NULL,
  priors = NULL,
  shiny_progress = FALSE,
  return_mcmc = FALSE,
  i18n = NULL,
  niters = formals(bboutools::bb_fit_survival)$niters,
  nthin = formals(bboutools::bb_fit_survival)$nthin,
  ...
)
```

## Arguments

- surv_data:

  dataframe. Survival data in bboudata format

- recruit_data:

  dataframe. Recruitment data in bboudata format

- N0:

  dataframe. Optional. Initial population estimates, required columns
  are PopulationName and N0

- disturbance:

  dataframe. Optional. If provided, fit a Beta model that includes
  disturbance covariates.

- priors:

  list. Optional. If disturbance is NA, this should be
  list(priors_survival=c(...),priors_recruitment=c(...)); see
  [`bboutools::bb_priors_survival`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_survival.html)
  and
  [`bboutools::bb_priors_recruitment`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_recruitment.html)
  for details. If disturbance is not NA, see `betaNationalPriors` for
  details.

- shiny_progress:

  logical. Should shiny progress bar be updated. Only set to TRUE if
  using in an app.

- return_mcmc:

  boolean. If TRUE return fitted survival and recruitment models.
  Default FALSE.

- niters:

  integer. The number of iterations per chain after thinning and
  burn-in.

- nthin:

  integer. The number of the thinning rate.

- ...:

  Other parameters passed on to
  [`bboutools::bb_fit_survival`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_survival.html)
  and
  [`bboutools::bb_fit_recruitment`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_recruitment.html).

## Value

If `return_mcmc` is TRUE then a list with results and fitted models, if
FALSE just the results summaries are returned.

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
s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)
estimateBayesianRates(s_data, r_data, N0 = 500)
#>   PopulationName     R_bar      R_sd R_iv_mean R_iv_shape R_bar_lower
#> 1              A 0.1893446 0.0822676 0.3148171   25.22859   0.1652187
#> 2              B 0.2038394 0.1014211 0.3148171   25.22859   0.1729125
#>   R_bar_upper     S_bar      S_sd S_iv_mean S_iv_shape S_bar_lower S_bar_upper
#> 1   0.2146941 0.8817129 0.2336958 0.4927553   13.06435   0.8273442   0.9230849
#> 2   0.2383758 0.9073243 0.2823650 0.4927553   13.06435   0.8524029   0.9456843
#>    N0 nCollarYears nSurvYears nCowsAllYears nRecruitYears
#> 1 500          900         31          2353            27
#> 2 500          519         18          2001            15
```
