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

  dataframe. Survival data in bboudata format.

- recruit_data:

  dataframe. Recruitment data in bboudata format.

- N0:

  number or dataframe. Optional. Initial population size(s). If NA
  (default) then population growth rate is \$\_t=S_t\*(1+cR_t)/s\$. If a
  data frame N0 column is required, and PopulationName column is
  required if there is more than one row. Additional (optional)
  variation columns will be used by `addN0Variation()`.

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
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
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
s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)
estimateBayesianRates(s_data, r_data, N0 = 500)
#>   PopulationName     R_bar       R_sd R_iv_mean R_iv_shape R_bar_lower
#> 1              A 0.1989716 0.07893634 0.2878181   16.94093   0.1755463
#> 2              B 0.2110102 0.10435515 0.2878181   16.94093   0.1780015
#>   R_bar_upper     S_bar      S_sd S_iv_mean S_iv_shape S_bar_lower S_bar_upper
#> 1   0.2244535 0.8834576 0.2371920 0.5044639   13.57145   0.8264123   0.9251399
#> 2   0.2455056 0.9054767 0.2964022 0.5044639   13.57145   0.8465901   0.9455474
#>    N0 nCollarYears nSurvYears nCowsAllYears nRecruitYears
#> 1 500          900         31          2047            27
#> 2 500          519         18          1645            15
```
