# Model demographic rates from survival and recruitment surveys

Fit Bayesian models for survival and recruitment based on boreal caribou
survey data. If disturbance data is provided a Beta model is fitted with
disturbance covariates and priors informed by national
demographic-disturbance relationships. If disturbance data is NULL
bboutools logistic models are fit with
[`bboutools::bb_fit_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_survival.html)
and
[`bboutools::bb_fit_recruitment()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_recruitment.html).
In both cases the fitted the models are used to estimate recruitment and
survival rates and return summary tables.

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
  (default) then population growth rate is
  \\\lambda_t=S_t\*(1+cR_t)/s\\. If a data frame N0 column is required,
  and PopulationName column is required if there is more than one row.
  Additional (optional) variation columns will be used by
  [`addN0Variation()`](https://landscitech.github.io/caribouMetrics/dev/reference/addN0Variation.md).

- disturbance:

  dataframe. Optional. If provided, fit a Beta model that includes
  disturbance covariates.

- priors:

  list. Optional. If disturbance is NULL, this should be
  list(priors_survival=c(...),priors_recruitment=c(...)); see
  [`bboutools::bb_priors_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_survival.html)
  and
  [`bboutools::bb_priors_recruitment()`](https://poissonconsulting.github.io/bboutools/reference/bb_priors_recruitment.html)
  for details. If disturbance is not NULL, see
  [`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md)
  for details.

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
  [`bboutools::bb_fit_survival()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_survival.html)
  and
  [`bboutools::bb_fit_recruitment()`](https://poissonconsulting.github.io/bboutools/reference/bb_fit_recruitment.html).

## Value

If `return_mcmc` FALSE a data.frame containing mean survival and
recruitment estimates, standard deviation, upper and lower credible
intervals, and interannual variability distribution parameters for each
population projected by the Bayesian model, as well as a summary of the
amount of input data.

If `return_mcmc` is TRUE then a list with elements:

- parTab: the data.frame returned when `return_mcmc` FALSE

- parList: a list of data.frames Rbar, Sbar, Siv and Riv. Sbar and Rbar
  contain yearly expected survival and recruitment estimates, standard
  deviation, and upper and lower credible intervals for each population.
  Estimates will only vary by year if disturbance was provided. Siv and
  Riv contains interannual variability distribution parameters for each
  population.

- surv_fit and recruit_fit: fitted model objects containing samples and
  input data

## Details

See
[`vignette("compare-bayesian-models")`](https://landscitech.github.io/caribouMetrics/dev/articles/compare-bayesian-models.md)
for more details on the differences between the two models.

## Interpretation of Years

Throughout the package a year is treated as the calendar year only if it
is in a table that also contains a month column. In all other contexts
year is treated as the caribou year which is the annual demographic
cycle of caribou, measured from one calving season to the next. By
default, April is set as the start of the caribou year.

## See also

Caribou demography functions:
[`addMissingYears()`](https://landscitech.github.io/caribouMetrics/dev/reference/addMissingYears.md),
[`addN0Variation()`](https://landscitech.github.io/caribouMetrics/dev/reference/addN0Variation.md),
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
[`getCaribouYear()`](https://landscitech.github.io/caribouMetrics/dev/reference/getCaribouYear.md),
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
#> 1              A 0.1989571 0.08037253 0.2857571   19.00788   0.1740515
#> 2              B 0.2130707 0.09865786 0.2857571   19.00788   0.1824531
#>   R_bar_upper     S_bar      S_sd S_iv_mean S_iv_shape S_bar_lower S_bar_upper
#> 1   0.2250004 0.8830194 0.2508802 0.5007712   14.57276   0.8222825   0.9258783
#> 2   0.2480182 0.9085359 0.3040960 0.5007712   14.57276   0.8486361   0.9504125
#>    N0 nCollarYears nSurvYears nCowsAllYears nRecruitYears
#> 1 500          900         31          2047            27
#> 2 500          519         18          1645            15
```
