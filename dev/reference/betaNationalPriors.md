# Get prior parameters for Bayesian Beta demographic rate models

Returns prior parameter values for Bayesian Beta demographic rate models
described by [Hughes et al.
(2025)](https://doi.org/10.1016/j.ecoinf.2025.103095). The starting
point is estimated coefficients from national demographic-disturbance
relationships in the table `popGrowthTableJohnsonECCC`.

## Usage

``` r
betaNationalPriors(
  modList = NULL,
  survivalModelNumber = "M1",
  recruitmentModelNumber = "M4",
  rAnthroSlopeSE = 0.006,
  rFireSlopeSE = 0.002,
  sAnthroSlopeSE = 5e-04,
  sIntSE = 0.06,
  sNuMin = 0.01,
  sNuMax = 0.13,
  rIntSE = 0.35,
  rNuMin = 0.01,
  rNuMax = 0.7,
  qMin = 0,
  qMax = 0,
  uMin = 0,
  uMax = 0,
  zMin = 0,
  zMax = 0,
  cowMult = 6,
  populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
  modelVersion = "Johnson",
  returnValues = TRUE
)
```

## Arguments

- modList:

  a named list of modifiers to use to change the priors. If a modifier
  is supplied here the corresponding argument below is ignored.

- survivalModelNumber, recruitmentModelNumber:

  character. Which model number to use see
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  for options.

- rAnthroSlopeSE:

  Standard deviation of effect of anthropogenic disturbance on
  recruitment.

- rFireSlopeSE:

  Standard deviation of effect of fire on recruitment.

- sAnthroSlopeSE:

  Standard deviation of effect of disturbance on survival.

- sIntSE:

  Standard deviation of survival intercept.

- sNuMin, sNuMax:

  Uniform prior for coefficient of variation among years for
  recruitment.

- rIntSE:

  Standard deviation of recruitment intercept.

- rNuMin, rNuMax:

  Uniform prior for coefficient of variation among years for
  recruitment.

- qMin:

  number in 0, 1. Minimum ratio of bulls to cows in composition survey
  groups.

- qMax:

  number in 0, 1. Maximum ratio of bulls to cows in composition survey
  groups.

- uMin:

  number in 0, 1. Minimum probability of misidentifying young bulls as
  adult females and vice versa in composition survey.

- uMax:

  number in 0, 1. Maximum probability of misidentifying young bulls as
  adult females and vice versa in composition survey.

- zMin:

  number in 0, 1. Minimum probability of missing calves in composition
  survey.

- zMax:

  number in 0, 1. Maximum probability of missing calves in composition
  survey.

- cowMult:

  number. The apparent number of adult females per collared animal in
  composition survey.

- populationGrowthTable:

  data.frame.[popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  is included in the package and should be used in most cases. A custom
  table of model coefficients and standard errors or confidence
  intervals can be provided but it must match the column names of
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md).
  If the table does not contain the standard error it is calculated from
  the confidence interval.

- modelVersion:

  character. Which model version to use. Currently the only option is
  "Johnson" for the model used in Johnson et. al. (2020), but additional
  options may be added in the future.

- returnValues:

  logical. Default is TRUE. If FALSE returns strings for some values
  showing the initial values and the modifier ie "0.9 \* 1.05"

## Value

a list with values:

- R_b0_mu: Recruitment intercept

- R_b0_sd: Recruitment intercept standard deviation,

- R_b1_mu: Recruitment anthropogenic disturbance slope,

- R_b1_sd: Recruitment anthropogenic disturbance standard deviation,

- R_b2_mu: Recruitment fire excluding anthropogenic disturbance slope,

- R_b2_sd: Recruitment fire excluding anthropogenic disturbance standard
  deviation,

- R_cv_min: Min of the prior distribution of the random effect of year
  on recruitment,

- R_cv_max: Max of the prior distribution of the random effect of year
  on recruitment,

- S_b0_mu: Adult female survival intercept,

- S_b0_sd: Adult female survival intercept standard error times
  modifier,

- S_b1_mu: Adult female survival anthropogenic disturbance slope,

- S_b1_sd: Adult female survival anthropogenic disturbance standard
  deviation,

- S_cv_min: Min of the prior distribution of the random effect of year
  on adult female survival,

- S_cv_max: Max of the prior distribution of the random effect of year
  on adult female survival,

- qMin,qMax,uMin,uMax,zMin,zMax,cowMult: Composition bias parameters.

## Details

Standard errors and random effects of year have been calibrated so that
the 95% prior prediction intervals for survival and recruitment from the
Bayesian model match the range between the 2.5% and 97.5% quantiles of
1000 survival and recruitment trajectories from the national demographic
model
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md).
A prior for the unknown composition survey bias correction term `c` is
set by specifying an apparent number of adult females per collared
animal(`cowMult`) and minimum and maximum values for each of the ratio
of bulls to cows (\\q\\), the probability of misidentifying young bulls
as adult females and vice versa (\\u\\), and the probability of missing
calves (\\z\\) in composition surveys. See
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md)
and [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095)
for additional details.

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
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

## Examples

``` r
betaNationalPriors()
#> $R_b0_mu
#> [1] -1.023
#> 
#> $R_b0_sd
#> [1] 0.35
#> 
#> $R_b1_mu
#> [1] -0.017
#> 
#> $R_b1_sd
#> [1] 0.006
#> 
#> $R_b2_mu
#> [1] -0.0081
#> 
#> $R_b2_sd
#> [1] 0.002
#> 
#> $R_cv_min
#> [1] 0.01
#> 
#> $R_cv_max
#> [1] 0.7
#> 
#> $S_b0_mu
#> [1] -0.142
#> 
#> $S_b0_sd
#> [1] 0.06
#> 
#> $S_b1_mu
#> [1] -8e-04
#> 
#> $S_b1_sd
#> [1] 5e-04
#> 
#> $S_cv_min
#> [1] 0.01
#> 
#> $S_cv_max
#> [1] 0.13
#> 
#> $cowMult
#> [1] 6
#> 
#> $qMin
#> [1] 0
#> 
#> $qMax
#> [1] 0
#> 
#> $uMin
#> [1] 0
#> 
#> $uMax
#> [1] 0
#> 
#> $zMin
#> [1] 0
#> 
#> $zMax
#> [1] 0
#> 
```
