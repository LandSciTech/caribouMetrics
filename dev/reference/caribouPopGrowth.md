# Caribou demographic model

A caribou demographic model with density dependence and interannual
variability following [Johnson et. al.
(2020)](doi:10.1111/1365-2664.13637) with modifications described in
[Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095) and
[Dyson et al. (in press)](https://doi.org/10.1101/2022.06.01.494350).
Default parameter values give the model in Dyson et al. (in press). Set
`probOption = "matchJohnson2020"` to reproduce the model used in Johnson
et al. 2020. Set `probOption = "continuous"`, `interannualVar = FALSE`,
and `K = FALSE` to reproduce the simpler 2-stage demographic model
without interannual variability, density dependence, or discrete numbers
of animals used by [Stewart et al.
(2023)](https://doi.org/10.1002/eap.2816).

## Usage

``` r
caribouPopGrowth(
  N0,
  numSteps,
  R_bar,
  S_bar,
  P_0 = 1,
  P_K = 0.6,
  a = 1,
  b = 4,
  K = 10000,
  r_max = 1.3,
  s = 0.5,
  l_R = 0,
  h_R = 0.82,
  l_S = 0.61,
  h_S = 1,
  c = 1,
  interannualVar = list(R_CV = 0.46, S_CV = 0.08696),
  probOption = "binomial",
  progress = interactive(),
  warn = T
)
```

## Arguments

- N0:

  Number or vector of numbers. Initial population size for one or more
  sample populations. If NA then population growth rate is
  \$\_t=S_t\*(1+cR_t)/s\$.

- numSteps:

  Number. Number of years to project.

- R_bar:

  Number or vector of numbers. Expected recruitment rate (calf:cow
  ratio) for one or more sample populations.

- S_bar:

  Number or vector of numbers. Expected adult female survival for one or
  more sample populations.

- P_0:

  Number. Maximum recruitment multiplier.

- P_K:

  Number. Recruitment multiplier at carrying capacity.

- a:

  Number. Density dependence shape parameter.

- b:

  Number. Allee effect parameter.

- K:

  Number. Carrying capacity.

- r_max:

  Number. Maximum population growth rate.

- s:

  Number. Sex ratio.

- l_R:

  Number. Minimum recruitment.

- h_R:

  Number. Maximum recruitment.

- l_S:

  Number. Minimum survival.

- h_S:

  Number. Maximum survival.

- c:

  Number. Bias correction term.

- interannualVar:

  list or logical. List containing interannual variability parameters.
  These can be either coefficients of variation (R_CV, S_CV), beta
  precision parameters (R_phi, S_phi), or random effects parameters from
  a logistic glmm (R_annual, S_annual). Set to `FALSE` to ignore
  interannual variability.

- probOption:

  Character. Choices are "binomial","continuous" or "matchJohnson2020".
  See description for details.

- progress:

  Logical. Should progress updates be shown?

## Value

A data.frame of population size (`N`), expected growth rate (`lambdaE`),
true growth rate (`lambda`), apparent annual reproduction rate (`R_t`),
adjusted reproduction (`X_t`), survival (`S_t`), number of recruits
(`n_recruits`), and surviving females (`surviving_adFemales`) for each
sample population projected for numSteps years.

## Details

If R_annual and S_annual are provided, interannual variation in survival
and recruitment is modelled as in a logistic glmm with random effect of
year.

If initial populations size N0 is NA then population growth rate is
\$\_t=S_t\*(1+cR_t/s)\$. In this case density dependence
(P_0,P_K,a,b,K,r_max) and demographic stochasticity (probOption) are
ignored.

See
[`vignette("caribouDemography")`](https://landscitech.github.io/caribouMetrics/dev/articles/caribouDemography.md)
and [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095)
for additional details and examples.

## References

Dyson, M., Endicott, S., Simpkins, C., Turner, J. W., Avery-Gomm, S.,
Johnson, C. A., Leblond, M., Neilson, E. W., Rempel, R., Wiebe, P. A.,
Baltzer, J. L., Stewart, F. E. C., & Hughes, J. (in press). Effective
conservation decisions require models designed for purpose: a case study
of boreal caribou in Ontario’s Ring of Fire. Ecology and Evolution In
press: <https://doi.org/10.1101/2022.06.01.494350>

Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
Integration of national demographic-disturbance relationships and local
data can improve caribou population viability projections and inform
monitoring decisions. Ecological Informatics, 87, p.103095.
<https://doi.org/10.1016/j.ecoinf.2025.103095>

Johnson, C.A., Sutherland, G.D., Neave, E., Leblond, M., Kirby, P.,
Superbie, C. and McLoughlin, P.D., 2020. Science to inform policy:
linking population dynamics to habitat for a threatened species in
Canada. Journal of Applied Ecology, 57(7), pp.1314-1327.
<https://doi.org/10.1111/1365-2664.13637>

Stewart, F.E., Micheletti, T., Cumming, S.G., Barros, C., Chubaty, A.M.,
Dookie, A.L., Duclos, I., Eddy, I., Haché, S., Hodson, J. and Hughes,
J., 2023. Climate-informed forecasts reveal dramatic local habitat
shifts and population uncertainty for northern boreal caribou.
Ecological Applications, 33(3), p.e2816.
<https://doi.org/10.1002/eap.2816>

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
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
caribouPopGrowth(100, 2, 0.5, 0.7)
#>    N0    lambda lambdaE  N       R_t      X_t       S_t n_recruits
#> 1 100 0.8602325   0.875 74 0.3259161 0.162958 0.6982619         16
#>   surviving_adFemales
#> 1                  58
```
