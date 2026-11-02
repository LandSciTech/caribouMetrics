# Sample demographic rates

Apply the sampled coefficients to the disturbance covariates to
calculate expected recruitment and survival according to the beta
regression models estimated by Johnson et al.
(2020).`estimateNationalRates` is a wrapper around
`estimateNationalRate` to sample both survival and recruitment rates
based on the result of
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md)
and using recommended defaults.

## Usage

``` r
estimateNationalRate(
  covTable,
  coefSamples,
  coefValues,
  modelVersion,
  resVar,
  ignorePrecision,
  returnSample,
  quantilesToUse = NULL,
  predInterval = c(0.025, 0.975),
  transformFn = function(y) {
     y
 }
)

estimateNationalRates(
  covTable,
  popGrowthPars,
  ignorePrecision = FALSE,
  returnSample = FALSE,
  useQuantiles = TRUE,
  predInterval = list(PI_R = c(0.025, 0.975), PI_S = c(0.025, 0.975)),
  transformFns = list(S_transform = function(y) {
(y * 46 - 0.5)/45
 }, R_transform
    = function(y) {
     y
 })
)
```

## Arguments

- covTable:

  data.frame. A table of covariate values to be used. Column names must
  match the coefficient names in
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md).
  Each row is a different scenario.

- coefSamples:

  matrix. Bootstrapped coefficients with one row per replicate and one
  column per coefficient

- coefValues:

  data.table. One row table with expected values for each coefficient

- modelVersion:

  character. Which model version to use. Currently the only option is
  "Johnson" for the model used in Johnson et. al. (2020), but additional
  options may be added in the future.

- resVar:

  character. Response variable, typically "femaleSurvival" or
  "recruitment"

- ignorePrecision:

  logical. Should the precision of the model be used if it is available?
  When precision is used variation among populations around the National
  mean responses is considered in addition to the uncertainty about the
  coefficient estimates.

- returnSample:

  logical. If TRUE the returned data.frame has replicates \* scenarios
  rows. If FALSE the returned data.frame has one row per scenario and
  additional columns summarizing the variation among replicates. See
  Value for details.

- quantilesToUse:

  numeric vector of length `coefSamples`. See `useQuantiles`.

- predInterval:

  numeric vector with length 2. The default 95% interval is
  (`c(0.025,0.975)`). Only relevant when `returnSample = TRUE` and
  `quantilesToUse = NULL`.

- transformFn:

  function used to transform demographic rates.

- popGrowthPars:

  list. Coefficient values and (optionally) quantiles returned by
  `getNationalCoefficients`.

- useQuantiles:

  logical or numeric. If it is a numeric vector it must be length 2 and
  give the low and high limits of the quantiles to use. Only relevant
  when `ignorePrecision = FALSE`. If `useQuantiles != FALSE`, each
  replicate population is assigned to a quantile of the distribution of
  variation around the expected values, and remains in that quantile as
  covariates change. If `useQuantiles != FALSE` and popGrowthPars
  contains quantiles, those quantiles will be used. If
  `useQuantiles = TRUE` and popGrowthPars does not contain quantiles,
  replicate populations will be assigned to quantiles in the default
  range of 0.025 and 0.975. If `useQuantiles = FALSE`, sampling is done
  independently for each combination of scenario and replicate, so the
  value for a particular replicate population in one scenario is
  unrelated to the values for that replicate in other scenarios. Useful
  for projecting impacts of changing disturbance on the trajectories of
  replicate populations.

- transformFns:

  list of functions used to transform demographic rates. The default is
  `list(S_transform = function(y){(y*46-0.5)/45},R_transform = function(y){y})`.
  The back transformation is applied to survival rates as in Johnson et
  al. 2020.

## Value

For `estimateNationalRate` a similar data frame for one response
variable

A data.frame of predictions. The data.frame includes all columns in
`covTable` with additional columns depending on `returnSample`.

If `returnSample = FALSE` the number of rows is the same as the number
of rows in `covTable`, additional columns are:

- "S_bar" and "R_bar": The mean estimated values of survival and
  recruitment (calves per cow)

- "S_stdErr" and "R_stdErr": Standard error of the estimated values

- "S_PIlow"/"S_PIhigh" and "R_PIlow"/"R_PIhigh": If not using quantiles,
  95\\ minimum values are returned.

If `returnSample = TRUE` the number of rows is
`nrow(covTable) * replicates` additional columns are:

- "scnID": A unique identifier for scenarios provided in `covTable`

- "replicate": A replicate identifier, unique within each scenario

- "S_bar" and "R_bar": The expected values of survival and recruitment
  (calves per cow)

## Details

Each population is optionally assigned to quantiles of the Beta error
distributions for survival and recruitment. Using quantiles means that
the population will stay in these quantiles as disturbance changes over
time, so there is persistent variation in recruitment and survival among
example populations.

A transformation function is also applied to survival to avoid survival
probabilities of 1.

A detailed description of the model is available in [Hughes et al.
(2025)](https://doi.org/10.1016/j.ecoinf.2025.103095)

## References

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
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
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
cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

cfSamps <- sampleNationalCoefs(cfs[[1]], 10)

# disturbance scenarios
distScen <- data.frame(Total_dist = 1:10/10)

# return summary across replicates
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = FALSE)
#>    Total_dist   average     stdErr     PIlow    PIhigh
#> 1         0.1 0.3838513 0.01116631 0.3825237 0.4147978
#> 2         0.2 0.3832760 0.01113087 0.3819419 0.4141116
#> 3         0.3 0.3827015 0.01109575 0.3813611 0.4134264
#> 4         0.4 0.3821279 0.01106095 0.3807811 0.4127425
#> 5         0.5 0.3815551 0.01102646 0.3802020 0.4120596
#> 6         0.6 0.3809832 0.01099228 0.3796238 0.4113779
#> 7         0.7 0.3804122 0.01095842 0.3790464 0.4106973
#> 8         0.8 0.3798420 0.01092487 0.3784700 0.4100178
#> 9         0.9 0.3792726 0.01089163 0.3778944 0.4093395
#> 10        1.0 0.3787041 0.01085871 0.3773197 0.4086623

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3941169
#> 2       1        0.1        V5 0.3833876
#> 3       1        0.1        V9 0.3882085
#> 4       1        0.1        V3 0.3825274
#> 5       1        0.1        V4 0.3989636
#> 6       1        0.1        V8 0.3886642
#> 7       1        0.1        V2 0.4193949
#> 8       1        0.1        V6 0.3889736
#> 9       1        0.1        V7 0.3825226
#> 10      1        0.1       V10 0.3860398
#> 11      2        0.2        V7 0.3819291
#> 12      2        0.2        V8 0.3881199
#> 13      2        0.2        V9 0.3876055
#> 14      2        0.2       V10 0.3854046
#> 15      2        0.2        V1 0.3935995
#> 16      2        0.2        V5 0.3828170
#> 17      2        0.2        V2 0.4187160
#> 18      2        0.2        V3 0.3819862
#> 19      2        0.2        V4 0.3982518
#> 20      2        0.2        V6 0.3883629
#> 21      3        0.3        V4 0.3975413
#> 22      3        0.3        V5 0.3822472
#> 23      3        0.3        V3 0.3814459
#> 24      3        0.3        V7 0.3813364
#> 25      3        0.3        V8 0.3875764
#> 26      3        0.3        V9 0.3870034
#> 27      3        0.3        V6 0.3877531
#> 28      3        0.3       V10 0.3847705
#> 29      3        0.3        V1 0.3930828
#> 30      3        0.3        V2 0.4180382
#> 31      4        0.4        V1 0.3925668
#> 32      4        0.4        V9 0.3864023
#> 33      4        0.4        V3 0.3809063
#> 34      4        0.4        V4 0.3968321
#> 35      4        0.4        V5 0.3816783
#> 36      4        0.4        V2 0.4173616
#> 37      4        0.4        V6 0.3871443
#> 38      4        0.4        V7 0.3807447
#> 39      4        0.4        V8 0.3870336
#> 40      4        0.4       V10 0.3841375
#> 41      5        0.5        V8 0.3864915
#> 42      5        0.5        V9 0.3858020
#> 43      5        0.5       V10 0.3835055
#> 44      5        0.5        V1 0.3920514
#> 45      5        0.5        V5 0.3811102
#> 46      5        0.5        V2 0.4166860
#> 47      5        0.5        V3 0.3803674
#> 48      5        0.5        V4 0.3961242
#> 49      5        0.5        V6 0.3865364
#> 50      5        0.5        V7 0.3801539
#> 51      6        0.6        V4 0.3954175
#> 52      6        0.6        V5 0.3805430
#> 53      6        0.6        V7 0.3795641
#> 54      6        0.6        V8 0.3859502
#> 55      6        0.6        V9 0.3852027
#> 56      6        0.6        V6 0.3859294
#> 57      6        0.6       V10 0.3828745
#> 58      6        0.6        V1 0.3915367
#> 59      6        0.6        V2 0.4160115
#> 60      6        0.6        V3 0.3798294
#> 61      7        0.7        V1 0.3910227
#> 62      7        0.7        V3 0.3792920
#> 63      7        0.7        V4 0.3947121
#> 64      7        0.7        V5 0.3799766
#> 65      7        0.7        V9 0.3846044
#> 66      7        0.7        V6 0.3853235
#> 67      7        0.7        V7 0.3789751
#> 68      7        0.7        V8 0.3854097
#> 69      7        0.7        V2 0.4153381
#> 70      7        0.7       V10 0.3822446
#> 71      8        0.8        V9 0.3840069
#> 72      8        0.8       V10 0.3816157
#> 73      8        0.8        V1 0.3905094
#> 74      8        0.8        V5 0.3794110
#> 75      8        0.8        V2 0.4146658
#> 76      8        0.8        V3 0.3787555
#> 77      8        0.8        V4 0.3940079
#> 78      8        0.8        V8 0.3848699
#> 79      8        0.8        V6 0.3847185
#> 80      8        0.8        V7 0.3783871
#> 81      9        0.9        V5 0.3788463
#> 82      9        0.9        V7 0.3777999
#> 83      9        0.9        V8 0.3843309
#> 84      9        0.9        V9 0.3834104
#> 85      9        0.9        V6 0.3841144
#> 86      9        0.9       V10 0.3809879
#> 87      9        0.9        V1 0.3899967
#> 88      9        0.9        V2 0.4139946
#> 89      9        0.9        V3 0.3782197
#> 90      9        0.9        V4 0.3933050
#> 91     10        1.0        V1 0.3894848
#> 92     10        1.0        V4 0.3926034
#> 93     10        1.0        V5 0.3782825
#> 94     10        1.0        V9 0.3828148
#> 95     10        1.0        V3 0.3776847
#> 96     10        1.0        V7 0.3772137
#> 97     10        1.0        V8 0.3837927
#> 98     10        1.0        V2 0.4133245
#> 99     10        1.0        V6 0.3835113
#> 100    10        1.0       V10 0.3803611

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3941169
#> 2       1        0.1        V5 0.3833876
#> 3       1        0.1        V9 0.3882085
#> 4       1        0.1        V3 0.3825274
#> 5       1        0.1        V4 0.3989636
#> 6       1        0.1        V8 0.3886642
#> 7       1        0.1        V2 0.4193949
#> 8       1        0.1        V6 0.3889736
#> 9       1        0.1        V7 0.3825226
#> 10      1        0.1       V10 0.3860398
#> 11      2        0.2        V7 0.3819291
#> 12      2        0.2        V8 0.3881199
#> 13      2        0.2        V9 0.3876055
#> 14      2        0.2       V10 0.3854046
#> 15      2        0.2        V1 0.3935995
#> 16      2        0.2        V5 0.3828170
#> 17      2        0.2        V2 0.4187160
#> 18      2        0.2        V3 0.3819862
#> 19      2        0.2        V4 0.3982518
#> 20      2        0.2        V6 0.3883629
#> 21      3        0.3        V4 0.3975413
#> 22      3        0.3        V5 0.3822472
#> 23      3        0.3        V3 0.3814459
#> 24      3        0.3        V7 0.3813364
#> 25      3        0.3        V8 0.3875764
#> 26      3        0.3        V9 0.3870034
#> 27      3        0.3        V6 0.3877531
#> 28      3        0.3       V10 0.3847705
#> 29      3        0.3        V1 0.3930828
#> 30      3        0.3        V2 0.4180382
#> 31      4        0.4        V1 0.3925668
#> 32      4        0.4        V9 0.3864023
#> 33      4        0.4        V3 0.3809063
#> 34      4        0.4        V4 0.3968321
#> 35      4        0.4        V5 0.3816783
#> 36      4        0.4        V2 0.4173616
#> 37      4        0.4        V6 0.3871443
#> 38      4        0.4        V7 0.3807447
#> 39      4        0.4        V8 0.3870336
#> 40      4        0.4       V10 0.3841375
#> 41      5        0.5        V8 0.3864915
#> 42      5        0.5        V9 0.3858020
#> 43      5        0.5       V10 0.3835055
#> 44      5        0.5        V1 0.3920514
#> 45      5        0.5        V5 0.3811102
#> 46      5        0.5        V2 0.4166860
#> 47      5        0.5        V3 0.3803674
#> 48      5        0.5        V4 0.3961242
#> 49      5        0.5        V6 0.3865364
#> 50      5        0.5        V7 0.3801539
#> 51      6        0.6        V4 0.3954175
#> 52      6        0.6        V5 0.3805430
#> 53      6        0.6        V7 0.3795641
#> 54      6        0.6        V8 0.3859502
#> 55      6        0.6        V9 0.3852027
#> 56      6        0.6        V6 0.3859294
#> 57      6        0.6       V10 0.3828745
#> 58      6        0.6        V1 0.3915367
#> 59      6        0.6        V2 0.4160115
#> 60      6        0.6        V3 0.3798294
#> 61      7        0.7        V1 0.3910227
#> 62      7        0.7        V3 0.3792920
#> 63      7        0.7        V4 0.3947121
#> 64      7        0.7        V5 0.3799766
#> 65      7        0.7        V9 0.3846044
#> 66      7        0.7        V6 0.3853235
#> 67      7        0.7        V7 0.3789751
#> 68      7        0.7        V8 0.3854097
#> 69      7        0.7        V2 0.4153381
#> 70      7        0.7       V10 0.3822446
#> 71      8        0.8        V9 0.3840069
#> 72      8        0.8       V10 0.3816157
#> 73      8        0.8        V1 0.3905094
#> 74      8        0.8        V5 0.3794110
#> 75      8        0.8        V2 0.4146658
#> 76      8        0.8        V3 0.3787555
#> 77      8        0.8        V4 0.3940079
#> 78      8        0.8        V8 0.3848699
#> 79      8        0.8        V6 0.3847185
#> 80      8        0.8        V7 0.3783871
#> 81      9        0.9        V5 0.3788463
#> 82      9        0.9        V7 0.3777999
#> 83      9        0.9        V8 0.3843309
#> 84      9        0.9        V9 0.3834104
#> 85      9        0.9        V6 0.3841144
#> 86      9        0.9       V10 0.3809879
#> 87      9        0.9        V1 0.3899967
#> 88      9        0.9        V2 0.4139946
#> 89      9        0.9        V3 0.3782197
#> 90      9        0.9        V4 0.3933050
#> 91     10        1.0        V1 0.3894848
#> 92     10        1.0        V4 0.3926034
#> 93     10        1.0        V5 0.3782825
#> 94     10        1.0        V9 0.3828148
#> 95     10        1.0        V3 0.3776847
#> 96     10        1.0        V7 0.3772137
#> 97     10        1.0        V8 0.3837927
#> 98     10        1.0        V2 0.4133245
#> 99     10        1.0        V6 0.3835113
#> 100    10        1.0       V10 0.3803611


# get coefficient samples
coefs <- getNationalCoefficients(10)

# table of different scenarios to test
covTableSim <- expand.grid(Anthro = seq(0, 90, by = 20),
                           Fire_excl_anthro = seq(0, 70, by = 20))
covTableSim$Total_dist = covTableSim$Anthro + covTableSim$Fire_excl_anthro

estimateNationalRates(covTableSim, coefs)
#> popGrowthPars contains quantiles so they are used instead of the defaults
#> popGrowthPars contains quantiles so they are used instead of the defaults
#>    Anthro Fire_excl_anthro Total_dist     S_bar   S_stdErr   S_PIlow  S_PIhigh
#> 1       0                0          0 0.8757906 0.04668115 0.7941926 0.9423547
#> 2       0               20         20 0.8757906 0.04668115 0.7941926 0.9423547
#> 3       0               40         40 0.8757906 0.04668115 0.7941926 0.9423547
#> 4       0               60         60 0.8757906 0.04668115 0.7941926 0.9423547
#> 5      20                0         20 0.8617131 0.04780900 0.7785798 0.9337019
#> 6      20               20         40 0.8617131 0.04780900 0.7785798 0.9337019
#> 7      20               40         60 0.8617131 0.04780900 0.7785798 0.9337019
#> 8      20               60         80 0.8617131 0.04780900 0.7785798 0.9337019
#> 9      40                0         40 0.8478591 0.04884037 0.7634280 0.9249523
#> 10     40               20         60 0.8478591 0.04884037 0.7634280 0.9249523
#> 11     40               40         80 0.8478591 0.04884037 0.7634280 0.9249523
#> 12     40               60        100 0.8478591 0.04884037 0.7634280 0.9249523
#> 13     60                0         60 0.8342249 0.04978961 0.7486983 0.9161344
#> 14     60               20         80 0.8342249 0.04978961 0.7486983 0.9161344
#> 15     60               40        100 0.8342249 0.04978961 0.7486983 0.9161344
#> 16     60               60        120 0.8342249 0.04978961 0.7486983 0.9161344
#> 17     80                0         80 0.8208071 0.05066783 0.7343592 0.9072708
#> 18     80               20        100 0.8208071 0.05066783 0.7343592 0.9072708
#> 19     80               40        120 0.8208071 0.05066783 0.7343592 0.9072708
#> 20     80               60        140 0.8208071 0.05066783 0.7343592 0.9072708
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.10713747 0.207439516 0.5671131
#> 2  0.30574618 0.10202780 0.174230180 0.4917608
#> 3  0.26001915 0.09863408 0.145829170 0.4270945
#> 4  0.22113100 0.09619349 0.121560054 0.3718042
#> 5  0.25589195 0.10054830 0.118607548 0.4587429
#> 6  0.21762106 0.09390219 0.098336498 0.3988482
#> 7  0.18507391 0.08880436 0.081083470 0.3476820
#> 8  0.15739448 0.08472413 0.066441505 0.3039881
#> 9  0.18213629 0.09112605 0.064668825 0.3727066
#> 10 0.15489621 0.08420332 0.052564403 0.3253606
#> 11 0.13173012 0.07858234 0.042379348 0.2849135
#> 12 0.11202872 0.07388967 0.033854602 0.2503001
#> 13 0.12963921 0.08080803 0.032831923 0.3047014
#> 14 0.11025053 0.07418588 0.025917059 0.2672441
#> 15 0.09376159 0.06862823 0.020214299 0.2351482
#> 16 0.07973872 0.06386740 0.015552971 0.2075523
#> 17 0.09227334 0.07062497 0.015002389 0.2508660
#> 18 0.07847305 0.06454468 0.011340446 0.2210800
#> 19 0.06673672 0.05932328 0.008420203 0.1954173
#> 20 0.05675565 0.05476577 0.006125794 0.1731973
```
