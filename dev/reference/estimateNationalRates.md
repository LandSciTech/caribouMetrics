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
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`demographyDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographyDefaults.md),
[`disturbanceDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceDefaults.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`monitoringDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/monitoringDefaults.md),
[`nationalTrajectoryDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/nationalTrajectoryDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`timeDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/timeDefaults.md),
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
#> 1         0.1 0.3838513 0.02142103 0.3585851 0.4159192
#> 2         0.2 0.3832760 0.02136499 0.3581126 0.4151887
#> 3         0.3 0.3827015 0.02130939 0.3576407 0.4144595
#> 4         0.4 0.3821279 0.02125422 0.3571694 0.4137316
#> 5         0.5 0.3815551 0.02119948 0.3566987 0.4130050
#> 6         0.6 0.3809832 0.02114516 0.3562287 0.4122796
#> 7         0.7 0.3804122 0.02109128 0.3557593 0.4115556
#> 8         0.8 0.3798420 0.02103782 0.3552905 0.4108329
#> 9         0.9 0.3792726 0.02098479 0.3548223 0.4101115
#> 10        1.0 0.3787041 0.02093218 0.3543548 0.4093913

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3843204
#> 2       1        0.1        V5 0.3739737
#> 3       1        0.1        V9 0.4182352
#> 4       1        0.1        V3 0.3655344
#> 5       1        0.1        V4 0.3608895
#> 6       1        0.1        V8 0.3866306
#> 7       1        0.1        V2 0.4079418
#> 8       1        0.1        V6 0.3981894
#> 9       1        0.1        V7 0.3579160
#> 10      1        0.1       V10 0.4077916
#> 11      2        0.2        V7 0.3574441
#> 12      2        0.2        V8 0.3860997
#> 13      2        0.2        V9 0.4174501
#> 14      2        0.2       V10 0.4071548
#> 15      2        0.2        V1 0.3836421
#> 16      2        0.2        V5 0.3733473
#> 17      2        0.2        V2 0.4073993
#> 18      2        0.2        V3 0.3649339
#> 19      2        0.2        V4 0.3604149
#> 20      2        0.2        V6 0.3976663
#> 21      3        0.3        V4 0.3599408
#> 22      3        0.3        V5 0.3727221
#> 23      3        0.3        V3 0.3643345
#> 24      3        0.3        V7 0.3569729
#> 25      3        0.3        V8 0.3855695
#> 26      3        0.3        V9 0.4166665
#> 27      3        0.3        V6 0.3971440
#> 28      3        0.3       V10 0.4065189
#> 29      3        0.3        V1 0.3829650
#> 30      3        0.3        V2 0.4068576
#> 31      4        0.4        V1 0.3822890
#> 32      4        0.4        V9 0.4158843
#> 33      4        0.4        V3 0.3637360
#> 34      4        0.4        V4 0.3594674
#> 35      4        0.4        V5 0.3720978
#> 36      4        0.4        V2 0.4063165
#> 37      4        0.4        V6 0.3966223
#> 38      4        0.4        V7 0.3565022
#> 39      4        0.4        V8 0.3850401
#> 40      4        0.4       V10 0.4058841
#> 41      5        0.5        V8 0.3845114
#> 42      5        0.5        V9 0.4151036
#> 43      5        0.5       V10 0.4052502
#> 44      5        0.5        V1 0.3816143
#> 45      5        0.5        V5 0.3714746
#> 46      5        0.5        V2 0.4057762
#> 47      5        0.5        V3 0.3631384
#> 48      5        0.5        V4 0.3589946
#> 49      5        0.5        V6 0.3961013
#> 50      5        0.5        V7 0.3560322
#> 51      6        0.6        V4 0.3585224
#> 52      6        0.6        V5 0.3708525
#> 53      6        0.6        V7 0.3555628
#> 54      6        0.6        V8 0.3839834
#> 55      6        0.6        V9 0.4143244
#> 56      6        0.6        V6 0.3955809
#> 57      6        0.6       V10 0.4046173
#> 58      6        0.6        V1 0.3809408
#> 59      6        0.6        V2 0.4052366
#> 60      6        0.6        V3 0.3625419
#> 61      7        0.7        V1 0.3802684
#> 62      7        0.7        V3 0.3619463
#> 63      7        0.7        V4 0.3580509
#> 64      7        0.7        V5 0.3702314
#> 65      7        0.7        V9 0.4135466
#> 66      7        0.7        V6 0.3950613
#> 67      7        0.7        V7 0.3550940
#> 68      7        0.7        V8 0.3834561
#> 69      7        0.7        V2 0.4046977
#> 70      7        0.7       V10 0.4039855
#> 71      8        0.8        V9 0.4127703
#> 72      8        0.8       V10 0.4033546
#> 73      8        0.8        V1 0.3795973
#> 74      8        0.8        V5 0.3696113
#> 75      8        0.8        V2 0.4041595
#> 76      8        0.8        V3 0.3613518
#> 77      8        0.8        V4 0.3575800
#> 78      8        0.8        V8 0.3829296
#> 79      8        0.8        V6 0.3945423
#> 80      8        0.8        V7 0.3546258
#> 81      9        0.9        V5 0.3689923
#> 82      9        0.9        V7 0.3541583
#> 83      9        0.9        V8 0.3824038
#> 84      9        0.9        V9 0.4119955
#> 85      9        0.9        V6 0.3940241
#> 86      9        0.9       V10 0.4027247
#> 87      9        0.9        V1 0.3789273
#> 88      9        0.9        V2 0.4036221
#> 89      9        0.9        V3 0.3607581
#> 90      9        0.9        V4 0.3571097
#> 91     10        1.0        V1 0.3782585
#> 92     10        1.0        V4 0.3566400
#> 93     10        1.0        V5 0.3683743
#> 94     10        1.0        V9 0.4112221
#> 95     10        1.0        V3 0.3601655
#> 96     10        1.0        V7 0.3536913
#> 97     10        1.0        V8 0.3818787
#> 98     10        1.0        V2 0.4030853
#> 99     10        1.0        V6 0.3935065
#> 100    10        1.0       V10 0.4020957

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3843204
#> 2       1        0.1        V5 0.3739737
#> 3       1        0.1        V9 0.4182352
#> 4       1        0.1        V3 0.3655344
#> 5       1        0.1        V4 0.3608895
#> 6       1        0.1        V8 0.3866306
#> 7       1        0.1        V2 0.4079418
#> 8       1        0.1        V6 0.3981894
#> 9       1        0.1        V7 0.3579160
#> 10      1        0.1       V10 0.4077916
#> 11      2        0.2        V7 0.3574441
#> 12      2        0.2        V8 0.3860997
#> 13      2        0.2        V9 0.4174501
#> 14      2        0.2       V10 0.4071548
#> 15      2        0.2        V1 0.3836421
#> 16      2        0.2        V5 0.3733473
#> 17      2        0.2        V2 0.4073993
#> 18      2        0.2        V3 0.3649339
#> 19      2        0.2        V4 0.3604149
#> 20      2        0.2        V6 0.3976663
#> 21      3        0.3        V4 0.3599408
#> 22      3        0.3        V5 0.3727221
#> 23      3        0.3        V3 0.3643345
#> 24      3        0.3        V7 0.3569729
#> 25      3        0.3        V8 0.3855695
#> 26      3        0.3        V9 0.4166665
#> 27      3        0.3        V6 0.3971440
#> 28      3        0.3       V10 0.4065189
#> 29      3        0.3        V1 0.3829650
#> 30      3        0.3        V2 0.4068576
#> 31      4        0.4        V1 0.3822890
#> 32      4        0.4        V9 0.4158843
#> 33      4        0.4        V3 0.3637360
#> 34      4        0.4        V4 0.3594674
#> 35      4        0.4        V5 0.3720978
#> 36      4        0.4        V2 0.4063165
#> 37      4        0.4        V6 0.3966223
#> 38      4        0.4        V7 0.3565022
#> 39      4        0.4        V8 0.3850401
#> 40      4        0.4       V10 0.4058841
#> 41      5        0.5        V8 0.3845114
#> 42      5        0.5        V9 0.4151036
#> 43      5        0.5       V10 0.4052502
#> 44      5        0.5        V1 0.3816143
#> 45      5        0.5        V5 0.3714746
#> 46      5        0.5        V2 0.4057762
#> 47      5        0.5        V3 0.3631384
#> 48      5        0.5        V4 0.3589946
#> 49      5        0.5        V6 0.3961013
#> 50      5        0.5        V7 0.3560322
#> 51      6        0.6        V4 0.3585224
#> 52      6        0.6        V5 0.3708525
#> 53      6        0.6        V7 0.3555628
#> 54      6        0.6        V8 0.3839834
#> 55      6        0.6        V9 0.4143244
#> 56      6        0.6        V6 0.3955809
#> 57      6        0.6       V10 0.4046173
#> 58      6        0.6        V1 0.3809408
#> 59      6        0.6        V2 0.4052366
#> 60      6        0.6        V3 0.3625419
#> 61      7        0.7        V1 0.3802684
#> 62      7        0.7        V3 0.3619463
#> 63      7        0.7        V4 0.3580509
#> 64      7        0.7        V5 0.3702314
#> 65      7        0.7        V9 0.4135466
#> 66      7        0.7        V6 0.3950613
#> 67      7        0.7        V7 0.3550940
#> 68      7        0.7        V8 0.3834561
#> 69      7        0.7        V2 0.4046977
#> 70      7        0.7       V10 0.4039855
#> 71      8        0.8        V9 0.4127703
#> 72      8        0.8       V10 0.4033546
#> 73      8        0.8        V1 0.3795973
#> 74      8        0.8        V5 0.3696113
#> 75      8        0.8        V2 0.4041595
#> 76      8        0.8        V3 0.3613518
#> 77      8        0.8        V4 0.3575800
#> 78      8        0.8        V8 0.3829296
#> 79      8        0.8        V6 0.3945423
#> 80      8        0.8        V7 0.3546258
#> 81      9        0.9        V5 0.3689923
#> 82      9        0.9        V7 0.3541583
#> 83      9        0.9        V8 0.3824038
#> 84      9        0.9        V9 0.4119955
#> 85      9        0.9        V6 0.3940241
#> 86      9        0.9       V10 0.4027247
#> 87      9        0.9        V1 0.3789273
#> 88      9        0.9        V2 0.4036221
#> 89      9        0.9        V3 0.3607581
#> 90      9        0.9        V4 0.3571097
#> 91     10        1.0        V1 0.3782585
#> 92     10        1.0        V4 0.3566400
#> 93     10        1.0        V5 0.3683743
#> 94     10        1.0        V9 0.4112221
#> 95     10        1.0        V3 0.3601655
#> 96     10        1.0        V7 0.3536913
#> 97     10        1.0        V8 0.3818787
#> 98     10        1.0        V2 0.4030853
#> 99     10        1.0        V6 0.3935065
#> 100    10        1.0       V10 0.4020957


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
#> 1       0                0          0 0.8757906 0.05114618 0.7783507 0.9443703
#> 2       0               20         20 0.8757906 0.05114618 0.7783507 0.9443703
#> 3       0               40         40 0.8757906 0.05114618 0.7783507 0.9443703
#> 4       0               60         60 0.8757906 0.05114618 0.7783507 0.9443703
#> 5      20                0         20 0.8617131 0.05253158 0.7645491 0.9364112
#> 6      20               20         40 0.8617131 0.05253158 0.7645491 0.9364112
#> 7      20               40         60 0.8617131 0.05253158 0.7645491 0.9364112
#> 8      20               60         80 0.8617131 0.05253158 0.7645491 0.9364112
#> 9      40                0         40 0.8478591 0.05386337 0.7511241 0.9283578
#> 10     40               20         60 0.8478591 0.05386337 0.7511241 0.9283578
#> 11     40               40         80 0.8478591 0.05386337 0.7511241 0.9283578
#> 12     40               60        100 0.8478591 0.05386337 0.7511241 0.9283578
#> 13     60                0         60 0.8342249 0.05515196 0.7380462 0.9202346
#> 14     60               20         80 0.8342249 0.05515196 0.7380462 0.9202346
#> 15     60               40        100 0.8342249 0.05515196 0.7380462 0.9202346
#> 16     60               60        120 0.8342249 0.05515196 0.7380462 0.9202346
#> 17     80                0         80 0.8208071 0.05640363 0.7252910 0.9120615
#> 18     80               20        100 0.8208071 0.05640363 0.7252910 0.9120615
#> 19     80               40        120 0.8208071 0.05640363 0.7252910 0.9120615
#> 20     80               60        140 0.8208071 0.05640363 0.7252910 0.9120615
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11705141 0.183115860 0.5745534
#> 2  0.30574618 0.11371140 0.127332631 0.5207941
#> 3  0.26001915 0.10913084 0.086876183 0.4722948
#> 4  0.22113100 0.10365861 0.057814288 0.4286524
#> 5  0.25589195 0.10947441 0.101035530 0.4721288
#> 6  0.21762106 0.10380705 0.067943129 0.4285032
#> 7  0.18507391 0.09776171 0.044377009 0.3893071
#> 8  0.15739448 0.09152569 0.027919199 0.3541199
#> 9  0.18213629 0.09909837 0.052561113 0.3891731
#> 10 0.15489621 0.09267963 0.033589995 0.3539997
#> 11 0.13173012 0.08632462 0.020548962 0.3224327
#> 12 0.11202872 0.08013398 0.011883270 0.2940981
#> 13 0.12963921 0.08782764 0.025014116 0.3223248
#> 14 0.11025053 0.08147858 0.014810206 0.2940013
#> 15 0.09376159 0.07541545 0.008210113 0.2685647
#> 16 0.07973872 0.06968874 0.004179145 0.2457016
#> 17 0.09227334 0.07675663 0.010415191 0.2684777
#> 18 0.07847305 0.07088850 0.005494647 0.2456234
#> 19 0.06673672 0.06540820 0.002622757 0.2250576
#> 20 0.05675565 0.06033423 0.001099626 0.2065247
```
