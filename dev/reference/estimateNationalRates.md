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
#> 1         0.1 0.3838513 0.02302574 0.3420799 0.4097783
#> 2         0.2 0.3832760 0.02298623 0.3416728 0.4092395
#> 3         0.3 0.3827015 0.02294699 0.3412661 0.4087015
#> 4         0.4 0.3821279 0.02290802 0.3408600 0.4081642
#> 5         0.5 0.3815551 0.02286932 0.3404544 0.4076276
#> 6         0.6 0.3809832 0.02283089 0.3400493 0.4070917
#> 7         0.7 0.3804122 0.02279273 0.3396447 0.4065566
#> 8         0.8 0.3798420 0.02275484 0.3392405 0.4060221
#> 9         0.9 0.3792726 0.02271721 0.3388369 0.4054884
#> 10        1.0 0.3787041 0.02267985 0.3384338 0.4049553

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4069748
#> 2       1        0.1        V5 0.3641097
#> 3       1        0.1        V9 0.3883336
#> 4       1        0.1        V3 0.3812044
#> 5       1        0.1        V4 0.3744925
#> 6       1        0.1        V8 0.3880021
#> 7       1        0.1        V2 0.4105922
#> 8       1        0.1        V6 0.3448941
#> 9       1        0.1        V7 0.3412628
#> 10      1        0.1       V10 0.3761226
#> 11      2        0.2        V7 0.3408927
#> 12      2        0.2        V8 0.3874544
#> 13      2        0.2        V9 0.3877462
#> 14      2        0.2       V10 0.3755125
#> 15      2        0.2        V1 0.4063325
#> 16      2        0.2        V5 0.3635502
#> 17      2        0.2        V2 0.4100835
#> 18      2        0.2        V3 0.3806649
#> 19      2        0.2        V4 0.3738430
#> 20      2        0.2        V6 0.3443598
#> 21      3        0.3        V4 0.3731947
#> 22      3        0.3        V5 0.3629917
#> 23      3        0.3        V3 0.3801261
#> 24      3        0.3        V7 0.3405229
#> 25      3        0.3        V8 0.3869075
#> 26      3        0.3        V9 0.3871596
#> 27      3        0.3        V6 0.3438262
#> 28      3        0.3       V10 0.3749033
#> 29      3        0.3        V1 0.4056911
#> 30      3        0.3        V2 0.4095755
#> 31      4        0.4        V1 0.4050508
#> 32      4        0.4        V9 0.3865739
#> 33      4        0.4        V3 0.3795881
#> 34      4        0.4        V4 0.3725474
#> 35      4        0.4        V5 0.3624339
#> 36      4        0.4        V2 0.4090681
#> 37      4        0.4        V6 0.3432935
#> 38      4        0.4        V7 0.3401536
#> 39      4        0.4        V8 0.3863613
#> 40      4        0.4       V10 0.3742952
#> 41      5        0.5        V8 0.3858160
#> 42      5        0.5        V9 0.3859891
#> 43      5        0.5       V10 0.3736880
#> 44      5        0.5        V1 0.4044115
#> 45      5        0.5        V5 0.3618771
#> 46      5        0.5        V2 0.4085613
#> 47      5        0.5        V3 0.3790508
#> 48      5        0.5        V4 0.3719013
#> 49      5        0.5        V6 0.3427616
#> 50      5        0.5        V7 0.3397846
#> 51      6        0.6        V4 0.3712563
#> 52      6        0.6        V5 0.3613211
#> 53      6        0.6        V7 0.3394161
#> 54      6        0.6        V8 0.3852714
#> 55      6        0.6        V9 0.3854052
#> 56      6        0.6        V6 0.3422305
#> 57      6        0.6       V10 0.3730819
#> 58      6        0.6        V1 0.4037732
#> 59      6        0.6        V2 0.4080552
#> 60      6        0.6        V3 0.3785143
#> 61      7        0.7        V1 0.4031359
#> 62      7        0.7        V3 0.3779786
#> 63      7        0.7        V4 0.3706125
#> 64      7        0.7        V5 0.3607659
#> 65      7        0.7        V9 0.3848221
#> 66      7        0.7        V6 0.3417002
#> 67      7        0.7        V7 0.3390479
#> 68      7        0.7        V8 0.3847275
#> 69      7        0.7        V2 0.4075496
#> 70      7        0.7       V10 0.3724767
#> 71      8        0.8        V9 0.3842400
#> 72      8        0.8       V10 0.3718725
#> 73      8        0.8        V1 0.4024996
#> 74      8        0.8        V5 0.3602116
#> 75      8        0.8        V2 0.4070447
#> 76      8        0.8        V3 0.3774436
#> 77      8        0.8        V4 0.3699697
#> 78      8        0.8        V8 0.3841845
#> 79      8        0.8        V6 0.3411708
#> 80      8        0.8        V7 0.3386801
#> 81      9        0.9        V5 0.3596581
#> 82      9        0.9        V7 0.3383128
#> 83      9        0.9        V8 0.3836422
#> 84      9        0.9        V9 0.3836587
#> 85      9        0.9        V6 0.3406421
#> 86      9        0.9       V10 0.3712693
#> 87      9        0.9        V1 0.4018644
#> 88      9        0.9        V2 0.4065405
#> 89      9        0.9        V3 0.3769094
#> 90      9        0.9        V4 0.3693281
#> 91     10        1.0        V1 0.4012301
#> 92     10        1.0        V4 0.3686876
#> 93     10        1.0        V5 0.3591055
#> 94     10        1.0        V9 0.3830783
#> 95     10        1.0        V3 0.3763760
#> 96     10        1.0        V7 0.3379458
#> 97     10        1.0        V8 0.3831006
#> 98     10        1.0        V2 0.4060368
#> 99     10        1.0        V6 0.3401143
#> 100    10        1.0       V10 0.3706670

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4069748
#> 2       1        0.1        V5 0.3641097
#> 3       1        0.1        V9 0.3883336
#> 4       1        0.1        V3 0.3812044
#> 5       1        0.1        V4 0.3744925
#> 6       1        0.1        V8 0.3880021
#> 7       1        0.1        V2 0.4105922
#> 8       1        0.1        V6 0.3448941
#> 9       1        0.1        V7 0.3412628
#> 10      1        0.1       V10 0.3761226
#> 11      2        0.2        V7 0.3408927
#> 12      2        0.2        V8 0.3874544
#> 13      2        0.2        V9 0.3877462
#> 14      2        0.2       V10 0.3755125
#> 15      2        0.2        V1 0.4063325
#> 16      2        0.2        V5 0.3635502
#> 17      2        0.2        V2 0.4100835
#> 18      2        0.2        V3 0.3806649
#> 19      2        0.2        V4 0.3738430
#> 20      2        0.2        V6 0.3443598
#> 21      3        0.3        V4 0.3731947
#> 22      3        0.3        V5 0.3629917
#> 23      3        0.3        V3 0.3801261
#> 24      3        0.3        V7 0.3405229
#> 25      3        0.3        V8 0.3869075
#> 26      3        0.3        V9 0.3871596
#> 27      3        0.3        V6 0.3438262
#> 28      3        0.3       V10 0.3749033
#> 29      3        0.3        V1 0.4056911
#> 30      3        0.3        V2 0.4095755
#> 31      4        0.4        V1 0.4050508
#> 32      4        0.4        V9 0.3865739
#> 33      4        0.4        V3 0.3795881
#> 34      4        0.4        V4 0.3725474
#> 35      4        0.4        V5 0.3624339
#> 36      4        0.4        V2 0.4090681
#> 37      4        0.4        V6 0.3432935
#> 38      4        0.4        V7 0.3401536
#> 39      4        0.4        V8 0.3863613
#> 40      4        0.4       V10 0.3742952
#> 41      5        0.5        V8 0.3858160
#> 42      5        0.5        V9 0.3859891
#> 43      5        0.5       V10 0.3736880
#> 44      5        0.5        V1 0.4044115
#> 45      5        0.5        V5 0.3618771
#> 46      5        0.5        V2 0.4085613
#> 47      5        0.5        V3 0.3790508
#> 48      5        0.5        V4 0.3719013
#> 49      5        0.5        V6 0.3427616
#> 50      5        0.5        V7 0.3397846
#> 51      6        0.6        V4 0.3712563
#> 52      6        0.6        V5 0.3613211
#> 53      6        0.6        V7 0.3394161
#> 54      6        0.6        V8 0.3852714
#> 55      6        0.6        V9 0.3854052
#> 56      6        0.6        V6 0.3422305
#> 57      6        0.6       V10 0.3730819
#> 58      6        0.6        V1 0.4037732
#> 59      6        0.6        V2 0.4080552
#> 60      6        0.6        V3 0.3785143
#> 61      7        0.7        V1 0.4031359
#> 62      7        0.7        V3 0.3779786
#> 63      7        0.7        V4 0.3706125
#> 64      7        0.7        V5 0.3607659
#> 65      7        0.7        V9 0.3848221
#> 66      7        0.7        V6 0.3417002
#> 67      7        0.7        V7 0.3390479
#> 68      7        0.7        V8 0.3847275
#> 69      7        0.7        V2 0.4075496
#> 70      7        0.7       V10 0.3724767
#> 71      8        0.8        V9 0.3842400
#> 72      8        0.8       V10 0.3718725
#> 73      8        0.8        V1 0.4024996
#> 74      8        0.8        V5 0.3602116
#> 75      8        0.8        V2 0.4070447
#> 76      8        0.8        V3 0.3774436
#> 77      8        0.8        V4 0.3699697
#> 78      8        0.8        V8 0.3841845
#> 79      8        0.8        V6 0.3411708
#> 80      8        0.8        V7 0.3386801
#> 81      9        0.9        V5 0.3596581
#> 82      9        0.9        V7 0.3383128
#> 83      9        0.9        V8 0.3836422
#> 84      9        0.9        V9 0.3836587
#> 85      9        0.9        V6 0.3406421
#> 86      9        0.9       V10 0.3712693
#> 87      9        0.9        V1 0.4018644
#> 88      9        0.9        V2 0.4065405
#> 89      9        0.9        V3 0.3769094
#> 90      9        0.9        V4 0.3693281
#> 91     10        1.0        V1 0.4012301
#> 92     10        1.0        V4 0.3686876
#> 93     10        1.0        V5 0.3591055
#> 94     10        1.0        V9 0.3830783
#> 95     10        1.0        V3 0.3763760
#> 96     10        1.0        V7 0.3379458
#> 97     10        1.0        V8 0.3831006
#> 98     10        1.0        V2 0.4060368
#> 99     10        1.0        V6 0.3401143
#> 100    10        1.0       V10 0.3706670


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
#> 1       0                0          0 0.8757906 0.04946405 0.7834086 0.9523261
#> 2       0               20         20 0.8757906 0.04946405 0.7834086 0.9523261
#> 3       0               40         40 0.8757906 0.04946405 0.7834086 0.9523261
#> 4       0               60         60 0.8757906 0.04946405 0.7834086 0.9523261
#> 5      20                0         20 0.8617131 0.05277263 0.7650697 0.9441833
#> 6      20               20         40 0.8617131 0.05277263 0.7650697 0.9441833
#> 7      20               40         60 0.8617131 0.05277263 0.7650697 0.9441833
#> 8      20               60         80 0.8617131 0.05277263 0.7650697 0.9441833
#> 9      40                0         40 0.8478591 0.05587881 0.7473631 0.9359011
#> 10     40               20         60 0.8478591 0.05587881 0.7473631 0.9359011
#> 11     40               40         80 0.8478591 0.05587881 0.7473631 0.9359011
#> 12     40               60        100 0.8478591 0.05587881 0.7473631 0.9359011
#> 13     60                0         60 0.8342249 0.05880321 0.7302298 0.9275146
#> 14     60               20         80 0.8342249 0.05880321 0.7302298 0.9275146
#> 15     60               40        100 0.8342249 0.05880321 0.7302298 0.9275146
#> 16     60               60        120 0.8342249 0.05880321 0.7302298 0.9275146
#> 17     80                0         80 0.8208071 0.06156263 0.7136226 0.9190513
#> 18     80               20        100 0.8208071 0.06156263 0.7136226 0.9190513
#> 19     80               40        120 0.8208071 0.06156263 0.7136226 0.9190513
#> 20     80               60        140 0.8208071 0.06156263 0.7136226 0.9190513
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13244719 0.166510010 0.6105063
#> 2  0.30574618 0.12589314 0.126014582 0.5448693
#> 3  0.26001915 0.11964093 0.094310459 0.4865799
#> 4  0.22113100 0.11327336 0.069620337 0.4350411
#> 5  0.25589195 0.11759259 0.094500542 0.4969796
#> 6  0.21762106 0.11006160 0.069767832 0.4442252
#> 7  0.18507391 0.10319903 0.050656306 0.3976715
#> 8  0.15739448 0.09663619 0.036044967 0.3566325
#> 9  0.18213629 0.10294774 0.050769907 0.4059634
#> 10 0.15489621 0.09549955 0.036131241 0.3639412
#> 11 0.13173012 0.08879470 0.025093900 0.3268968
#> 12 0.11202872 0.08254174 0.016920963 0.2942129
#> 13 0.12963921 0.08933309 0.025158514 0.3334954
#> 14 0.11025053 0.08243304 0.016968276 0.3000381
#> 15 0.09376159 0.07623753 0.011041785 0.2704811
#> 16 0.07973872 0.07052449 0.006878981 0.2443112
#> 17 0.09227334 0.07705374 0.011075607 0.2757526
#> 18 0.07847305 0.07087904 0.006902302 0.2489838
#> 19 0.06673672 0.06532868 0.004079513 0.2252272
#> 20 0.05675565 0.06023519 0.002260411 0.2040706
```
