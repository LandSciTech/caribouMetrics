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
#> 1         0.1 0.3838513 0.02828098 0.3377376 0.4292900
#> 2         0.2 0.3832760 0.02825374 0.3372524 0.4286995
#> 3         0.3 0.3827015 0.02822655 0.3367680 0.4281097
#> 4         0.4 0.3821279 0.02819940 0.3362842 0.4275207
#> 5         0.5 0.3815551 0.02817229 0.3358011 0.4269326
#> 6         0.6 0.3809832 0.02814522 0.3353187 0.4263452
#> 7         0.7 0.3804122 0.02811820 0.3348370 0.4257587
#> 8         0.8 0.3798420 0.02809122 0.3343560 0.4251730
#> 9         0.9 0.3792726 0.02806429 0.3338757 0.4245880
#> 10        1.0 0.3787041 0.02803739 0.3333961 0.4240039

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3874884
#> 2       1        0.1        V5 0.3788389
#> 3       1        0.1        V9 0.3551936
#> 4       1        0.1        V3 0.3592282
#> 5       1        0.1        V4 0.4412645
#> 6       1        0.1        V8 0.3724569
#> 7       1        0.1        V2 0.3690937
#> 8       1        0.1        V6 0.3326698
#> 9       1        0.1        V7 0.3880446
#> 10      1        0.1       V10 0.3746163
#> 11      2        0.2        V7 0.3874895
#> 12      2        0.2        V8 0.3718921
#> 13      2        0.2        V9 0.3546603
#> 14      2        0.2       V10 0.3741004
#> 15      2        0.2        V1 0.3869349
#> 16      2        0.2        V5 0.3782802
#> 17      2        0.2        V2 0.3685823
#> 18      2        0.2        V3 0.3586466
#> 19      2        0.2        V4 0.4406636
#> 20      2        0.2        V6 0.3321986
#> 21      3        0.3        V4 0.4400636
#> 22      3        0.3        V5 0.3777222
#> 23      3        0.3        V3 0.3580660
#> 24      3        0.3        V7 0.3869352
#> 25      3        0.3        V8 0.3713282
#> 26      3        0.3        V9 0.3541278
#> 27      3        0.3        V6 0.3317280
#> 28      3        0.3       V10 0.3735853
#> 29      3        0.3        V1 0.3863822
#> 30      3        0.3        V2 0.3680717
#> 31      4        0.4        V1 0.3858302
#> 32      4        0.4        V9 0.3535961
#> 33      4        0.4        V3 0.3574863
#> 34      4        0.4        V4 0.4394643
#> 35      4        0.4        V5 0.3771651
#> 36      4        0.4        V2 0.3675617
#> 37      4        0.4        V6 0.3312582
#> 38      4        0.4        V7 0.3863816
#> 39      4        0.4        V8 0.3707652
#> 40      4        0.4       V10 0.3730709
#> 41      5        0.5        V8 0.3702030
#> 42      5        0.5        V9 0.3530651
#> 43      5        0.5       V10 0.3725572
#> 44      5        0.5        V1 0.3852791
#> 45      5        0.5        V5 0.3766089
#> 46      5        0.5        V2 0.3670524
#> 47      5        0.5        V3 0.3569076
#> 48      5        0.5        V4 0.4388659
#> 49      5        0.5        V6 0.3307890
#> 50      5        0.5        V7 0.3858289
#> 51      6        0.6        V4 0.4382683
#> 52      6        0.6        V5 0.3760534
#> 53      6        0.6        V7 0.3852770
#> 54      6        0.6        V8 0.3696416
#> 55      6        0.6        V9 0.3525350
#> 56      6        0.6        V6 0.3303204
#> 57      6        0.6       V10 0.3720442
#> 58      6        0.6        V1 0.3847287
#> 59      6        0.6        V2 0.3665439
#> 60      6        0.6        V3 0.3563297
#> 61      7        0.7        V1 0.3841792
#> 62      7        0.7        V3 0.3557529
#> 63      7        0.7        V4 0.4376714
#> 64      7        0.7        V5 0.3754987
#> 65      7        0.7        V9 0.3520057
#> 66      7        0.7        V6 0.3298526
#> 67      7        0.7        V7 0.3847258
#> 68      7        0.7        V8 0.3690811
#> 69      7        0.7        V2 0.3660361
#> 70      7        0.7       V10 0.3715319
#> 71      8        0.8        V9 0.3514772
#> 72      8        0.8       V10 0.3710204
#> 73      8        0.8        V1 0.3836304
#> 74      8        0.8        V5 0.3749449
#> 75      8        0.8        V2 0.3655289
#> 76      8        0.8        V3 0.3551769
#> 77      8        0.8        V4 0.4370755
#> 78      8        0.8        V8 0.3685215
#> 79      8        0.8        V6 0.3293854
#> 80      8        0.8        V7 0.3841755
#> 81      9        0.9        V5 0.3743919
#> 82      9        0.9        V7 0.3836259
#> 83      9        0.9        V8 0.3679627
#> 84      9        0.9        V9 0.3509495
#> 85      9        0.9        V6 0.3289188
#> 86      9        0.9       V10 0.3705095
#> 87      9        0.9        V1 0.3830824
#> 88      9        0.9        V2 0.3650225
#> 89      9        0.9        V3 0.3546019
#> 90      9        0.9        V4 0.4364803
#> 91     10        1.0        V1 0.3825352
#> 92     10        1.0        V4 0.4358859
#> 93     10        1.0        V5 0.3738397
#> 94     10        1.0        V9 0.3504225
#> 95     10        1.0        V3 0.3540278
#> 96     10        1.0        V7 0.3830771
#> 97     10        1.0        V8 0.3674047
#> 98     10        1.0        V2 0.3645168
#> 99     10        1.0        V6 0.3284529
#> 100    10        1.0       V10 0.3699993

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3874884
#> 2       1        0.1        V5 0.3788389
#> 3       1        0.1        V9 0.3551936
#> 4       1        0.1        V3 0.3592282
#> 5       1        0.1        V4 0.4412645
#> 6       1        0.1        V8 0.3724569
#> 7       1        0.1        V2 0.3690937
#> 8       1        0.1        V6 0.3326698
#> 9       1        0.1        V7 0.3880446
#> 10      1        0.1       V10 0.3746163
#> 11      2        0.2        V7 0.3874895
#> 12      2        0.2        V8 0.3718921
#> 13      2        0.2        V9 0.3546603
#> 14      2        0.2       V10 0.3741004
#> 15      2        0.2        V1 0.3869349
#> 16      2        0.2        V5 0.3782802
#> 17      2        0.2        V2 0.3685823
#> 18      2        0.2        V3 0.3586466
#> 19      2        0.2        V4 0.4406636
#> 20      2        0.2        V6 0.3321986
#> 21      3        0.3        V4 0.4400636
#> 22      3        0.3        V5 0.3777222
#> 23      3        0.3        V3 0.3580660
#> 24      3        0.3        V7 0.3869352
#> 25      3        0.3        V8 0.3713282
#> 26      3        0.3        V9 0.3541278
#> 27      3        0.3        V6 0.3317280
#> 28      3        0.3       V10 0.3735853
#> 29      3        0.3        V1 0.3863822
#> 30      3        0.3        V2 0.3680717
#> 31      4        0.4        V1 0.3858302
#> 32      4        0.4        V9 0.3535961
#> 33      4        0.4        V3 0.3574863
#> 34      4        0.4        V4 0.4394643
#> 35      4        0.4        V5 0.3771651
#> 36      4        0.4        V2 0.3675617
#> 37      4        0.4        V6 0.3312582
#> 38      4        0.4        V7 0.3863816
#> 39      4        0.4        V8 0.3707652
#> 40      4        0.4       V10 0.3730709
#> 41      5        0.5        V8 0.3702030
#> 42      5        0.5        V9 0.3530651
#> 43      5        0.5       V10 0.3725572
#> 44      5        0.5        V1 0.3852791
#> 45      5        0.5        V5 0.3766089
#> 46      5        0.5        V2 0.3670524
#> 47      5        0.5        V3 0.3569076
#> 48      5        0.5        V4 0.4388659
#> 49      5        0.5        V6 0.3307890
#> 50      5        0.5        V7 0.3858289
#> 51      6        0.6        V4 0.4382683
#> 52      6        0.6        V5 0.3760534
#> 53      6        0.6        V7 0.3852770
#> 54      6        0.6        V8 0.3696416
#> 55      6        0.6        V9 0.3525350
#> 56      6        0.6        V6 0.3303204
#> 57      6        0.6       V10 0.3720442
#> 58      6        0.6        V1 0.3847287
#> 59      6        0.6        V2 0.3665439
#> 60      6        0.6        V3 0.3563297
#> 61      7        0.7        V1 0.3841792
#> 62      7        0.7        V3 0.3557529
#> 63      7        0.7        V4 0.4376714
#> 64      7        0.7        V5 0.3754987
#> 65      7        0.7        V9 0.3520057
#> 66      7        0.7        V6 0.3298526
#> 67      7        0.7        V7 0.3847258
#> 68      7        0.7        V8 0.3690811
#> 69      7        0.7        V2 0.3660361
#> 70      7        0.7       V10 0.3715319
#> 71      8        0.8        V9 0.3514772
#> 72      8        0.8       V10 0.3710204
#> 73      8        0.8        V1 0.3836304
#> 74      8        0.8        V5 0.3749449
#> 75      8        0.8        V2 0.3655289
#> 76      8        0.8        V3 0.3551769
#> 77      8        0.8        V4 0.4370755
#> 78      8        0.8        V8 0.3685215
#> 79      8        0.8        V6 0.3293854
#> 80      8        0.8        V7 0.3841755
#> 81      9        0.9        V5 0.3743919
#> 82      9        0.9        V7 0.3836259
#> 83      9        0.9        V8 0.3679627
#> 84      9        0.9        V9 0.3509495
#> 85      9        0.9        V6 0.3289188
#> 86      9        0.9       V10 0.3705095
#> 87      9        0.9        V1 0.3830824
#> 88      9        0.9        V2 0.3650225
#> 89      9        0.9        V3 0.3546019
#> 90      9        0.9        V4 0.4364803
#> 91     10        1.0        V1 0.3825352
#> 92     10        1.0        V4 0.4358859
#> 93     10        1.0        V5 0.3738397
#> 94     10        1.0        V9 0.3504225
#> 95     10        1.0        V3 0.3540278
#> 96     10        1.0        V7 0.3830771
#> 97     10        1.0        V8 0.3674047
#> 98     10        1.0        V2 0.3645168
#> 99     10        1.0        V6 0.3284529
#> 100    10        1.0       V10 0.3699993


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
#> 1       0                0          0 0.8757906 0.04464314 0.8050178 0.9561196
#> 2       0               20         20 0.8757906 0.04464314 0.8050178 0.9561196
#> 3       0               40         40 0.8757906 0.04464314 0.8050178 0.9561196
#> 4       0               60         60 0.8757906 0.04464314 0.8050178 0.9561196
#> 5      20                0         20 0.8617131 0.04615282 0.7873118 0.9435103
#> 6      20               20         40 0.8617131 0.04615282 0.7873118 0.9435103
#> 7      20               40         60 0.8617131 0.04615282 0.7873118 0.9435103
#> 8      20               60         80 0.8617131 0.04615282 0.7873118 0.9435103
#> 9      40                0         40 0.8478591 0.04754140 0.7702018 0.9305427
#> 10     40               20         60 0.8478591 0.04754140 0.7702018 0.9305427
#> 11     40               40         80 0.8478591 0.04754140 0.7702018 0.9305427
#> 12     40               60        100 0.8478591 0.04754140 0.7702018 0.9305427
#> 13     60                0         60 0.8342249 0.04883100 0.7536294 0.9173439
#> 14     60               20         80 0.8342249 0.04883100 0.7536294 0.9173439
#> 15     60               40        100 0.8342249 0.04883100 0.7536294 0.9173439
#> 16     60               60        120 0.8342249 0.04883100 0.7536294 0.9173439
#> 17     80                0         80 0.8208071 0.05003750 0.7375488 0.9040042
#> 18     80               20        100 0.8208071 0.05003750 0.7375488 0.9040042
#> 19     80               40        120 0.8208071 0.05003750 0.7375488 0.9040042
#> 20     80               60        140 0.8208071 0.05003750 0.7375488 0.9040042
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12101443 0.182019638 0.5692486
#> 2  0.30574618 0.11830656 0.139262735 0.5080497
#> 3  0.26001915 0.11509190 0.105521721 0.4538430
#> 4  0.22113100 0.11110512 0.079005296 0.4059731
#> 5  0.25589195 0.11217625 0.106448733 0.4771664
#> 6  0.21762106 0.10721989 0.079731741 0.4265578
#> 7  0.18507391 0.10221937 0.058864788 0.3819084
#> 8  0.15739448 0.09700398 0.042709325 0.3425484
#> 9  0.18213629 0.10193273 0.059434233 0.4011050
#> 10 0.15489621 0.09601196 0.043147844 0.3594695
#> 11 0.13173012 0.09028146 0.030678729 0.3227674
#> 12 0.11202872 0.08462793 0.021272716 0.2903953
#> 13 0.12963921 0.09132159 0.031014921 0.3385469
#> 14 0.11025053 0.08517407 0.021524077 0.3043170
#> 15 0.09376159 0.07932756 0.014494755 0.2741065
#> 16 0.07973872 0.07370279 0.009412254 0.2473997
#> 17 0.09227334 0.08099557 0.014680548 0.2871017
#> 18 0.07847305 0.07502910 0.009544660 0.2588942
#> 19 0.06673672 0.06940568 0.005939389 0.2339288
#> 20 0.05675565 0.06406719 0.003504524 0.2117763
```
