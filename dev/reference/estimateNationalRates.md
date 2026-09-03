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
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
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
#> 1         0.1 0.3838513 0.03007032 0.3366512 0.4255391
#> 2         0.2 0.3832760 0.03001784 0.3361571 0.4248774
#> 3         0.3 0.3827015 0.02996548 0.3356638 0.4242167
#> 4         0.4 0.3821279 0.02991323 0.3351713 0.4235570
#> 5         0.5 0.3815551 0.02986111 0.3346794 0.4228983
#> 6         0.6 0.3809832 0.02980910 0.3341883 0.4222407
#> 7         0.7 0.3804122 0.02975721 0.3336979 0.4215841
#> 8         0.8 0.3798420 0.02970544 0.3332082 0.4209285
#> 9         0.9 0.3792726 0.02965379 0.3327192 0.4202739
#> 10        1.0 0.3787041 0.02960225 0.3322310 0.4196204

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3835802
#> 2       1        0.1        V5 0.3917792
#> 3       1        0.1        V9 0.3352095
#> 4       1        0.1        V3 0.3954149
#> 5       1        0.1        V4 0.3847256
#> 6       1        0.1        V8 0.3776426
#> 7       1        0.1        V2 0.3501786
#> 8       1        0.1        V6 0.4305793
#> 9       1        0.1        V7 0.4081786
#> 10      1        0.1       V10 0.3416170
#> 11      2        0.2        V7 0.4074955
#> 12      2        0.2        V8 0.3770827
#> 13      2        0.2        V9 0.3347264
#> 14      2        0.2       V10 0.3410852
#> 15      2        0.2        V1 0.3830506
#> 16      2        0.2        V5 0.3911962
#> 17      2        0.2        V2 0.3496359
#> 18      2        0.2        V3 0.3948165
#> 19      2        0.2        V4 0.3841466
#> 20      2        0.2        V6 0.4299238
#> 21      3        0.3        V4 0.3835686
#> 22      3        0.3        V5 0.3906140
#> 23      3        0.3        V3 0.3942190
#> 24      3        0.3        V7 0.4068134
#> 25      3        0.3        V8 0.3765236
#> 26      3        0.3        V9 0.3342441
#> 27      3        0.3        V6 0.4292692
#> 28      3        0.3       V10 0.3405542
#> 29      3        0.3        V1 0.3825218
#> 30      3        0.3        V2 0.3490940
#> 31      4        0.4        V1 0.3819936
#> 32      4        0.4        V9 0.3337624
#> 33      4        0.4        V3 0.3936223
#> 34      4        0.4        V4 0.3829914
#> 35      4        0.4        V5 0.3900326
#> 36      4        0.4        V2 0.3485530
#> 37      4        0.4        V6 0.4286157
#> 38      4        0.4        V7 0.4061325
#> 39      4        0.4        V8 0.3759653
#> 40      4        0.4       V10 0.3400241
#> 41      5        0.5        V8 0.3754079
#> 42      5        0.5        V9 0.3332814
#> 43      5        0.5       V10 0.3394947
#> 44      5        0.5        V1 0.3814663
#> 45      5        0.5        V5 0.3894522
#> 46      5        0.5        V2 0.3480128
#> 47      5        0.5        V3 0.3930266
#> 48      5        0.5        V4 0.3824150
#> 49      5        0.5        V6 0.4279632
#> 50      5        0.5        V7 0.4054528
#> 51      6        0.6        V4 0.3818396
#> 52      6        0.6        V5 0.3888726
#> 53      6        0.6        V7 0.4047742
#> 54      6        0.6        V8 0.3748513
#> 55      6        0.6        V9 0.3328011
#> 56      6        0.6        V6 0.4273116
#> 57      6        0.6       V10 0.3389663
#> 58      6        0.6        V1 0.3809396
#> 59      6        0.6        V2 0.3474734
#> 60      6        0.6        V3 0.3924318
#> 61      7        0.7        V1 0.3804137
#> 62      7        0.7        V3 0.3918379
#> 63      7        0.7        V4 0.3812650
#> 64      7        0.7        V5 0.3882938
#> 65      7        0.7        V9 0.3323215
#> 66      7        0.7        V6 0.4266611
#> 67      7        0.7        V7 0.4040967
#> 68      7        0.7        V8 0.3742955
#> 69      7        0.7        V2 0.3469348
#> 70      7        0.7       V10 0.3384386
#> 71      8        0.8        V9 0.3318426
#> 72      8        0.8       V10 0.3379117
#> 73      8        0.8        V1 0.3798885
#> 74      8        0.8        V5 0.3877159
#> 75      8        0.8        V2 0.3463971
#> 76      8        0.8        V3 0.3912449
#> 77      8        0.8        V4 0.3806912
#> 78      8        0.8        V8 0.3737406
#> 79      8        0.8        V6 0.4260115
#> 80      8        0.8        V7 0.4034204
#> 81      9        0.9        V5 0.3871389
#> 82      9        0.9        V7 0.4027452
#> 83      9        0.9        V8 0.3731864
#> 84      9        0.9        V9 0.3313644
#> 85      9        0.9        V6 0.4253629
#> 86      9        0.9       V10 0.3373857
#> 87      9        0.9        V1 0.3793640
#> 88      9        0.9        V2 0.3458603
#> 89      9        0.9        V3 0.3906528
#> 90      9        0.9        V4 0.3801184
#> 91     10        1.0        V1 0.3788402
#> 92     10        1.0        V4 0.3795464
#> 93     10        1.0        V5 0.3865627
#> 94     10        1.0        V9 0.3308869
#> 95     10        1.0        V3 0.3900616
#> 96     10        1.0        V7 0.4020711
#> 97     10        1.0        V8 0.3726331
#> 98     10        1.0        V2 0.3453242
#> 99     10        1.0        V6 0.4247153
#> 100    10        1.0       V10 0.3368605

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3835802
#> 2       1        0.1        V5 0.3917792
#> 3       1        0.1        V9 0.3352095
#> 4       1        0.1        V3 0.3954149
#> 5       1        0.1        V4 0.3847256
#> 6       1        0.1        V8 0.3776426
#> 7       1        0.1        V2 0.3501786
#> 8       1        0.1        V6 0.4305793
#> 9       1        0.1        V7 0.4081786
#> 10      1        0.1       V10 0.3416170
#> 11      2        0.2        V7 0.4074955
#> 12      2        0.2        V8 0.3770827
#> 13      2        0.2        V9 0.3347264
#> 14      2        0.2       V10 0.3410852
#> 15      2        0.2        V1 0.3830506
#> 16      2        0.2        V5 0.3911962
#> 17      2        0.2        V2 0.3496359
#> 18      2        0.2        V3 0.3948165
#> 19      2        0.2        V4 0.3841466
#> 20      2        0.2        V6 0.4299238
#> 21      3        0.3        V4 0.3835686
#> 22      3        0.3        V5 0.3906140
#> 23      3        0.3        V3 0.3942190
#> 24      3        0.3        V7 0.4068134
#> 25      3        0.3        V8 0.3765236
#> 26      3        0.3        V9 0.3342441
#> 27      3        0.3        V6 0.4292692
#> 28      3        0.3       V10 0.3405542
#> 29      3        0.3        V1 0.3825218
#> 30      3        0.3        V2 0.3490940
#> 31      4        0.4        V1 0.3819936
#> 32      4        0.4        V9 0.3337624
#> 33      4        0.4        V3 0.3936223
#> 34      4        0.4        V4 0.3829914
#> 35      4        0.4        V5 0.3900326
#> 36      4        0.4        V2 0.3485530
#> 37      4        0.4        V6 0.4286157
#> 38      4        0.4        V7 0.4061325
#> 39      4        0.4        V8 0.3759653
#> 40      4        0.4       V10 0.3400241
#> 41      5        0.5        V8 0.3754079
#> 42      5        0.5        V9 0.3332814
#> 43      5        0.5       V10 0.3394947
#> 44      5        0.5        V1 0.3814663
#> 45      5        0.5        V5 0.3894522
#> 46      5        0.5        V2 0.3480128
#> 47      5        0.5        V3 0.3930266
#> 48      5        0.5        V4 0.3824150
#> 49      5        0.5        V6 0.4279632
#> 50      5        0.5        V7 0.4054528
#> 51      6        0.6        V4 0.3818396
#> 52      6        0.6        V5 0.3888726
#> 53      6        0.6        V7 0.4047742
#> 54      6        0.6        V8 0.3748513
#> 55      6        0.6        V9 0.3328011
#> 56      6        0.6        V6 0.4273116
#> 57      6        0.6       V10 0.3389663
#> 58      6        0.6        V1 0.3809396
#> 59      6        0.6        V2 0.3474734
#> 60      6        0.6        V3 0.3924318
#> 61      7        0.7        V1 0.3804137
#> 62      7        0.7        V3 0.3918379
#> 63      7        0.7        V4 0.3812650
#> 64      7        0.7        V5 0.3882938
#> 65      7        0.7        V9 0.3323215
#> 66      7        0.7        V6 0.4266611
#> 67      7        0.7        V7 0.4040967
#> 68      7        0.7        V8 0.3742955
#> 69      7        0.7        V2 0.3469348
#> 70      7        0.7       V10 0.3384386
#> 71      8        0.8        V9 0.3318426
#> 72      8        0.8       V10 0.3379117
#> 73      8        0.8        V1 0.3798885
#> 74      8        0.8        V5 0.3877159
#> 75      8        0.8        V2 0.3463971
#> 76      8        0.8        V3 0.3912449
#> 77      8        0.8        V4 0.3806912
#> 78      8        0.8        V8 0.3737406
#> 79      8        0.8        V6 0.4260115
#> 80      8        0.8        V7 0.4034204
#> 81      9        0.9        V5 0.3871389
#> 82      9        0.9        V7 0.4027452
#> 83      9        0.9        V8 0.3731864
#> 84      9        0.9        V9 0.3313644
#> 85      9        0.9        V6 0.4253629
#> 86      9        0.9       V10 0.3373857
#> 87      9        0.9        V1 0.3793640
#> 88      9        0.9        V2 0.3458603
#> 89      9        0.9        V3 0.3906528
#> 90      9        0.9        V4 0.3801184
#> 91     10        1.0        V1 0.3788402
#> 92     10        1.0        V4 0.3795464
#> 93     10        1.0        V5 0.3865627
#> 94     10        1.0        V9 0.3308869
#> 95     10        1.0        V3 0.3900616
#> 96     10        1.0        V7 0.4020711
#> 97     10        1.0        V8 0.3726331
#> 98     10        1.0        V2 0.3453242
#> 99     10        1.0        V6 0.4247153
#> 100    10        1.0       V10 0.3368605


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
#> 1       0                0          0 0.8757906 0.04960049 0.7932878 0.9508908
#> 2       0               20         20 0.8757906 0.04960049 0.7932878 0.9508908
#> 3       0               40         40 0.8757906 0.04960049 0.7932878 0.9508908
#> 4       0               60         60 0.8757906 0.04960049 0.7932878 0.9508908
#> 5      20                0         20 0.8617131 0.05171050 0.7795096 0.9407489
#> 6      20               20         40 0.8617131 0.05171050 0.7795096 0.9407489
#> 7      20               40         60 0.8617131 0.05171050 0.7795096 0.9407489
#> 8      20               60         80 0.8617131 0.05171050 0.7795096 0.9407489
#> 9      40                0         40 0.8478591 0.05369424 0.7660997 0.9304097
#> 10     40               20         60 0.8478591 0.05369424 0.7660997 0.9304097
#> 11     40               40         80 0.8478591 0.05369424 0.7660997 0.9304097
#> 12     40               60        100 0.8478591 0.05369424 0.7660997 0.9304097
#> 13     60                0         60 0.8342249 0.05556146 0.7530291 0.9199328
#> 14     60               20         80 0.8342249 0.05556146 0.7530291 0.9199328
#> 15     60               40        100 0.8342249 0.05556146 0.7530291 0.9199328
#> 16     60               60        120 0.8342249 0.05556146 0.7530291 0.9199328
#> 17     80                0         80 0.8208071 0.05732067 0.7402741 0.9093638
#> 18     80               20        100 0.8208071 0.05732067 0.7402741 0.9093638
#> 19     80               40        120 0.8208071 0.05732067 0.7402741 0.9093638
#> 20     80               60        140 0.8208071 0.05732067 0.7402741 0.9093638
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12038873 0.1565326562 0.5792469
#> 2  0.30574618 0.11174771 0.1227777197 0.5179666
#> 3  0.26001915 0.10417157 0.0954941718 0.4635262
#> 4  0.22113100 0.09732887 0.0735290846 0.4153179
#> 5  0.25589195 0.10503915 0.0781093527 0.4669586
#> 6  0.21762106 0.09729505 0.0596012236 0.4183543
#> 7  0.18507391 0.09039503 0.0448648361 0.3753825
#> 8  0.15739448 0.08412894 0.0332346350 0.3374149
#> 9  0.18213629 0.09073929 0.0356372799 0.3780880
#> 10 0.15489621 0.08395241 0.0260228844 0.3398053
#> 11 0.13173012 0.07785043 0.0185967978 0.3059752
#> 12 0.11202872 0.07229384 0.0129544523 0.2760544
#> 13 0.12963921 0.07790277 0.0140996952 0.3081057
#> 14 0.11025053 0.07206156 0.0095986145 0.2779398
#> 15 0.09376159 0.06677963 0.0063114897 0.2512240
#> 16 0.07973872 0.06195849 0.0039820304 0.2275166
#> 17 0.09227334 0.06672109 0.0044396173 0.2529088
#> 18 0.07847305 0.06175738 0.0026985698 0.2290133
#> 19 0.06673672 0.05724521 0.0015504628 0.2077590
#> 20 0.05675565 0.05310985 0.0008331458 0.1887970
```
