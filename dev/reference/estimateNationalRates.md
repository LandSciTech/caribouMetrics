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
#> 1         0.1 0.3838513 0.02633998 0.3420246 0.4189988
#> 2         0.2 0.3832760 0.02629701 0.3415114 0.4183493
#> 3         0.3 0.3827015 0.02625412 0.3409990 0.4177008
#> 4         0.4 0.3821279 0.02621132 0.3404874 0.4170533
#> 5         0.5 0.3815551 0.02616861 0.3399765 0.4164068
#> 6         0.6 0.3809832 0.02612599 0.3394664 0.4157613
#> 7         0.7 0.3804122 0.02608345 0.3389571 0.4151168
#> 8         0.8 0.3798420 0.02604100 0.3384485 0.4144733
#> 9         0.9 0.3792726 0.02599864 0.3379407 0.4138308
#> 10        1.0 0.3787041 0.02595636 0.3374336 0.4131893

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3773750
#> 2       1        0.1        V5 0.3540246
#> 3       1        0.1        V9 0.3900188
#> 4       1        0.1        V3 0.3780412
#> 5       1        0.1        V4 0.3839406
#> 6       1        0.1        V8 0.4156926
#> 7       1        0.1        V2 0.4199587
#> 8       1        0.1        V6 0.3393590
#> 9       1        0.1        V7 0.3895618
#> 10      1        0.1       V10 0.3512062
#> 11      2        0.2        V7 0.3889902
#> 12      2        0.2        V8 0.4150650
#> 13      2        0.2        V9 0.3894185
#> 14      2        0.2       V10 0.3506667
#> 15      2        0.2        V1 0.3767883
#> 16      2        0.2        V5 0.3534807
#> 17      2        0.2        V2 0.4193029
#> 18      2        0.2        V3 0.3775168
#> 19      2        0.2        V4 0.3833317
#> 20      2        0.2        V6 0.3388535
#> 21      3        0.3        V4 0.3827237
#> 22      3        0.3        V5 0.3529378
#> 23      3        0.3        V3 0.3769931
#> 24      3        0.3        V7 0.3884194
#> 25      3        0.3        V8 0.4144383
#> 26      3        0.3        V9 0.3888191
#> 27      3        0.3        V6 0.3383487
#> 28      3        0.3       V10 0.3501280
#> 29      3        0.3        V1 0.3762026
#> 30      3        0.3        V2 0.4186480
#> 31      4        0.4        V1 0.3756178
#> 32      4        0.4        V9 0.3882207
#> 33      4        0.4        V3 0.3764702
#> 34      4        0.4        V4 0.3821167
#> 35      4        0.4        V5 0.3523956
#> 36      4        0.4        V2 0.4179942
#> 37      4        0.4        V6 0.3378446
#> 38      4        0.4        V7 0.3878494
#> 39      4        0.4        V8 0.4138126
#> 40      4        0.4       V10 0.3495902
#> 41      5        0.5        V8 0.4131878
#> 42      5        0.5        V9 0.3876231
#> 43      5        0.5       V10 0.3490532
#> 44      5        0.5        V1 0.3750339
#> 45      5        0.5        V5 0.3518543
#> 46      5        0.5        V2 0.4173414
#> 47      5        0.5        V3 0.3759479
#> 48      5        0.5        V4 0.3815106
#> 49      5        0.5        V6 0.3373413
#> 50      5        0.5        V7 0.3872803
#> 51      6        0.6        V4 0.3809055
#> 52      6        0.6        V5 0.3513138
#> 53      6        0.6        V7 0.3867120
#> 54      6        0.6        V8 0.4125640
#> 55      6        0.6        V9 0.3870265
#> 56      6        0.6        V6 0.3368388
#> 57      6        0.6       V10 0.3485170
#> 58      6        0.6        V1 0.3744509
#> 59      6        0.6        V2 0.4166896
#> 60      6        0.6        V3 0.3754264
#> 61      7        0.7        V1 0.3738688
#> 62      7        0.7        V3 0.3749057
#> 63      7        0.7        V4 0.3803013
#> 64      7        0.7        V5 0.3507742
#> 65      7        0.7        V9 0.3864309
#> 66      7        0.7        V6 0.3363370
#> 67      7        0.7        V7 0.3861446
#> 68      7        0.7        V8 0.4119411
#> 69      7        0.7        V2 0.4160388
#> 70      7        0.7       V10 0.3479816
#> 71      8        0.8        V9 0.3858361
#> 72      8        0.8       V10 0.3474471
#> 73      8        0.8        V1 0.3732876
#> 74      8        0.8        V5 0.3502354
#> 75      8        0.8        V2 0.4153890
#> 76      8        0.8        V3 0.3743856
#> 77      8        0.8        V4 0.3796982
#> 78      8        0.8        V8 0.4113192
#> 79      8        0.8        V6 0.3358360
#> 80      8        0.8        V7 0.3855779
#> 81      9        0.9        V5 0.3496974
#> 82      9        0.9        V7 0.3850122
#> 83      9        0.9        V8 0.4106982
#> 84      9        0.9        V9 0.3852422
#> 85      9        0.9        V6 0.3353357
#> 86      9        0.9       V10 0.3469134
#> 87      9        0.9        V1 0.3727073
#> 88      9        0.9        V2 0.4147403
#> 89      9        0.9        V3 0.3738663
#> 90      9        0.9        V4 0.3790959
#> 91     10        1.0        V1 0.3721279
#> 92     10        1.0        V4 0.3784947
#> 93     10        1.0        V5 0.3491602
#> 94     10        1.0        V9 0.3846493
#> 95     10        1.0        V3 0.3733477
#> 96     10        1.0        V7 0.3844472
#> 97     10        1.0        V8 0.4100781
#> 98     10        1.0        V2 0.4140926
#> 99     10        1.0        V6 0.3348361
#> 100    10        1.0       V10 0.3463805

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3773750
#> 2       1        0.1        V5 0.3540246
#> 3       1        0.1        V9 0.3900188
#> 4       1        0.1        V3 0.3780412
#> 5       1        0.1        V4 0.3839406
#> 6       1        0.1        V8 0.4156926
#> 7       1        0.1        V2 0.4199587
#> 8       1        0.1        V6 0.3393590
#> 9       1        0.1        V7 0.3895618
#> 10      1        0.1       V10 0.3512062
#> 11      2        0.2        V7 0.3889902
#> 12      2        0.2        V8 0.4150650
#> 13      2        0.2        V9 0.3894185
#> 14      2        0.2       V10 0.3506667
#> 15      2        0.2        V1 0.3767883
#> 16      2        0.2        V5 0.3534807
#> 17      2        0.2        V2 0.4193029
#> 18      2        0.2        V3 0.3775168
#> 19      2        0.2        V4 0.3833317
#> 20      2        0.2        V6 0.3388535
#> 21      3        0.3        V4 0.3827237
#> 22      3        0.3        V5 0.3529378
#> 23      3        0.3        V3 0.3769931
#> 24      3        0.3        V7 0.3884194
#> 25      3        0.3        V8 0.4144383
#> 26      3        0.3        V9 0.3888191
#> 27      3        0.3        V6 0.3383487
#> 28      3        0.3       V10 0.3501280
#> 29      3        0.3        V1 0.3762026
#> 30      3        0.3        V2 0.4186480
#> 31      4        0.4        V1 0.3756178
#> 32      4        0.4        V9 0.3882207
#> 33      4        0.4        V3 0.3764702
#> 34      4        0.4        V4 0.3821167
#> 35      4        0.4        V5 0.3523956
#> 36      4        0.4        V2 0.4179942
#> 37      4        0.4        V6 0.3378446
#> 38      4        0.4        V7 0.3878494
#> 39      4        0.4        V8 0.4138126
#> 40      4        0.4       V10 0.3495902
#> 41      5        0.5        V8 0.4131878
#> 42      5        0.5        V9 0.3876231
#> 43      5        0.5       V10 0.3490532
#> 44      5        0.5        V1 0.3750339
#> 45      5        0.5        V5 0.3518543
#> 46      5        0.5        V2 0.4173414
#> 47      5        0.5        V3 0.3759479
#> 48      5        0.5        V4 0.3815106
#> 49      5        0.5        V6 0.3373413
#> 50      5        0.5        V7 0.3872803
#> 51      6        0.6        V4 0.3809055
#> 52      6        0.6        V5 0.3513138
#> 53      6        0.6        V7 0.3867120
#> 54      6        0.6        V8 0.4125640
#> 55      6        0.6        V9 0.3870265
#> 56      6        0.6        V6 0.3368388
#> 57      6        0.6       V10 0.3485170
#> 58      6        0.6        V1 0.3744509
#> 59      6        0.6        V2 0.4166896
#> 60      6        0.6        V3 0.3754264
#> 61      7        0.7        V1 0.3738688
#> 62      7        0.7        V3 0.3749057
#> 63      7        0.7        V4 0.3803013
#> 64      7        0.7        V5 0.3507742
#> 65      7        0.7        V9 0.3864309
#> 66      7        0.7        V6 0.3363370
#> 67      7        0.7        V7 0.3861446
#> 68      7        0.7        V8 0.4119411
#> 69      7        0.7        V2 0.4160388
#> 70      7        0.7       V10 0.3479816
#> 71      8        0.8        V9 0.3858361
#> 72      8        0.8       V10 0.3474471
#> 73      8        0.8        V1 0.3732876
#> 74      8        0.8        V5 0.3502354
#> 75      8        0.8        V2 0.4153890
#> 76      8        0.8        V3 0.3743856
#> 77      8        0.8        V4 0.3796982
#> 78      8        0.8        V8 0.4113192
#> 79      8        0.8        V6 0.3358360
#> 80      8        0.8        V7 0.3855779
#> 81      9        0.9        V5 0.3496974
#> 82      9        0.9        V7 0.3850122
#> 83      9        0.9        V8 0.4106982
#> 84      9        0.9        V9 0.3852422
#> 85      9        0.9        V6 0.3353357
#> 86      9        0.9       V10 0.3469134
#> 87      9        0.9        V1 0.3727073
#> 88      9        0.9        V2 0.4147403
#> 89      9        0.9        V3 0.3738663
#> 90      9        0.9        V4 0.3790959
#> 91     10        1.0        V1 0.3721279
#> 92     10        1.0        V4 0.3784947
#> 93     10        1.0        V5 0.3491602
#> 94     10        1.0        V9 0.3846493
#> 95     10        1.0        V3 0.3733477
#> 96     10        1.0        V7 0.3844472
#> 97     10        1.0        V8 0.4100781
#> 98     10        1.0        V2 0.4140926
#> 99     10        1.0        V6 0.3348361
#> 100    10        1.0       V10 0.3463805


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
#> 1       0                0          0 0.8757906 0.04319503 0.7946990 0.9409968
#> 2       0               20         20 0.8757906 0.04319503 0.7946990 0.9409968
#> 3       0               40         40 0.8757906 0.04319503 0.7946990 0.9409968
#> 4       0               60         60 0.8757906 0.04319503 0.7946990 0.9409968
#> 5      20                0         20 0.8617131 0.04635155 0.7716287 0.9298100
#> 6      20               20         40 0.8617131 0.04635155 0.7716287 0.9298100
#> 7      20               40         60 0.8617131 0.04635155 0.7716287 0.9298100
#> 8      20               60         80 0.8617131 0.04635155 0.7716287 0.9298100
#> 9      40                0         40 0.8478591 0.04930878 0.7495411 0.9184784
#> 10     40               20         60 0.8478591 0.04930878 0.7495411 0.9184784
#> 11     40               40         80 0.8478591 0.04930878 0.7495411 0.9184784
#> 12     40               60        100 0.8478591 0.04930878 0.7495411 0.9184784
#> 13     60                0         60 0.8342249 0.05208319 0.7283240 0.9070562
#> 14     60               20         80 0.8342249 0.05208319 0.7283240 0.9070562
#> 15     60               40        100 0.8342249 0.05208319 0.7283240 0.9070562
#> 16     60               60        120 0.8342249 0.05208319 0.7283240 0.9070562
#> 17     80                0         80 0.8208071 0.05468946 0.7078929 0.8955848
#> 18     80               20        100 0.8208071 0.05468946 0.7078929 0.8955848
#> 19     80               40        120 0.8208071 0.05468946 0.7078929 0.8955848
#> 20     80               60        140 0.8208071 0.05468946 0.7078929 0.8955848
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11855613 0.168173462 0.5775846
#> 2  0.30574618 0.11594768 0.124850981 0.5236567
#> 3  0.26001915 0.11246886 0.091495991 0.4749828
#> 4  0.22113100 0.10821351 0.065976219 0.4311654
#> 5  0.25589195 0.11327999 0.090387547 0.4737063
#> 6  0.21762106 0.10799246 0.065132048 0.4300174
#> 7  0.18507391 0.10256379 0.046000618 0.3907512
#> 8  0.15739448 0.09700780 0.031699700 0.3554899
#> 9  0.18213629 0.10377909 0.045372039 0.3897230
#> 10 0.15489621 0.09750290 0.031233954 0.3545668
#> 11 0.13173012 0.09144305 0.020858952 0.3230062
#> 12 0.11202872 0.08557742 0.013421670 0.2946695
#> 13 0.12963921 0.09272764 0.020525140 0.3221800
#> 14 0.11025053 0.08635503 0.013186033 0.2939276
#> 15 0.09376159 0.08036324 0.008089867 0.2685482
#> 16 0.07973872 0.07471322 0.004689836 0.2457307
#> 17 0.09227334 0.08157233 0.007931758 0.2678834
#> 18 0.07847305 0.07556275 0.004587108 0.2451327
#> 19 0.06673672 0.06999238 0.002472731 0.2246549
#> 20 0.05675565 0.06481359 0.001222636 0.2061970
```
