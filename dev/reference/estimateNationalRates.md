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
#> 1         0.1 0.3838513 0.03031820 0.3481473 0.4335226
#> 2         0.2 0.3832760 0.03029509 0.3476357 0.4328229
#> 3         0.3 0.3827015 0.03027212 0.3471248 0.4321242
#> 4         0.4 0.3821279 0.03024931 0.3466147 0.4314267
#> 5         0.5 0.3815551 0.03022666 0.3461054 0.4307303
#> 6         0.6 0.3809832 0.03020415 0.3455968 0.4300351
#> 7         0.7 0.3804122 0.03018179 0.3450889 0.4293410
#> 8         0.8 0.3798420 0.03015958 0.3445818 0.4286479
#> 9         0.9 0.3792726 0.03013752 0.3440755 0.4279561
#> 10        1.0 0.3787041 0.03011561 0.3435699 0.4272653

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4260353
#> 2       1        0.1        V5 0.3914379
#> 3       1        0.1        V9 0.4114308
#> 4       1        0.1        V3 0.3697281
#> 5       1        0.1        V4 0.4033079
#> 6       1        0.1        V8 0.3444331
#> 7       1        0.1        V2 0.4239104
#> 8       1        0.1        V6 0.4081531
#> 9       1        0.1        V7 0.3609406
#> 10      1        0.1       V10 0.4356964
#> 11      2        0.2        V7 0.3603335
#> 12      2        0.2        V8 0.3439492
#> 13      2        0.2        V9 0.4108929
#> 14      2        0.2       V10 0.4349717
#> 15      2        0.2        V1 0.4254213
#> 16      2        0.2        V5 0.3907915
#> 17      2        0.2        V2 0.4233915
#> 18      2        0.2        V3 0.3691036
#> 19      2        0.2        V4 0.4027135
#> 20      2        0.2        V6 0.4076225
#> 21      3        0.3        V4 0.4021199
#> 22      3        0.3        V5 0.3901463
#> 23      3        0.3        V3 0.3684803
#> 24      3        0.3        V7 0.3597275
#> 25      3        0.3        V8 0.3434660
#> 26      3        0.3        V9 0.4103558
#> 27      3        0.3        V6 0.4070925
#> 28      3        0.3       V10 0.4342482
#> 29      3        0.3        V1 0.4248083
#> 30      3        0.3        V2 0.4228732
#> 31      4        0.4        V1 0.4241961
#> 32      4        0.4        V9 0.4098194
#> 33      4        0.4        V3 0.3678580
#> 34      4        0.4        V4 0.4015273
#> 35      4        0.4        V5 0.3895021
#> 36      4        0.4        V2 0.4223555
#> 37      4        0.4        V6 0.4065633
#> 38      4        0.4        V7 0.3591225
#> 39      4        0.4        V8 0.3429834
#> 40      4        0.4       V10 0.4335259
#> 41      5        0.5        V8 0.3425015
#> 42      5        0.5        V9 0.4092837
#> 43      5        0.5       V10 0.4328049
#> 44      5        0.5        V1 0.4235848
#> 45      5        0.5        V5 0.3888589
#> 46      5        0.5        V2 0.4218384
#> 47      5        0.5        V3 0.3672367
#> 48      5        0.5        V4 0.4009355
#> 49      5        0.5        V6 0.4060347
#> 50      5        0.5        V7 0.3585186
#> 51      6        0.6        V4 0.4003445
#> 52      6        0.6        V5 0.3882169
#> 53      6        0.6        V7 0.3579156
#> 54      6        0.6        V8 0.3420203
#> 55      6        0.6        V9 0.4087486
#> 56      6        0.6        V6 0.4055068
#> 57      6        0.6       V10 0.4320850
#> 58      6        0.6        V1 0.4229743
#> 59      6        0.6        V2 0.4213220
#> 60      6        0.6        V3 0.3666165
#> 61      7        0.7        V1 0.4223648
#> 62      7        0.7        V3 0.3659973
#> 63      7        0.7        V4 0.3997545
#> 64      7        0.7        V5 0.3875759
#> 65      7        0.7        V9 0.4082143
#> 66      7        0.7        V6 0.4049796
#> 67      7        0.7        V7 0.3573137
#> 68      7        0.7        V8 0.3415398
#> 69      7        0.7        V2 0.4208062
#> 70      7        0.7       V10 0.4313663
#> 71      8        0.8        V9 0.4076807
#> 72      8        0.8       V10 0.4306488
#> 73      8        0.8        V1 0.4217561
#> 74      8        0.8        V5 0.3869359
#> 75      8        0.8        V2 0.4202911
#> 76      8        0.8        V3 0.3653792
#> 77      8        0.8        V4 0.3991653
#> 78      8        0.8        V8 0.3410600
#> 79      8        0.8        V6 0.4044531
#> 80      8        0.8        V7 0.3567127
#> 81      9        0.9        V5 0.3862970
#> 82      9        0.9        V7 0.3561128
#> 83      9        0.9        V8 0.3405808
#> 84      9        0.9        V9 0.4071478
#> 85      9        0.9        V6 0.4039272
#> 86      9        0.9       V10 0.4299325
#> 87      9        0.9        V1 0.4211483
#> 88      9        0.9        V2 0.4197766
#> 89      9        0.9        V3 0.3647621
#> 90      9        0.9        V4 0.3985770
#> 91     10        1.0        V1 0.4205414
#> 92     10        1.0        V4 0.3979895
#> 93     10        1.0        V5 0.3856592
#> 94     10        1.0        V9 0.4066155
#> 95     10        1.0        V3 0.3641461
#> 96     10        1.0        V7 0.3555139
#> 97     10        1.0        V8 0.3401023
#> 98     10        1.0        V2 0.4192627
#> 99     10        1.0        V6 0.4034021
#> 100    10        1.0       V10 0.4292174

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4260353
#> 2       1        0.1        V5 0.3914379
#> 3       1        0.1        V9 0.4114308
#> 4       1        0.1        V3 0.3697281
#> 5       1        0.1        V4 0.4033079
#> 6       1        0.1        V8 0.3444331
#> 7       1        0.1        V2 0.4239104
#> 8       1        0.1        V6 0.4081531
#> 9       1        0.1        V7 0.3609406
#> 10      1        0.1       V10 0.4356964
#> 11      2        0.2        V7 0.3603335
#> 12      2        0.2        V8 0.3439492
#> 13      2        0.2        V9 0.4108929
#> 14      2        0.2       V10 0.4349717
#> 15      2        0.2        V1 0.4254213
#> 16      2        0.2        V5 0.3907915
#> 17      2        0.2        V2 0.4233915
#> 18      2        0.2        V3 0.3691036
#> 19      2        0.2        V4 0.4027135
#> 20      2        0.2        V6 0.4076225
#> 21      3        0.3        V4 0.4021199
#> 22      3        0.3        V5 0.3901463
#> 23      3        0.3        V3 0.3684803
#> 24      3        0.3        V7 0.3597275
#> 25      3        0.3        V8 0.3434660
#> 26      3        0.3        V9 0.4103558
#> 27      3        0.3        V6 0.4070925
#> 28      3        0.3       V10 0.4342482
#> 29      3        0.3        V1 0.4248083
#> 30      3        0.3        V2 0.4228732
#> 31      4        0.4        V1 0.4241961
#> 32      4        0.4        V9 0.4098194
#> 33      4        0.4        V3 0.3678580
#> 34      4        0.4        V4 0.4015273
#> 35      4        0.4        V5 0.3895021
#> 36      4        0.4        V2 0.4223555
#> 37      4        0.4        V6 0.4065633
#> 38      4        0.4        V7 0.3591225
#> 39      4        0.4        V8 0.3429834
#> 40      4        0.4       V10 0.4335259
#> 41      5        0.5        V8 0.3425015
#> 42      5        0.5        V9 0.4092837
#> 43      5        0.5       V10 0.4328049
#> 44      5        0.5        V1 0.4235848
#> 45      5        0.5        V5 0.3888589
#> 46      5        0.5        V2 0.4218384
#> 47      5        0.5        V3 0.3672367
#> 48      5        0.5        V4 0.4009355
#> 49      5        0.5        V6 0.4060347
#> 50      5        0.5        V7 0.3585186
#> 51      6        0.6        V4 0.4003445
#> 52      6        0.6        V5 0.3882169
#> 53      6        0.6        V7 0.3579156
#> 54      6        0.6        V8 0.3420203
#> 55      6        0.6        V9 0.4087486
#> 56      6        0.6        V6 0.4055068
#> 57      6        0.6       V10 0.4320850
#> 58      6        0.6        V1 0.4229743
#> 59      6        0.6        V2 0.4213220
#> 60      6        0.6        V3 0.3666165
#> 61      7        0.7        V1 0.4223648
#> 62      7        0.7        V3 0.3659973
#> 63      7        0.7        V4 0.3997545
#> 64      7        0.7        V5 0.3875759
#> 65      7        0.7        V9 0.4082143
#> 66      7        0.7        V6 0.4049796
#> 67      7        0.7        V7 0.3573137
#> 68      7        0.7        V8 0.3415398
#> 69      7        0.7        V2 0.4208062
#> 70      7        0.7       V10 0.4313663
#> 71      8        0.8        V9 0.4076807
#> 72      8        0.8       V10 0.4306488
#> 73      8        0.8        V1 0.4217561
#> 74      8        0.8        V5 0.3869359
#> 75      8        0.8        V2 0.4202911
#> 76      8        0.8        V3 0.3653792
#> 77      8        0.8        V4 0.3991653
#> 78      8        0.8        V8 0.3410600
#> 79      8        0.8        V6 0.4044531
#> 80      8        0.8        V7 0.3567127
#> 81      9        0.9        V5 0.3862970
#> 82      9        0.9        V7 0.3561128
#> 83      9        0.9        V8 0.3405808
#> 84      9        0.9        V9 0.4071478
#> 85      9        0.9        V6 0.4039272
#> 86      9        0.9       V10 0.4299325
#> 87      9        0.9        V1 0.4211483
#> 88      9        0.9        V2 0.4197766
#> 89      9        0.9        V3 0.3647621
#> 90      9        0.9        V4 0.3985770
#> 91     10        1.0        V1 0.4205414
#> 92     10        1.0        V4 0.3979895
#> 93     10        1.0        V5 0.3856592
#> 94     10        1.0        V9 0.4066155
#> 95     10        1.0        V3 0.3641461
#> 96     10        1.0        V7 0.3555139
#> 97     10        1.0        V8 0.3401023
#> 98     10        1.0        V2 0.4192627
#> 99     10        1.0        V6 0.4034021
#> 100    10        1.0       V10 0.4292174


# get coefficient samples
coefs <- getNationalCoefficients(10)

# table of different scenarios to test
covTableSim <- expand.grid(Anthro = seq(0, 90, by = 20),
                           fire_excl_anthro = seq(0, 70, by = 20))
covTableSim$Total_dist = covTableSim$Anthro + covTableSim$fire_excl_anthro

estimateNationalRates(covTableSim, coefs)
#> popGrowthPars contains quantiles so they are used instead of the defaults
#> popGrowthPars contains quantiles so they are used instead of the defaults
#>    Anthro fire_excl_anthro Total_dist     S_bar   S_stdErr   S_PIlow  S_PIhigh
#> 1       0                0          0 0.8757906 0.04645088 0.7887908 0.9496866
#> 2       0               20         20 0.8757906 0.04645088 0.7887908 0.9496866
#> 3       0               40         40 0.8757906 0.04645088 0.7887908 0.9496866
#> 4       0               60         60 0.8757906 0.04645088 0.7887908 0.9496866
#> 5      20                0         20 0.8617131 0.04877288 0.7704674 0.9381235
#> 6      20               20         40 0.8617131 0.04877288 0.7704674 0.9381235
#> 7      20               40         60 0.8617131 0.04877288 0.7704674 0.9381235
#> 8      20               60         80 0.8617131 0.04877288 0.7704674 0.9381235
#> 9      40                0         40 0.8478591 0.05091374 0.7527691 0.9263462
#> 10     40               20         60 0.8478591 0.05091374 0.7527691 0.9263462
#> 11     40               40         80 0.8478591 0.05091374 0.7527691 0.9263462
#> 12     40               60        100 0.8478591 0.05091374 0.7527691 0.9263462
#> 13     60                0         60 0.8342249 0.05290062 0.7356375 0.9144317
#> 14     60               20         80 0.8342249 0.05290062 0.7356375 0.9144317
#> 15     60               40        100 0.8342249 0.05290062 0.7356375 0.9144317
#> 16     60               60        120 0.8342249 0.05290062 0.7356375 0.9144317
#> 17     80                0         80 0.8208071 0.05475400 0.7190265 0.9024365
#> 18     80               20        100 0.8208071 0.05475400 0.7190265 0.9024365
#> 19     80               40        120 0.8208071 0.05475400 0.7190265 0.9024365
#> 20     80               60        140 0.8208071 0.05475400 0.7190265 0.9024365
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11041328 0.181638405 0.5455950
#> 2  0.30574618 0.10414466 0.137149336 0.4878056
#> 3  0.26001915 0.09841776 0.102455094 0.4365454
#> 4  0.22113100 0.09297178 0.075526428 0.3911830
#> 5  0.25589195 0.09971022 0.104825489 0.4376795
#> 6  0.21762106 0.09316064 0.077360598 0.3921858
#> 7  0.18507391 0.08714509 0.056184104 0.3519730
#> 8  0.15739448 0.08150262 0.040016479 0.3164403
#> 9  0.18213629 0.08797835 0.057620730 0.3528617
#> 10 0.15489621 0.08167361 0.041107000 0.3172257
#> 11 0.13173012 0.07587617 0.028649349 0.2857278
#> 12 0.11202872 0.07048037 0.019407846 0.2578639
#> 13 0.12963921 0.07642107 0.029483839 0.2864242
#> 14 0.11025053 0.07062278 0.020020948 0.2584803
#> 15 0.09376159 0.06529370 0.013139082 0.2337284
#> 16 0.07973872 0.06036008 0.008270226 0.2117646
#> 17 0.09227334 0.06565375 0.013590437 0.2342764
#> 18 0.07847305 0.06047272 0.008584495 0.2122514
#> 19 0.06673672 0.05571709 0.005156047 0.1926652
#> 20 0.05675565 0.05133252 0.002911139 0.1752016
```
