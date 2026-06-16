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
#> 1         0.1 0.3838513 0.01653373 0.3674743 0.4121371
#> 2         0.2 0.3832760 0.01652487 0.3668945 0.4115770
#> 3         0.3 0.3827015 0.01651619 0.3663156 0.4110176
#> 4         0.4 0.3821279 0.01650768 0.3657375 0.4104589
#> 5         0.5 0.3815551 0.01649934 0.3651604 0.4099010
#> 6         0.6 0.3809832 0.01649117 0.3645843 0.4093439
#> 7         0.7 0.3804122 0.01648317 0.3640090 0.4087875
#> 8         0.8 0.3798420 0.01647534 0.3634346 0.4082319
#> 9         0.9 0.3792726 0.01646768 0.3628611 0.4076771
#> 10        1.0 0.3787041 0.01646018 0.3622886 0.4071230

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4042522
#> 2       1        0.1        V5 0.3713235
#> 3       1        0.1        V9 0.3807732
#> 4       1        0.1        V3 0.4124558
#> 5       1        0.1        V4 0.3960519
#> 6       1        0.1        V8 0.4025183
#> 7       1        0.1        V2 0.4110394
#> 8       1        0.1        V6 0.3663569
#> 9       1        0.1        V7 0.4066073
#> 10      1        0.1       V10 0.3880663
#> 11      2        0.2        V7 0.4059107
#> 12      2        0.2        V8 0.4018896
#> 13      2        0.2        V9 0.3802304
#> 14      2        0.2       V10 0.3875164
#> 15      2        0.2        V1 0.4037237
#> 16      2        0.2        V5 0.3707269
#> 17      2        0.2        V2 0.4104449
#> 18      2        0.2        V3 0.4119056
#> 19      2        0.2        V4 0.3955470
#> 20      2        0.2        V6 0.3657819
#> 21      3        0.3        V4 0.3950428
#> 22      3        0.3        V5 0.3701313
#> 23      3        0.3        V3 0.4113562
#> 24      3        0.3        V7 0.4052154
#> 25      3        0.3        V8 0.4012619
#> 26      3        0.3        V9 0.3796884
#> 27      3        0.3        V6 0.3652078
#> 28      3        0.3       V10 0.3869673
#> 29      3        0.3        V1 0.4031959
#> 30      3        0.3        V2 0.4098512
#> 31      4        0.4        V1 0.4026687
#> 32      4        0.4        V9 0.3791471
#> 33      4        0.4        V3 0.4108074
#> 34      4        0.4        V4 0.3945393
#> 35      4        0.4        V5 0.3695367
#> 36      4        0.4        V2 0.4092584
#> 37      4        0.4        V6 0.3646346
#> 38      4        0.4        V7 0.4045212
#> 39      4        0.4        V8 0.4006352
#> 40      4        0.4       V10 0.3864190
#> 41      5        0.5        V8 0.4000095
#> 42      5        0.5        V9 0.3786067
#> 43      5        0.5       V10 0.3858714
#> 44      5        0.5        V1 0.4021423
#> 45      5        0.5        V5 0.3689430
#> 46      5        0.5        V2 0.4086665
#> 47      5        0.5        V3 0.4102594
#> 48      5        0.5        V4 0.3940364
#> 49      5        0.5        V6 0.3640623
#> 50      5        0.5        V7 0.4038283
#> 51      6        0.6        V4 0.3935341
#> 52      6        0.6        V5 0.3683503
#> 53      6        0.6        V7 0.4031365
#> 54      6        0.6        V8 0.3993848
#> 55      6        0.6        V9 0.3780670
#> 56      6        0.6        V6 0.3634909
#> 57      6        0.6       V10 0.3853246
#> 58      6        0.6        V1 0.4016165
#> 59      6        0.6        V2 0.4080754
#> 60      6        0.6        V3 0.4097122
#> 61      7        0.7        V1 0.4010914
#> 62      7        0.7        V3 0.4091656
#> 63      7        0.7        V4 0.3930325
#> 64      7        0.7        V5 0.3677586
#> 65      7        0.7        V9 0.3775281
#> 66      7        0.7        V6 0.3629204
#> 67      7        0.7        V7 0.4024459
#> 68      7        0.7        V8 0.3987610
#> 69      7        0.7        V2 0.4074851
#> 70      7        0.7       V10 0.3847786
#> 71      8        0.8        V9 0.3769899
#> 72      8        0.8       V10 0.3842333
#> 73      8        0.8        V1 0.4005670
#> 74      8        0.8        V5 0.3671677
#> 75      8        0.8        V2 0.4068957
#> 76      8        0.8        V3 0.4086198
#> 77      8        0.8        V4 0.3925315
#> 78      8        0.8        V8 0.3981382
#> 79      8        0.8        V6 0.3623508
#> 80      8        0.8        V7 0.4017565
#> 81      9        0.9        V5 0.3665779
#> 82      9        0.9        V7 0.4010683
#> 83      9        0.9        V8 0.3975164
#> 84      9        0.9        V9 0.3764525
#> 85      9        0.9        V6 0.3617821
#> 86      9        0.9       V10 0.3836889
#> 87      9        0.9        V1 0.4000433
#> 88      9        0.9        V2 0.4063072
#> 89      9        0.9        V3 0.4080748
#> 90      9        0.9        V4 0.3920311
#> 91     10        1.0        V1 0.3995203
#> 92     10        1.0        V4 0.3915314
#> 93     10        1.0        V5 0.3659890
#> 94     10        1.0        V9 0.3759159
#> 95     10        1.0        V3 0.4075304
#> 96     10        1.0        V7 0.4003812
#> 97     10        1.0        V8 0.3968955
#> 98     10        1.0        V2 0.4057195
#> 99     10        1.0        V6 0.3612143
#> 100    10        1.0       V10 0.3831452

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4042522
#> 2       1        0.1        V5 0.3713235
#> 3       1        0.1        V9 0.3807732
#> 4       1        0.1        V3 0.4124558
#> 5       1        0.1        V4 0.3960519
#> 6       1        0.1        V8 0.4025183
#> 7       1        0.1        V2 0.4110394
#> 8       1        0.1        V6 0.3663569
#> 9       1        0.1        V7 0.4066073
#> 10      1        0.1       V10 0.3880663
#> 11      2        0.2        V7 0.4059107
#> 12      2        0.2        V8 0.4018896
#> 13      2        0.2        V9 0.3802304
#> 14      2        0.2       V10 0.3875164
#> 15      2        0.2        V1 0.4037237
#> 16      2        0.2        V5 0.3707269
#> 17      2        0.2        V2 0.4104449
#> 18      2        0.2        V3 0.4119056
#> 19      2        0.2        V4 0.3955470
#> 20      2        0.2        V6 0.3657819
#> 21      3        0.3        V4 0.3950428
#> 22      3        0.3        V5 0.3701313
#> 23      3        0.3        V3 0.4113562
#> 24      3        0.3        V7 0.4052154
#> 25      3        0.3        V8 0.4012619
#> 26      3        0.3        V9 0.3796884
#> 27      3        0.3        V6 0.3652078
#> 28      3        0.3       V10 0.3869673
#> 29      3        0.3        V1 0.4031959
#> 30      3        0.3        V2 0.4098512
#> 31      4        0.4        V1 0.4026687
#> 32      4        0.4        V9 0.3791471
#> 33      4        0.4        V3 0.4108074
#> 34      4        0.4        V4 0.3945393
#> 35      4        0.4        V5 0.3695367
#> 36      4        0.4        V2 0.4092584
#> 37      4        0.4        V6 0.3646346
#> 38      4        0.4        V7 0.4045212
#> 39      4        0.4        V8 0.4006352
#> 40      4        0.4       V10 0.3864190
#> 41      5        0.5        V8 0.4000095
#> 42      5        0.5        V9 0.3786067
#> 43      5        0.5       V10 0.3858714
#> 44      5        0.5        V1 0.4021423
#> 45      5        0.5        V5 0.3689430
#> 46      5        0.5        V2 0.4086665
#> 47      5        0.5        V3 0.4102594
#> 48      5        0.5        V4 0.3940364
#> 49      5        0.5        V6 0.3640623
#> 50      5        0.5        V7 0.4038283
#> 51      6        0.6        V4 0.3935341
#> 52      6        0.6        V5 0.3683503
#> 53      6        0.6        V7 0.4031365
#> 54      6        0.6        V8 0.3993848
#> 55      6        0.6        V9 0.3780670
#> 56      6        0.6        V6 0.3634909
#> 57      6        0.6       V10 0.3853246
#> 58      6        0.6        V1 0.4016165
#> 59      6        0.6        V2 0.4080754
#> 60      6        0.6        V3 0.4097122
#> 61      7        0.7        V1 0.4010914
#> 62      7        0.7        V3 0.4091656
#> 63      7        0.7        V4 0.3930325
#> 64      7        0.7        V5 0.3677586
#> 65      7        0.7        V9 0.3775281
#> 66      7        0.7        V6 0.3629204
#> 67      7        0.7        V7 0.4024459
#> 68      7        0.7        V8 0.3987610
#> 69      7        0.7        V2 0.4074851
#> 70      7        0.7       V10 0.3847786
#> 71      8        0.8        V9 0.3769899
#> 72      8        0.8       V10 0.3842333
#> 73      8        0.8        V1 0.4005670
#> 74      8        0.8        V5 0.3671677
#> 75      8        0.8        V2 0.4068957
#> 76      8        0.8        V3 0.4086198
#> 77      8        0.8        V4 0.3925315
#> 78      8        0.8        V8 0.3981382
#> 79      8        0.8        V6 0.3623508
#> 80      8        0.8        V7 0.4017565
#> 81      9        0.9        V5 0.3665779
#> 82      9        0.9        V7 0.4010683
#> 83      9        0.9        V8 0.3975164
#> 84      9        0.9        V9 0.3764525
#> 85      9        0.9        V6 0.3617821
#> 86      9        0.9       V10 0.3836889
#> 87      9        0.9        V1 0.4000433
#> 88      9        0.9        V2 0.4063072
#> 89      9        0.9        V3 0.4080748
#> 90      9        0.9        V4 0.3920311
#> 91     10        1.0        V1 0.3995203
#> 92     10        1.0        V4 0.3915314
#> 93     10        1.0        V5 0.3659890
#> 94     10        1.0        V9 0.3759159
#> 95     10        1.0        V3 0.4075304
#> 96     10        1.0        V7 0.4003812
#> 97     10        1.0        V8 0.3968955
#> 98     10        1.0        V2 0.4057195
#> 99     10        1.0        V6 0.3612143
#> 100    10        1.0       V10 0.3831452


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
#> 1       0                0          0 0.8757906 0.05100525 0.7649523 0.9489344
#> 2       0               20         20 0.8757906 0.05100525 0.7649523 0.9489344
#> 3       0               40         40 0.8757906 0.05100525 0.7649523 0.9489344
#> 4       0               60         60 0.8757906 0.05100525 0.7649523 0.9489344
#> 5      20                0         20 0.8617131 0.05387618 0.7456745 0.9406895
#> 6      20               20         40 0.8617131 0.05387618 0.7456745 0.9406895
#> 7      20               40         60 0.8617131 0.05387618 0.7456745 0.9406895
#> 8      20               60         80 0.8617131 0.05387618 0.7456745 0.9406895
#> 9      40                0         40 0.8478591 0.05654157 0.7270950 0.9323113
#> 10     40               20         60 0.8478591 0.05654157 0.7270950 0.9323113
#> 11     40               40         80 0.8478591 0.05654157 0.7270950 0.9323113
#> 12     40               60        100 0.8478591 0.05654157 0.7270950 0.9323113
#> 13     60                0         60 0.8342249 0.05902482 0.7091487 0.9238329
#> 14     60               20         80 0.8342249 0.05902482 0.7091487 0.9238329
#> 15     60               40        100 0.8342249 0.05902482 0.7091487 0.9238329
#> 16     60               60        120 0.8342249 0.05902482 0.7091487 0.9238329
#> 17     80                0         80 0.8208071 0.06134468 0.6917835 0.9152807
#> 18     80               20        100 0.8208071 0.06134468 0.6917835 0.9152807
#> 19     80               40        120 0.8208071 0.06134468 0.6917835 0.9152807
#> 20     80               60        140 0.8208071 0.06134468 0.6917835 0.9152807
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12982617 0.177024659 0.6012647
#> 2  0.30574618 0.12694280 0.132210227 0.5442181
#> 3  0.26001915 0.12299357 0.097577740 0.4928332
#> 4  0.22113100 0.11816194 0.070959565 0.4467068
#> 5  0.25589195 0.12094901 0.102647607 0.4993210
#> 6  0.21762106 0.11566658 0.074843798 0.4525246
#> 7  0.18507391 0.11004708 0.053618739 0.4105890
#> 8  0.15739448 0.10417335 0.037593618 0.3730436
#> 9  0.18213629 0.10952479 0.056704553 0.4158759
#> 10 0.15489621 0.10337867 0.039909851 0.3777763
#> 11 0.13173012 0.09725332 0.027381699 0.3436679
#> 12 0.11202872 0.09118499 0.018208169 0.3131145
#> 13 0.12963921 0.09744870 0.029181065 0.3479679
#> 14 0.11025053 0.09121883 0.019512966 0.3169681
#> 15 0.09376159 0.08518300 0.012572695 0.2891709
#> 16 0.07973872 0.07935287 0.007739115 0.2642029
#> 17 0.09227334 0.08569532 0.013549609 0.2926789
#> 18 0.07847305 0.07976810 0.008408700 0.2673567
#> 19 0.06673672 0.07411091 0.004941725 0.2445682
#> 20 0.05675565 0.06872200 0.002715308 0.2240050
```
