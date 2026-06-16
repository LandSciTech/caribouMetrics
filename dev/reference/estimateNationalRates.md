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
#> 1         0.1 0.3838513 0.02392626 0.3477371 0.4154090
#> 2         0.2 0.3832760 0.02395601 0.3470815 0.4148972
#> 3         0.3 0.3827015 0.02398571 0.3464271 0.4143861
#> 4         0.4 0.3821279 0.02401537 0.3457740 0.4138755
#> 5         0.5 0.3815551 0.02404499 0.3451221 0.4133656
#> 6         0.6 0.3809832 0.02407456 0.3444714 0.4128564
#> 7         0.7 0.3804122 0.02410409 0.3438220 0.4123477
#> 8         0.8 0.3798420 0.02413357 0.3431739 0.4118397
#> 9         0.9 0.3792726 0.02416301 0.3425269 0.4113324
#> 10        1.0 0.3787041 0.02419240 0.3418812 0.4108256

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3627122
#> 2       1        0.1        V5 0.3821052
#> 3       1        0.1        V9 0.3743900
#> 4       1        0.1        V3 0.3758358
#> 5       1        0.1        V4 0.4161135
#> 6       1        0.1        V8 0.3984481
#> 7       1        0.1        V2 0.4129825
#> 8       1        0.1        V6 0.3433894
#> 9       1        0.1        V7 0.4005134
#> 10      1        0.1       V10 0.4098128
#> 11      2        0.2        V7 0.3999187
#> 12      2        0.2        V8 0.3978444
#> 13      2        0.2        V9 0.3737938
#> 14      2        0.2       V10 0.4091640
#> 15      2        0.2        V1 0.3621253
#> 16      2        0.2        V5 0.3816043
#> 17      2        0.2        V2 0.4124678
#> 18      2        0.2        V3 0.3752079
#> 19      2        0.2        V4 0.4156025
#> 20      2        0.2        V6 0.3427139
#> 21      3        0.3        V4 0.4150922
#> 22      3        0.3        V5 0.3811040
#> 23      3        0.3        V3 0.3745811
#> 24      3        0.3        V7 0.3993249
#> 25      3        0.3        V8 0.3972416
#> 26      3        0.3        V9 0.3731986
#> 27      3        0.3        V6 0.3420397
#> 28      3        0.3       V10 0.4085163
#> 29      3        0.3        V1 0.3615394
#> 30      3        0.3        V2 0.4119538
#> 31      4        0.4        V1 0.3609544
#> 32      4        0.4        V9 0.3726044
#> 33      4        0.4        V3 0.3739553
#> 34      4        0.4        V4 0.4145825
#> 35      4        0.4        V5 0.3806043
#> 36      4        0.4        V2 0.4114403
#> 37      4        0.4        V6 0.3413667
#> 38      4        0.4        V7 0.3987320
#> 39      4        0.4        V8 0.3966397
#> 40      4        0.4       V10 0.4078696
#> 41      5        0.5        V8 0.3960387
#> 42      5        0.5        V9 0.3720111
#> 43      5        0.5       V10 0.4072240
#> 44      5        0.5        V1 0.3603704
#> 45      5        0.5        V5 0.3801054
#> 46      5        0.5        V2 0.4109276
#> 47      5        0.5        V3 0.3733305
#> 48      5        0.5        V4 0.4140735
#> 49      5        0.5        V6 0.3406952
#> 50      5        0.5        V7 0.3981400
#> 51      6        0.6        V4 0.4135651
#> 52      6        0.6        V5 0.3796070
#> 53      6        0.6        V7 0.3975488
#> 54      6        0.6        V8 0.3954387
#> 55      6        0.6        V9 0.3714188
#> 56      6        0.6        V6 0.3400249
#> 57      6        0.6       V10 0.4065793
#> 58      6        0.6        V1 0.3597873
#> 59      6        0.6        V2 0.4104154
#> 60      6        0.6        V3 0.3727068
#> 61      7        0.7        V1 0.3592052
#> 62      7        0.7        V3 0.3720841
#> 63      7        0.7        V4 0.4130572
#> 64      7        0.7        V5 0.3791093
#> 65      7        0.7        V9 0.3708274
#> 66      7        0.7        V6 0.3393560
#> 67      7        0.7        V7 0.3969585
#> 68      7        0.7        V8 0.3948395
#> 69      7        0.7        V2 0.4099039
#> 70      7        0.7       V10 0.4059357
#> 71      8        0.8        V9 0.3702369
#> 72      8        0.8       V10 0.4052931
#> 73      8        0.8        V1 0.3586240
#> 74      8        0.8        V5 0.3786123
#> 75      8        0.8        V2 0.4093930
#> 76      8        0.8        V3 0.3714625
#> 77      8        0.8        V4 0.4125501
#> 78      8        0.8        V8 0.3942412
#> 79      8        0.8        V6 0.3386883
#> 80      8        0.8        V7 0.3963691
#> 81      9        0.9        V5 0.3781159
#> 82      9        0.9        V7 0.3957806
#> 83      9        0.9        V8 0.3936439
#> 84      9        0.9        V9 0.3696474
#> 85      9        0.9        V6 0.3380220
#> 86      9        0.9       V10 0.4046515
#> 87      9        0.9        V1 0.3580437
#> 88      9        0.9        V2 0.4088828
#> 89      9        0.9        V3 0.3708419
#> 90      9        0.9        V4 0.4120435
#> 91     10        1.0        V1 0.3574644
#> 92     10        1.0        V4 0.4115376
#> 93     10        1.0        V5 0.3776202
#> 94     10        1.0        V9 0.3690588
#> 95     10        1.0        V3 0.3702223
#> 96     10        1.0        V7 0.3951929
#> 97     10        1.0        V8 0.3930475
#> 98     10        1.0        V2 0.4083732
#> 99     10        1.0        V6 0.3373570
#> 100    10        1.0       V10 0.4040109

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3627122
#> 2       1        0.1        V5 0.3821052
#> 3       1        0.1        V9 0.3743900
#> 4       1        0.1        V3 0.3758358
#> 5       1        0.1        V4 0.4161135
#> 6       1        0.1        V8 0.3984481
#> 7       1        0.1        V2 0.4129825
#> 8       1        0.1        V6 0.3433894
#> 9       1        0.1        V7 0.4005134
#> 10      1        0.1       V10 0.4098128
#> 11      2        0.2        V7 0.3999187
#> 12      2        0.2        V8 0.3978444
#> 13      2        0.2        V9 0.3737938
#> 14      2        0.2       V10 0.4091640
#> 15      2        0.2        V1 0.3621253
#> 16      2        0.2        V5 0.3816043
#> 17      2        0.2        V2 0.4124678
#> 18      2        0.2        V3 0.3752079
#> 19      2        0.2        V4 0.4156025
#> 20      2        0.2        V6 0.3427139
#> 21      3        0.3        V4 0.4150922
#> 22      3        0.3        V5 0.3811040
#> 23      3        0.3        V3 0.3745811
#> 24      3        0.3        V7 0.3993249
#> 25      3        0.3        V8 0.3972416
#> 26      3        0.3        V9 0.3731986
#> 27      3        0.3        V6 0.3420397
#> 28      3        0.3       V10 0.4085163
#> 29      3        0.3        V1 0.3615394
#> 30      3        0.3        V2 0.4119538
#> 31      4        0.4        V1 0.3609544
#> 32      4        0.4        V9 0.3726044
#> 33      4        0.4        V3 0.3739553
#> 34      4        0.4        V4 0.4145825
#> 35      4        0.4        V5 0.3806043
#> 36      4        0.4        V2 0.4114403
#> 37      4        0.4        V6 0.3413667
#> 38      4        0.4        V7 0.3987320
#> 39      4        0.4        V8 0.3966397
#> 40      4        0.4       V10 0.4078696
#> 41      5        0.5        V8 0.3960387
#> 42      5        0.5        V9 0.3720111
#> 43      5        0.5       V10 0.4072240
#> 44      5        0.5        V1 0.3603704
#> 45      5        0.5        V5 0.3801054
#> 46      5        0.5        V2 0.4109276
#> 47      5        0.5        V3 0.3733305
#> 48      5        0.5        V4 0.4140735
#> 49      5        0.5        V6 0.3406952
#> 50      5        0.5        V7 0.3981400
#> 51      6        0.6        V4 0.4135651
#> 52      6        0.6        V5 0.3796070
#> 53      6        0.6        V7 0.3975488
#> 54      6        0.6        V8 0.3954387
#> 55      6        0.6        V9 0.3714188
#> 56      6        0.6        V6 0.3400249
#> 57      6        0.6       V10 0.4065793
#> 58      6        0.6        V1 0.3597873
#> 59      6        0.6        V2 0.4104154
#> 60      6        0.6        V3 0.3727068
#> 61      7        0.7        V1 0.3592052
#> 62      7        0.7        V3 0.3720841
#> 63      7        0.7        V4 0.4130572
#> 64      7        0.7        V5 0.3791093
#> 65      7        0.7        V9 0.3708274
#> 66      7        0.7        V6 0.3393560
#> 67      7        0.7        V7 0.3969585
#> 68      7        0.7        V8 0.3948395
#> 69      7        0.7        V2 0.4099039
#> 70      7        0.7       V10 0.4059357
#> 71      8        0.8        V9 0.3702369
#> 72      8        0.8       V10 0.4052931
#> 73      8        0.8        V1 0.3586240
#> 74      8        0.8        V5 0.3786123
#> 75      8        0.8        V2 0.4093930
#> 76      8        0.8        V3 0.3714625
#> 77      8        0.8        V4 0.4125501
#> 78      8        0.8        V8 0.3942412
#> 79      8        0.8        V6 0.3386883
#> 80      8        0.8        V7 0.3963691
#> 81      9        0.9        V5 0.3781159
#> 82      9        0.9        V7 0.3957806
#> 83      9        0.9        V8 0.3936439
#> 84      9        0.9        V9 0.3696474
#> 85      9        0.9        V6 0.3380220
#> 86      9        0.9       V10 0.4046515
#> 87      9        0.9        V1 0.3580437
#> 88      9        0.9        V2 0.4088828
#> 89      9        0.9        V3 0.3708419
#> 90      9        0.9        V4 0.4120435
#> 91     10        1.0        V1 0.3574644
#> 92     10        1.0        V4 0.4115376
#> 93     10        1.0        V5 0.3776202
#> 94     10        1.0        V9 0.3690588
#> 95     10        1.0        V3 0.3702223
#> 96     10        1.0        V7 0.3951929
#> 97     10        1.0        V8 0.3930475
#> 98     10        1.0        V2 0.4083732
#> 99     10        1.0        V6 0.3373570
#> 100    10        1.0       V10 0.4040109


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
#> 1       0                0          0 0.8757906 0.04467556 0.7886374 0.9377437
#> 2       0               20         20 0.8757906 0.04467556 0.7886374 0.9377437
#> 3       0               40         40 0.8757906 0.04467556 0.7886374 0.9377437
#> 4       0               60         60 0.8757906 0.04467556 0.7886374 0.9377437
#> 5      20                0         20 0.8617131 0.04673028 0.7741115 0.9310478
#> 6      20               20         40 0.8617131 0.04673028 0.7741115 0.9310478
#> 7      20               40         60 0.8617131 0.04673028 0.7741115 0.9310478
#> 8      20               60         80 0.8617131 0.04673028 0.7741115 0.9310478
#> 9      40                0         40 0.8478591 0.04871985 0.7599813 0.9243022
#> 10     40               20         60 0.8478591 0.04871985 0.7599813 0.9243022
#> 11     40               40         80 0.8478591 0.04871985 0.7599813 0.9243022
#> 12     40               60        100 0.8478591 0.04871985 0.7599813 0.9243022
#> 13     60                0         60 0.8342249 0.05064534 0.7462166 0.9175192
#> 14     60               20         80 0.8342249 0.05064534 0.7462166 0.9175192
#> 15     60               40        100 0.8342249 0.05064534 0.7462166 0.9175192
#> 16     60               60        120 0.8342249 0.05064534 0.7462166 0.9175192
#> 17     80                0         80 0.8208071 0.05250806 0.7327922 0.9107090
#> 18     80               20        100 0.8208071 0.05250806 0.7327922 0.9107090
#> 19     80               40        120 0.8208071 0.05250806 0.7327922 0.9107090
#> 20     80               60        140 0.8208071 0.05250806 0.7327922 0.9107090
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11699767 0.173181646 0.5602325
#> 2  0.30574618 0.11246774 0.127952179 0.5085488
#> 3  0.26001915 0.10798791 0.093317997 0.4619131
#> 4  0.22113100 0.10322761 0.066961093 0.4199267
#> 5  0.25589195 0.11088738 0.091029989 0.4561239
#> 6  0.21762106 0.10429503 0.065228132 0.4147193
#> 7  0.18507391 0.09827075 0.045793843 0.3774956
#> 8  0.15739448 0.09250193 0.031350290 0.3440488
#> 9  0.18213629 0.10083390 0.044525337 0.3728811
#> 10 0.15489621 0.09371522 0.030416221 0.3399029
#> 11 0.13173012 0.08733815 0.020132694 0.3102705
#> 12 0.11202872 0.08143948 0.012816272 0.2836323
#> 13 0.12963921 0.08942293 0.019476495 0.3065966
#> 14 0.11025053 0.08254444 0.012357093 0.2803283
#> 15 0.09376159 0.07642311 0.007463934 0.2566916
#> 16 0.07973872 0.07084634 0.004241009 0.2353971
#> 17 0.09227334 0.07805205 0.007164087 0.2537578
#> 18 0.07847305 0.07176646 0.004049273 0.2327518
#> 19 0.06673672 0.06618171 0.002118342 0.2137943
#> 20 0.05675565 0.06112858 0.001006619 0.1966534
```
