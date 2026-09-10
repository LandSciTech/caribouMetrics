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
#> 1         0.1 0.3838513 0.01997574 0.3486924 0.4106717
#> 2         0.2 0.3832760 0.01994409 0.3482327 0.4101114
#> 3         0.3 0.3827015 0.01991258 0.3477737 0.4095518
#> 4         0.4 0.3821279 0.01988119 0.3473153 0.4089931
#> 5         0.5 0.3815551 0.01984994 0.3468575 0.4084350
#> 6         0.6 0.3809832 0.01981883 0.3464003 0.4078777
#> 7         0.7 0.3804122 0.01978784 0.3459437 0.4073212
#> 8         0.8 0.3798420 0.01975698 0.3454878 0.4067655
#> 9         0.9 0.3792726 0.01972626 0.3450324 0.4062105
#> 10        1.0 0.3787041 0.01969567 0.3445776 0.4056563

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4032941
#> 2       1        0.1        V5 0.3763029
#> 3       1        0.1        V9 0.3455169
#> 4       1        0.1        V3 0.4128136
#> 5       1        0.1        V4 0.3703960
#> 6       1        0.1        V8 0.3746600
#> 7       1        0.1        V2 0.3765711
#> 8       1        0.1        V6 0.3891441
#> 9       1        0.1        V7 0.3596302
#> 10      1        0.1       V10 0.3917523
#> 11      2        0.2        V7 0.3590854
#> 12      2        0.2        V8 0.3741470
#> 13      2        0.2        V9 0.3450820
#> 14      2        0.2       V10 0.3911149
#> 15      2        0.2        V1 0.4027447
#> 16      2        0.2        V5 0.3757989
#> 17      2        0.2        V2 0.3760192
#> 18      2        0.2        V3 0.4122501
#> 19      2        0.2        V4 0.3698668
#> 20      2        0.2        V6 0.3886374
#> 21      3        0.3        V4 0.3693383
#> 22      3        0.3        V5 0.3752955
#> 23      3        0.3        V3 0.4116874
#> 24      3        0.3        V7 0.3585415
#> 25      3        0.3        V8 0.3736347
#> 26      3        0.3        V9 0.3446476
#> 27      3        0.3        V6 0.3881313
#> 28      3        0.3       V10 0.3904786
#> 29      3        0.3        V1 0.4021960
#> 30      3        0.3        V2 0.3754682
#> 31      4        0.4        V1 0.4016481
#> 32      4        0.4        V9 0.3442138
#> 33      4        0.4        V3 0.4111255
#> 34      4        0.4        V4 0.3688106
#> 35      4        0.4        V5 0.3747928
#> 36      4        0.4        V2 0.3749180
#> 37      4        0.4        V6 0.3876259
#> 38      4        0.4        V7 0.3579984
#> 39      4        0.4        V8 0.3731231
#> 40      4        0.4       V10 0.3898433
#> 41      5        0.5        V8 0.3726122
#> 42      5        0.5        V9 0.3437805
#> 43      5        0.5       V10 0.3892090
#> 44      5        0.5        V1 0.4011009
#> 45      5        0.5        V5 0.3742908
#> 46      5        0.5        V2 0.3743685
#> 47      5        0.5        V3 0.4105643
#> 48      5        0.5        V4 0.3682836
#> 49      5        0.5        V6 0.3871212
#> 50      5        0.5        V7 0.3574561
#> 51      6        0.6        V4 0.3677574
#> 52      6        0.6        V5 0.3737894
#> 53      6        0.6        V7 0.3569146
#> 54      6        0.6        V8 0.3721020
#> 55      6        0.6        V9 0.3433478
#> 56      6        0.6        V6 0.3866171
#> 57      6        0.6       V10 0.3885758
#> 58      6        0.6        V1 0.4005544
#> 59      6        0.6        V2 0.3738199
#> 60      6        0.6        V3 0.4100039
#> 61      7        0.7        V1 0.4000087
#> 62      7        0.7        V3 0.4094442
#> 63      7        0.7        V4 0.3672320
#> 64      7        0.7        V5 0.3732887
#> 65      7        0.7        V9 0.3429156
#> 66      7        0.7        V6 0.3861137
#> 67      7        0.7        V7 0.3563740
#> 68      7        0.7        V8 0.3715925
#> 69      7        0.7        V2 0.3732721
#> 70      7        0.7       V10 0.3879436
#> 71      8        0.8        V9 0.3424840
#> 72      8        0.8       V10 0.3873124
#> 73      8        0.8        V1 0.3994638
#> 74      8        0.8        V5 0.3727887
#> 75      8        0.8        V2 0.3727251
#> 76      8        0.8        V3 0.4088853
#> 77      8        0.8        V4 0.3667073
#> 78      8        0.8        V8 0.3710837
#> 79      8        0.8        V6 0.3856110
#> 80      8        0.8        V7 0.3558341
#> 81      9        0.9        V5 0.3722894
#> 82      9        0.9        V7 0.3552951
#> 83      9        0.9        V8 0.3705756
#> 84      9        0.9        V9 0.3420529
#> 85      9        0.9        V6 0.3851089
#> 86      9        0.9       V10 0.3866823
#> 87      9        0.9        V1 0.3989196
#> 88      9        0.9        V2 0.3721789
#> 89      9        0.9        V3 0.4083272
#> 90      9        0.9        V4 0.3661833
#> 91     10        1.0        V1 0.3983761
#> 92     10        1.0        V4 0.3656601
#> 93     10        1.0        V5 0.3717907
#> 94     10        1.0        V9 0.3416223
#> 95     10        1.0        V3 0.4077699
#> 96     10        1.0        V7 0.3547569
#> 97     10        1.0        V8 0.3700682
#> 98     10        1.0        V2 0.3716334
#> 99     10        1.0        V6 0.3846074
#> 100    10        1.0       V10 0.3860532

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4032941
#> 2       1        0.1        V5 0.3763029
#> 3       1        0.1        V9 0.3455169
#> 4       1        0.1        V3 0.4128136
#> 5       1        0.1        V4 0.3703960
#> 6       1        0.1        V8 0.3746600
#> 7       1        0.1        V2 0.3765711
#> 8       1        0.1        V6 0.3891441
#> 9       1        0.1        V7 0.3596302
#> 10      1        0.1       V10 0.3917523
#> 11      2        0.2        V7 0.3590854
#> 12      2        0.2        V8 0.3741470
#> 13      2        0.2        V9 0.3450820
#> 14      2        0.2       V10 0.3911149
#> 15      2        0.2        V1 0.4027447
#> 16      2        0.2        V5 0.3757989
#> 17      2        0.2        V2 0.3760192
#> 18      2        0.2        V3 0.4122501
#> 19      2        0.2        V4 0.3698668
#> 20      2        0.2        V6 0.3886374
#> 21      3        0.3        V4 0.3693383
#> 22      3        0.3        V5 0.3752955
#> 23      3        0.3        V3 0.4116874
#> 24      3        0.3        V7 0.3585415
#> 25      3        0.3        V8 0.3736347
#> 26      3        0.3        V9 0.3446476
#> 27      3        0.3        V6 0.3881313
#> 28      3        0.3       V10 0.3904786
#> 29      3        0.3        V1 0.4021960
#> 30      3        0.3        V2 0.3754682
#> 31      4        0.4        V1 0.4016481
#> 32      4        0.4        V9 0.3442138
#> 33      4        0.4        V3 0.4111255
#> 34      4        0.4        V4 0.3688106
#> 35      4        0.4        V5 0.3747928
#> 36      4        0.4        V2 0.3749180
#> 37      4        0.4        V6 0.3876259
#> 38      4        0.4        V7 0.3579984
#> 39      4        0.4        V8 0.3731231
#> 40      4        0.4       V10 0.3898433
#> 41      5        0.5        V8 0.3726122
#> 42      5        0.5        V9 0.3437805
#> 43      5        0.5       V10 0.3892090
#> 44      5        0.5        V1 0.4011009
#> 45      5        0.5        V5 0.3742908
#> 46      5        0.5        V2 0.3743685
#> 47      5        0.5        V3 0.4105643
#> 48      5        0.5        V4 0.3682836
#> 49      5        0.5        V6 0.3871212
#> 50      5        0.5        V7 0.3574561
#> 51      6        0.6        V4 0.3677574
#> 52      6        0.6        V5 0.3737894
#> 53      6        0.6        V7 0.3569146
#> 54      6        0.6        V8 0.3721020
#> 55      6        0.6        V9 0.3433478
#> 56      6        0.6        V6 0.3866171
#> 57      6        0.6       V10 0.3885758
#> 58      6        0.6        V1 0.4005544
#> 59      6        0.6        V2 0.3738199
#> 60      6        0.6        V3 0.4100039
#> 61      7        0.7        V1 0.4000087
#> 62      7        0.7        V3 0.4094442
#> 63      7        0.7        V4 0.3672320
#> 64      7        0.7        V5 0.3732887
#> 65      7        0.7        V9 0.3429156
#> 66      7        0.7        V6 0.3861137
#> 67      7        0.7        V7 0.3563740
#> 68      7        0.7        V8 0.3715925
#> 69      7        0.7        V2 0.3732721
#> 70      7        0.7       V10 0.3879436
#> 71      8        0.8        V9 0.3424840
#> 72      8        0.8       V10 0.3873124
#> 73      8        0.8        V1 0.3994638
#> 74      8        0.8        V5 0.3727887
#> 75      8        0.8        V2 0.3727251
#> 76      8        0.8        V3 0.4088853
#> 77      8        0.8        V4 0.3667073
#> 78      8        0.8        V8 0.3710837
#> 79      8        0.8        V6 0.3856110
#> 80      8        0.8        V7 0.3558341
#> 81      9        0.9        V5 0.3722894
#> 82      9        0.9        V7 0.3552951
#> 83      9        0.9        V8 0.3705756
#> 84      9        0.9        V9 0.3420529
#> 85      9        0.9        V6 0.3851089
#> 86      9        0.9       V10 0.3866823
#> 87      9        0.9        V1 0.3989196
#> 88      9        0.9        V2 0.3721789
#> 89      9        0.9        V3 0.4083272
#> 90      9        0.9        V4 0.3661833
#> 91     10        1.0        V1 0.3983761
#> 92     10        1.0        V4 0.3656601
#> 93     10        1.0        V5 0.3717907
#> 94     10        1.0        V9 0.3416223
#> 95     10        1.0        V3 0.4077699
#> 96     10        1.0        V7 0.3547569
#> 97     10        1.0        V8 0.3700682
#> 98     10        1.0        V2 0.3716334
#> 99     10        1.0        V6 0.3846074
#> 100    10        1.0       V10 0.3860532


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
#> 1       0                0          0 0.8757906 0.04541497 0.7992461 0.9482765
#> 2       0               20         20 0.8757906 0.04541497 0.7992461 0.9482765
#> 3       0               40         40 0.8757906 0.04541497 0.7992461 0.9482765
#> 4       0               60         60 0.8757906 0.04541497 0.7992461 0.9482765
#> 5      20                0         20 0.8617131 0.04829298 0.7794644 0.9387651
#> 6      20               20         40 0.8617131 0.04829298 0.7794644 0.9387651
#> 7      20               40         60 0.8617131 0.04829298 0.7794644 0.9387651
#> 8      20               60         80 0.8617131 0.04829298 0.7794644 0.9387651
#> 9      40                0         40 0.8478591 0.05099457 0.7604321 0.9290915
#> 10     40               20         60 0.8478591 0.05099457 0.7604321 0.9290915
#> 11     40               40         80 0.8478591 0.05099457 0.7604321 0.9290915
#> 12     40               60        100 0.8478591 0.05099457 0.7604321 0.9290915
#> 13     60                0         60 0.8342249 0.05353891 0.7420686 0.9193024
#> 14     60               20         80 0.8342249 0.05353891 0.7420686 0.9193024
#> 15     60               40        100 0.8342249 0.05353891 0.7420686 0.9193024
#> 16     60               60        120 0.8342249 0.05353891 0.7420686 0.9193024
#> 17     80                0         80 0.8208071 0.05594114 0.7243121 0.9094339
#> 18     80               20        100 0.8208071 0.05594114 0.7243121 0.9094339
#> 19     80               40        120 0.8208071 0.05594114 0.7243121 0.9094339
#> 20     80               60        140 0.8208071 0.05594114 0.7243121 0.9094339
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12730456 0.148376873 0.5715158
#> 2  0.30574618 0.12104717 0.115182878 0.5110390
#> 3  0.26001915 0.11471901 0.088557702 0.4573305
#> 4  0.22113100 0.10825849 0.067304197 0.4097754
#> 5  0.25589195 0.11091820 0.080360721 0.4554033
#> 6  0.21762106 0.10410785 0.060789644 0.4080707
#> 7  0.18507391 0.09753744 0.045316589 0.3662274
#> 8  0.15739448 0.09114002 0.033200665 0.3292554
#> 9  0.18213629 0.09526640 0.040612967 0.3647282
#> 10 0.15489621 0.08870497 0.029547346 0.3279307
#> 11 0.13173012 0.08249134 0.021031504 0.2954073
#> 12 0.11202872 0.07656927 0.014586078 0.2666346
#> 13 0.12963921 0.08111631 0.018500710 0.2942414
#> 14 0.11025053 0.07514013 0.012697046 0.2656025
#> 15 0.09376159 0.06953135 0.008431936 0.2402265
#> 16 0.07973872 0.06424540 0.005384248 0.2176957
#> 17 0.09227334 0.06870485 0.007212367 0.2393152
#> 18 0.07847305 0.06342462 0.004533067 0.2168857
#> 19 0.06673672 0.05849215 0.002710750 0.1969181
#> 20 0.05675565 0.05387276 0.001527729 0.1790879
```
