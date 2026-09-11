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
#> 1         0.1 0.3838513 0.01826857 0.3542525 0.4075959
#> 2         0.2 0.3832760 0.01829582 0.3536437 0.4070814
#> 3         0.3 0.3827015 0.01832314 0.3530360 0.4065676
#> 4         0.4 0.3821279 0.01835053 0.3524293 0.4060544
#> 5         0.5 0.3815551 0.01837799 0.3518237 0.4055419
#> 6         0.6 0.3809832 0.01840552 0.3512191 0.4050300
#> 7         0.7 0.3804122 0.01843311 0.3506156 0.4045188
#> 8         0.8 0.3798420 0.01846076 0.3500132 0.4040082
#> 9         0.9 0.3792726 0.01848846 0.3494118 0.4034983
#> 10        1.0 0.3787041 0.01851623 0.3488114 0.4029890

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3637927
#> 2       1        0.1        V5 0.3776920
#> 3       1        0.1        V9 0.4027052
#> 4       1        0.1        V3 0.3731877
#> 5       1        0.1        V4 0.3688440
#> 6       1        0.1        V8 0.3524545
#> 7       1        0.1        V2 0.4090157
#> 8       1        0.1        V6 0.3604455
#> 9       1        0.1        V7 0.3618967
#> 10      1        0.1       V10 0.3702457
#> 11      2        0.2        V7 0.3613647
#> 12      2        0.2        V8 0.3518125
#> 13      2        0.2        V9 0.4022720
#> 14      2        0.2       V10 0.3697438
#> 15      2        0.2        V1 0.3632404
#> 16      2        0.2        V5 0.3771733
#> 17      2        0.2        V2 0.4084777
#> 18      2        0.2        V3 0.3727802
#> 19      2        0.2        V4 0.3683226
#> 20      2        0.2        V6 0.3599510
#> 21      3        0.3        V4 0.3678020
#> 22      3        0.3        V5 0.3766553
#> 23      3        0.3        V3 0.3723730
#> 24      3        0.3        V7 0.3608335
#> 25      3        0.3        V8 0.3511717
#> 26      3        0.3        V9 0.4018392
#> 27      3        0.3        V6 0.3594572
#> 28      3        0.3       V10 0.3692426
#> 29      3        0.3        V1 0.3626891
#> 30      3        0.3        V2 0.4079403
#> 31      4        0.4        V1 0.3621385
#> 32      4        0.4        V9 0.4014069
#> 33      4        0.4        V3 0.3719664
#> 34      4        0.4        V4 0.3672820
#> 35      4        0.4        V5 0.3761381
#> 36      4        0.4        V2 0.4074037
#> 37      4        0.4        V6 0.3589640
#> 38      4        0.4        V7 0.3603031
#> 39      4        0.4        V8 0.3505321
#> 40      4        0.4       V10 0.3687421
#> 41      5        0.5        V8 0.3498937
#> 42      5        0.5        V9 0.4009751
#> 43      5        0.5       V10 0.3682423
#> 44      5        0.5        V1 0.3615888
#> 45      5        0.5        V5 0.3756215
#> 46      5        0.5        V2 0.4068677
#> 47      5        0.5        V3 0.3715602
#> 48      5        0.5        V4 0.3667628
#> 49      5        0.5        V6 0.3584716
#> 50      5        0.5        V7 0.3597735
#> 51      6        0.6        V4 0.3662444
#> 52      6        0.6        V5 0.3751057
#> 53      6        0.6        V7 0.3592446
#> 54      6        0.6        V8 0.3492564
#> 55      6        0.6        V9 0.4005438
#> 56      6        0.6        V6 0.3579798
#> 57      6        0.6       V10 0.3677431
#> 58      6        0.6        V1 0.3610399
#> 59      6        0.6        V2 0.4063325
#> 60      6        0.6        V3 0.3711544
#> 61      7        0.7        V1 0.3604919
#> 62      7        0.7        V3 0.3707491
#> 63      7        0.7        V4 0.3657266
#> 64      7        0.7        V5 0.3745906
#> 65      7        0.7        V9 0.4001129
#> 66      7        0.7        V6 0.3574887
#> 67      7        0.7        V7 0.3587165
#> 68      7        0.7        V8 0.3486203
#> 69      7        0.7        V2 0.4057979
#> 70      7        0.7       V10 0.3672447
#> 71      8        0.8        V9 0.3996825
#> 72      8        0.8       V10 0.3667469
#> 73      8        0.8        V1 0.3599447
#> 74      8        0.8        V5 0.3740762
#> 75      8        0.8        V2 0.4052641
#> 76      8        0.8        V3 0.3703442
#> 77      8        0.8        V4 0.3652096
#> 78      8        0.8        V8 0.3479853
#> 79      8        0.8        V6 0.3569982
#> 80      8        0.8        V7 0.3581892
#> 81      9        0.9        V5 0.3735625
#> 82      9        0.9        V7 0.3576627
#> 83      9        0.9        V8 0.3473515
#> 84      9        0.9        V9 0.3992525
#> 85      9        0.9        V6 0.3565085
#> 86      9        0.9       V10 0.3662497
#> 87      9        0.9        V1 0.3593983
#> 88      9        0.9        V2 0.4047310
#> 89      9        0.9        V3 0.3699397
#> 90      9        0.9        V4 0.3646934
#> 91     10        1.0        V1 0.3588528
#> 92     10        1.0        V4 0.3641778
#> 93     10        1.0        V5 0.3730495
#> 94     10        1.0        V9 0.3988230
#> 95     10        1.0        V3 0.3695357
#> 96     10        1.0        V7 0.3571369
#> 97     10        1.0        V8 0.3467188
#> 98     10        1.0        V2 0.4041985
#> 99     10        1.0        V6 0.3560194
#> 100    10        1.0       V10 0.3657533

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3637927
#> 2       1        0.1        V5 0.3776920
#> 3       1        0.1        V9 0.4027052
#> 4       1        0.1        V3 0.3731877
#> 5       1        0.1        V4 0.3688440
#> 6       1        0.1        V8 0.3524545
#> 7       1        0.1        V2 0.4090157
#> 8       1        0.1        V6 0.3604455
#> 9       1        0.1        V7 0.3618967
#> 10      1        0.1       V10 0.3702457
#> 11      2        0.2        V7 0.3613647
#> 12      2        0.2        V8 0.3518125
#> 13      2        0.2        V9 0.4022720
#> 14      2        0.2       V10 0.3697438
#> 15      2        0.2        V1 0.3632404
#> 16      2        0.2        V5 0.3771733
#> 17      2        0.2        V2 0.4084777
#> 18      2        0.2        V3 0.3727802
#> 19      2        0.2        V4 0.3683226
#> 20      2        0.2        V6 0.3599510
#> 21      3        0.3        V4 0.3678020
#> 22      3        0.3        V5 0.3766553
#> 23      3        0.3        V3 0.3723730
#> 24      3        0.3        V7 0.3608335
#> 25      3        0.3        V8 0.3511717
#> 26      3        0.3        V9 0.4018392
#> 27      3        0.3        V6 0.3594572
#> 28      3        0.3       V10 0.3692426
#> 29      3        0.3        V1 0.3626891
#> 30      3        0.3        V2 0.4079403
#> 31      4        0.4        V1 0.3621385
#> 32      4        0.4        V9 0.4014069
#> 33      4        0.4        V3 0.3719664
#> 34      4        0.4        V4 0.3672820
#> 35      4        0.4        V5 0.3761381
#> 36      4        0.4        V2 0.4074037
#> 37      4        0.4        V6 0.3589640
#> 38      4        0.4        V7 0.3603031
#> 39      4        0.4        V8 0.3505321
#> 40      4        0.4       V10 0.3687421
#> 41      5        0.5        V8 0.3498937
#> 42      5        0.5        V9 0.4009751
#> 43      5        0.5       V10 0.3682423
#> 44      5        0.5        V1 0.3615888
#> 45      5        0.5        V5 0.3756215
#> 46      5        0.5        V2 0.4068677
#> 47      5        0.5        V3 0.3715602
#> 48      5        0.5        V4 0.3667628
#> 49      5        0.5        V6 0.3584716
#> 50      5        0.5        V7 0.3597735
#> 51      6        0.6        V4 0.3662444
#> 52      6        0.6        V5 0.3751057
#> 53      6        0.6        V7 0.3592446
#> 54      6        0.6        V8 0.3492564
#> 55      6        0.6        V9 0.4005438
#> 56      6        0.6        V6 0.3579798
#> 57      6        0.6       V10 0.3677431
#> 58      6        0.6        V1 0.3610399
#> 59      6        0.6        V2 0.4063325
#> 60      6        0.6        V3 0.3711544
#> 61      7        0.7        V1 0.3604919
#> 62      7        0.7        V3 0.3707491
#> 63      7        0.7        V4 0.3657266
#> 64      7        0.7        V5 0.3745906
#> 65      7        0.7        V9 0.4001129
#> 66      7        0.7        V6 0.3574887
#> 67      7        0.7        V7 0.3587165
#> 68      7        0.7        V8 0.3486203
#> 69      7        0.7        V2 0.4057979
#> 70      7        0.7       V10 0.3672447
#> 71      8        0.8        V9 0.3996825
#> 72      8        0.8       V10 0.3667469
#> 73      8        0.8        V1 0.3599447
#> 74      8        0.8        V5 0.3740762
#> 75      8        0.8        V2 0.4052641
#> 76      8        0.8        V3 0.3703442
#> 77      8        0.8        V4 0.3652096
#> 78      8        0.8        V8 0.3479853
#> 79      8        0.8        V6 0.3569982
#> 80      8        0.8        V7 0.3581892
#> 81      9        0.9        V5 0.3735625
#> 82      9        0.9        V7 0.3576627
#> 83      9        0.9        V8 0.3473515
#> 84      9        0.9        V9 0.3992525
#> 85      9        0.9        V6 0.3565085
#> 86      9        0.9       V10 0.3662497
#> 87      9        0.9        V1 0.3593983
#> 88      9        0.9        V2 0.4047310
#> 89      9        0.9        V3 0.3699397
#> 90      9        0.9        V4 0.3646934
#> 91     10        1.0        V1 0.3588528
#> 92     10        1.0        V4 0.3641778
#> 93     10        1.0        V5 0.3730495
#> 94     10        1.0        V9 0.3988230
#> 95     10        1.0        V3 0.3695357
#> 96     10        1.0        V7 0.3571369
#> 97     10        1.0        V8 0.3467188
#> 98     10        1.0        V2 0.4041985
#> 99     10        1.0        V6 0.3560194
#> 100    10        1.0       V10 0.3657533


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
#> 1       0                0          0 0.8757906 0.04248407 0.7978637 0.9400070
#> 2       0               20         20 0.8757906 0.04248407 0.7978637 0.9400070
#> 3       0               40         40 0.8757906 0.04248407 0.7978637 0.9400070
#> 4       0               60         60 0.8757906 0.04248407 0.7978637 0.9400070
#> 5      20                0         20 0.8617131 0.04563763 0.7783529 0.9291569
#> 6      20               20         40 0.8617131 0.04563763 0.7783529 0.9291569
#> 7      20               40         60 0.8617131 0.04563763 0.7783529 0.9291569
#> 8      20               60         80 0.8617131 0.04563763 0.7783529 0.9291569
#> 9      40                0         40 0.8478591 0.04861966 0.7595477 0.9181893
#> 10     40               20         60 0.8478591 0.04861966 0.7595477 0.9181893
#> 11     40               40         80 0.8478591 0.04861966 0.7595477 0.9181893
#> 12     40               60        100 0.8478591 0.04861966 0.7595477 0.9181893
#> 13     60                0         60 0.8342249 0.05144079 0.7413774 0.9071504
#> 14     60               20         80 0.8342249 0.05144079 0.7413774 0.9071504
#> 15     60               40        100 0.8342249 0.05144079 0.7413774 0.9071504
#> 16     60               60        120 0.8342249 0.05144079 0.7413774 0.9071504
#> 17     80                0         80 0.8208071 0.05411129 0.7237869 0.8960758
#> 18     80               20        100 0.8208071 0.05411129 0.7237869 0.8960758
#> 19     80               40        120 0.8208071 0.05411129 0.7237869 0.8960758
#> 20     80               60        140 0.8208071 0.05411129 0.7237869 0.8960758
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11411804 0.160453035 0.5447216
#> 2  0.30574618 0.11363210 0.121956439 0.4870525
#> 3  0.26001915 0.11324439 0.091690129 0.4359169
#> 4  0.22113100 0.11236508 0.068020676 0.3906797
#> 5  0.25589195 0.10213864 0.084649277 0.4281654
#> 6  0.21762106 0.09896825 0.062540073 0.3838279
#> 7  0.18507391 0.09628602 0.045425737 0.3446521
#> 8  0.15739448 0.09364060 0.032324582 0.3100415
#> 9  0.18213629 0.08914165 0.041496121 0.3387200
#> 10 0.15489621 0.08497325 0.029343207 0.3047995
#> 11 0.13173012 0.08136909 0.020215758 0.2748113
#> 12 0.11202872 0.07800581 0.013495406 0.2482679
#> 13 0.12963921 0.07659768 0.018169997 0.2702662
#> 14 0.11025053 0.07226518 0.012012817 0.2442413
#> 15 0.09376159 0.06845978 0.007628587 0.2211594
#> 16 0.07973872 0.06495336 0.004613068 0.2006417
#> 17 0.09227334 0.06516877 0.006686911 0.2176528
#> 18 0.07847305 0.06107027 0.003983137 0.1975198
#> 19 0.06673672 0.05742420 0.002228289 0.1795659
#> 20 0.05675565 0.05406822 0.001155113 0.1635029
```
