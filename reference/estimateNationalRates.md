# Sample demographic rates

Apply the sampled coefficients to the disturbance covariates to
calculate expected recruitment and survival according to the beta
regression models estimated by Johnson et al.
(2020).`estimateNationalRates` is a wrapper around
`estimateNationalRate` to sample both survival and recruitment rates
based on the result of
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md)
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
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md).
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
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateBayesianRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummaryForApp.md)

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
#> 1         0.1 0.3838513 0.02541271 0.3435217 0.4178856
#> 2         0.2 0.3832760 0.02537358 0.3429980 0.4172900
#> 3         0.3 0.3827015 0.02533460 0.3424751 0.4166952
#> 4         0.4 0.3821279 0.02529578 0.3419530 0.4161012
#> 5         0.5 0.3815551 0.02525712 0.3414317 0.4155081
#> 6         0.6 0.3809832 0.02521862 0.3409111 0.4149159
#> 7         0.7 0.3804122 0.02518026 0.3403914 0.4143245
#> 8         0.8 0.3798420 0.02514207 0.3398725 0.4137340
#> 9         0.9 0.3792726 0.02510402 0.3393544 0.4131443
#> 10        1.0 0.3787041 0.02506613 0.3388370 0.4125554

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3622781
#> 2       1        0.1        V5 0.3380762
#> 3       1        0.1        V9 0.4139746
#> 4       1        0.1        V3 0.3954644
#> 5       1        0.1        V4 0.3658831
#> 6       1        0.1        V8 0.4190211
#> 7       1        0.1        V2 0.3727711
#> 8       1        0.1        V6 0.3956350
#> 9       1        0.1        V7 0.3667731
#> 10      1        0.1       V10 0.3946132
#> 11      2        0.2        V7 0.3662089
#> 12      2        0.2        V8 0.4184495
#> 13      2        0.2        V9 0.4132961
#> 14      2        0.2       V10 0.3940569
#> 15      2        0.2        V1 0.3617408
#> 16      2        0.2        V5 0.3375565
#> 17      2        0.2        V2 0.3722533
#> 18      2        0.2        V3 0.3947634
#> 19      2        0.2        V4 0.3653173
#> 20      2        0.2        V6 0.3950358
#> 21      3        0.3        V4 0.3647525
#> 22      3        0.3        V5 0.3370375
#> 23      3        0.3        V3 0.3940636
#> 24      3        0.3        V7 0.3656456
#> 25      3        0.3        V8 0.4178787
#> 26      3        0.3        V9 0.4126187
#> 27      3        0.3        V6 0.3944374
#> 28      3        0.3       V10 0.3935013
#> 29      3        0.3        V1 0.3612044
#> 30      3        0.3        V2 0.3717362
#> 31      4        0.4        V1 0.3606687
#> 32      4        0.4        V9 0.4119424
#> 33      4        0.4        V3 0.3933651
#> 34      4        0.4        V4 0.3641885
#> 35      4        0.4        V5 0.3365194
#> 36      4        0.4        V2 0.3712198
#> 37      4        0.4        V6 0.3938399
#> 38      4        0.4        V7 0.3650831
#> 39      4        0.4        V8 0.4173086
#> 40      4        0.4       V10 0.3929465
#> 41      5        0.5        V8 0.4167394
#> 42      5        0.5        V9 0.4112672
#> 43      5        0.5       V10 0.3923924
#> 44      5        0.5        V1 0.3601338
#> 45      5        0.5        V5 0.3360020
#> 46      5        0.5        V2 0.3707041
#> 47      5        0.5        V3 0.3926678
#> 48      5        0.5        V4 0.3636253
#> 49      5        0.5        V6 0.3932434
#> 50      5        0.5        V7 0.3645215
#> 51      6        0.6        V4 0.3630631
#> 52      6        0.6        V5 0.3354854
#> 53      6        0.6        V7 0.3639608
#> 54      6        0.6        V8 0.4161709
#> 55      6        0.6        V9 0.4105931
#> 56      6        0.6        V6 0.3926477
#> 57      6        0.6       V10 0.3918392
#> 58      6        0.6        V1 0.3595997
#> 59      6        0.6        V2 0.3701891
#> 60      6        0.6        V3 0.3919717
#> 61      7        0.7        V1 0.3590664
#> 62      7        0.7        V3 0.3912769
#> 63      7        0.7        V4 0.3625017
#> 64      7        0.7        V5 0.3349696
#> 65      7        0.7        V9 0.4099201
#> 66      7        0.7        V6 0.3920530
#> 67      7        0.7        V7 0.3634009
#> 68      7        0.7        V8 0.4156032
#> 69      7        0.7        V2 0.3696749
#> 70      7        0.7       V10 0.3912867
#> 71      8        0.8        V9 0.4092482
#> 72      8        0.8       V10 0.3907350
#> 73      8        0.8        V1 0.3585339
#> 74      8        0.8        V5 0.3344547
#> 75      8        0.8        V2 0.3691613
#> 76      8        0.8        V3 0.3905833
#> 77      8        0.8        V4 0.3619412
#> 78      8        0.8        V8 0.4150363
#> 79      8        0.8        V6 0.3914591
#> 80      8        0.8        V7 0.3628419
#> 81      9        0.9        V5 0.3339405
#> 82      9        0.9        V7 0.3622837
#> 83      9        0.9        V8 0.4144701
#> 84      9        0.9        V9 0.4085774
#> 85      9        0.9        V6 0.3908662
#> 86      9        0.9       V10 0.3901841
#> 87      9        0.9        V1 0.3580022
#> 88      9        0.9        V2 0.3686485
#> 89      9        0.9        V3 0.3898909
#> 90      9        0.9        V4 0.3613816
#> 91     10        1.0        V1 0.3574713
#> 92     10        1.0        V4 0.3608228
#> 93     10        1.0        V5 0.3334271
#> 94     10        1.0        V9 0.4079078
#> 95     10        1.0        V3 0.3891998
#> 96     10        1.0        V7 0.3617264
#> 97     10        1.0        V8 0.4139048
#> 98     10        1.0        V2 0.3681364
#> 99     10        1.0        V6 0.3902742
#> 100    10        1.0       V10 0.3896340

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3622781
#> 2       1        0.1        V5 0.3380762
#> 3       1        0.1        V9 0.4139746
#> 4       1        0.1        V3 0.3954644
#> 5       1        0.1        V4 0.3658831
#> 6       1        0.1        V8 0.4190211
#> 7       1        0.1        V2 0.3727711
#> 8       1        0.1        V6 0.3956350
#> 9       1        0.1        V7 0.3667731
#> 10      1        0.1       V10 0.3946132
#> 11      2        0.2        V7 0.3662089
#> 12      2        0.2        V8 0.4184495
#> 13      2        0.2        V9 0.4132961
#> 14      2        0.2       V10 0.3940569
#> 15      2        0.2        V1 0.3617408
#> 16      2        0.2        V5 0.3375565
#> 17      2        0.2        V2 0.3722533
#> 18      2        0.2        V3 0.3947634
#> 19      2        0.2        V4 0.3653173
#> 20      2        0.2        V6 0.3950358
#> 21      3        0.3        V4 0.3647525
#> 22      3        0.3        V5 0.3370375
#> 23      3        0.3        V3 0.3940636
#> 24      3        0.3        V7 0.3656456
#> 25      3        0.3        V8 0.4178787
#> 26      3        0.3        V9 0.4126187
#> 27      3        0.3        V6 0.3944374
#> 28      3        0.3       V10 0.3935013
#> 29      3        0.3        V1 0.3612044
#> 30      3        0.3        V2 0.3717362
#> 31      4        0.4        V1 0.3606687
#> 32      4        0.4        V9 0.4119424
#> 33      4        0.4        V3 0.3933651
#> 34      4        0.4        V4 0.3641885
#> 35      4        0.4        V5 0.3365194
#> 36      4        0.4        V2 0.3712198
#> 37      4        0.4        V6 0.3938399
#> 38      4        0.4        V7 0.3650831
#> 39      4        0.4        V8 0.4173086
#> 40      4        0.4       V10 0.3929465
#> 41      5        0.5        V8 0.4167394
#> 42      5        0.5        V9 0.4112672
#> 43      5        0.5       V10 0.3923924
#> 44      5        0.5        V1 0.3601338
#> 45      5        0.5        V5 0.3360020
#> 46      5        0.5        V2 0.3707041
#> 47      5        0.5        V3 0.3926678
#> 48      5        0.5        V4 0.3636253
#> 49      5        0.5        V6 0.3932434
#> 50      5        0.5        V7 0.3645215
#> 51      6        0.6        V4 0.3630631
#> 52      6        0.6        V5 0.3354854
#> 53      6        0.6        V7 0.3639608
#> 54      6        0.6        V8 0.4161709
#> 55      6        0.6        V9 0.4105931
#> 56      6        0.6        V6 0.3926477
#> 57      6        0.6       V10 0.3918392
#> 58      6        0.6        V1 0.3595997
#> 59      6        0.6        V2 0.3701891
#> 60      6        0.6        V3 0.3919717
#> 61      7        0.7        V1 0.3590664
#> 62      7        0.7        V3 0.3912769
#> 63      7        0.7        V4 0.3625017
#> 64      7        0.7        V5 0.3349696
#> 65      7        0.7        V9 0.4099201
#> 66      7        0.7        V6 0.3920530
#> 67      7        0.7        V7 0.3634009
#> 68      7        0.7        V8 0.4156032
#> 69      7        0.7        V2 0.3696749
#> 70      7        0.7       V10 0.3912867
#> 71      8        0.8        V9 0.4092482
#> 72      8        0.8       V10 0.3907350
#> 73      8        0.8        V1 0.3585339
#> 74      8        0.8        V5 0.3344547
#> 75      8        0.8        V2 0.3691613
#> 76      8        0.8        V3 0.3905833
#> 77      8        0.8        V4 0.3619412
#> 78      8        0.8        V8 0.4150363
#> 79      8        0.8        V6 0.3914591
#> 80      8        0.8        V7 0.3628419
#> 81      9        0.9        V5 0.3339405
#> 82      9        0.9        V7 0.3622837
#> 83      9        0.9        V8 0.4144701
#> 84      9        0.9        V9 0.4085774
#> 85      9        0.9        V6 0.3908662
#> 86      9        0.9       V10 0.3901841
#> 87      9        0.9        V1 0.3580022
#> 88      9        0.9        V2 0.3686485
#> 89      9        0.9        V3 0.3898909
#> 90      9        0.9        V4 0.3613816
#> 91     10        1.0        V1 0.3574713
#> 92     10        1.0        V4 0.3608228
#> 93     10        1.0        V5 0.3334271
#> 94     10        1.0        V9 0.4079078
#> 95     10        1.0        V3 0.3891998
#> 96     10        1.0        V7 0.3617264
#> 97     10        1.0        V8 0.4139048
#> 98     10        1.0        V2 0.3681364
#> 99     10        1.0        V6 0.3902742
#> 100    10        1.0       V10 0.3896340


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
#> 1       0                0          0 0.8757906 0.04898661 0.7798998 0.9396401
#> 2       0               20         20 0.8757906 0.04898661 0.7798998 0.9396401
#> 3       0               40         40 0.8757906 0.04898661 0.7798998 0.9396401
#> 4       0               60         60 0.8757906 0.04898661 0.7798998 0.9396401
#> 5      20                0         20 0.8617131 0.05123290 0.7602785 0.9307419
#> 6      20               20         40 0.8617131 0.05123290 0.7602785 0.9307419
#> 7      20               40         60 0.8617131 0.05123290 0.7602785 0.9307419
#> 8      20               60         80 0.8617131 0.05123290 0.7602785 0.9307419
#> 9      40                0         40 0.8478591 0.05332441 0.7414081 0.9217556
#> 10     40               20         60 0.8478591 0.05332441 0.7414081 0.9217556
#> 11     40               40         80 0.8478591 0.05332441 0.7414081 0.9217556
#> 12     40               60        100 0.8478591 0.05332441 0.7414081 0.9217556
#> 13     60                0         60 0.8342249 0.05527907 0.7232098 0.9127092
#> 14     60               20         80 0.8342249 0.05527907 0.7232098 0.9127092
#> 15     60               40        100 0.8342249 0.05527907 0.7232098 0.9127092
#> 16     60               60        120 0.8342249 0.05527907 0.7232098 0.9127092
#> 17     80                0         80 0.8208071 0.05711100 0.7056223 0.9036248
#> 18     80               20        100 0.8208071 0.05711100 0.7056223 0.9036248
#> 19     80               40        120 0.8208071 0.05711100 0.7056223 0.9036248
#> 20     80               60        140 0.8208071 0.05711100 0.7056223 0.9036248
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13526746 0.172861570 0.6278812
#> 2  0.30574618 0.12907480 0.137197194 0.5601922
#> 3  0.26001915 0.12264076 0.108116564 0.4999717
#> 4  0.22113100 0.11601523 0.084475654 0.4466607
#> 5  0.25589195 0.12210615 0.098824969 0.5160747
#> 6  0.21762106 0.11484522 0.076943789 0.4608977
#> 7  0.18507391 0.10780835 0.059268573 0.4121564
#> 8  0.15739448 0.10096895 0.045082184 0.3691711
#> 9  0.18213629 0.10817949 0.053668631 0.4251661
#> 10 0.15489621 0.10085708 0.040611925 0.3806405
#> 11 0.13173012 0.09394609 0.030253415 0.3413913
#> 12 0.11202872 0.08739498 0.022125316 0.3067887
#> 13 0.12963921 0.09465696 0.027020868 0.3518636
#> 14 0.11025053 0.08775952 0.019611890 0.3160232
#> 15 0.09376159 0.08132685 0.013909885 0.2844079
#> 16 0.07973872 0.07530276 0.009600641 0.2564814
#> 17 0.09227334 0.08214321 0.012174907 0.2928479
#> 18 0.07847305 0.07587922 0.008309186 0.2639414
#> 19 0.06673672 0.07006982 0.005479729 0.2383729
#> 20 0.05675565 0.06466302 0.003469873 0.2157028
```
