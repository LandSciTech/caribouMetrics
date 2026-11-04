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
#> 1         0.1 0.3838513 0.03818345 0.3277876 0.4485008
#> 2         0.2 0.3832760 0.03815823 0.3272520 0.4478837
#> 3         0.3 0.3827015 0.03813303 0.3267173 0.4472673
#> 4         0.4 0.3821279 0.03810786 0.3261834 0.4466519
#> 5         0.5 0.3815551 0.03808270 0.3256504 0.4460372
#> 6         0.6 0.3809832 0.03805756 0.3251183 0.4454235
#> 7         0.7 0.3804122 0.03803245 0.3245871 0.4448105
#> 8         0.8 0.3798420 0.03800735 0.3240567 0.4441984
#> 9         0.9 0.3792726 0.03798228 0.3235272 0.4435872
#> 10        1.0 0.3787041 0.03795722 0.3229986 0.4429767

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3954915
#> 2       1        0.1        V5 0.3707448
#> 3       1        0.1        V9 0.3549671
#> 4       1        0.1        V3 0.3871387
#> 5       1        0.1        V4 0.3860854
#> 6       1        0.1        V8 0.4024568
#> 7       1        0.1        V2 0.4618685
#> 8       1        0.1        V6 0.3206759
#> 9       1        0.1        V7 0.3554291
#> 10      1        0.1       V10 0.3522836
#> 11      2        0.2        V7 0.3548169
#> 12      2        0.2        V8 0.4018959
#> 13      2        0.2        V9 0.3544153
#> 14      2        0.2       V10 0.3517898
#> 15      2        0.2        V1 0.3949098
#> 16      2        0.2        V5 0.3701746
#> 17      2        0.2        V2 0.4612349
#> 18      2        0.2        V3 0.3865504
#> 19      2        0.2        V4 0.3854471
#> 20      2        0.2        V6 0.3201281
#> 21      3        0.3        V4 0.3848098
#> 22      3        0.3        V5 0.3696053
#> 23      3        0.3        V3 0.3859630
#> 24      3        0.3        V7 0.3542057
#> 25      3        0.3        V8 0.4013359
#> 26      3        0.3        V9 0.3538644
#> 27      3        0.3        V6 0.3195813
#> 28      3        0.3       V10 0.3512967
#> 29      3        0.3        V1 0.3943290
#> 30      3        0.3        V2 0.4606023
#> 31      4        0.4        V1 0.3937490
#> 32      4        0.4        V9 0.3533143
#> 33      4        0.4        V3 0.3853765
#> 34      4        0.4        V4 0.3841736
#> 35      4        0.4        V5 0.3690369
#> 36      4        0.4        V2 0.4599705
#> 37      4        0.4        V6 0.3190354
#> 38      4        0.4        V7 0.3535956
#> 39      4        0.4        V8 0.4007766
#> 40      4        0.4       V10 0.3508043
#> 41      5        0.5        V8 0.4002181
#> 42      5        0.5        V9 0.3527651
#> 43      5        0.5       V10 0.3503126
#> 44      5        0.5        V1 0.3931699
#> 45      5        0.5        V5 0.3684693
#> 46      5        0.5        V2 0.4593396
#> 47      5        0.5        V3 0.3847909
#> 48      5        0.5        V4 0.3835385
#> 49      5        0.5        V6 0.3184904
#> 50      5        0.5        V7 0.3529866
#> 51      6        0.6        V4 0.3829044
#> 52      6        0.6        V5 0.3679027
#> 53      6        0.6        V7 0.3523786
#> 54      6        0.6        V8 0.3996604
#> 55      6        0.6        V9 0.3522168
#> 56      6        0.6        V6 0.3179464
#> 57      6        0.6       V10 0.3498216
#> 58      6        0.6        V1 0.3925916
#> 59      6        0.6        V2 0.4587095
#> 60      6        0.6        V3 0.3842062
#> 61      7        0.7        V1 0.3920141
#> 62      7        0.7        V3 0.3836224
#> 63      7        0.7        V4 0.3822714
#> 64      7        0.7        V5 0.3673369
#> 65      7        0.7        V9 0.3516693
#> 66      7        0.7        V6 0.3174033
#> 67      7        0.7        V7 0.3517716
#> 68      7        0.7        V8 0.3991034
#> 69      7        0.7        V2 0.4580803
#> 70      7        0.7       V10 0.3493313
#> 71      8        0.8        V9 0.3511226
#> 72      8        0.8       V10 0.3488417
#> 73      8        0.8        V1 0.3914376
#> 74      8        0.8        V5 0.3667720
#> 75      8        0.8        V2 0.4574520
#> 76      8        0.8        V3 0.3830394
#> 77      8        0.8        V4 0.3816394
#> 78      8        0.8        V8 0.3985473
#> 79      8        0.8        V6 0.3168611
#> 80      8        0.8        V7 0.3511657
#> 81      9        0.9        V5 0.3662079
#> 82      9        0.9        V7 0.3505608
#> 83      9        0.9        V8 0.3979919
#> 84      9        0.9        V9 0.3505768
#> 85      9        0.9        V6 0.3163198
#> 86      9        0.9       V10 0.3483527
#> 87      9        0.9        V1 0.3908618
#> 88      9        0.9        V2 0.4568245
#> 89      9        0.9        V3 0.3824574
#> 90      9        0.9        V4 0.3810084
#> 91     10        1.0        V1 0.3902869
#> 92     10        1.0        V4 0.3803785
#> 93     10        1.0        V5 0.3656447
#> 94     10        1.0        V9 0.3500319
#> 95     10        1.0        V3 0.3818762
#> 96     10        1.0        V7 0.3499570
#> 97     10        1.0        V8 0.3974372
#> 98     10        1.0        V2 0.4561979
#> 99     10        1.0        V6 0.3157795
#> 100    10        1.0       V10 0.3478645

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3954915
#> 2       1        0.1        V5 0.3707448
#> 3       1        0.1        V9 0.3549671
#> 4       1        0.1        V3 0.3871387
#> 5       1        0.1        V4 0.3860854
#> 6       1        0.1        V8 0.4024568
#> 7       1        0.1        V2 0.4618685
#> 8       1        0.1        V6 0.3206759
#> 9       1        0.1        V7 0.3554291
#> 10      1        0.1       V10 0.3522836
#> 11      2        0.2        V7 0.3548169
#> 12      2        0.2        V8 0.4018959
#> 13      2        0.2        V9 0.3544153
#> 14      2        0.2       V10 0.3517898
#> 15      2        0.2        V1 0.3949098
#> 16      2        0.2        V5 0.3701746
#> 17      2        0.2        V2 0.4612349
#> 18      2        0.2        V3 0.3865504
#> 19      2        0.2        V4 0.3854471
#> 20      2        0.2        V6 0.3201281
#> 21      3        0.3        V4 0.3848098
#> 22      3        0.3        V5 0.3696053
#> 23      3        0.3        V3 0.3859630
#> 24      3        0.3        V7 0.3542057
#> 25      3        0.3        V8 0.4013359
#> 26      3        0.3        V9 0.3538644
#> 27      3        0.3        V6 0.3195813
#> 28      3        0.3       V10 0.3512967
#> 29      3        0.3        V1 0.3943290
#> 30      3        0.3        V2 0.4606023
#> 31      4        0.4        V1 0.3937490
#> 32      4        0.4        V9 0.3533143
#> 33      4        0.4        V3 0.3853765
#> 34      4        0.4        V4 0.3841736
#> 35      4        0.4        V5 0.3690369
#> 36      4        0.4        V2 0.4599705
#> 37      4        0.4        V6 0.3190354
#> 38      4        0.4        V7 0.3535956
#> 39      4        0.4        V8 0.4007766
#> 40      4        0.4       V10 0.3508043
#> 41      5        0.5        V8 0.4002181
#> 42      5        0.5        V9 0.3527651
#> 43      5        0.5       V10 0.3503126
#> 44      5        0.5        V1 0.3931699
#> 45      5        0.5        V5 0.3684693
#> 46      5        0.5        V2 0.4593396
#> 47      5        0.5        V3 0.3847909
#> 48      5        0.5        V4 0.3835385
#> 49      5        0.5        V6 0.3184904
#> 50      5        0.5        V7 0.3529866
#> 51      6        0.6        V4 0.3829044
#> 52      6        0.6        V5 0.3679027
#> 53      6        0.6        V7 0.3523786
#> 54      6        0.6        V8 0.3996604
#> 55      6        0.6        V9 0.3522168
#> 56      6        0.6        V6 0.3179464
#> 57      6        0.6       V10 0.3498216
#> 58      6        0.6        V1 0.3925916
#> 59      6        0.6        V2 0.4587095
#> 60      6        0.6        V3 0.3842062
#> 61      7        0.7        V1 0.3920141
#> 62      7        0.7        V3 0.3836224
#> 63      7        0.7        V4 0.3822714
#> 64      7        0.7        V5 0.3673369
#> 65      7        0.7        V9 0.3516693
#> 66      7        0.7        V6 0.3174033
#> 67      7        0.7        V7 0.3517716
#> 68      7        0.7        V8 0.3991034
#> 69      7        0.7        V2 0.4580803
#> 70      7        0.7       V10 0.3493313
#> 71      8        0.8        V9 0.3511226
#> 72      8        0.8       V10 0.3488417
#> 73      8        0.8        V1 0.3914376
#> 74      8        0.8        V5 0.3667720
#> 75      8        0.8        V2 0.4574520
#> 76      8        0.8        V3 0.3830394
#> 77      8        0.8        V4 0.3816394
#> 78      8        0.8        V8 0.3985473
#> 79      8        0.8        V6 0.3168611
#> 80      8        0.8        V7 0.3511657
#> 81      9        0.9        V5 0.3662079
#> 82      9        0.9        V7 0.3505608
#> 83      9        0.9        V8 0.3979919
#> 84      9        0.9        V9 0.3505768
#> 85      9        0.9        V6 0.3163198
#> 86      9        0.9       V10 0.3483527
#> 87      9        0.9        V1 0.3908618
#> 88      9        0.9        V2 0.4568245
#> 89      9        0.9        V3 0.3824574
#> 90      9        0.9        V4 0.3810084
#> 91     10        1.0        V1 0.3902869
#> 92     10        1.0        V4 0.3803785
#> 93     10        1.0        V5 0.3656447
#> 94     10        1.0        V9 0.3500319
#> 95     10        1.0        V3 0.3818762
#> 96     10        1.0        V7 0.3499570
#> 97     10        1.0        V8 0.3974372
#> 98     10        1.0        V2 0.4561979
#> 99     10        1.0        V6 0.3157795
#> 100    10        1.0       V10 0.3478645


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
#> 1       0                0          0 0.8757906 0.04971413 0.7917709 0.9495354
#> 2       0               20         20 0.8757906 0.04971413 0.7917709 0.9495354
#> 3       0               40         40 0.8757906 0.04971413 0.7917709 0.9495354
#> 4       0               60         60 0.8757906 0.04971413 0.7917709 0.9495354
#> 5      20                0         20 0.8617131 0.05218458 0.7733166 0.9412604
#> 6      20               20         40 0.8617131 0.05218458 0.7733166 0.9412604
#> 7      20               40         60 0.8617131 0.05218458 0.7733166 0.9412604
#> 8      20               60         80 0.8617131 0.05218458 0.7733166 0.9412604
#> 9      40                0         40 0.8478591 0.05447778 0.7554887 0.9328579
#> 10     40               20         60 0.8478591 0.05447778 0.7554887 0.9328579
#> 11     40               40         80 0.8478591 0.05447778 0.7554887 0.9328579
#> 12     40               60        100 0.8478591 0.05447778 0.7554887 0.9328579
#> 13     60                0         60 0.8342249 0.05661438 0.7382291 0.9243607
#> 14     60               20         80 0.8342249 0.05661438 0.7382291 0.9243607
#> 15     60               40        100 0.8342249 0.05661438 0.7382291 0.9243607
#> 16     60               60        120 0.8342249 0.05661438 0.7382291 0.9243607
#> 17     80                0         80 0.8208071 0.05861088 0.7214921 0.9157950
#> 18     80               20        100 0.8208071 0.05861088 0.7214921 0.9157950
#> 19     80               40        120 0.8208071 0.05861088 0.7214921 0.9157950
#> 20     80               60        140 0.8208071 0.05861088 0.7214921 0.9157950
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12511307 0.174321612 0.6103338
#> 2  0.30574618 0.12333156 0.128071769 0.5576788
#> 3  0.26001915 0.12002599 0.092816852 0.5096932
#> 4  0.22113100 0.11567779 0.066118335 0.4661020
#> 5  0.25589195 0.11984382 0.093839641 0.5045444
#> 6  0.21762106 0.11488213 0.066889406 0.4614308
#> 7  0.18507391 0.10949677 0.046677435 0.4223505
#> 8  0.15739448 0.10387244 0.031729581 0.3869652
#> 9  0.18213629 0.10971563 0.047257561 0.4181662
#> 10 0.15489621 0.10354294 0.032154903 0.3831779
#> 11 0.13173012 0.09748043 0.021186229 0.3515140
#> 12 0.11202872 0.09157846 0.013414645 0.3228575
#> 13 0.12963921 0.09764698 0.021494863 0.3481252
#> 14 0.11025053 0.09132131 0.013630011 0.3197901
#> 15 0.09376159 0.08534699 0.008229649 0.2941330
#> 16 0.07973872 0.07971307 0.004675549 0.2708819
#> 17 0.09227334 0.08534649 0.008376426 0.2913854
#> 18 0.07847305 0.07940521 0.004769627 0.2683904
#> 19 0.06673672 0.07390327 0.002519347 0.2475247
#> 20 0.05675565 0.06880229 0.001212262 0.2285630
```
