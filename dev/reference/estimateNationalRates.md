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
#> 1         0.1 0.3838513 0.02151957 0.3563248 0.4117289
#> 2         0.2 0.3832760 0.02147204 0.3557733 0.4110540
#> 3         0.3 0.3827015 0.02142472 0.3552226 0.4103803
#> 4         0.4 0.3821279 0.02137761 0.3546729 0.4097076
#> 5         0.5 0.3815551 0.02133070 0.3541239 0.4090361
#> 6         0.6 0.3809832 0.02128401 0.3535758 0.4083657
#> 7         0.7 0.3804122 0.02123753 0.3530286 0.4076964
#> 8         0.8 0.3798420 0.02119125 0.3524822 0.4070282
#> 9         0.9 0.3792726 0.02114518 0.3519367 0.4063611
#> 10        1.0 0.3787041 0.02109932 0.3513920 0.4056951

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3998696
#> 2       1        0.1        V5 0.3836235
#> 3       1        0.1        V9 0.3561557
#> 4       1        0.1        V3 0.4076408
#> 5       1        0.1        V4 0.3569070
#> 6       1        0.1        V8 0.3810604
#> 7       1        0.1        V2 0.4107237
#> 8       1        0.1        V6 0.3870568
#> 9       1        0.1        V7 0.4120207
#> 10      1        0.1       V10 0.3647738
#> 11      2        0.2        V7 0.4113699
#> 12      2        0.2        V8 0.3804268
#> 13      2        0.2        V9 0.3555996
#> 14      2        0.2       V10 0.3642169
#> 15      2        0.2        V1 0.3992441
#> 16      2        0.2        V5 0.3829353
#> 17      2        0.2        V2 0.4099660
#> 18      2        0.2        V3 0.4070476
#> 19      2        0.2        V4 0.3563713
#> 20      2        0.2        V6 0.3863983
#> 21      3        0.3        V4 0.3558365
#> 22      3        0.3        V5 0.3822483
#> 23      3        0.3        V3 0.4064553
#> 24      3        0.3        V7 0.4107201
#> 25      3        0.3        V8 0.3797943
#> 26      3        0.3        V9 0.3550444
#> 27      3        0.3        V6 0.3857410
#> 28      3        0.3       V10 0.3636609
#> 29      3        0.3        V1 0.3986195
#> 30      3        0.3        V2 0.4092097
#> 31      4        0.4        V1 0.3979959
#> 32      4        0.4        V9 0.3544901
#> 33      4        0.4        V3 0.4058638
#> 34      4        0.4        V4 0.3553025
#> 35      4        0.4        V5 0.3815626
#> 36      4        0.4        V2 0.4084548
#> 37      4        0.4        V6 0.3850848
#> 38      4        0.4        V7 0.4100714
#> 39      4        0.4        V8 0.3791628
#> 40      4        0.4       V10 0.3631057
#> 41      5        0.5        V8 0.3785324
#> 42      5        0.5        V9 0.3539366
#> 43      5        0.5       V10 0.3625514
#> 44      5        0.5        V1 0.3973733
#> 45      5        0.5        V5 0.3808780
#> 46      5        0.5        V2 0.4077013
#> 47      5        0.5        V3 0.4052732
#> 48      5        0.5        V4 0.3547692
#> 49      5        0.5        V6 0.3844298
#> 50      5        0.5        V7 0.4094237
#> 51      6        0.6        V4 0.3542368
#> 52      6        0.6        V5 0.3801947
#> 53      6        0.6        V7 0.4087770
#> 54      6        0.6        V8 0.3779030
#> 55      6        0.6        V9 0.3533840
#> 56      6        0.6        V6 0.3837758
#> 57      6        0.6       V10 0.3619979
#> 58      6        0.6        V1 0.3967517
#> 59      6        0.6        V2 0.4069491
#> 60      6        0.6        V3 0.4046835
#> 61      7        0.7        V1 0.3961310
#> 62      7        0.7        V3 0.4040946
#> 63      7        0.7        V4 0.3537052
#> 64      7        0.7        V5 0.3795127
#> 65      7        0.7        V9 0.3528322
#> 66      7        0.7        V6 0.3831229
#> 67      7        0.7        V7 0.4081313
#> 68      7        0.7        V8 0.3772747
#> 69      7        0.7        V2 0.4061984
#> 70      7        0.7       V10 0.3614453
#> 71      8        0.8        V9 0.3522813
#> 72      8        0.8       V10 0.3608935
#> 73      8        0.8        V1 0.3955113
#> 74      8        0.8        V5 0.3788318
#> 75      8        0.8        V2 0.4054491
#> 76      8        0.8        V3 0.4035066
#> 77      8        0.8        V4 0.3531744
#> 78      8        0.8        V8 0.3766475
#> 79      8        0.8        V6 0.3824712
#> 80      8        0.8        V7 0.4074866
#> 81      9        0.9        V5 0.3781522
#> 82      9        0.9        V7 0.4068430
#> 83      9        0.9        V8 0.3760212
#> 84      9        0.9        V9 0.3517313
#> 85      9        0.9        V6 0.3818206
#> 86      9        0.9       V10 0.3603426
#> 87      9        0.9        V1 0.3948926
#> 88      9        0.9        V2 0.4047011
#> 89      9        0.9        V3 0.4029194
#> 90      9        0.9        V4 0.3526443
#> 91     10        1.0        V1 0.3942748
#> 92     10        1.0        V4 0.3521151
#> 93     10        1.0        V5 0.3774738
#> 94     10        1.0        V9 0.3511821
#> 95     10        1.0        V3 0.4023331
#> 96     10        1.0        V7 0.4062004
#> 97     10        1.0        V8 0.3753961
#> 98     10        1.0        V2 0.4039545
#> 99     10        1.0        V6 0.3811711
#> 100    10        1.0       V10 0.3597925

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3998696
#> 2       1        0.1        V5 0.3836235
#> 3       1        0.1        V9 0.3561557
#> 4       1        0.1        V3 0.4076408
#> 5       1        0.1        V4 0.3569070
#> 6       1        0.1        V8 0.3810604
#> 7       1        0.1        V2 0.4107237
#> 8       1        0.1        V6 0.3870568
#> 9       1        0.1        V7 0.4120207
#> 10      1        0.1       V10 0.3647738
#> 11      2        0.2        V7 0.4113699
#> 12      2        0.2        V8 0.3804268
#> 13      2        0.2        V9 0.3555996
#> 14      2        0.2       V10 0.3642169
#> 15      2        0.2        V1 0.3992441
#> 16      2        0.2        V5 0.3829353
#> 17      2        0.2        V2 0.4099660
#> 18      2        0.2        V3 0.4070476
#> 19      2        0.2        V4 0.3563713
#> 20      2        0.2        V6 0.3863983
#> 21      3        0.3        V4 0.3558365
#> 22      3        0.3        V5 0.3822483
#> 23      3        0.3        V3 0.4064553
#> 24      3        0.3        V7 0.4107201
#> 25      3        0.3        V8 0.3797943
#> 26      3        0.3        V9 0.3550444
#> 27      3        0.3        V6 0.3857410
#> 28      3        0.3       V10 0.3636609
#> 29      3        0.3        V1 0.3986195
#> 30      3        0.3        V2 0.4092097
#> 31      4        0.4        V1 0.3979959
#> 32      4        0.4        V9 0.3544901
#> 33      4        0.4        V3 0.4058638
#> 34      4        0.4        V4 0.3553025
#> 35      4        0.4        V5 0.3815626
#> 36      4        0.4        V2 0.4084548
#> 37      4        0.4        V6 0.3850848
#> 38      4        0.4        V7 0.4100714
#> 39      4        0.4        V8 0.3791628
#> 40      4        0.4       V10 0.3631057
#> 41      5        0.5        V8 0.3785324
#> 42      5        0.5        V9 0.3539366
#> 43      5        0.5       V10 0.3625514
#> 44      5        0.5        V1 0.3973733
#> 45      5        0.5        V5 0.3808780
#> 46      5        0.5        V2 0.4077013
#> 47      5        0.5        V3 0.4052732
#> 48      5        0.5        V4 0.3547692
#> 49      5        0.5        V6 0.3844298
#> 50      5        0.5        V7 0.4094237
#> 51      6        0.6        V4 0.3542368
#> 52      6        0.6        V5 0.3801947
#> 53      6        0.6        V7 0.4087770
#> 54      6        0.6        V8 0.3779030
#> 55      6        0.6        V9 0.3533840
#> 56      6        0.6        V6 0.3837758
#> 57      6        0.6       V10 0.3619979
#> 58      6        0.6        V1 0.3967517
#> 59      6        0.6        V2 0.4069491
#> 60      6        0.6        V3 0.4046835
#> 61      7        0.7        V1 0.3961310
#> 62      7        0.7        V3 0.4040946
#> 63      7        0.7        V4 0.3537052
#> 64      7        0.7        V5 0.3795127
#> 65      7        0.7        V9 0.3528322
#> 66      7        0.7        V6 0.3831229
#> 67      7        0.7        V7 0.4081313
#> 68      7        0.7        V8 0.3772747
#> 69      7        0.7        V2 0.4061984
#> 70      7        0.7       V10 0.3614453
#> 71      8        0.8        V9 0.3522813
#> 72      8        0.8       V10 0.3608935
#> 73      8        0.8        V1 0.3955113
#> 74      8        0.8        V5 0.3788318
#> 75      8        0.8        V2 0.4054491
#> 76      8        0.8        V3 0.4035066
#> 77      8        0.8        V4 0.3531744
#> 78      8        0.8        V8 0.3766475
#> 79      8        0.8        V6 0.3824712
#> 80      8        0.8        V7 0.4074866
#> 81      9        0.9        V5 0.3781522
#> 82      9        0.9        V7 0.4068430
#> 83      9        0.9        V8 0.3760212
#> 84      9        0.9        V9 0.3517313
#> 85      9        0.9        V6 0.3818206
#> 86      9        0.9       V10 0.3603426
#> 87      9        0.9        V1 0.3948926
#> 88      9        0.9        V2 0.4047011
#> 89      9        0.9        V3 0.4029194
#> 90      9        0.9        V4 0.3526443
#> 91     10        1.0        V1 0.3942748
#> 92     10        1.0        V4 0.3521151
#> 93     10        1.0        V5 0.3774738
#> 94     10        1.0        V9 0.3511821
#> 95     10        1.0        V3 0.4023331
#> 96     10        1.0        V7 0.4062004
#> 97     10        1.0        V8 0.3753961
#> 98     10        1.0        V2 0.4039545
#> 99     10        1.0        V6 0.3811711
#> 100    10        1.0       V10 0.3597925


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
#> 1       0                0          0 0.8757906 0.05273441 0.7620860 0.9451029
#> 2       0               20         20 0.8757906 0.05273441 0.7620860 0.9451029
#> 3       0               40         40 0.8757906 0.05273441 0.7620860 0.9451029
#> 4       0               60         60 0.8757906 0.05273441 0.7620860 0.9451029
#> 5      20                0         20 0.8617131 0.05516920 0.7422561 0.9331362
#> 6      20               20         40 0.8617131 0.05516920 0.7422561 0.9331362
#> 7      20               40         60 0.8617131 0.05516920 0.7422561 0.9331362
#> 8      20               60         80 0.8617131 0.05516920 0.7422561 0.9331362
#> 9      40                0         40 0.8478591 0.05739016 0.7231515 0.9209829
#> 10     40               20         60 0.8478591 0.05739016 0.7231515 0.9209829
#> 11     40               40         80 0.8478591 0.05739016 0.7231515 0.9209829
#> 12     40               60        100 0.8478591 0.05739016 0.7231515 0.9209829
#> 13     60                0         60 0.8342249 0.05942657 0.7047056 0.9087154
#> 14     60               20         80 0.8342249 0.05942657 0.7047056 0.9087154
#> 15     60               40        100 0.8342249 0.05942657 0.7047056 0.9087154
#> 16     60               60        120 0.8342249 0.05942657 0.7047056 0.9087154
#> 17     80                0         80 0.8208071 0.06130115 0.6868652 0.8963872
#> 18     80               20        100 0.8208071 0.06130115 0.6868652 0.8963872
#> 19     80               40        120 0.8208071 0.06130115 0.6868652 0.8963872
#> 20     80               60        140 0.8208071 0.06130115 0.6868652 0.8963872
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13680035 0.148713658 0.6094232
#> 2  0.30574618 0.12970516 0.114099028 0.5363775
#> 3  0.26001915 0.12333628 0.086622541 0.4724374
#> 4  0.22113100 0.11733529 0.064928766 0.4167394
#> 5  0.25589195 0.12422900 0.080897002 0.4936294
#> 6  0.21762106 0.11594900 0.060428296 0.4351794
#> 7  0.18507391 0.10860858 0.044422774 0.3843540
#> 8  0.15739448 0.10190735 0.032038458 0.3402142
#> 9  0.18213629 0.10976639 0.041127564 0.4011748
#> 10 0.15489621 0.10156453 0.029509492 0.3548201
#> 11 0.13173012 0.09428599 0.020674080 0.3145666
#> 12 0.11202872 0.08769682 0.014075771 0.2795851
#> 13 0.12963921 0.09534770 0.018893663 0.3278880
#> 14 0.11025053 0.08776977 0.012764590 0.2911668
#> 15 0.09376159 0.08102098 0.008319753 0.2592233
#> 16 0.07973872 0.07492184 0.005193530 0.2313743
#> 17 0.09227334 0.08190334 0.007456256 0.2698039
#> 18 0.07847305 0.07515057 0.004600416 0.2406069
#> 19 0.06673672 0.06911119 0.002688241 0.2151047
#> 20 0.05675565 0.06364663 0.001471469 0.1927515
```
