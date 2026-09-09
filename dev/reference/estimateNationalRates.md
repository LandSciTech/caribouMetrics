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
#> 1         0.1 0.3838513 0.02003159 0.3589187 0.4152430
#> 2         0.2 0.3832760 0.02003211 0.3583177 0.4146614
#> 3         0.3 0.3827015 0.02003264 0.3577177 0.4140807
#> 4         0.4 0.3821279 0.02003318 0.3571187 0.4135007
#> 5         0.5 0.3815551 0.02003373 0.3565207 0.4129216
#> 6         0.6 0.3809832 0.02003428 0.3559237 0.4123433
#> 7         0.7 0.3804122 0.02003483 0.3553277 0.4117658
#> 8         0.8 0.3798420 0.02003539 0.3547327 0.4111891
#> 9         0.9 0.3792726 0.02003596 0.3541387 0.4106133
#> 10        1.0 0.3787041 0.02003653 0.3535457 0.4100382

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3883621
#> 2       1        0.1        V5 0.3789057
#> 3       1        0.1        V9 0.4179214
#> 4       1        0.1        V3 0.4060172
#> 5       1        0.1        V4 0.3946402
#> 6       1        0.1        V8 0.3899764
#> 7       1        0.1        V2 0.3622496
#> 8       1        0.1        V6 0.3617010
#> 9       1        0.1        V7 0.3716539
#> 10      1        0.1       V10 0.3581110
#> 11      2        0.2        V7 0.3710941
#> 12      2        0.2        V8 0.3894793
#> 13      2        0.2        V9 0.4173446
#> 14      2        0.2       V10 0.3575056
#> 15      2        0.2        V1 0.3878009
#> 16      2        0.2        V5 0.3783786
#> 17      2        0.2        V2 0.3617061
#> 18      2        0.2        V3 0.4054193
#> 19      2        0.2        V4 0.3940728
#> 20      2        0.2        V6 0.3611150
#> 21      3        0.3        V4 0.3935063
#> 22      3        0.3        V5 0.3778522
#> 23      3        0.3        V3 0.4048223
#> 24      3        0.3        V7 0.3705352
#> 25      3        0.3        V8 0.3889829
#> 26      3        0.3        V9 0.4167686
#> 27      3        0.3        V6 0.3605299
#> 28      3        0.3       V10 0.3569012
#> 29      3        0.3        V1 0.3872404
#> 30      3        0.3        V2 0.3611634
#> 31      4        0.4        V1 0.3866809
#> 32      4        0.4        V9 0.4161934
#> 33      4        0.4        V3 0.4042261
#> 34      4        0.4        V4 0.3929405
#> 35      4        0.4        V5 0.3773266
#> 36      4        0.4        V2 0.3606215
#> 37      4        0.4        V6 0.3599459
#> 38      4        0.4        V7 0.3699771
#> 39      4        0.4        V8 0.3884871
#> 40      4        0.4       V10 0.3562979
#> 41      5        0.5        V8 0.3879919
#> 42      5        0.5        V9 0.4156189
#> 43      5        0.5       V10 0.3556956
#> 44      5        0.5        V1 0.3861221
#> 45      5        0.5        V5 0.3768017
#> 46      5        0.5        V2 0.3600804
#> 47      5        0.5        V3 0.4036309
#> 48      5        0.5        V4 0.3923756
#> 49      5        0.5        V6 0.3593627
#> 50      5        0.5        V7 0.3694198
#> 51      6        0.6        V4 0.3918115
#> 52      6        0.6        V5 0.3762775
#> 53      6        0.6        V7 0.3688634
#> 54      6        0.6        V8 0.3874973
#> 55      6        0.6        V9 0.4150453
#> 56      6        0.6        V6 0.3587805
#> 57      6        0.6       V10 0.3550943
#> 58      6        0.6        V1 0.3855641
#> 59      6        0.6        V2 0.3595402
#> 60      6        0.6        V3 0.4030365
#> 61      7        0.7        V1 0.3850069
#> 62      7        0.7        V3 0.4024430
#> 63      7        0.7        V4 0.3912482
#> 64      7        0.7        V5 0.3757540
#> 65      7        0.7        V9 0.4144725
#> 66      7        0.7        V6 0.3581992
#> 67      7        0.7        V7 0.3683079
#> 68      7        0.7        V8 0.3870034
#> 69      7        0.7        V2 0.3590007
#> 70      7        0.7       V10 0.3544941
#> 71      8        0.8        V9 0.4139004
#> 72      8        0.8       V10 0.3538948
#> 73      8        0.8        V1 0.3844506
#> 74      8        0.8        V5 0.3752313
#> 75      8        0.8        V2 0.3584621
#> 76      8        0.8        V3 0.4018504
#> 77      8        0.8        V4 0.3906857
#> 78      8        0.8        V8 0.3865101
#> 79      8        0.8        V6 0.3576189
#> 80      8        0.8        V7 0.3677531
#> 81      9        0.9        V5 0.3747093
#> 82      9        0.9        V7 0.3671992
#> 83      9        0.9        V8 0.3860175
#> 84      9        0.9        V9 0.4133291
#> 85      9        0.9        V6 0.3570395
#> 86      9        0.9       V10 0.3532966
#> 87      9        0.9        V1 0.3838950
#> 88      9        0.9        V2 0.3579242
#> 89      9        0.9        V3 0.4012586
#> 90      9        0.9        V4 0.3901240
#> 91     10        1.0        V1 0.3833402
#> 92     10        1.0        V4 0.3895631
#> 93     10        1.0        V5 0.3741880
#> 94     10        1.0        V9 0.4127587
#> 95     10        1.0        V3 0.4006677
#> 96     10        1.0        V7 0.3666462
#> 97     10        1.0        V8 0.3855255
#> 98     10        1.0        V2 0.3573872
#> 99     10        1.0        V6 0.3564611
#> 100    10        1.0       V10 0.3526993

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3883621
#> 2       1        0.1        V5 0.3789057
#> 3       1        0.1        V9 0.4179214
#> 4       1        0.1        V3 0.4060172
#> 5       1        0.1        V4 0.3946402
#> 6       1        0.1        V8 0.3899764
#> 7       1        0.1        V2 0.3622496
#> 8       1        0.1        V6 0.3617010
#> 9       1        0.1        V7 0.3716539
#> 10      1        0.1       V10 0.3581110
#> 11      2        0.2        V7 0.3710941
#> 12      2        0.2        V8 0.3894793
#> 13      2        0.2        V9 0.4173446
#> 14      2        0.2       V10 0.3575056
#> 15      2        0.2        V1 0.3878009
#> 16      2        0.2        V5 0.3783786
#> 17      2        0.2        V2 0.3617061
#> 18      2        0.2        V3 0.4054193
#> 19      2        0.2        V4 0.3940728
#> 20      2        0.2        V6 0.3611150
#> 21      3        0.3        V4 0.3935063
#> 22      3        0.3        V5 0.3778522
#> 23      3        0.3        V3 0.4048223
#> 24      3        0.3        V7 0.3705352
#> 25      3        0.3        V8 0.3889829
#> 26      3        0.3        V9 0.4167686
#> 27      3        0.3        V6 0.3605299
#> 28      3        0.3       V10 0.3569012
#> 29      3        0.3        V1 0.3872404
#> 30      3        0.3        V2 0.3611634
#> 31      4        0.4        V1 0.3866809
#> 32      4        0.4        V9 0.4161934
#> 33      4        0.4        V3 0.4042261
#> 34      4        0.4        V4 0.3929405
#> 35      4        0.4        V5 0.3773266
#> 36      4        0.4        V2 0.3606215
#> 37      4        0.4        V6 0.3599459
#> 38      4        0.4        V7 0.3699771
#> 39      4        0.4        V8 0.3884871
#> 40      4        0.4       V10 0.3562979
#> 41      5        0.5        V8 0.3879919
#> 42      5        0.5        V9 0.4156189
#> 43      5        0.5       V10 0.3556956
#> 44      5        0.5        V1 0.3861221
#> 45      5        0.5        V5 0.3768017
#> 46      5        0.5        V2 0.3600804
#> 47      5        0.5        V3 0.4036309
#> 48      5        0.5        V4 0.3923756
#> 49      5        0.5        V6 0.3593627
#> 50      5        0.5        V7 0.3694198
#> 51      6        0.6        V4 0.3918115
#> 52      6        0.6        V5 0.3762775
#> 53      6        0.6        V7 0.3688634
#> 54      6        0.6        V8 0.3874973
#> 55      6        0.6        V9 0.4150453
#> 56      6        0.6        V6 0.3587805
#> 57      6        0.6       V10 0.3550943
#> 58      6        0.6        V1 0.3855641
#> 59      6        0.6        V2 0.3595402
#> 60      6        0.6        V3 0.4030365
#> 61      7        0.7        V1 0.3850069
#> 62      7        0.7        V3 0.4024430
#> 63      7        0.7        V4 0.3912482
#> 64      7        0.7        V5 0.3757540
#> 65      7        0.7        V9 0.4144725
#> 66      7        0.7        V6 0.3581992
#> 67      7        0.7        V7 0.3683079
#> 68      7        0.7        V8 0.3870034
#> 69      7        0.7        V2 0.3590007
#> 70      7        0.7       V10 0.3544941
#> 71      8        0.8        V9 0.4139004
#> 72      8        0.8       V10 0.3538948
#> 73      8        0.8        V1 0.3844506
#> 74      8        0.8        V5 0.3752313
#> 75      8        0.8        V2 0.3584621
#> 76      8        0.8        V3 0.4018504
#> 77      8        0.8        V4 0.3906857
#> 78      8        0.8        V8 0.3865101
#> 79      8        0.8        V6 0.3576189
#> 80      8        0.8        V7 0.3677531
#> 81      9        0.9        V5 0.3747093
#> 82      9        0.9        V7 0.3671992
#> 83      9        0.9        V8 0.3860175
#> 84      9        0.9        V9 0.4133291
#> 85      9        0.9        V6 0.3570395
#> 86      9        0.9       V10 0.3532966
#> 87      9        0.9        V1 0.3838950
#> 88      9        0.9        V2 0.3579242
#> 89      9        0.9        V3 0.4012586
#> 90      9        0.9        V4 0.3901240
#> 91     10        1.0        V1 0.3833402
#> 92     10        1.0        V4 0.3895631
#> 93     10        1.0        V5 0.3741880
#> 94     10        1.0        V9 0.4127587
#> 95     10        1.0        V3 0.4006677
#> 96     10        1.0        V7 0.3666462
#> 97     10        1.0        V8 0.3855255
#> 98     10        1.0        V2 0.3573872
#> 99     10        1.0        V6 0.3564611
#> 100    10        1.0       V10 0.3526993


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
#> 1       0                0          0 0.8757906 0.05119748 0.7640166 0.9408379
#> 2       0               20         20 0.8757906 0.05119748 0.7640166 0.9408379
#> 3       0               40         40 0.8757906 0.05119748 0.7640166 0.9408379
#> 4       0               60         60 0.8757906 0.05119748 0.7640166 0.9408379
#> 5      20                0         20 0.8617131 0.05394881 0.7470336 0.9328925
#> 6      20               20         40 0.8617131 0.05394881 0.7470336 0.9328925
#> 7      20               40         60 0.8617131 0.05394881 0.7470336 0.9328925
#> 8      20               60         80 0.8617131 0.05394881 0.7470336 0.9328925
#> 9      40                0         40 0.8478591 0.05656293 0.7305970 0.9248689
#> 10     40               20         60 0.8478591 0.05656293 0.7305970 0.9248689
#> 11     40               40         80 0.8478591 0.05656293 0.7305970 0.9248689
#> 12     40               60        100 0.8478591 0.05656293 0.7305970 0.9248689
#> 13     60                0         60 0.8342249 0.05904733 0.7146606 0.9167883
#> 14     60               20         80 0.8342249 0.05904733 0.7146606 0.9167883
#> 15     60               40        100 0.8342249 0.05904733 0.7146606 0.9167883
#> 16     60               60        120 0.8342249 0.05904733 0.7146606 0.9167883
#> 17     80                0         80 0.8208071 0.06140924 0.6991870 0.9086685
#> 18     80               20        100 0.8208071 0.06140924 0.6991870 0.9086685
#> 19     80               40        120 0.8208071 0.06140924 0.6991870 0.9086685
#> 20     80               60        140 0.8208071 0.06140924 0.6991870 0.9086685
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12213067 0.180820067 0.5913706
#> 2  0.30574618 0.11847498 0.139238170 0.5299872
#> 3  0.26001915 0.11351001 0.106231648 0.4752682
#> 4  0.22113100 0.10765890 0.080134382 0.4266596
#> 5  0.25589195 0.11023596 0.100414731 0.4634868
#> 6  0.21762106 0.10463619 0.075550788 0.4162083
#> 7  0.18507391 0.09860974 0.056039428 0.3743010
#> 8  0.15739448 0.09235703 0.040862641 0.3371770
#> 9  0.18213629 0.09663402 0.052632844 0.3652971
#> 10 0.15489621 0.09053536 0.038230502 0.3292008
#> 11 0.13173012 0.08443098 0.027185642 0.2972150
#> 12 0.11202872 0.07841550 0.018845528 0.2688448
#> 13 0.12963921 0.08321735 0.025290362 0.2903389
#> 14 0.11025053 0.07729665 0.017430791 0.2627408
#> 15 0.09376159 0.07156683 0.011638049 0.2382156
#> 16 0.07973872 0.06606838 0.007478644 0.2163750
#> 17 0.09227334 0.07080780 0.010673052 0.2329310
#> 18 0.07847305 0.06536580 0.006799162 0.2116616
#> 19 0.06673672 0.06019624 0.004128514 0.1926585
#> 20 0.05675565 0.05531206 0.002365922 0.1756263
```
