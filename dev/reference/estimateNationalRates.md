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
#> 1         0.1 0.3838513 0.01848770 0.3512042 0.4053856
#> 2         0.2 0.3832760 0.01845546 0.3507187 0.4047827
#> 3         0.3 0.3827015 0.01842350 0.3502338 0.4041806
#> 4         0.4 0.3821279 0.01839181 0.3497497 0.4035794
#> 5         0.5 0.3815551 0.01836040 0.3492662 0.4029792
#> 6         0.6 0.3809832 0.01832926 0.3487833 0.4023798
#> 7         0.7 0.3804122 0.01829840 0.3483012 0.4017813
#> 8         0.8 0.3798420 0.01826780 0.3478197 0.4011837
#> 9         0.9 0.3792726 0.01823748 0.3473390 0.4005870
#> 10        1.0 0.3787041 0.01820743 0.3468588 0.3999912

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3673344
#> 2       1        0.1        V5 0.3604310
#> 3       1        0.1        V9 0.3538942
#> 4       1        0.1        V3 0.4083409
#> 5       1        0.1        V4 0.3504233
#> 6       1        0.1        V8 0.3720030
#> 7       1        0.1        V2 0.3734789
#> 8       1        0.1        V6 0.3816146
#> 9       1        0.1        V7 0.3884478
#> 10      1        0.1       V10 0.3952064
#> 11      2        0.2        V7 0.3878427
#> 12      2        0.2        V8 0.3715132
#> 13      2        0.2        V9 0.3532968
#> 14      2        0.2       V10 0.3945537
#> 15      2        0.2        V1 0.3667797
#> 16      2        0.2        V5 0.3598969
#> 17      2        0.2        V2 0.3727963
#> 18      2        0.2        V3 0.4077524
#> 19      2        0.2        V4 0.3499702
#> 20      2        0.2        V6 0.3811030
#> 21      3        0.3        V4 0.3495177
#> 22      3        0.3        V5 0.3593637
#> 23      3        0.3        V3 0.4071647
#> 24      3        0.3        V7 0.3872387
#> 25      3        0.3        V8 0.3710239
#> 26      3        0.3        V9 0.3527004
#> 27      3        0.3        V6 0.3805920
#> 28      3        0.3       V10 0.3939020
#> 29      3        0.3        V1 0.3662258
#> 30      3        0.3        V2 0.3721149
#> 31      4        0.4        V1 0.3656728
#> 32      4        0.4        V9 0.3521050
#> 33      4        0.4        V3 0.4065779
#> 34      4        0.4        V4 0.3490658
#> 35      4        0.4        V5 0.3588312
#> 36      4        0.4        V2 0.3714347
#> 37      4        0.4        V6 0.3800818
#> 38      4        0.4        V7 0.3866355
#> 39      4        0.4        V8 0.3705353
#> 40      4        0.4       V10 0.3932514
#> 41      5        0.5        V8 0.3700474
#> 42      5        0.5        V9 0.3515106
#> 43      5        0.5       V10 0.3926019
#> 44      5        0.5        V1 0.3651206
#> 45      5        0.5        V5 0.3582996
#> 46      5        0.5        V2 0.3707558
#> 47      5        0.5        V3 0.4059919
#> 48      5        0.5        V4 0.3486145
#> 49      5        0.5        V6 0.3795723
#> 50      5        0.5        V7 0.3860333
#> 51      6        0.6        V4 0.3481638
#> 52      6        0.6        V5 0.3577687
#> 53      6        0.6        V7 0.3854320
#> 54      6        0.6        V8 0.3695601
#> 55      6        0.6        V9 0.3509172
#> 56      6        0.6        V6 0.3790634
#> 57      6        0.6       V10 0.3919535
#> 58      6        0.6        V1 0.3645693
#> 59      6        0.6        V2 0.3700781
#> 60      6        0.6        V3 0.4054068
#> 61      7        0.7        V1 0.3640188
#> 62      7        0.7        V3 0.4048225
#> 63      7        0.7        V4 0.3477137
#> 64      7        0.7        V5 0.3572386
#> 65      7        0.7        V9 0.3503248
#> 66      7        0.7        V6 0.3785552
#> 67      7        0.7        V7 0.3848317
#> 68      7        0.7        V8 0.3690734
#> 69      7        0.7        V2 0.3694017
#> 70      7        0.7       V10 0.3913061
#> 71      8        0.8        V9 0.3497334
#> 72      8        0.8       V10 0.3906598
#> 73      8        0.8        V1 0.3634691
#> 74      8        0.8        V5 0.3567093
#> 75      8        0.8        V2 0.3687265
#> 76      8        0.8        V3 0.4042390
#> 77      8        0.8        V4 0.3472642
#> 78      8        0.8        V8 0.3685874
#> 79      8        0.8        V6 0.3780477
#> 80      8        0.8        V7 0.3842323
#> 81      9        0.9        V5 0.3561807
#> 82      9        0.9        V7 0.3836338
#> 83      9        0.9        V8 0.3681020
#> 84      9        0.9        V9 0.3491430
#> 85      9        0.9        V6 0.3775409
#> 86      9        0.9       V10 0.3900146
#> 87      9        0.9        V1 0.3629202
#> 88      9        0.9        V2 0.3680525
#> 89      9        0.9        V3 0.4036564
#> 90      9        0.9        V4 0.3468152
#> 91     10        1.0        V1 0.3623722
#> 92     10        1.0        V4 0.3463668
#> 93     10        1.0        V5 0.3556530
#> 94     10        1.0        V9 0.3485537
#> 95     10        1.0        V3 0.4030747
#> 96     10        1.0        V7 0.3830363
#> 97     10        1.0        V8 0.3676172
#> 98     10        1.0        V2 0.3673798
#> 99     10        1.0        V6 0.3770347
#> 100    10        1.0       V10 0.3893704

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3673344
#> 2       1        0.1        V5 0.3604310
#> 3       1        0.1        V9 0.3538942
#> 4       1        0.1        V3 0.4083409
#> 5       1        0.1        V4 0.3504233
#> 6       1        0.1        V8 0.3720030
#> 7       1        0.1        V2 0.3734789
#> 8       1        0.1        V6 0.3816146
#> 9       1        0.1        V7 0.3884478
#> 10      1        0.1       V10 0.3952064
#> 11      2        0.2        V7 0.3878427
#> 12      2        0.2        V8 0.3715132
#> 13      2        0.2        V9 0.3532968
#> 14      2        0.2       V10 0.3945537
#> 15      2        0.2        V1 0.3667797
#> 16      2        0.2        V5 0.3598969
#> 17      2        0.2        V2 0.3727963
#> 18      2        0.2        V3 0.4077524
#> 19      2        0.2        V4 0.3499702
#> 20      2        0.2        V6 0.3811030
#> 21      3        0.3        V4 0.3495177
#> 22      3        0.3        V5 0.3593637
#> 23      3        0.3        V3 0.4071647
#> 24      3        0.3        V7 0.3872387
#> 25      3        0.3        V8 0.3710239
#> 26      3        0.3        V9 0.3527004
#> 27      3        0.3        V6 0.3805920
#> 28      3        0.3       V10 0.3939020
#> 29      3        0.3        V1 0.3662258
#> 30      3        0.3        V2 0.3721149
#> 31      4        0.4        V1 0.3656728
#> 32      4        0.4        V9 0.3521050
#> 33      4        0.4        V3 0.4065779
#> 34      4        0.4        V4 0.3490658
#> 35      4        0.4        V5 0.3588312
#> 36      4        0.4        V2 0.3714347
#> 37      4        0.4        V6 0.3800818
#> 38      4        0.4        V7 0.3866355
#> 39      4        0.4        V8 0.3705353
#> 40      4        0.4       V10 0.3932514
#> 41      5        0.5        V8 0.3700474
#> 42      5        0.5        V9 0.3515106
#> 43      5        0.5       V10 0.3926019
#> 44      5        0.5        V1 0.3651206
#> 45      5        0.5        V5 0.3582996
#> 46      5        0.5        V2 0.3707558
#> 47      5        0.5        V3 0.4059919
#> 48      5        0.5        V4 0.3486145
#> 49      5        0.5        V6 0.3795723
#> 50      5        0.5        V7 0.3860333
#> 51      6        0.6        V4 0.3481638
#> 52      6        0.6        V5 0.3577687
#> 53      6        0.6        V7 0.3854320
#> 54      6        0.6        V8 0.3695601
#> 55      6        0.6        V9 0.3509172
#> 56      6        0.6        V6 0.3790634
#> 57      6        0.6       V10 0.3919535
#> 58      6        0.6        V1 0.3645693
#> 59      6        0.6        V2 0.3700781
#> 60      6        0.6        V3 0.4054068
#> 61      7        0.7        V1 0.3640188
#> 62      7        0.7        V3 0.4048225
#> 63      7        0.7        V4 0.3477137
#> 64      7        0.7        V5 0.3572386
#> 65      7        0.7        V9 0.3503248
#> 66      7        0.7        V6 0.3785552
#> 67      7        0.7        V7 0.3848317
#> 68      7        0.7        V8 0.3690734
#> 69      7        0.7        V2 0.3694017
#> 70      7        0.7       V10 0.3913061
#> 71      8        0.8        V9 0.3497334
#> 72      8        0.8       V10 0.3906598
#> 73      8        0.8        V1 0.3634691
#> 74      8        0.8        V5 0.3567093
#> 75      8        0.8        V2 0.3687265
#> 76      8        0.8        V3 0.4042390
#> 77      8        0.8        V4 0.3472642
#> 78      8        0.8        V8 0.3685874
#> 79      8        0.8        V6 0.3780477
#> 80      8        0.8        V7 0.3842323
#> 81      9        0.9        V5 0.3561807
#> 82      9        0.9        V7 0.3836338
#> 83      9        0.9        V8 0.3681020
#> 84      9        0.9        V9 0.3491430
#> 85      9        0.9        V6 0.3775409
#> 86      9        0.9       V10 0.3900146
#> 87      9        0.9        V1 0.3629202
#> 88      9        0.9        V2 0.3680525
#> 89      9        0.9        V3 0.4036564
#> 90      9        0.9        V4 0.3468152
#> 91     10        1.0        V1 0.3623722
#> 92     10        1.0        V4 0.3463668
#> 93     10        1.0        V5 0.3556530
#> 94     10        1.0        V9 0.3485537
#> 95     10        1.0        V3 0.4030747
#> 96     10        1.0        V7 0.3830363
#> 97     10        1.0        V8 0.3676172
#> 98     10        1.0        V2 0.3673798
#> 99     10        1.0        V6 0.3770347
#> 100    10        1.0       V10 0.3893704


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
#> 1       0                0          0 0.8757906 0.05084446 0.7741489 0.9508605
#> 2       0               20         20 0.8757906 0.05084446 0.7741489 0.9508605
#> 3       0               40         40 0.8757906 0.05084446 0.7741489 0.9508605
#> 4       0               60         60 0.8757906 0.05084446 0.7741489 0.9508605
#> 5      20                0         20 0.8617131 0.05258205 0.7549047 0.9395294
#> 6      20               20         40 0.8617131 0.05258205 0.7549047 0.9395294
#> 7      20               40         60 0.8617131 0.05258205 0.7549047 0.9395294
#> 8      20               60         80 0.8617131 0.05258205 0.7549047 0.9395294
#> 9      40                0         40 0.8478591 0.05413945 0.7363167 0.9279576
#> 10     40               20         60 0.8478591 0.05413945 0.7363167 0.9279576
#> 11     40               40         80 0.8478591 0.05413945 0.7363167 0.9279576
#> 12     40               60        100 0.8478591 0.05413945 0.7363167 0.9279576
#> 13     60                0         60 0.8342249 0.05554587 0.7183291 0.9162251
#> 14     60               20         80 0.8342249 0.05554587 0.7183291 0.9162251
#> 15     60               40        100 0.8342249 0.05554587 0.7183291 0.9162251
#> 16     60               60        120 0.8342249 0.05554587 0.7183291 0.9162251
#> 17     80                0         80 0.8208071 0.05682330 0.7008968 0.9043914
#> 18     80               20        100 0.8208071 0.05682330 0.7008968 0.9043914
#> 19     80               40        120 0.8208071 0.05682330 0.7008968 0.9043914
#> 20     80               60        140 0.8208071 0.05682330 0.7008968 0.9043914
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12213599 0.168182060 0.5861407
#> 2  0.30574618 0.11612806 0.134544092 0.5225299
#> 3  0.26001915 0.11037148 0.106889223 0.4661060
#> 4  0.22113100 0.10467972 0.084219527 0.4162271
#> 5  0.25589195 0.11143593 0.094026902 0.4657150
#> 6  0.21762106 0.10397401 0.073709517 0.4158819
#> 7  0.18507391 0.09707847 0.057174193 0.3719129
#> 8  0.15739448 0.09059760 0.043804564 0.3331516
#> 9  0.18213629 0.09847835 0.049562660 0.3716087
#> 10 0.15489621 0.09093215 0.037688466 0.3328835
#> 11 0.13173012 0.08402866 0.028215867 0.2987490
#> 12 0.11202872 0.07764734 0.020744352 0.2686465
#> 13 0.12963921 0.08519361 0.023936861 0.2985129
#> 14 0.11025053 0.07817779 0.017404797 0.2684382
#> 15 0.09376159 0.07178320 0.012367905 0.2418897
#> 16 0.07973872 0.06591390 0.008556302 0.2184205
#> 17 0.09227334 0.07261120 0.010163494 0.2417059
#> 18 0.07847305 0.06636853 0.006917095 0.2182578
#> 19 0.06673672 0.06068815 0.004550162 0.1974903
#> 20 0.05675565 0.05549453 0.002876603 0.1790549
```
