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
#> 1         0.1 0.3838513 0.02249828 0.3478785 0.4196239
#> 2         0.2 0.3832760 0.02246845 0.3473445 0.4189916
#> 3         0.3 0.3827015 0.02243874 0.3468113 0.4183603
#> 4         0.4 0.3821279 0.02240915 0.3462789 0.4177299
#> 5         0.5 0.3815551 0.02237966 0.3457473 0.4171005
#> 6         0.6 0.3809832 0.02235029 0.3452165 0.4164720
#> 7         0.7 0.3804122 0.02232104 0.3446866 0.4158445
#> 8         0.8 0.3798420 0.02229190 0.3441575 0.4152179
#> 9         0.9 0.3792726 0.02226287 0.3436291 0.4145923
#> 10        1.0 0.3787041 0.02223395 0.3431016 0.4139676

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3821575
#> 2       1        0.1        V5 0.3597539
#> 3       1        0.1        V9 0.3959668
#> 4       1        0.1        V3 0.3847881
#> 5       1        0.1        V4 0.3949819
#> 6       1        0.1        V8 0.4231491
#> 7       1        0.1        V2 0.3967499
#> 8       1        0.1        V6 0.3840757
#> 9       1        0.1        V7 0.3444308
#> 10      1        0.1       V10 0.4074815
#> 11      2        0.2        V7 0.3439142
#> 12      2        0.2        V8 0.4225138
#> 13      2        0.2        V9 0.3953827
#> 14      2        0.2       V10 0.4068597
#> 15      2        0.2        V1 0.3816790
#> 16      2        0.2        V5 0.3591598
#> 17      2        0.2        V2 0.3961154
#> 18      2        0.2        V3 0.3841802
#> 19      2        0.2        V4 0.3944232
#> 20      2        0.2        V6 0.3835121
#> 21      3        0.3        V4 0.3938653
#> 22      3        0.3        V5 0.3585666
#> 23      3        0.3        V3 0.3835732
#> 24      3        0.3        V7 0.3433984
#> 25      3        0.3        V8 0.4218794
#> 26      3        0.3        V9 0.3947994
#> 27      3        0.3        V6 0.3829494
#> 28      3        0.3       V10 0.4062389
#> 29      3        0.3        V1 0.3812011
#> 30      3        0.3        V2 0.3954819
#> 31      4        0.4        V1 0.3807238
#> 32      4        0.4        V9 0.3942170
#> 33      4        0.4        V3 0.3829671
#> 34      4        0.4        V4 0.3933082
#> 35      4        0.4        V5 0.3579744
#> 36      4        0.4        V2 0.3948495
#> 37      4        0.4        V6 0.3823874
#> 38      4        0.4        V7 0.3428834
#> 39      4        0.4        V8 0.4212460
#> 40      4        0.4       V10 0.4056190
#> 41      5        0.5        V8 0.4206135
#> 42      5        0.5        V9 0.3936354
#> 43      5        0.5       V10 0.4050001
#> 44      5        0.5        V1 0.3802471
#> 45      5        0.5        V5 0.3573832
#> 46      5        0.5        V2 0.3942180
#> 47      5        0.5        V3 0.3823621
#> 48      5        0.5        V4 0.3927518
#> 49      5        0.5        V6 0.3818263
#> 50      5        0.5        V7 0.3423691
#> 51      6        0.6        V4 0.3921963
#> 52      6        0.6        V5 0.3567929
#> 53      6        0.6        V7 0.3418556
#> 54      6        0.6        V8 0.4199820
#> 55      6        0.6        V9 0.3930547
#> 56      6        0.6        V6 0.3812660
#> 57      6        0.6       V10 0.4043821
#> 58      6        0.6        V1 0.3797710
#> 59      6        0.6        V2 0.3935876
#> 60      6        0.6        V3 0.3817580
#> 61      7        0.7        V1 0.3792955
#> 62      7        0.7        V3 0.3811548
#> 63      7        0.7        V4 0.3916415
#> 64      7        0.7        V5 0.3562037
#> 65      7        0.7        V9 0.3924749
#> 66      7        0.7        V6 0.3807066
#> 67      7        0.7        V7 0.3413429
#> 68      7        0.7        V8 0.4193514
#> 69      7        0.7        V2 0.3929581
#> 70      7        0.7       V10 0.4037651
#> 71      8        0.8        V9 0.3918959
#> 72      8        0.8       V10 0.4031490
#> 73      8        0.8        V1 0.3788206
#> 74      8        0.8        V5 0.3556154
#> 75      8        0.8        V2 0.3923297
#> 76      8        0.8        V3 0.3805526
#> 77      8        0.8        V4 0.3910875
#> 78      8        0.8        V8 0.4187218
#> 79      8        0.8        V6 0.3801479
#> 80      8        0.8        V7 0.3408310
#> 81      9        0.9        V5 0.3550281
#> 82      9        0.9        V7 0.3403198
#> 83      9        0.9        V8 0.4180931
#> 84      9        0.9        V9 0.3913177
#> 85      9        0.9        V6 0.3795901
#> 86      9        0.9       V10 0.4025339
#> 87      9        0.9        V1 0.3783463
#> 88      9        0.9        V2 0.3917023
#> 89      9        0.9        V3 0.3799513
#> 90      9        0.9        V4 0.3905343
#> 91     10        1.0        V1 0.3778726
#> 92     10        1.0        V4 0.3899819
#> 93     10        1.0        V5 0.3544417
#> 94     10        1.0        V9 0.3907404
#> 95     10        1.0        V3 0.3793510
#> 96     10        1.0        V7 0.3398093
#> 97     10        1.0        V8 0.4174653
#> 98     10        1.0        V2 0.3910758
#> 99     10        1.0        V6 0.3790331
#> 100    10        1.0       V10 0.4019197

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3821575
#> 2       1        0.1        V5 0.3597539
#> 3       1        0.1        V9 0.3959668
#> 4       1        0.1        V3 0.3847881
#> 5       1        0.1        V4 0.3949819
#> 6       1        0.1        V8 0.4231491
#> 7       1        0.1        V2 0.3967499
#> 8       1        0.1        V6 0.3840757
#> 9       1        0.1        V7 0.3444308
#> 10      1        0.1       V10 0.4074815
#> 11      2        0.2        V7 0.3439142
#> 12      2        0.2        V8 0.4225138
#> 13      2        0.2        V9 0.3953827
#> 14      2        0.2       V10 0.4068597
#> 15      2        0.2        V1 0.3816790
#> 16      2        0.2        V5 0.3591598
#> 17      2        0.2        V2 0.3961154
#> 18      2        0.2        V3 0.3841802
#> 19      2        0.2        V4 0.3944232
#> 20      2        0.2        V6 0.3835121
#> 21      3        0.3        V4 0.3938653
#> 22      3        0.3        V5 0.3585666
#> 23      3        0.3        V3 0.3835732
#> 24      3        0.3        V7 0.3433984
#> 25      3        0.3        V8 0.4218794
#> 26      3        0.3        V9 0.3947994
#> 27      3        0.3        V6 0.3829494
#> 28      3        0.3       V10 0.4062389
#> 29      3        0.3        V1 0.3812011
#> 30      3        0.3        V2 0.3954819
#> 31      4        0.4        V1 0.3807238
#> 32      4        0.4        V9 0.3942170
#> 33      4        0.4        V3 0.3829671
#> 34      4        0.4        V4 0.3933082
#> 35      4        0.4        V5 0.3579744
#> 36      4        0.4        V2 0.3948495
#> 37      4        0.4        V6 0.3823874
#> 38      4        0.4        V7 0.3428834
#> 39      4        0.4        V8 0.4212460
#> 40      4        0.4       V10 0.4056190
#> 41      5        0.5        V8 0.4206135
#> 42      5        0.5        V9 0.3936354
#> 43      5        0.5       V10 0.4050001
#> 44      5        0.5        V1 0.3802471
#> 45      5        0.5        V5 0.3573832
#> 46      5        0.5        V2 0.3942180
#> 47      5        0.5        V3 0.3823621
#> 48      5        0.5        V4 0.3927518
#> 49      5        0.5        V6 0.3818263
#> 50      5        0.5        V7 0.3423691
#> 51      6        0.6        V4 0.3921963
#> 52      6        0.6        V5 0.3567929
#> 53      6        0.6        V7 0.3418556
#> 54      6        0.6        V8 0.4199820
#> 55      6        0.6        V9 0.3930547
#> 56      6        0.6        V6 0.3812660
#> 57      6        0.6       V10 0.4043821
#> 58      6        0.6        V1 0.3797710
#> 59      6        0.6        V2 0.3935876
#> 60      6        0.6        V3 0.3817580
#> 61      7        0.7        V1 0.3792955
#> 62      7        0.7        V3 0.3811548
#> 63      7        0.7        V4 0.3916415
#> 64      7        0.7        V5 0.3562037
#> 65      7        0.7        V9 0.3924749
#> 66      7        0.7        V6 0.3807066
#> 67      7        0.7        V7 0.3413429
#> 68      7        0.7        V8 0.4193514
#> 69      7        0.7        V2 0.3929581
#> 70      7        0.7       V10 0.4037651
#> 71      8        0.8        V9 0.3918959
#> 72      8        0.8       V10 0.4031490
#> 73      8        0.8        V1 0.3788206
#> 74      8        0.8        V5 0.3556154
#> 75      8        0.8        V2 0.3923297
#> 76      8        0.8        V3 0.3805526
#> 77      8        0.8        V4 0.3910875
#> 78      8        0.8        V8 0.4187218
#> 79      8        0.8        V6 0.3801479
#> 80      8        0.8        V7 0.3408310
#> 81      9        0.9        V5 0.3550281
#> 82      9        0.9        V7 0.3403198
#> 83      9        0.9        V8 0.4180931
#> 84      9        0.9        V9 0.3913177
#> 85      9        0.9        V6 0.3795901
#> 86      9        0.9       V10 0.4025339
#> 87      9        0.9        V1 0.3783463
#> 88      9        0.9        V2 0.3917023
#> 89      9        0.9        V3 0.3799513
#> 90      9        0.9        V4 0.3905343
#> 91     10        1.0        V1 0.3778726
#> 92     10        1.0        V4 0.3899819
#> 93     10        1.0        V5 0.3544417
#> 94     10        1.0        V9 0.3907404
#> 95     10        1.0        V3 0.3793510
#> 96     10        1.0        V7 0.3398093
#> 97     10        1.0        V8 0.4174653
#> 98     10        1.0        V2 0.3910758
#> 99     10        1.0        V6 0.3790331
#> 100    10        1.0       V10 0.4019197


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
#> 1       0                0          0 0.8757906 0.04712009 0.7788986 0.9400293
#> 2       0               20         20 0.8757906 0.04712009 0.7788986 0.9400293
#> 3       0               40         40 0.8757906 0.04712009 0.7788986 0.9400293
#> 4       0               60         60 0.8757906 0.04712009 0.7788986 0.9400293
#> 5      20                0         20 0.8617131 0.05008992 0.7573298 0.9298567
#> 6      20               20         40 0.8617131 0.05008992 0.7573298 0.9298567
#> 7      20               40         60 0.8617131 0.05008992 0.7573298 0.9298567
#> 8      20               60         80 0.8617131 0.05008992 0.7573298 0.9298567
#> 9      40                0         40 0.8478591 0.05283184 0.7366247 0.9195689
#> 10     40               20         60 0.8478591 0.05283184 0.7366247 0.9195689
#> 11     40               40         80 0.8478591 0.05283184 0.7366247 0.9195689
#> 12     40               60        100 0.8478591 0.05283184 0.7366247 0.9195689
#> 13     60                0         60 0.8342249 0.05537348 0.7166932 0.9092064
#> 14     60               20         80 0.8342249 0.05537348 0.7166932 0.9092064
#> 15     60               40        100 0.8342249 0.05537348 0.7166932 0.9092064
#> 16     60               60        120 0.8342249 0.05537348 0.7166932 0.9092064
#> 17     80                0         80 0.8208071 0.05773675 0.6974654 0.8988009
#> 18     80               20        100 0.8208071 0.05773675 0.6974654 0.8988009
#> 19     80               40        120 0.8208071 0.05773675 0.6974654 0.8988009
#> 20     80               60        140 0.8208071 0.05773675 0.6974654 0.8988009
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12530264 0.176628739 0.5742897
#> 2  0.30574618 0.11478838 0.133390322 0.4957822
#> 3  0.26001915 0.10524221 0.099637695 0.4287492
#> 4  0.22113100 0.09653459 0.073420058 0.3717445
#> 5  0.25589195 0.11133245 0.105420328 0.4632116
#> 6  0.21762106 0.10113337 0.077899239 0.4010344
#> 7  0.18507391 0.09204984 0.056648575 0.3482025
#> 8  0.15739448 0.08389078 0.040401189 0.3033169
#> 9  0.18213629 0.09729221 0.060267879 0.3753467
#> 10 0.15489621 0.08798493 0.043154444 0.3263832
#> 11 0.13173012 0.07973624 0.030205289 0.2847602
#> 12 0.11202872 0.07236036 0.020565754 0.2492933
#> 13 0.12963921 0.08412660 0.032388188 0.3061547
#> 14 0.11025053 0.07588103 0.022177463 0.2675372
#> 15 0.09376159 0.06856869 0.014702194 0.2345789
#> 16 0.07973872 0.06202825 0.009369671 0.2063267
#> 17 0.09227334 0.07221497 0.015941598 0.2515403
#> 18 0.07847305 0.06502116 0.010242298 0.2208839
#> 19 0.06673672 0.05862098 0.006282445 0.1945394
#> 20 0.05675565 0.05288017 0.003641234 0.1717587
```
