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
#> 1         0.1 0.3838513 0.01788040 0.3669151 0.4179322
#> 2         0.2 0.3832760 0.01785940 0.3663631 0.4173533
#> 3         0.3 0.3827015 0.01783862 0.3658120 0.4167752
#> 4         0.4 0.3821279 0.01781808 0.3652616 0.4161979
#> 5         0.5 0.3815551 0.01779778 0.3647121 0.4156215
#> 6         0.6 0.3809832 0.01777770 0.3641634 0.4150458
#> 7         0.7 0.3804122 0.01775785 0.3636156 0.4144709
#> 8         0.8 0.3798420 0.01773824 0.3630686 0.4138969
#> 9         0.9 0.3792726 0.01771885 0.3625223 0.4133236
#> 10        1.0 0.3787041 0.01769968 0.3619770 0.4127511

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3769221
#> 2       1        0.1        V5 0.3780271
#> 3       1        0.1        V9 0.4154263
#> 4       1        0.1        V3 0.3814589
#> 5       1        0.1        V4 0.3775088
#> 6       1        0.1        V8 0.3640099
#> 7       1        0.1        V2 0.4186596
#> 8       1        0.1        V6 0.3849456
#> 9       1        0.1        V7 0.3929474
#> 10      1        0.1       V10 0.4027224
#> 11      2        0.2        V7 0.3923435
#> 12      2        0.2        V8 0.3634661
#> 13      2        0.2        V9 0.4148600
#> 14      2        0.2       V10 0.4020904
#> 15      2        0.2        V1 0.3763420
#> 16      2        0.2        V5 0.3775679
#> 17      2        0.2        V2 0.4180772
#> 18      2        0.2        V3 0.3809881
#> 19      2        0.2        V4 0.3769430
#> 20      2        0.2        V6 0.3842739
#> 21      3        0.3        V4 0.3763781
#> 22      3        0.3        V5 0.3771092
#> 23      3        0.3        V3 0.3805179
#> 24      3        0.3        V7 0.3917406
#> 25      3        0.3        V8 0.3629230
#> 26      3        0.3        V9 0.4142943
#> 27      3        0.3        V6 0.3836034
#> 28      3        0.3       V10 0.4014593
#> 29      3        0.3        V1 0.3757628
#> 30      3        0.3        V2 0.4174955
#> 31      4        0.4        V1 0.3751844
#> 32      4        0.4        V9 0.4137295
#> 33      4        0.4        V3 0.3800483
#> 34      4        0.4        V4 0.3758140
#> 35      4        0.4        V5 0.3766510
#> 36      4        0.4        V2 0.4169146
#> 37      4        0.4        V6 0.3829340
#> 38      4        0.4        V7 0.3911386
#> 39      4        0.4        V8 0.3623808
#> 40      4        0.4       V10 0.4008292
#> 41      5        0.5        V8 0.3618394
#> 42      5        0.5        V9 0.4131654
#> 43      5        0.5       V10 0.4002001
#> 44      5        0.5        V1 0.3746070
#> 45      5        0.5        V5 0.3761935
#> 46      5        0.5        V2 0.4163345
#> 47      5        0.5        V3 0.3795793
#> 48      5        0.5        V4 0.3752507
#> 49      5        0.5        V6 0.3822658
#> 50      5        0.5        V7 0.3905375
#> 51      6        0.6        V4 0.3746883
#> 52      6        0.6        V5 0.3757364
#> 53      6        0.6        V7 0.3899373
#> 54      6        0.6        V8 0.3612988
#> 55      6        0.6        V9 0.4126021
#> 56      6        0.6        V6 0.3815988
#> 57      6        0.6       V10 0.3995720
#> 58      6        0.6        V1 0.3740304
#> 59      6        0.6        V2 0.4157553
#> 60      6        0.6        V3 0.3791108
#> 61      7        0.7        V1 0.3734547
#> 62      7        0.7        V3 0.3786430
#> 63      7        0.7        V4 0.3741267
#> 64      7        0.7        V5 0.3752800
#> 65      7        0.7        V9 0.4120396
#> 66      7        0.7        V6 0.3809330
#> 67      7        0.7        V7 0.3893381
#> 68      7        0.7        V8 0.3607591
#> 69      7        0.7        V2 0.4151768
#> 70      7        0.7       V10 0.3989449
#> 71      8        0.8        V9 0.4114778
#> 72      8        0.8       V10 0.3983187
#> 73      8        0.8        V1 0.3728800
#> 74      8        0.8        V5 0.3748240
#> 75      8        0.8        V2 0.4145992
#> 76      8        0.8        V3 0.3781757
#> 77      8        0.8        V4 0.3735659
#> 78      8        0.8        V8 0.3602201
#> 79      8        0.8        V6 0.3802683
#> 80      8        0.8        V7 0.3887398
#> 81      9        0.9        V5 0.3743687
#> 82      9        0.9        V7 0.3881424
#> 83      9        0.9        V8 0.3596819
#> 84      9        0.9        V9 0.4109168
#> 85      9        0.9        V6 0.3796047
#> 86      9        0.9       V10 0.3976936
#> 87      9        0.9        V1 0.3723061
#> 88      9        0.9        V2 0.4140223
#> 89      9        0.9        V3 0.3777089
#> 90      9        0.9        V4 0.3730060
#> 91     10        1.0        V1 0.3717330
#> 92     10        1.0        V4 0.3724470
#> 93     10        1.0        V5 0.3739139
#> 94     10        1.0        V9 0.4103565
#> 95     10        1.0        V3 0.3772428
#> 96     10        1.0        V7 0.3875459
#> 97     10        1.0        V8 0.3591445
#> 98     10        1.0        V2 0.4134463
#> 99     10        1.0        V6 0.3789424
#> 100    10        1.0       V10 0.3970694

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3769221
#> 2       1        0.1        V5 0.3780271
#> 3       1        0.1        V9 0.4154263
#> 4       1        0.1        V3 0.3814589
#> 5       1        0.1        V4 0.3775088
#> 6       1        0.1        V8 0.3640099
#> 7       1        0.1        V2 0.4186596
#> 8       1        0.1        V6 0.3849456
#> 9       1        0.1        V7 0.3929474
#> 10      1        0.1       V10 0.4027224
#> 11      2        0.2        V7 0.3923435
#> 12      2        0.2        V8 0.3634661
#> 13      2        0.2        V9 0.4148600
#> 14      2        0.2       V10 0.4020904
#> 15      2        0.2        V1 0.3763420
#> 16      2        0.2        V5 0.3775679
#> 17      2        0.2        V2 0.4180772
#> 18      2        0.2        V3 0.3809881
#> 19      2        0.2        V4 0.3769430
#> 20      2        0.2        V6 0.3842739
#> 21      3        0.3        V4 0.3763781
#> 22      3        0.3        V5 0.3771092
#> 23      3        0.3        V3 0.3805179
#> 24      3        0.3        V7 0.3917406
#> 25      3        0.3        V8 0.3629230
#> 26      3        0.3        V9 0.4142943
#> 27      3        0.3        V6 0.3836034
#> 28      3        0.3       V10 0.4014593
#> 29      3        0.3        V1 0.3757628
#> 30      3        0.3        V2 0.4174955
#> 31      4        0.4        V1 0.3751844
#> 32      4        0.4        V9 0.4137295
#> 33      4        0.4        V3 0.3800483
#> 34      4        0.4        V4 0.3758140
#> 35      4        0.4        V5 0.3766510
#> 36      4        0.4        V2 0.4169146
#> 37      4        0.4        V6 0.3829340
#> 38      4        0.4        V7 0.3911386
#> 39      4        0.4        V8 0.3623808
#> 40      4        0.4       V10 0.4008292
#> 41      5        0.5        V8 0.3618394
#> 42      5        0.5        V9 0.4131654
#> 43      5        0.5       V10 0.4002001
#> 44      5        0.5        V1 0.3746070
#> 45      5        0.5        V5 0.3761935
#> 46      5        0.5        V2 0.4163345
#> 47      5        0.5        V3 0.3795793
#> 48      5        0.5        V4 0.3752507
#> 49      5        0.5        V6 0.3822658
#> 50      5        0.5        V7 0.3905375
#> 51      6        0.6        V4 0.3746883
#> 52      6        0.6        V5 0.3757364
#> 53      6        0.6        V7 0.3899373
#> 54      6        0.6        V8 0.3612988
#> 55      6        0.6        V9 0.4126021
#> 56      6        0.6        V6 0.3815988
#> 57      6        0.6       V10 0.3995720
#> 58      6        0.6        V1 0.3740304
#> 59      6        0.6        V2 0.4157553
#> 60      6        0.6        V3 0.3791108
#> 61      7        0.7        V1 0.3734547
#> 62      7        0.7        V3 0.3786430
#> 63      7        0.7        V4 0.3741267
#> 64      7        0.7        V5 0.3752800
#> 65      7        0.7        V9 0.4120396
#> 66      7        0.7        V6 0.3809330
#> 67      7        0.7        V7 0.3893381
#> 68      7        0.7        V8 0.3607591
#> 69      7        0.7        V2 0.4151768
#> 70      7        0.7       V10 0.3989449
#> 71      8        0.8        V9 0.4114778
#> 72      8        0.8       V10 0.3983187
#> 73      8        0.8        V1 0.3728800
#> 74      8        0.8        V5 0.3748240
#> 75      8        0.8        V2 0.4145992
#> 76      8        0.8        V3 0.3781757
#> 77      8        0.8        V4 0.3735659
#> 78      8        0.8        V8 0.3602201
#> 79      8        0.8        V6 0.3802683
#> 80      8        0.8        V7 0.3887398
#> 81      9        0.9        V5 0.3743687
#> 82      9        0.9        V7 0.3881424
#> 83      9        0.9        V8 0.3596819
#> 84      9        0.9        V9 0.4109168
#> 85      9        0.9        V6 0.3796047
#> 86      9        0.9       V10 0.3976936
#> 87      9        0.9        V1 0.3723061
#> 88      9        0.9        V2 0.4140223
#> 89      9        0.9        V3 0.3777089
#> 90      9        0.9        V4 0.3730060
#> 91     10        1.0        V1 0.3717330
#> 92     10        1.0        V4 0.3724470
#> 93     10        1.0        V5 0.3739139
#> 94     10        1.0        V9 0.4103565
#> 95     10        1.0        V3 0.3772428
#> 96     10        1.0        V7 0.3875459
#> 97     10        1.0        V8 0.3591445
#> 98     10        1.0        V2 0.4134463
#> 99     10        1.0        V6 0.3789424
#> 100    10        1.0       V10 0.3970694


# get coefficient samples
coefs <- getNationalCoefficients(10)

# table of different scenarios to test
covTableSim <- expand.grid(Anthro = seq(0, 90, by = 20),
                           fire_excl_anthro = seq(0, 70, by = 20))
covTableSim$Total_dist = covTableSim$Anthro + covTableSim$fire_excl_anthro

estimateNationalRates(covTableSim, coefs)
#> popGrowthPars contains quantiles so they are used instead of the defaults
#> popGrowthPars contains quantiles so they are used instead of the defaults
#>    Anthro fire_excl_anthro Total_dist     S_bar   S_stdErr   S_PIlow  S_PIhigh
#> 1       0                0          0 0.8757906 0.04485075 0.7890023 0.9480141
#> 2       0               20         20 0.8757906 0.04485075 0.7890023 0.9480141
#> 3       0               40         40 0.8757906 0.04485075 0.7890023 0.9480141
#> 4       0               60         60 0.8757906 0.04485075 0.7890023 0.9480141
#> 5      20                0         20 0.8617131 0.04804631 0.7704946 0.9411687
#> 6      20               20         40 0.8617131 0.04804631 0.7704946 0.9411687
#> 7      20               40         60 0.8617131 0.04804631 0.7704946 0.9411687
#> 8      20               60         80 0.8617131 0.04804631 0.7704946 0.9411687
#> 9      40                0         40 0.8478591 0.05111076 0.7526235 0.9342362
#> 10     40               20         60 0.8478591 0.05111076 0.7526235 0.9342362
#> 11     40               40         80 0.8478591 0.05111076 0.7526235 0.9342362
#> 12     40               60        100 0.8478591 0.05111076 0.7526235 0.9342362
#> 13     60                0         60 0.8342249 0.05404813 0.7353292 0.9272356
#> 14     60               20         80 0.8342249 0.05404813 0.7353292 0.9272356
#> 15     60               40        100 0.8342249 0.05404813 0.7353292 0.9272356
#> 16     60               60        120 0.8342249 0.05404813 0.7353292 0.9272356
#> 17     80                0         80 0.8208071 0.05686330 0.7185643 0.9201824
#> 18     80               20        100 0.8208071 0.05686330 0.7185643 0.9201824
#> 19     80               40        120 0.8208071 0.05686330 0.7185643 0.9201824
#> 20     80               60        140 0.8208071 0.05686330 0.7185643 0.9201824
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12331206 0.157816151 0.5673626
#> 2  0.30574618 0.11707535 0.129119939 0.5234490
#> 3  0.26001915 0.11186714 0.105050889 0.4831025
#> 4  0.22113100 0.10719614 0.084912290 0.4461022
#> 5  0.25589195 0.11048816 0.089591413 0.4546858
#> 6  0.21762106 0.10460050 0.072014743 0.4200713
#> 7  0.18507391 0.09942287 0.057403589 0.3883857
#> 8  0.15739448 0.09467747 0.045319096 0.3593945
#> 9  0.18213629 0.09721622 0.048112435 0.3661158
#> 10 0.15489621 0.09185617 0.037675970 0.3390214
#> 11 0.13173012 0.08701803 0.029144836 0.3142330
#> 12 0.11202872 0.08253746 0.022231175 0.2915475
#> 13 0.12963921 0.08438128 0.023815113 0.2968082
#> 14 0.11025053 0.07962840 0.017950762 0.2755939
#> 15 0.09376159 0.07527402 0.013290011 0.2561594
#> 16 0.07973872 0.07121631 0.009638174 0.2383401
#> 17 0.09227334 0.07249714 0.010462544 0.2424758
#> 18 0.07847305 0.06835901 0.007455438 0.2257825
#> 19 0.06673672 0.06453029 0.005173463 0.2104469
#> 20 0.05675565 0.06094480 0.003481606 0.1963397
```
