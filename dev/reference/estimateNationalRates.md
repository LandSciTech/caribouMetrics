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
#> 1         0.1 0.3838513 0.01609206 0.3517118 0.4036301
#> 2         0.2 0.3832760 0.01606029 0.3511706 0.4029841
#> 3         0.3 0.3827015 0.01602865 0.3506303 0.4023391
#> 4         0.4 0.3821279 0.01599713 0.3500908 0.4016951
#> 5         0.5 0.3815551 0.01596574 0.3495522 0.4010522
#> 6         0.6 0.3809832 0.01593447 0.3490143 0.4004103
#> 7         0.7 0.3804122 0.01590333 0.3484773 0.3997695
#> 8         0.8 0.3798420 0.01587231 0.3479412 0.3991296
#> 9         0.9 0.3792726 0.01584142 0.3474058 0.3984908
#> 10        1.0 0.3787041 0.01581066 0.3468713 0.3978530

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3819224
#> 2       1        0.1        V5 0.3719531
#> 3       1        0.1        V9 0.3894035
#> 4       1        0.1        V3 0.3729329
#> 5       1        0.1        V4 0.3678971
#> 6       1        0.1        V8 0.3736089
#> 7       1        0.1        V2 0.3852026
#> 8       1        0.1        V6 0.3633122
#> 9       1        0.1        V7 0.3483439
#> 10      1        0.1       V10 0.4077604
#> 11      2        0.2        V7 0.3478268
#> 12      2        0.2        V8 0.3730604
#> 13      2        0.2        V9 0.3888140
#> 14      2        0.2       V10 0.4070980
#> 15      2        0.2        V1 0.3813558
#> 16      2        0.2        V5 0.3713867
#> 17      2        0.2        V2 0.3846270
#> 18      2        0.2        V3 0.3723881
#> 19      2        0.2        V4 0.3673837
#> 20      2        0.2        V6 0.3626882
#> 21      3        0.3        V4 0.3668710
#> 22      3        0.3        V5 0.3708211
#> 23      3        0.3        V3 0.3718441
#> 24      3        0.3        V7 0.3473104
#> 25      3        0.3        V8 0.3725127
#> 26      3        0.3        V9 0.3882253
#> 27      3        0.3        V6 0.3620654
#> 28      3        0.3       V10 0.4064366
#> 29      3        0.3        V1 0.3807900
#> 30      3        0.3        V2 0.3840522
#> 31      4        0.4        V1 0.3802251
#> 32      4        0.4        V9 0.3876376
#> 33      4        0.4        V3 0.3713009
#> 34      4        0.4        V4 0.3663591
#> 35      4        0.4        V5 0.3702564
#> 36      4        0.4        V2 0.3834782
#> 37      4        0.4        V6 0.3614436
#> 38      4        0.4        V7 0.3467948
#> 39      4        0.4        V8 0.3719658
#> 40      4        0.4       V10 0.4057764
#> 41      5        0.5        V8 0.3714197
#> 42      5        0.5        V9 0.3870508
#> 43      5        0.5       V10 0.4051172
#> 44      5        0.5        V1 0.3796610
#> 45      5        0.5        V5 0.3696925
#> 46      5        0.5        V2 0.3829052
#> 47      5        0.5        V3 0.3707585
#> 48      5        0.5        V4 0.3658478
#> 49      5        0.5        V6 0.3608228
#> 50      5        0.5        V7 0.3462800
#> 51      6        0.6        V4 0.3653372
#> 52      6        0.6        V5 0.3691295
#> 53      6        0.6        V7 0.3457660
#> 54      6        0.6        V8 0.3708744
#> 55      6        0.6        V9 0.3864648
#> 56      6        0.6        V6 0.3602032
#> 57      6        0.6       V10 0.4044590
#> 58      6        0.6        V1 0.3790978
#> 59      6        0.6        V2 0.3823329
#> 60      6        0.6        V3 0.3702169
#> 61      7        0.7        V1 0.3785354
#> 62      7        0.7        V3 0.3696761
#> 63      7        0.7        V4 0.3648274
#> 64      7        0.7        V5 0.3685674
#> 65      7        0.7        V9 0.3858798
#> 66      7        0.7        V6 0.3595846
#> 67      7        0.7        V7 0.3452527
#> 68      7        0.7        V8 0.3703299
#> 69      7        0.7        V2 0.3817616
#> 70      7        0.7       V10 0.4038020
#> 71      8        0.8        V9 0.3852956
#> 72      8        0.8       V10 0.4031460
#> 73      8        0.8        V1 0.3779738
#> 74      8        0.8        V5 0.3680061
#> 75      8        0.8        V2 0.3811911
#> 76      8        0.8        V3 0.3691361
#> 77      8        0.8        V4 0.3643183
#> 78      8        0.8        V8 0.3697862
#> 79      8        0.8        V6 0.3589670
#> 80      8        0.8        V7 0.3447401
#> 81      9        0.9        V5 0.3674456
#> 82      9        0.9        V7 0.3442283
#> 83      9        0.9        V8 0.3692433
#> 84      9        0.9        V9 0.3847123
#> 85      9        0.9        V6 0.3583506
#> 86      9        0.9       V10 0.4024911
#> 87      9        0.9        V1 0.3774130
#> 88      9        0.9        V2 0.3806214
#> 89      9        0.9        V3 0.3685969
#> 90      9        0.9        V4 0.3638099
#> 91     10        1.0        V1 0.3768531
#> 92     10        1.0        V4 0.3633022
#> 93     10        1.0        V5 0.3668861
#> 94     10        1.0        V9 0.3841299
#> 95     10        1.0        V3 0.3680584
#> 96     10        1.0        V7 0.3437173
#> 97     10        1.0        V8 0.3687012
#> 98     10        1.0        V2 0.3800526
#> 99     10        1.0        V6 0.3577351
#> 100    10        1.0       V10 0.4018372

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3819224
#> 2       1        0.1        V5 0.3719531
#> 3       1        0.1        V9 0.3894035
#> 4       1        0.1        V3 0.3729329
#> 5       1        0.1        V4 0.3678971
#> 6       1        0.1        V8 0.3736089
#> 7       1        0.1        V2 0.3852026
#> 8       1        0.1        V6 0.3633122
#> 9       1        0.1        V7 0.3483439
#> 10      1        0.1       V10 0.4077604
#> 11      2        0.2        V7 0.3478268
#> 12      2        0.2        V8 0.3730604
#> 13      2        0.2        V9 0.3888140
#> 14      2        0.2       V10 0.4070980
#> 15      2        0.2        V1 0.3813558
#> 16      2        0.2        V5 0.3713867
#> 17      2        0.2        V2 0.3846270
#> 18      2        0.2        V3 0.3723881
#> 19      2        0.2        V4 0.3673837
#> 20      2        0.2        V6 0.3626882
#> 21      3        0.3        V4 0.3668710
#> 22      3        0.3        V5 0.3708211
#> 23      3        0.3        V3 0.3718441
#> 24      3        0.3        V7 0.3473104
#> 25      3        0.3        V8 0.3725127
#> 26      3        0.3        V9 0.3882253
#> 27      3        0.3        V6 0.3620654
#> 28      3        0.3       V10 0.4064366
#> 29      3        0.3        V1 0.3807900
#> 30      3        0.3        V2 0.3840522
#> 31      4        0.4        V1 0.3802251
#> 32      4        0.4        V9 0.3876376
#> 33      4        0.4        V3 0.3713009
#> 34      4        0.4        V4 0.3663591
#> 35      4        0.4        V5 0.3702564
#> 36      4        0.4        V2 0.3834782
#> 37      4        0.4        V6 0.3614436
#> 38      4        0.4        V7 0.3467948
#> 39      4        0.4        V8 0.3719658
#> 40      4        0.4       V10 0.4057764
#> 41      5        0.5        V8 0.3714197
#> 42      5        0.5        V9 0.3870508
#> 43      5        0.5       V10 0.4051172
#> 44      5        0.5        V1 0.3796610
#> 45      5        0.5        V5 0.3696925
#> 46      5        0.5        V2 0.3829052
#> 47      5        0.5        V3 0.3707585
#> 48      5        0.5        V4 0.3658478
#> 49      5        0.5        V6 0.3608228
#> 50      5        0.5        V7 0.3462800
#> 51      6        0.6        V4 0.3653372
#> 52      6        0.6        V5 0.3691295
#> 53      6        0.6        V7 0.3457660
#> 54      6        0.6        V8 0.3708744
#> 55      6        0.6        V9 0.3864648
#> 56      6        0.6        V6 0.3602032
#> 57      6        0.6       V10 0.4044590
#> 58      6        0.6        V1 0.3790978
#> 59      6        0.6        V2 0.3823329
#> 60      6        0.6        V3 0.3702169
#> 61      7        0.7        V1 0.3785354
#> 62      7        0.7        V3 0.3696761
#> 63      7        0.7        V4 0.3648274
#> 64      7        0.7        V5 0.3685674
#> 65      7        0.7        V9 0.3858798
#> 66      7        0.7        V6 0.3595846
#> 67      7        0.7        V7 0.3452527
#> 68      7        0.7        V8 0.3703299
#> 69      7        0.7        V2 0.3817616
#> 70      7        0.7       V10 0.4038020
#> 71      8        0.8        V9 0.3852956
#> 72      8        0.8       V10 0.4031460
#> 73      8        0.8        V1 0.3779738
#> 74      8        0.8        V5 0.3680061
#> 75      8        0.8        V2 0.3811911
#> 76      8        0.8        V3 0.3691361
#> 77      8        0.8        V4 0.3643183
#> 78      8        0.8        V8 0.3697862
#> 79      8        0.8        V6 0.3589670
#> 80      8        0.8        V7 0.3447401
#> 81      9        0.9        V5 0.3674456
#> 82      9        0.9        V7 0.3442283
#> 83      9        0.9        V8 0.3692433
#> 84      9        0.9        V9 0.3847123
#> 85      9        0.9        V6 0.3583506
#> 86      9        0.9       V10 0.4024911
#> 87      9        0.9        V1 0.3774130
#> 88      9        0.9        V2 0.3806214
#> 89      9        0.9        V3 0.3685969
#> 90      9        0.9        V4 0.3638099
#> 91     10        1.0        V1 0.3768531
#> 92     10        1.0        V4 0.3633022
#> 93     10        1.0        V5 0.3668861
#> 94     10        1.0        V9 0.3841299
#> 95     10        1.0        V3 0.3680584
#> 96     10        1.0        V7 0.3437173
#> 97     10        1.0        V8 0.3687012
#> 98     10        1.0        V2 0.3800526
#> 99     10        1.0        V6 0.3577351
#> 100    10        1.0       V10 0.4018372


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
#> 1       0                0          0 0.8757906 0.04773220 0.7788377 0.9425038
#> 2       0               20         20 0.8757906 0.04773220 0.7788377 0.9425038
#> 3       0               40         40 0.8757906 0.04773220 0.7788377 0.9425038
#> 4       0               60         60 0.8757906 0.04773220 0.7788377 0.9425038
#> 5      20                0         20 0.8617131 0.04988550 0.7605282 0.9330522
#> 6      20               20         40 0.8617131 0.04988550 0.7605282 0.9330522
#> 7      20               40         60 0.8617131 0.04988550 0.7605282 0.9330522
#> 8      20               60         80 0.8617131 0.04988550 0.7605282 0.9330522
#> 9      40                0         40 0.8478591 0.05188544 0.7428134 0.9234888
#> 10     40               20         60 0.8478591 0.05188544 0.7428134 0.9234888
#> 11     40               40         80 0.8478591 0.05188544 0.7428134 0.9234888
#> 12     40               60        100 0.8478591 0.05188544 0.7428134 0.9234888
#> 13     60                0         60 0.8342249 0.05374923 0.7256441 0.9138492
#> 14     60               20         80 0.8342249 0.05374923 0.7256441 0.9138492
#> 15     60               40        100 0.8342249 0.05374923 0.7256441 0.9138492
#> 16     60               60        120 0.8342249 0.05374923 0.7256441 0.9138492
#> 17     80                0         80 0.8208071 0.05549062 0.7089804 0.9041618
#> 18     80               20        100 0.8208071 0.05549062 0.7089804 0.9041618
#> 19     80               40        120 0.8208071 0.05549062 0.7089804 0.9041618
#> 20     80               60        140 0.8208071 0.05549062 0.7089804 0.9041618
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12449422 0.1716884127 0.5722469
#> 2  0.30574618 0.11855317 0.1324809913 0.5093957
#> 3  0.26001915 0.11198500 0.1012665862 0.4538513
#> 4  0.22113100 0.10509956 0.0765198268 0.4049189
#> 5  0.25589195 0.10756997 0.0871145417 0.4397460
#> 6  0.21762106 0.10091817 0.0653528196 0.3925077
#> 7  0.18507391 0.09422694 0.0482832464 0.3509672
#> 8  0.15739448 0.08762138 0.0350264354 0.3144444
#> 9  0.18213629 0.09063282 0.0406631913 0.3404356
#> 10 0.15489621 0.08427303 0.0291655114 0.3051823
#> 11 0.13173012 0.07810609 0.0204225512 0.2741596
#> 12 0.11202872 0.07217984 0.0138957099 0.2468198
#> 13 0.12963921 0.07520640 0.0166351433 0.2662843
#> 14 0.11025053 0.06951034 0.0111179853 0.2398709
#> 15 0.09376159 0.06408831 0.0071538926 0.2165296
#> 16 0.07973872 0.05895386 0.0043982604 0.1958456
#> 17 0.09227334 0.06177595 0.0055278080 0.2105843
#> 18 0.07847305 0.05684701 0.0033032876 0.1905662
#> 19 0.06673672 0.05220612 0.0018584663 0.1727496
#> 20 0.05675565 0.04785185 0.0009727971 0.1568285
```
