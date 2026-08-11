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
#> 1         0.1 0.3838513 0.03122636 0.3232526 0.4171631
#> 2         0.2 0.3832760 0.03116684 0.3227824 0.4165240
#> 3         0.3 0.3827015 0.03110757 0.3223130 0.4158858
#> 4         0.4 0.3821279 0.03104853 0.3218442 0.4152486
#> 5         0.5 0.3815551 0.03098974 0.3213761 0.4146123
#> 6         0.6 0.3809832 0.03093118 0.3209086 0.4139771
#> 7         0.7 0.3804122 0.03087286 0.3204419 0.4133428
#> 8         0.8 0.3798420 0.03081477 0.3199758 0.4127095
#> 9         0.9 0.3792726 0.03075692 0.3195104 0.4120772
#> 10        1.0 0.3787041 0.03069931 0.3190457 0.4114458

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3665041
#> 2       1        0.1        V5 0.3778057
#> 3       1        0.1        V9 0.3991960
#> 4       1        0.1        V3 0.4210201
#> 5       1        0.1        V4 0.3404738
#> 6       1        0.1        V8 0.4013361
#> 7       1        0.1        V2 0.3182529
#> 8       1        0.1        V6 0.3718323
#> 9       1        0.1        V7 0.4038782
#> 10      1        0.1       V10 0.3890496
#> 11      2        0.2        V7 0.4033483
#> 12      2        0.2        V8 0.4007197
#> 13      2        0.2        V9 0.3984659
#> 14      2        0.2       V10 0.3885413
#> 15      2        0.2        V1 0.3659562
#> 16      2        0.2        V5 0.3772267
#> 17      2        0.2        V2 0.3177815
#> 18      2        0.2        V3 0.4203492
#> 19      2        0.2        V4 0.3400080
#> 20      2        0.2        V6 0.3712040
#> 21      3        0.3        V4 0.3395428
#> 22      3        0.3        V5 0.3766486
#> 23      3        0.3        V3 0.4196793
#> 24      3        0.3        V7 0.4028191
#> 25      3        0.3        V8 0.4001044
#> 26      3        0.3        V9 0.3977371
#> 27      3        0.3        V6 0.3705769
#> 28      3        0.3       V10 0.3880336
#> 29      3        0.3        V1 0.3654092
#> 30      3        0.3        V2 0.3173107
#> 31      4        0.4        V1 0.3648629
#> 32      4        0.4        V9 0.3970096
#> 33      4        0.4        V3 0.4190105
#> 34      4        0.4        V4 0.3390782
#> 35      4        0.4        V5 0.3760714
#> 36      4        0.4        V2 0.3168407
#> 37      4        0.4        V6 0.3699508
#> 38      4        0.4        V7 0.4022906
#> 39      4        0.4        V8 0.3994899
#> 40      4        0.4       V10 0.3875267
#> 41      5        0.5        V8 0.3988765
#> 42      5        0.5        V9 0.3962835
#> 43      5        0.5       V10 0.3870204
#> 44      5        0.5        V1 0.3643174
#> 45      5        0.5        V5 0.3754951
#> 46      5        0.5        V2 0.3163714
#> 47      5        0.5        V3 0.4183428
#> 48      5        0.5        V4 0.3386143
#> 49      5        0.5        V6 0.3693257
#> 50      5        0.5        V7 0.4017628
#> 51      6        0.6        V4 0.3381510
#> 52      6        0.6        V5 0.3749197
#> 53      6        0.6        V7 0.4012357
#> 54      6        0.6        V8 0.3982639
#> 55      6        0.6        V9 0.3955587
#> 56      6        0.6        V6 0.3687017
#> 57      6        0.6       V10 0.3865147
#> 58      6        0.6        V1 0.3637728
#> 59      6        0.6        V2 0.3159028
#> 60      6        0.6        V3 0.4176762
#> 61      7        0.7        V1 0.3632290
#> 62      7        0.7        V3 0.4170106
#> 63      7        0.7        V4 0.3376883
#> 64      7        0.7        V5 0.3743451
#> 65      7        0.7        V9 0.3948353
#> 66      7        0.7        V6 0.3680788
#> 67      7        0.7        V7 0.4007092
#> 68      7        0.7        V8 0.3976523
#> 69      7        0.7        V2 0.3154349
#> 70      7        0.7       V10 0.3860097
#> 71      8        0.8        V9 0.3941131
#> 72      8        0.8       V10 0.3855054
#> 73      8        0.8        V1 0.3626860
#> 74      8        0.8        V5 0.3737714
#> 75      8        0.8        V2 0.3149676
#> 76      8        0.8        V3 0.4163461
#> 77      8        0.8        V4 0.3372263
#> 78      8        0.8        V8 0.3970417
#> 79      8        0.8        V6 0.3674569
#> 80      8        0.8        V7 0.4001835
#> 81      9        0.9        V5 0.3731986
#> 82      9        0.9        V7 0.3996584
#> 83      9        0.9        V8 0.3964319
#> 84      9        0.9        V9 0.3933923
#> 85      9        0.9        V6 0.3668361
#> 86      9        0.9       V10 0.3850018
#> 87      9        0.9        V1 0.3621438
#> 88      9        0.9        V2 0.3145011
#> 89      9        0.9        V3 0.4156826
#> 90      9        0.9        V4 0.3367649
#> 91     10        1.0        V1 0.3616024
#> 92     10        1.0        V4 0.3363041
#> 93     10        1.0        V5 0.3726267
#> 94     10        1.0        V9 0.3926728
#> 95     10        1.0        V3 0.4150202
#> 96     10        1.0        V7 0.3991341
#> 97     10        1.0        V8 0.3958232
#> 98     10        1.0        V2 0.3140352
#> 99     10        1.0        V6 0.3662163
#> 100    10        1.0       V10 0.3844988

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3665041
#> 2       1        0.1        V5 0.3778057
#> 3       1        0.1        V9 0.3991960
#> 4       1        0.1        V3 0.4210201
#> 5       1        0.1        V4 0.3404738
#> 6       1        0.1        V8 0.4013361
#> 7       1        0.1        V2 0.3182529
#> 8       1        0.1        V6 0.3718323
#> 9       1        0.1        V7 0.4038782
#> 10      1        0.1       V10 0.3890496
#> 11      2        0.2        V7 0.4033483
#> 12      2        0.2        V8 0.4007197
#> 13      2        0.2        V9 0.3984659
#> 14      2        0.2       V10 0.3885413
#> 15      2        0.2        V1 0.3659562
#> 16      2        0.2        V5 0.3772267
#> 17      2        0.2        V2 0.3177815
#> 18      2        0.2        V3 0.4203492
#> 19      2        0.2        V4 0.3400080
#> 20      2        0.2        V6 0.3712040
#> 21      3        0.3        V4 0.3395428
#> 22      3        0.3        V5 0.3766486
#> 23      3        0.3        V3 0.4196793
#> 24      3        0.3        V7 0.4028191
#> 25      3        0.3        V8 0.4001044
#> 26      3        0.3        V9 0.3977371
#> 27      3        0.3        V6 0.3705769
#> 28      3        0.3       V10 0.3880336
#> 29      3        0.3        V1 0.3654092
#> 30      3        0.3        V2 0.3173107
#> 31      4        0.4        V1 0.3648629
#> 32      4        0.4        V9 0.3970096
#> 33      4        0.4        V3 0.4190105
#> 34      4        0.4        V4 0.3390782
#> 35      4        0.4        V5 0.3760714
#> 36      4        0.4        V2 0.3168407
#> 37      4        0.4        V6 0.3699508
#> 38      4        0.4        V7 0.4022906
#> 39      4        0.4        V8 0.3994899
#> 40      4        0.4       V10 0.3875267
#> 41      5        0.5        V8 0.3988765
#> 42      5        0.5        V9 0.3962835
#> 43      5        0.5       V10 0.3870204
#> 44      5        0.5        V1 0.3643174
#> 45      5        0.5        V5 0.3754951
#> 46      5        0.5        V2 0.3163714
#> 47      5        0.5        V3 0.4183428
#> 48      5        0.5        V4 0.3386143
#> 49      5        0.5        V6 0.3693257
#> 50      5        0.5        V7 0.4017628
#> 51      6        0.6        V4 0.3381510
#> 52      6        0.6        V5 0.3749197
#> 53      6        0.6        V7 0.4012357
#> 54      6        0.6        V8 0.3982639
#> 55      6        0.6        V9 0.3955587
#> 56      6        0.6        V6 0.3687017
#> 57      6        0.6       V10 0.3865147
#> 58      6        0.6        V1 0.3637728
#> 59      6        0.6        V2 0.3159028
#> 60      6        0.6        V3 0.4176762
#> 61      7        0.7        V1 0.3632290
#> 62      7        0.7        V3 0.4170106
#> 63      7        0.7        V4 0.3376883
#> 64      7        0.7        V5 0.3743451
#> 65      7        0.7        V9 0.3948353
#> 66      7        0.7        V6 0.3680788
#> 67      7        0.7        V7 0.4007092
#> 68      7        0.7        V8 0.3976523
#> 69      7        0.7        V2 0.3154349
#> 70      7        0.7       V10 0.3860097
#> 71      8        0.8        V9 0.3941131
#> 72      8        0.8       V10 0.3855054
#> 73      8        0.8        V1 0.3626860
#> 74      8        0.8        V5 0.3737714
#> 75      8        0.8        V2 0.3149676
#> 76      8        0.8        V3 0.4163461
#> 77      8        0.8        V4 0.3372263
#> 78      8        0.8        V8 0.3970417
#> 79      8        0.8        V6 0.3674569
#> 80      8        0.8        V7 0.4001835
#> 81      9        0.9        V5 0.3731986
#> 82      9        0.9        V7 0.3996584
#> 83      9        0.9        V8 0.3964319
#> 84      9        0.9        V9 0.3933923
#> 85      9        0.9        V6 0.3668361
#> 86      9        0.9       V10 0.3850018
#> 87      9        0.9        V1 0.3621438
#> 88      9        0.9        V2 0.3145011
#> 89      9        0.9        V3 0.4156826
#> 90      9        0.9        V4 0.3367649
#> 91     10        1.0        V1 0.3616024
#> 92     10        1.0        V4 0.3363041
#> 93     10        1.0        V5 0.3726267
#> 94     10        1.0        V9 0.3926728
#> 95     10        1.0        V3 0.4150202
#> 96     10        1.0        V7 0.3991341
#> 97     10        1.0        V8 0.3958232
#> 98     10        1.0        V2 0.3140352
#> 99     10        1.0        V6 0.3662163
#> 100    10        1.0       V10 0.3844988


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
#> 1       0                0          0 0.8757906 0.04438655 0.7911087 0.9524123
#> 2       0               20         20 0.8757906 0.04438655 0.7911087 0.9524123
#> 3       0               40         40 0.8757906 0.04438655 0.7911087 0.9524123
#> 4       0               60         60 0.8757906 0.04438655 0.7911087 0.9524123
#> 5      20                0         20 0.8617131 0.04614491 0.7733481 0.9393243
#> 6      20               20         40 0.8617131 0.04614491 0.7733481 0.9393243
#> 7      20               40         60 0.8617131 0.04614491 0.7733481 0.9393243
#> 8      20               60         80 0.8617131 0.04614491 0.7733481 0.9393243
#> 9      40                0         40 0.8478591 0.04779893 0.7561795 0.9259257
#> 10     40               20         60 0.8478591 0.04779893 0.7561795 0.9259257
#> 11     40               40         80 0.8478591 0.04779893 0.7561795 0.9259257
#> 12     40               60        100 0.8478591 0.04779893 0.7561795 0.9259257
#> 13     60                0         60 0.8342249 0.04936465 0.7395481 0.9123352
#> 14     60               20         80 0.8342249 0.04936465 0.7395481 0.9123352
#> 15     60               40        100 0.8342249 0.04936465 0.7395481 0.9123352
#> 16     60               60        120 0.8342249 0.04936465 0.7395481 0.9123352
#> 17     80                0         80 0.8208071 0.05085349 0.7234102 0.8986374
#> 18     80               20        100 0.8208071 0.05085349 0.7234102 0.8986374
#> 19     80               40        120 0.8208071 0.05085349 0.7234102 0.8986374
#> 20     80               60        140 0.8208071 0.05085349 0.7234102 0.8986374
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13266094 0.188081335 0.6375684
#> 2  0.30574618 0.12827157 0.147209506 0.5806481
#> 3  0.26001915 0.12399177 0.114325227 0.5288676
#> 4  0.22113100 0.11947814 0.087947518 0.4819632
#> 5  0.25589195 0.12338954 0.100947490 0.5169364
#> 6  0.21762106 0.11730323 0.077252951 0.4711756
#> 7  0.18507391 0.11160140 0.058392311 0.4298514
#> 8  0.15739448 0.10604990 0.043491718 0.3925810
#> 9  0.18213629 0.11037938 0.050802630 0.4203592
#> 10 0.15489621 0.10397096 0.037537449 0.3840235
#> 11 0.13173012 0.09802163 0.027215422 0.3512691
#> 12 0.11202872 0.09237133 0.019292588 0.3217355
#> 13 0.12963921 0.09652255 0.023147875 0.3437482
#> 14 0.11025053 0.09042219 0.016209249 0.3149515
#> 15 0.09376159 0.08477412 0.011024431 0.2889601
#> 16 0.07973872 0.07946686 0.007241625 0.2654703
#> 17 0.09227334 0.08321016 0.009055577 0.2829843
#> 18 0.07847305 0.07768513 0.005836287 0.2600642
#> 19 0.06673672 0.07257255 0.003594724 0.2393046
#> 20 0.05675565 0.06779278 0.002097207 0.2204608
```
