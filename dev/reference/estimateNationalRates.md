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
#> 1         0.1 0.3838513 0.01912314 0.3401919 0.3925434
#> 2         0.2 0.3832760 0.01909284 0.3396532 0.3919561
#> 3         0.3 0.3827015 0.01906265 0.3391155 0.3913697
#> 4         0.4 0.3821279 0.01903259 0.3385785 0.3907842
#> 5         0.5 0.3815551 0.01900263 0.3380424 0.3901996
#> 6         0.6 0.3809832 0.01897280 0.3375072 0.3896158
#> 7         0.7 0.3804122 0.01894307 0.3369728 0.3890329
#> 8         0.8 0.3798420 0.01891347 0.3364393 0.3884509
#> 9         0.9 0.3792726 0.01888397 0.3359066 0.3878697
#> 10        1.0 0.3787041 0.01885459 0.3353747 0.3872894

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3461393
#> 2       1        0.1        V5 0.3755669
#> 3       1        0.1        V9 0.3400546
#> 4       1        0.1        V3 0.3711769
#> 5       1        0.1        V4 0.3640544
#> 6       1        0.1        V8 0.3934812
#> 7       1        0.1        V2 0.3643627
#> 8       1        0.1        V6 0.3406648
#> 9       1        0.1        V7 0.3786267
#> 10      1        0.1       V10 0.3893133
#> 11      2        0.2        V7 0.3779687
#> 12      2        0.2        V8 0.3928977
#> 13      2        0.2        V9 0.3395057
#> 14      2        0.2       V10 0.3887129
#> 15      2        0.2        V1 0.3455855
#> 16      2        0.2        V5 0.3749844
#> 17      2        0.2        V2 0.3638293
#> 18      2        0.2        V3 0.3705733
#> 19      2        0.2        V4 0.3635422
#> 20      2        0.2        V6 0.3401614
#> 21      3        0.3        V4 0.3630308
#> 22      3        0.3        V5 0.3744027
#> 23      3        0.3        V3 0.3699708
#> 24      3        0.3        V7 0.3773118
#> 25      3        0.3        V8 0.3923151
#> 26      3        0.3        V9 0.3389577
#> 27      3        0.3        V6 0.3396587
#> 28      3        0.3       V10 0.3881133
#> 29      3        0.3        V1 0.3450326
#> 30      3        0.3        V2 0.3632967
#> 31      4        0.4        V1 0.3444806
#> 32      4        0.4        V9 0.3384106
#> 33      4        0.4        V3 0.3693692
#> 34      4        0.4        V4 0.3625201
#> 35      4        0.4        V5 0.3738219
#> 36      4        0.4        V2 0.3627648
#> 37      4        0.4        V6 0.3391568
#> 38      4        0.4        V7 0.3766561
#> 39      4        0.4        V8 0.3917334
#> 40      4        0.4       V10 0.3875147
#> 41      5        0.5        V8 0.3911526
#> 42      5        0.5        V9 0.3378644
#> 43      5        0.5       V10 0.3869171
#> 44      5        0.5        V1 0.3439294
#> 45      5        0.5        V5 0.3732421
#> 46      5        0.5        V2 0.3622338
#> 47      5        0.5        V3 0.3687686
#> 48      5        0.5        V4 0.3620101
#> 49      5        0.5        V6 0.3386555
#> 50      5        0.5        V7 0.3760015
#> 51      6        0.6        V4 0.3615008
#> 52      6        0.6        V5 0.3726631
#> 53      6        0.6        V7 0.3753481
#> 54      6        0.6        V8 0.3905726
#> 55      6        0.6        V9 0.3373191
#> 56      6        0.6        V6 0.3381551
#> 57      6        0.6       V10 0.3863203
#> 58      6        0.6        V1 0.3433792
#> 59      6        0.6        V2 0.3617035
#> 60      6        0.6        V3 0.3681690
#> 61      7        0.7        V1 0.3428298
#> 62      7        0.7        V3 0.3675704
#> 63      7        0.7        V4 0.3609922
#> 64      7        0.7        V5 0.3720850
#> 65      7        0.7        V9 0.3367747
#> 66      7        0.7        V6 0.3376553
#> 67      7        0.7        V7 0.3746958
#> 68      7        0.7        V8 0.3899934
#> 69      7        0.7        V2 0.3611740
#> 70      7        0.7       V10 0.3857245
#> 71      8        0.8        V9 0.3362311
#> 72      8        0.8       V10 0.3851295
#> 73      8        0.8        V1 0.3422813
#> 74      8        0.8        V5 0.3715079
#> 75      8        0.8        V2 0.3606452
#> 76      8        0.8        V3 0.3669727
#> 77      8        0.8        V4 0.3604844
#> 78      8        0.8        V8 0.3894151
#> 79      8        0.8        V6 0.3371563
#> 80      8        0.8        V7 0.3740446
#> 81      9        0.9        V5 0.3709316
#> 82      9        0.9        V7 0.3733945
#> 83      9        0.9        V8 0.3888377
#> 84      9        0.9        V9 0.3356884
#> 85      9        0.9        V6 0.3366581
#> 86      9        0.9       V10 0.3845356
#> 87      9        0.9        V1 0.3417337
#> 88      9        0.9        V2 0.3601173
#> 89      9        0.9        V3 0.3663760
#> 90      9        0.9        V4 0.3599773
#> 91     10        1.0        V1 0.3411869
#> 92     10        1.0        V4 0.3594708
#> 93     10        1.0        V5 0.3703562
#> 94     10        1.0        V9 0.3351466
#> 95     10        1.0        V3 0.3657803
#> 96     10        1.0        V7 0.3727456
#> 97     10        1.0        V8 0.3882611
#> 98     10        1.0        V2 0.3595901
#> 99     10        1.0        V6 0.3361606
#> 100    10        1.0       V10 0.3839425

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3461393
#> 2       1        0.1        V5 0.3755669
#> 3       1        0.1        V9 0.3400546
#> 4       1        0.1        V3 0.3711769
#> 5       1        0.1        V4 0.3640544
#> 6       1        0.1        V8 0.3934812
#> 7       1        0.1        V2 0.3643627
#> 8       1        0.1        V6 0.3406648
#> 9       1        0.1        V7 0.3786267
#> 10      1        0.1       V10 0.3893133
#> 11      2        0.2        V7 0.3779687
#> 12      2        0.2        V8 0.3928977
#> 13      2        0.2        V9 0.3395057
#> 14      2        0.2       V10 0.3887129
#> 15      2        0.2        V1 0.3455855
#> 16      2        0.2        V5 0.3749844
#> 17      2        0.2        V2 0.3638293
#> 18      2        0.2        V3 0.3705733
#> 19      2        0.2        V4 0.3635422
#> 20      2        0.2        V6 0.3401614
#> 21      3        0.3        V4 0.3630308
#> 22      3        0.3        V5 0.3744027
#> 23      3        0.3        V3 0.3699708
#> 24      3        0.3        V7 0.3773118
#> 25      3        0.3        V8 0.3923151
#> 26      3        0.3        V9 0.3389577
#> 27      3        0.3        V6 0.3396587
#> 28      3        0.3       V10 0.3881133
#> 29      3        0.3        V1 0.3450326
#> 30      3        0.3        V2 0.3632967
#> 31      4        0.4        V1 0.3444806
#> 32      4        0.4        V9 0.3384106
#> 33      4        0.4        V3 0.3693692
#> 34      4        0.4        V4 0.3625201
#> 35      4        0.4        V5 0.3738219
#> 36      4        0.4        V2 0.3627648
#> 37      4        0.4        V6 0.3391568
#> 38      4        0.4        V7 0.3766561
#> 39      4        0.4        V8 0.3917334
#> 40      4        0.4       V10 0.3875147
#> 41      5        0.5        V8 0.3911526
#> 42      5        0.5        V9 0.3378644
#> 43      5        0.5       V10 0.3869171
#> 44      5        0.5        V1 0.3439294
#> 45      5        0.5        V5 0.3732421
#> 46      5        0.5        V2 0.3622338
#> 47      5        0.5        V3 0.3687686
#> 48      5        0.5        V4 0.3620101
#> 49      5        0.5        V6 0.3386555
#> 50      5        0.5        V7 0.3760015
#> 51      6        0.6        V4 0.3615008
#> 52      6        0.6        V5 0.3726631
#> 53      6        0.6        V7 0.3753481
#> 54      6        0.6        V8 0.3905726
#> 55      6        0.6        V9 0.3373191
#> 56      6        0.6        V6 0.3381551
#> 57      6        0.6       V10 0.3863203
#> 58      6        0.6        V1 0.3433792
#> 59      6        0.6        V2 0.3617035
#> 60      6        0.6        V3 0.3681690
#> 61      7        0.7        V1 0.3428298
#> 62      7        0.7        V3 0.3675704
#> 63      7        0.7        V4 0.3609922
#> 64      7        0.7        V5 0.3720850
#> 65      7        0.7        V9 0.3367747
#> 66      7        0.7        V6 0.3376553
#> 67      7        0.7        V7 0.3746958
#> 68      7        0.7        V8 0.3899934
#> 69      7        0.7        V2 0.3611740
#> 70      7        0.7       V10 0.3857245
#> 71      8        0.8        V9 0.3362311
#> 72      8        0.8       V10 0.3851295
#> 73      8        0.8        V1 0.3422813
#> 74      8        0.8        V5 0.3715079
#> 75      8        0.8        V2 0.3606452
#> 76      8        0.8        V3 0.3669727
#> 77      8        0.8        V4 0.3604844
#> 78      8        0.8        V8 0.3894151
#> 79      8        0.8        V6 0.3371563
#> 80      8        0.8        V7 0.3740446
#> 81      9        0.9        V5 0.3709316
#> 82      9        0.9        V7 0.3733945
#> 83      9        0.9        V8 0.3888377
#> 84      9        0.9        V9 0.3356884
#> 85      9        0.9        V6 0.3366581
#> 86      9        0.9       V10 0.3845356
#> 87      9        0.9        V1 0.3417337
#> 88      9        0.9        V2 0.3601173
#> 89      9        0.9        V3 0.3663760
#> 90      9        0.9        V4 0.3599773
#> 91     10        1.0        V1 0.3411869
#> 92     10        1.0        V4 0.3594708
#> 93     10        1.0        V5 0.3703562
#> 94     10        1.0        V9 0.3351466
#> 95     10        1.0        V3 0.3657803
#> 96     10        1.0        V7 0.3727456
#> 97     10        1.0        V8 0.3882611
#> 98     10        1.0        V2 0.3595901
#> 99     10        1.0        V6 0.3361606
#> 100    10        1.0       V10 0.3839425


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
#> 1       0                0          0 0.8757906 0.05050990 0.7929711 0.9554551
#> 2       0               20         20 0.8757906 0.05050990 0.7929711 0.9554551
#> 3       0               40         40 0.8757906 0.05050990 0.7929711 0.9554551
#> 4       0               60         60 0.8757906 0.05050990 0.7929711 0.9554551
#> 5      20                0         20 0.8617131 0.05091099 0.7789719 0.9422848
#> 6      20               20         40 0.8617131 0.05091099 0.7789719 0.9422848
#> 7      20               40         60 0.8617131 0.05091099 0.7789719 0.9422848
#> 8      20               60         80 0.8617131 0.05091099 0.7789719 0.9422848
#> 9      40                0         40 0.8478591 0.05114598 0.7653480 0.9287910
#> 10     40               20         60 0.8478591 0.05114598 0.7653480 0.9287910
#> 11     40               40         80 0.8478591 0.05114598 0.7653480 0.9287910
#> 12     40               60        100 0.8478591 0.05114598 0.7653480 0.9287910
#> 13     60                0         60 0.8342249 0.05125587 0.7520702 0.9151019
#> 14     60               20         80 0.8342249 0.05125587 0.7520702 0.9151019
#> 15     60               40        100 0.8342249 0.05125587 0.7520702 0.9151019
#> 16     60               60        120 0.8342249 0.05125587 0.7520702 0.9151019
#> 17     80                0         80 0.8208071 0.05127014 0.7391142 0.9013066
#> 18     80               20        100 0.8208071 0.05127014 0.7391142 0.9013066
#> 19     80               40        120 0.8208071 0.05127014 0.7391142 0.9013066
#> 20     80               60        140 0.8208071 0.05127014 0.7391142 0.9013066
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11772084 0.144342445 0.5440717
#> 2  0.30574618 0.11652729 0.103476409 0.4948492
#> 3  0.26001915 0.11454075 0.072867094 0.4503885
#> 4  0.22113100 0.11169759 0.050170260 0.4103044
#> 5  0.25589195 0.10948098 0.084695272 0.4478497
#> 6  0.21762106 0.10607344 0.058907021 0.4080171
#> 7  0.18507391 0.10239909 0.039935307 0.3721450
#> 8  0.15739448 0.09841131 0.026223865 0.3398526
#> 9  0.18213629 0.09943433 0.047216769 0.3700987
#> 10 0.15489621 0.09512019 0.031452222 0.3380106
#> 11 0.13173012 0.09079420 0.020206141 0.3091217
#> 12 0.11202872 0.08643070 0.012408348 0.2830998
#> 13 0.12963921 0.08869291 0.024474029 0.3074735
#> 14 0.11025053 0.08417442 0.015336908 0.2816145
#> 15 0.09376159 0.07977032 0.009130304 0.2582999
#> 16 0.07973872 0.07546636 0.005097750 0.2372537
#> 17 0.09227334 0.07806601 0.011443682 0.2569683
#> 18 0.07847305 0.07370721 0.006576179 0.2360508
#> 19 0.06673672 0.06952118 0.003513072 0.2171375
#> 20 0.05675565 0.06549668 0.001713207 0.2000046
```
