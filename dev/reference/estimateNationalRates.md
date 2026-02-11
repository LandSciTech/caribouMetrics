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
#> 1         0.1 0.3838513 0.01944594 0.3700670 0.4278760
#> 2         0.2 0.3832760 0.01942336 0.3694178 0.4271929
#> 3         0.3 0.3827015 0.01940095 0.3687697 0.4265109
#> 4         0.4 0.3821279 0.01937873 0.3681228 0.4258300
#> 5         0.5 0.3815551 0.01935669 0.3674750 0.4251502
#> 6         0.6 0.3809832 0.01933483 0.3668273 0.4244714
#> 7         0.7 0.3804122 0.01931314 0.3661807 0.4237938
#> 8         0.8 0.3798420 0.01929163 0.3655354 0.4231172
#> 9         0.9 0.3792726 0.01927030 0.3648911 0.4224417
#> 10        1.0 0.3787041 0.01924914 0.3642480 0.4217673

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3800414
#> 2       1        0.1        V5 0.3944177
#> 3       1        0.1        V9 0.3789580
#> 4       1        0.1        V3 0.3674857
#> 5       1        0.1        V4 0.4051277
#> 6       1        0.1        V8 0.4043673
#> 7       1        0.1        V2 0.3832472
#> 8       1        0.1        V6 0.3790039
#> 9       1        0.1        V7 0.4344804
#> 10      1        0.1       V10 0.3809434
#> 11      2        0.2        V7 0.4337697
#> 12      2        0.2        V8 0.4038099
#> 13      2        0.2        V9 0.3784125
#> 14      2        0.2       V10 0.3803939
#> 15      2        0.2        V1 0.3794698
#> 16      2        0.2        V5 0.3938945
#> 17      2        0.2        V2 0.3826576
#> 18      2        0.2        V3 0.3668064
#> 19      2        0.2        V4 0.4045398
#> 20      2        0.2        V6 0.3784446
#> 21      3        0.3        V4 0.4039527
#> 22      3        0.3        V5 0.3933719
#> 23      3        0.3        V3 0.3661283
#> 24      3        0.3        V7 0.4330601
#> 25      3        0.3        V8 0.4032533
#> 26      3        0.3        V9 0.3778678
#> 27      3        0.3        V6 0.3778861
#> 28      3        0.3       V10 0.3798452
#> 29      3        0.3        V1 0.3788990
#> 30      3        0.3        V2 0.3820690
#> 31      4        0.4        V1 0.3783291
#> 32      4        0.4        V9 0.3773239
#> 33      4        0.4        V3 0.3654515
#> 34      4        0.4        V4 0.4033665
#> 35      4        0.4        V5 0.3928501
#> 36      4        0.4        V2 0.3814813
#> 37      4        0.4        V6 0.3773285
#> 38      4        0.4        V7 0.4323517
#> 39      4        0.4        V8 0.4026975
#> 40      4        0.4       V10 0.3792973
#> 41      5        0.5        V8 0.4021425
#> 42      5        0.5        V9 0.3767807
#> 43      5        0.5       V10 0.3787502
#> 44      5        0.5        V1 0.3777600
#> 45      5        0.5        V5 0.3923289
#> 46      5        0.5        V2 0.3808944
#> 47      5        0.5        V3 0.3647759
#> 48      5        0.5        V4 0.4027812
#> 49      5        0.5        V6 0.3767717
#> 50      5        0.5        V7 0.4316444
#> 51      6        0.6        V4 0.4021967
#> 52      6        0.6        V5 0.3918085
#> 53      6        0.6        V7 0.4309383
#> 54      6        0.6        V8 0.4015882
#> 55      6        0.6        V9 0.3762383
#> 56      6        0.6        V6 0.3762157
#> 57      6        0.6       V10 0.3782039
#> 58      6        0.6        V1 0.3771918
#> 59      6        0.6        V2 0.3803085
#> 60      6        0.6        V3 0.3641016
#> 61      7        0.7        V1 0.3766245
#> 62      7        0.7        V3 0.3634286
#> 63      7        0.7        V4 0.4016130
#> 64      7        0.7        V5 0.3912887
#> 65      7        0.7        V9 0.3756967
#> 66      7        0.7        V6 0.3756605
#> 67      7        0.7        V7 0.4302334
#> 68      7        0.7        V8 0.4010347
#> 69      7        0.7        V2 0.3797235
#> 70      7        0.7       V10 0.3776584
#> 71      8        0.8        V9 0.3751559
#> 72      8        0.8       V10 0.3771137
#> 73      8        0.8        V1 0.3760580
#> 74      8        0.8        V5 0.3907696
#> 75      8        0.8        V2 0.3791394
#> 76      8        0.8        V3 0.3627568
#> 77      8        0.8        V4 0.4010302
#> 78      8        0.8        V8 0.4004819
#> 79      8        0.8        V6 0.3751061
#> 80      8        0.8        V7 0.4295296
#> 81      9        0.9        V5 0.3902512
#> 82      9        0.9        V7 0.4288269
#> 83      9        0.9        V8 0.3999299
#> 84      9        0.9        V9 0.3746159
#> 85      9        0.9        V6 0.3745526
#> 86      9        0.9       V10 0.3765697
#> 87      9        0.9        V1 0.3754924
#> 88      9        0.9        V2 0.3785561
#> 89      9        0.9        V3 0.3620862
#> 90      9        0.9        V4 0.4004483
#> 91     10        1.0        V1 0.3749276
#> 92     10        1.0        V4 0.3998671
#> 93     10        1.0        V5 0.3897335
#> 94     10        1.0        V9 0.3740766
#> 95     10        1.0        V3 0.3614168
#> 96     10        1.0        V7 0.4281254
#> 97     10        1.0        V8 0.3993787
#> 98     10        1.0        V2 0.3779738
#> 99     10        1.0        V6 0.3739998
#> 100    10        1.0       V10 0.3760265

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3800414
#> 2       1        0.1        V5 0.3944177
#> 3       1        0.1        V9 0.3789580
#> 4       1        0.1        V3 0.3674857
#> 5       1        0.1        V4 0.4051277
#> 6       1        0.1        V8 0.4043673
#> 7       1        0.1        V2 0.3832472
#> 8       1        0.1        V6 0.3790039
#> 9       1        0.1        V7 0.4344804
#> 10      1        0.1       V10 0.3809434
#> 11      2        0.2        V7 0.4337697
#> 12      2        0.2        V8 0.4038099
#> 13      2        0.2        V9 0.3784125
#> 14      2        0.2       V10 0.3803939
#> 15      2        0.2        V1 0.3794698
#> 16      2        0.2        V5 0.3938945
#> 17      2        0.2        V2 0.3826576
#> 18      2        0.2        V3 0.3668064
#> 19      2        0.2        V4 0.4045398
#> 20      2        0.2        V6 0.3784446
#> 21      3        0.3        V4 0.4039527
#> 22      3        0.3        V5 0.3933719
#> 23      3        0.3        V3 0.3661283
#> 24      3        0.3        V7 0.4330601
#> 25      3        0.3        V8 0.4032533
#> 26      3        0.3        V9 0.3778678
#> 27      3        0.3        V6 0.3778861
#> 28      3        0.3       V10 0.3798452
#> 29      3        0.3        V1 0.3788990
#> 30      3        0.3        V2 0.3820690
#> 31      4        0.4        V1 0.3783291
#> 32      4        0.4        V9 0.3773239
#> 33      4        0.4        V3 0.3654515
#> 34      4        0.4        V4 0.4033665
#> 35      4        0.4        V5 0.3928501
#> 36      4        0.4        V2 0.3814813
#> 37      4        0.4        V6 0.3773285
#> 38      4        0.4        V7 0.4323517
#> 39      4        0.4        V8 0.4026975
#> 40      4        0.4       V10 0.3792973
#> 41      5        0.5        V8 0.4021425
#> 42      5        0.5        V9 0.3767807
#> 43      5        0.5       V10 0.3787502
#> 44      5        0.5        V1 0.3777600
#> 45      5        0.5        V5 0.3923289
#> 46      5        0.5        V2 0.3808944
#> 47      5        0.5        V3 0.3647759
#> 48      5        0.5        V4 0.4027812
#> 49      5        0.5        V6 0.3767717
#> 50      5        0.5        V7 0.4316444
#> 51      6        0.6        V4 0.4021967
#> 52      6        0.6        V5 0.3918085
#> 53      6        0.6        V7 0.4309383
#> 54      6        0.6        V8 0.4015882
#> 55      6        0.6        V9 0.3762383
#> 56      6        0.6        V6 0.3762157
#> 57      6        0.6       V10 0.3782039
#> 58      6        0.6        V1 0.3771918
#> 59      6        0.6        V2 0.3803085
#> 60      6        0.6        V3 0.3641016
#> 61      7        0.7        V1 0.3766245
#> 62      7        0.7        V3 0.3634286
#> 63      7        0.7        V4 0.4016130
#> 64      7        0.7        V5 0.3912887
#> 65      7        0.7        V9 0.3756967
#> 66      7        0.7        V6 0.3756605
#> 67      7        0.7        V7 0.4302334
#> 68      7        0.7        V8 0.4010347
#> 69      7        0.7        V2 0.3797235
#> 70      7        0.7       V10 0.3776584
#> 71      8        0.8        V9 0.3751559
#> 72      8        0.8       V10 0.3771137
#> 73      8        0.8        V1 0.3760580
#> 74      8        0.8        V5 0.3907696
#> 75      8        0.8        V2 0.3791394
#> 76      8        0.8        V3 0.3627568
#> 77      8        0.8        V4 0.4010302
#> 78      8        0.8        V8 0.4004819
#> 79      8        0.8        V6 0.3751061
#> 80      8        0.8        V7 0.4295296
#> 81      9        0.9        V5 0.3902512
#> 82      9        0.9        V7 0.4288269
#> 83      9        0.9        V8 0.3999299
#> 84      9        0.9        V9 0.3746159
#> 85      9        0.9        V6 0.3745526
#> 86      9        0.9       V10 0.3765697
#> 87      9        0.9        V1 0.3754924
#> 88      9        0.9        V2 0.3785561
#> 89      9        0.9        V3 0.3620862
#> 90      9        0.9        V4 0.4004483
#> 91     10        1.0        V1 0.3749276
#> 92     10        1.0        V4 0.3998671
#> 93     10        1.0        V5 0.3897335
#> 94     10        1.0        V9 0.3740766
#> 95     10        1.0        V3 0.3614168
#> 96     10        1.0        V7 0.4281254
#> 97     10        1.0        V8 0.3993787
#> 98     10        1.0        V2 0.3779738
#> 99     10        1.0        V6 0.3739998
#> 100    10        1.0       V10 0.3760265


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
#> 1       0                0          0 0.8757906 0.04787619 0.7803809 0.9452473
#> 2       0               20         20 0.8757906 0.04787619 0.7803809 0.9452473
#> 3       0               40         40 0.8757906 0.04787619 0.7803809 0.9452473
#> 4       0               60         60 0.8757906 0.04787619 0.7803809 0.9452473
#> 5      20                0         20 0.8617131 0.04953776 0.7646692 0.9336509
#> 6      20               20         40 0.8617131 0.04953776 0.7646692 0.9336509
#> 7      20               40         60 0.8617131 0.04953776 0.7646692 0.9336509
#> 8      20               60         80 0.8617131 0.04953776 0.7646692 0.9336509
#> 9      40                0         40 0.8478591 0.05105123 0.7494294 0.9218947
#> 10     40               20         60 0.8478591 0.05105123 0.7494294 0.9218947
#> 11     40               40         80 0.8478591 0.05105123 0.7494294 0.9218947
#> 12     40               60        100 0.8478591 0.05105123 0.7494294 0.9218947
#> 13     60                0         60 0.8342249 0.05243728 0.7346222 0.9100425
#> 14     60               20         80 0.8342249 0.05243728 0.7346222 0.9100425
#> 15     60               40        100 0.8342249 0.05243728 0.7346222 0.9100425
#> 16     60               60        120 0.8342249 0.05243728 0.7346222 0.9100425
#> 17     80                0         80 0.8208071 0.05371206 0.7202156 0.8981417
#> 18     80               20        100 0.8208071 0.05371206 0.7202156 0.8981417
#> 19     80               40        120 0.8208071 0.05371206 0.7202156 0.8981417
#> 20     80               60        140 0.8208071 0.05371206 0.7202156 0.8981417
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12512664 0.148093237 0.5849012
#> 2  0.30574618 0.11848385 0.111602952 0.5345798
#> 3  0.26001915 0.11267109 0.083052631 0.4887818
#> 4  0.22113100 0.10726130 0.060859713 0.4472072
#> 5  0.25589195 0.11607953 0.084073556 0.4863983
#> 6  0.21762106 0.10848996 0.061649924 0.4450455
#> 7  0.18507391 0.10181981 0.044377539 0.4075680
#> 8  0.15739448 0.09575393 0.031241514 0.3736315
#> 9  0.18213629 0.10501455 0.044989109 0.4056205
#> 10 0.15489621 0.09754899 0.031703083 0.3718684
#> 11 0.13173012 0.09094256 0.021757400 0.3413152
#> 12 0.11202872 0.08496770 0.014466642 0.3136540
#> 13 0.12963921 0.09337496 0.022103599 0.3397280
#> 14 0.11025053 0.08650272 0.014717297 0.3122167
#> 15 0.09376159 0.08038888 0.009437913 0.2872961
#> 16 0.07973872 0.07486691 0.005786611 0.2647032
#> 17 0.09227334 0.08205102 0.009616721 0.2860006
#> 18 0.07847305 0.07595058 0.005907887 0.2635281
#> 19 0.06673672 0.07050130 0.003441341 0.2431299
#> 20 0.05675565 0.06557764 0.001881455 0.2245878
```
