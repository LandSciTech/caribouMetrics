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
[`bbouNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/bbouNationalPriors.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/dev/reference/trajectoriesFromSummary.md)

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
#> 1         0.1 0.3838513 0.02595484 0.3477194 0.4225883
#> 2         0.2 0.3832760 0.02590916 0.3472202 0.4219611
#> 3         0.3 0.3827015 0.02586362 0.3467218 0.4213348
#> 4         0.4 0.3821279 0.02581823 0.3462240 0.4207094
#> 5         0.5 0.3815551 0.02577298 0.3457269 0.4200850
#> 6         0.6 0.3809832 0.02572787 0.3452306 0.4194614
#> 7         0.7 0.3804122 0.02568291 0.3447350 0.4188389
#> 8         0.8 0.3798420 0.02563809 0.3442400 0.4182172
#> 9         0.9 0.3792726 0.02559341 0.3437458 0.4175964
#> 10        1.0 0.3787041 0.02554887 0.3432523 0.4169766

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3850991
#> 2       1        0.1        V5 0.3784574
#> 3       1        0.1        V9 0.3479168
#> 4       1        0.1        V3 0.4209873
#> 5       1        0.1        V4 0.3749622
#> 6       1        0.1        V8 0.3920063
#> 7       1        0.1        V2 0.3476621
#> 8       1        0.1        V6 0.3941711
#> 9       1        0.1        V7 0.4230532
#> 10      1        0.1       V10 0.4023645
#> 11      2        0.2        V7 0.4224269
#> 12      2        0.2        V8 0.3914388
#> 13      2        0.2        V9 0.3474058
#> 14      2        0.2       V10 0.4018263
#> 15      2        0.2        V1 0.3845717
#> 16      2        0.2        V5 0.3780214
#> 17      2        0.2        V2 0.3471664
#> 18      2        0.2        V3 0.4203566
#> 19      2        0.2        V4 0.3744290
#> 20      2        0.2        V6 0.3935632
#> 21      3        0.3        V4 0.3738967
#> 22      3        0.3        V5 0.3775859
#> 23      3        0.3        V3 0.4197269
#> 24      3        0.3        V7 0.4218016
#> 25      3        0.3        V8 0.3908721
#> 26      3        0.3        V9 0.3468955
#> 27      3        0.3        V6 0.3929563
#> 28      3        0.3       V10 0.4012889
#> 29      3        0.3        V1 0.3840451
#> 30      3        0.3        V2 0.3466713
#> 31      4        0.4        V1 0.3835191
#> 32      4        0.4        V9 0.3463859
#> 33      4        0.4        V3 0.4190981
#> 34      4        0.4        V4 0.3733651
#> 35      4        0.4        V5 0.3771509
#> 36      4        0.4        V2 0.3461770
#> 37      4        0.4        V6 0.3923503
#> 38      4        0.4        V7 0.4211772
#> 39      4        0.4        V8 0.3903063
#> 40      4        0.4       V10 0.4007522
#> 41      5        0.5        V8 0.3897413
#> 42      5        0.5        V9 0.3458771
#> 43      5        0.5       V10 0.4002162
#> 44      5        0.5        V1 0.3829939
#> 45      5        0.5        V5 0.3767164
#> 46      5        0.5        V2 0.3456833
#> 47      5        0.5        V3 0.4184702
#> 48      5        0.5        V4 0.3728342
#> 49      5        0.5        V6 0.3917452
#> 50      5        0.5        V7 0.4205538
#> 51      6        0.6        V4 0.3723041
#> 52      6        0.6        V5 0.3762825
#> 53      6        0.6        V7 0.4199312
#> 54      6        0.6        V8 0.3891771
#> 55      6        0.6        V9 0.3453690
#> 56      6        0.6        V6 0.3911411
#> 57      6        0.6       V10 0.3996810
#> 58      6        0.6        V1 0.3824695
#> 59      6        0.6        V2 0.3451904
#> 60      6        0.6        V3 0.4178433
#> 61      7        0.7        V1 0.3819457
#> 62      7        0.7        V3 0.4172174
#> 63      7        0.7        V4 0.3717748
#> 64      7        0.7        V5 0.3758490
#> 65      7        0.7        V9 0.3448617
#> 66      7        0.7        V6 0.3905379
#> 67      7        0.7        V7 0.4193096
#> 68      7        0.7        V8 0.3886138
#> 69      7        0.7        V2 0.3446982
#> 70      7        0.7       V10 0.3991464
#> 71      8        0.8        V9 0.3443552
#> 72      8        0.8       V10 0.3986126
#> 73      8        0.8        V1 0.3814227
#> 74      8        0.8        V5 0.3754160
#> 75      8        0.8        V2 0.3442066
#> 76      8        0.8        V3 0.4165923
#> 77      8        0.8        V4 0.3712462
#> 78      8        0.8        V8 0.3880512
#> 79      8        0.8        V6 0.3899356
#> 80      8        0.8        V7 0.4186889
#> 81      9        0.9        V5 0.3749835
#> 82      9        0.9        V7 0.4180691
#> 83      9        0.9        V8 0.3874895
#> 84      9        0.9        V9 0.3438493
#> 85      9        0.9        V6 0.3893342
#> 86      9        0.9       V10 0.3980795
#> 87      9        0.9        V1 0.3809003
#> 88      9        0.9        V2 0.3437158
#> 89      9        0.9        V3 0.4159682
#> 90      9        0.9        V4 0.3707184
#> 91     10        1.0        V1 0.3803787
#> 92     10        1.0        V4 0.3701913
#> 93     10        1.0        V5 0.3745515
#> 94     10        1.0        V9 0.3433442
#> 95     10        1.0        V3 0.4153451
#> 96     10        1.0        V7 0.4174503
#> 97     10        1.0        V8 0.3869285
#> 98     10        1.0        V2 0.3432256
#> 99     10        1.0        V6 0.3887338
#> 100    10        1.0       V10 0.3975471

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3850991
#> 2       1        0.1        V5 0.3784574
#> 3       1        0.1        V9 0.3479168
#> 4       1        0.1        V3 0.4209873
#> 5       1        0.1        V4 0.3749622
#> 6       1        0.1        V8 0.3920063
#> 7       1        0.1        V2 0.3476621
#> 8       1        0.1        V6 0.3941711
#> 9       1        0.1        V7 0.4230532
#> 10      1        0.1       V10 0.4023645
#> 11      2        0.2        V7 0.4224269
#> 12      2        0.2        V8 0.3914388
#> 13      2        0.2        V9 0.3474058
#> 14      2        0.2       V10 0.4018263
#> 15      2        0.2        V1 0.3845717
#> 16      2        0.2        V5 0.3780214
#> 17      2        0.2        V2 0.3471664
#> 18      2        0.2        V3 0.4203566
#> 19      2        0.2        V4 0.3744290
#> 20      2        0.2        V6 0.3935632
#> 21      3        0.3        V4 0.3738967
#> 22      3        0.3        V5 0.3775859
#> 23      3        0.3        V3 0.4197269
#> 24      3        0.3        V7 0.4218016
#> 25      3        0.3        V8 0.3908721
#> 26      3        0.3        V9 0.3468955
#> 27      3        0.3        V6 0.3929563
#> 28      3        0.3       V10 0.4012889
#> 29      3        0.3        V1 0.3840451
#> 30      3        0.3        V2 0.3466713
#> 31      4        0.4        V1 0.3835191
#> 32      4        0.4        V9 0.3463859
#> 33      4        0.4        V3 0.4190981
#> 34      4        0.4        V4 0.3733651
#> 35      4        0.4        V5 0.3771509
#> 36      4        0.4        V2 0.3461770
#> 37      4        0.4        V6 0.3923503
#> 38      4        0.4        V7 0.4211772
#> 39      4        0.4        V8 0.3903063
#> 40      4        0.4       V10 0.4007522
#> 41      5        0.5        V8 0.3897413
#> 42      5        0.5        V9 0.3458771
#> 43      5        0.5       V10 0.4002162
#> 44      5        0.5        V1 0.3829939
#> 45      5        0.5        V5 0.3767164
#> 46      5        0.5        V2 0.3456833
#> 47      5        0.5        V3 0.4184702
#> 48      5        0.5        V4 0.3728342
#> 49      5        0.5        V6 0.3917452
#> 50      5        0.5        V7 0.4205538
#> 51      6        0.6        V4 0.3723041
#> 52      6        0.6        V5 0.3762825
#> 53      6        0.6        V7 0.4199312
#> 54      6        0.6        V8 0.3891771
#> 55      6        0.6        V9 0.3453690
#> 56      6        0.6        V6 0.3911411
#> 57      6        0.6       V10 0.3996810
#> 58      6        0.6        V1 0.3824695
#> 59      6        0.6        V2 0.3451904
#> 60      6        0.6        V3 0.4178433
#> 61      7        0.7        V1 0.3819457
#> 62      7        0.7        V3 0.4172174
#> 63      7        0.7        V4 0.3717748
#> 64      7        0.7        V5 0.3758490
#> 65      7        0.7        V9 0.3448617
#> 66      7        0.7        V6 0.3905379
#> 67      7        0.7        V7 0.4193096
#> 68      7        0.7        V8 0.3886138
#> 69      7        0.7        V2 0.3446982
#> 70      7        0.7       V10 0.3991464
#> 71      8        0.8        V9 0.3443552
#> 72      8        0.8       V10 0.3986126
#> 73      8        0.8        V1 0.3814227
#> 74      8        0.8        V5 0.3754160
#> 75      8        0.8        V2 0.3442066
#> 76      8        0.8        V3 0.4165923
#> 77      8        0.8        V4 0.3712462
#> 78      8        0.8        V8 0.3880512
#> 79      8        0.8        V6 0.3899356
#> 80      8        0.8        V7 0.4186889
#> 81      9        0.9        V5 0.3749835
#> 82      9        0.9        V7 0.4180691
#> 83      9        0.9        V8 0.3874895
#> 84      9        0.9        V9 0.3438493
#> 85      9        0.9        V6 0.3893342
#> 86      9        0.9       V10 0.3980795
#> 87      9        0.9        V1 0.3809003
#> 88      9        0.9        V2 0.3437158
#> 89      9        0.9        V3 0.4159682
#> 90      9        0.9        V4 0.3707184
#> 91     10        1.0        V1 0.3803787
#> 92     10        1.0        V4 0.3701913
#> 93     10        1.0        V5 0.3745515
#> 94     10        1.0        V9 0.3433442
#> 95     10        1.0        V3 0.4153451
#> 96     10        1.0        V7 0.4174503
#> 97     10        1.0        V8 0.3869285
#> 98     10        1.0        V2 0.3432256
#> 99     10        1.0        V6 0.3887338
#> 100    10        1.0       V10 0.3975471


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
#> 1       0                0          0 0.8757906 0.04106652 0.8018306 0.9398397
#> 2       0               20         20 0.8757906 0.04106652 0.8018306 0.9398397
#> 3       0               40         40 0.8757906 0.04106652 0.8018306 0.9398397
#> 4       0               60         60 0.8757906 0.04106652 0.8018306 0.9398397
#> 5      20                0         20 0.8617131 0.04364739 0.7818989 0.9302220
#> 6      20               20         40 0.8617131 0.04364739 0.7818989 0.9302220
#> 7      20               40         60 0.8617131 0.04364739 0.7818989 0.9302220
#> 8      20               60         80 0.8617131 0.04364739 0.7818989 0.9302220
#> 9      40                0         40 0.8478591 0.04604521 0.7627256 0.9204931
#> 10     40               20         60 0.8478591 0.04604521 0.7627256 0.9204931
#> 11     40               40         80 0.8478591 0.04604521 0.7627256 0.9204931
#> 12     40               60        100 0.8478591 0.04604521 0.7627256 0.9204931
#> 13     60                0         60 0.8342249 0.04827972 0.7442285 0.9106886
#> 14     60               20         80 0.8342249 0.04827972 0.7442285 0.9106886
#> 15     60               40        100 0.8342249 0.04827972 0.7442285 0.9106886
#> 16     60               60        120 0.8342249 0.04827972 0.7442285 0.9106886
#> 17     80                0         80 0.8208071 0.05036716 0.7263444 0.9008369
#> 18     80               20        100 0.8208071 0.05036716 0.7263444 0.9008369
#> 19     80               40        120 0.8208071 0.05036716 0.7263444 0.9008369
#> 20     80               60        140 0.8208071 0.05036716 0.7263444 0.9008369
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13896094 0.148566874 0.5951615
#> 2  0.30574618 0.13585154 0.117295420 0.5503457
#> 3  0.26001915 0.13216999 0.091834421 0.5090091
#> 4  0.22113100 0.12794468 0.071188690 0.4709677
#> 5  0.25589195 0.12333681 0.088262272 0.4808906
#> 6  0.21762106 0.11899319 0.068301111 0.4451262
#> 7  0.18507391 0.11447575 0.052223112 0.4122887
#> 8  0.15739448 0.10978260 0.039371392 0.3821592
#> 9  0.18213629 0.10813414 0.049985261 0.3900117
#> 10 0.15489621 0.10336225 0.037592237 0.3617260
#> 11 0.13173012 0.09863244 0.027797787 0.3357838
#> 12 0.11202872 0.09392846 0.020152079 0.3119880
#> 13 0.12963921 0.09387510 0.026452608 0.3181901
#> 14 0.11025053 0.08914365 0.019111033 0.2958455
#> 15 0.09376159 0.08456028 0.013482022 0.2753349
#> 16 0.07973872 0.08010342 0.009248718 0.2564949
#> 17 0.09227334 0.08086670 0.012725119 0.2614085
#> 18 0.07847305 0.07642224 0.008687134 0.2436935
#> 19 0.06673672 0.07216887 0.005733928 0.2273955
#> 20 0.05675565 0.06808423 0.003637677 0.2123838
```
