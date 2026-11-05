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
[`addMissingYears()`](https://landscitech.github.io/caribouMetrics/dev/reference/addMissingYears.md),
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
[`getCaribouYear()`](https://landscitech.github.io/caribouMetrics/dev/reference/getCaribouYear.md),
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
#> 1         0.1 0.3838513 0.02151409 0.3556476 0.4133959
#> 2         0.2 0.3832760 0.02149922 0.3550959 0.4128648
#> 3         0.3 0.3827015 0.02148452 0.3545449 0.4123344
#> 4         0.4 0.3821279 0.02146998 0.3539949 0.4118047
#> 5         0.5 0.3815551 0.02145560 0.3534456 0.4112757
#> 6         0.6 0.3809832 0.02144138 0.3528973 0.4107473
#> 7         0.7 0.3804122 0.02142733 0.3523498 0.4102197
#> 8         0.8 0.3798420 0.02141343 0.3518031 0.4096927
#> 9         0.9 0.3792726 0.02139969 0.3512573 0.4091664
#> 10        1.0 0.3787041 0.02138610 0.3507124 0.4086408

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4149447
#> 2       1        0.1        V5 0.3529868
#> 3       1        0.1        V9 0.4080610
#> 4       1        0.1        V3 0.3768586
#> 5       1        0.1        V4 0.4023708
#> 6       1        0.1        V8 0.3656943
#> 7       1        0.1        V2 0.3787004
#> 8       1        0.1        V6 0.3745220
#> 9       1        0.1        V7 0.3648128
#> 10      1        0.1       V10 0.4054395
#> 11      2        0.2        V7 0.3643088
#> 12      2        0.2        V8 0.3651300
#> 13      2        0.2        V9 0.4074364
#> 14      2        0.2       V10 0.4047466
#> 15      2        0.2        V1 0.4144408
#> 16      2        0.2        V5 0.3524211
#> 17      2        0.2        V2 0.3780739
#> 18      2        0.2        V3 0.3763369
#> 19      2        0.2        V4 0.4018312
#> 20      2        0.2        V6 0.3739524
#> 21      3        0.3        V4 0.4012923
#> 22      3        0.3        V5 0.3518564
#> 23      3        0.3        V3 0.3758159
#> 24      3        0.3        V7 0.3638054
#> 25      3        0.3        V8 0.3645665
#> 26      3        0.3        V9 0.4068129
#> 27      3        0.3        V6 0.3733838
#> 28      3        0.3       V10 0.4040549
#> 29      3        0.3        V1 0.4139374
#> 30      3        0.3        V2 0.3774485
#> 31      4        0.4        V1 0.4134347
#> 32      4        0.4        V9 0.4061902
#> 33      4        0.4        V3 0.3752956
#> 34      4        0.4        V4 0.4007542
#> 35      4        0.4        V5 0.3512926
#> 36      4        0.4        V2 0.3768241
#> 37      4        0.4        V6 0.3728159
#> 38      4        0.4        V7 0.3633028
#> 39      4        0.4        V8 0.3640040
#> 40      4        0.4       V10 0.4033643
#> 41      5        0.5        V8 0.3634422
#> 42      5        0.5        V9 0.4055685
#> 43      5        0.5       V10 0.4026750
#> 44      5        0.5        V1 0.4129326
#> 45      5        0.5        V5 0.3507296
#> 46      5        0.5        V2 0.3762008
#> 47      5        0.5        V3 0.3747761
#> 48      5        0.5        V4 0.4002168
#> 49      5        0.5        V6 0.3722490
#> 50      5        0.5        V7 0.3628008
#> 51      6        0.6        V4 0.3996800
#> 52      6        0.6        V5 0.3501676
#> 53      6        0.6        V7 0.3622996
#> 54      6        0.6        V8 0.3628814
#> 55      6        0.6        V9 0.4049478
#> 56      6        0.6        V6 0.3716829
#> 57      6        0.6       V10 0.4019868
#> 58      6        0.6        V1 0.4124311
#> 59      6        0.6        V2 0.3755785
#> 60      6        0.6        V3 0.3742572
#> 61      7        0.7        V1 0.4119302
#> 62      7        0.7        V3 0.3737391
#> 63      7        0.7        V4 0.3991441
#> 64      7        0.7        V5 0.3496065
#> 65      7        0.7        V9 0.4043280
#> 66      7        0.7        V6 0.3711176
#> 67      7        0.7        V7 0.3617990
#> 68      7        0.7        V8 0.3623214
#> 69      7        0.7        V2 0.3749572
#> 70      7        0.7       V10 0.4012998
#> 71      8        0.8        V9 0.4037092
#> 72      8        0.8       V10 0.4006139
#> 73      8        0.8        V1 0.4114299
#> 74      8        0.8        V5 0.3490462
#> 75      8        0.8        V2 0.3743369
#> 76      8        0.8        V3 0.3732217
#> 77      8        0.8        V4 0.3986088
#> 78      8        0.8        V8 0.3617623
#> 79      8        0.8        V6 0.3705532
#> 80      8        0.8        V7 0.3612991
#> 81      9        0.9        V5 0.3484869
#> 82      9        0.9        V7 0.3607999
#> 83      9        0.9        V8 0.3612041
#> 84      9        0.9        V9 0.4030913
#> 85      9        0.9        V6 0.3699897
#> 86      9        0.9       V10 0.3999292
#> 87      9        0.9        V1 0.4109302
#> 88      9        0.9        V2 0.3737177
#> 89      9        0.9        V3 0.3727050
#> 90      9        0.9        V4 0.3980742
#> 91     10        1.0        V1 0.4104311
#> 92     10        1.0        V4 0.3975404
#> 93     10        1.0        V5 0.3479285
#> 94     10        1.0        V9 0.4024744
#> 95     10        1.0        V3 0.3721891
#> 96     10        1.0        V7 0.3603014
#> 97     10        1.0        V8 0.3606467
#> 98     10        1.0        V2 0.3730995
#> 99     10        1.0        V6 0.3694270
#> 100    10        1.0       V10 0.3992457

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4149447
#> 2       1        0.1        V5 0.3529868
#> 3       1        0.1        V9 0.4080610
#> 4       1        0.1        V3 0.3768586
#> 5       1        0.1        V4 0.4023708
#> 6       1        0.1        V8 0.3656943
#> 7       1        0.1        V2 0.3787004
#> 8       1        0.1        V6 0.3745220
#> 9       1        0.1        V7 0.3648128
#> 10      1        0.1       V10 0.4054395
#> 11      2        0.2        V7 0.3643088
#> 12      2        0.2        V8 0.3651300
#> 13      2        0.2        V9 0.4074364
#> 14      2        0.2       V10 0.4047466
#> 15      2        0.2        V1 0.4144408
#> 16      2        0.2        V5 0.3524211
#> 17      2        0.2        V2 0.3780739
#> 18      2        0.2        V3 0.3763369
#> 19      2        0.2        V4 0.4018312
#> 20      2        0.2        V6 0.3739524
#> 21      3        0.3        V4 0.4012923
#> 22      3        0.3        V5 0.3518564
#> 23      3        0.3        V3 0.3758159
#> 24      3        0.3        V7 0.3638054
#> 25      3        0.3        V8 0.3645665
#> 26      3        0.3        V9 0.4068129
#> 27      3        0.3        V6 0.3733838
#> 28      3        0.3       V10 0.4040549
#> 29      3        0.3        V1 0.4139374
#> 30      3        0.3        V2 0.3774485
#> 31      4        0.4        V1 0.4134347
#> 32      4        0.4        V9 0.4061902
#> 33      4        0.4        V3 0.3752956
#> 34      4        0.4        V4 0.4007542
#> 35      4        0.4        V5 0.3512926
#> 36      4        0.4        V2 0.3768241
#> 37      4        0.4        V6 0.3728159
#> 38      4        0.4        V7 0.3633028
#> 39      4        0.4        V8 0.3640040
#> 40      4        0.4       V10 0.4033643
#> 41      5        0.5        V8 0.3634422
#> 42      5        0.5        V9 0.4055685
#> 43      5        0.5       V10 0.4026750
#> 44      5        0.5        V1 0.4129326
#> 45      5        0.5        V5 0.3507296
#> 46      5        0.5        V2 0.3762008
#> 47      5        0.5        V3 0.3747761
#> 48      5        0.5        V4 0.4002168
#> 49      5        0.5        V6 0.3722490
#> 50      5        0.5        V7 0.3628008
#> 51      6        0.6        V4 0.3996800
#> 52      6        0.6        V5 0.3501676
#> 53      6        0.6        V7 0.3622996
#> 54      6        0.6        V8 0.3628814
#> 55      6        0.6        V9 0.4049478
#> 56      6        0.6        V6 0.3716829
#> 57      6        0.6       V10 0.4019868
#> 58      6        0.6        V1 0.4124311
#> 59      6        0.6        V2 0.3755785
#> 60      6        0.6        V3 0.3742572
#> 61      7        0.7        V1 0.4119302
#> 62      7        0.7        V3 0.3737391
#> 63      7        0.7        V4 0.3991441
#> 64      7        0.7        V5 0.3496065
#> 65      7        0.7        V9 0.4043280
#> 66      7        0.7        V6 0.3711176
#> 67      7        0.7        V7 0.3617990
#> 68      7        0.7        V8 0.3623214
#> 69      7        0.7        V2 0.3749572
#> 70      7        0.7       V10 0.4012998
#> 71      8        0.8        V9 0.4037092
#> 72      8        0.8       V10 0.4006139
#> 73      8        0.8        V1 0.4114299
#> 74      8        0.8        V5 0.3490462
#> 75      8        0.8        V2 0.3743369
#> 76      8        0.8        V3 0.3732217
#> 77      8        0.8        V4 0.3986088
#> 78      8        0.8        V8 0.3617623
#> 79      8        0.8        V6 0.3705532
#> 80      8        0.8        V7 0.3612991
#> 81      9        0.9        V5 0.3484869
#> 82      9        0.9        V7 0.3607999
#> 83      9        0.9        V8 0.3612041
#> 84      9        0.9        V9 0.4030913
#> 85      9        0.9        V6 0.3699897
#> 86      9        0.9       V10 0.3999292
#> 87      9        0.9        V1 0.4109302
#> 88      9        0.9        V2 0.3737177
#> 89      9        0.9        V3 0.3727050
#> 90      9        0.9        V4 0.3980742
#> 91     10        1.0        V1 0.4104311
#> 92     10        1.0        V4 0.3975404
#> 93     10        1.0        V5 0.3479285
#> 94     10        1.0        V9 0.4024744
#> 95     10        1.0        V3 0.3721891
#> 96     10        1.0        V7 0.3603014
#> 97     10        1.0        V8 0.3606467
#> 98     10        1.0        V2 0.3730995
#> 99     10        1.0        V6 0.3694270
#> 100    10        1.0       V10 0.3992457


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
#> 1       0                0          0 0.8757906 0.05152284 0.7655996 0.9516038
#> 2       0               20         20 0.8757906 0.05152284 0.7655996 0.9516038
#> 3       0               40         40 0.8757906 0.05152284 0.7655996 0.9516038
#> 4       0               60         60 0.8757906 0.05152284 0.7655996 0.9516038
#> 5      20                0         20 0.8617131 0.05374186 0.7489903 0.9397324
#> 6      20               20         40 0.8617131 0.05374186 0.7489903 0.9397324
#> 7      20               40         60 0.8617131 0.05374186 0.7489903 0.9397324
#> 8      20               60         80 0.8617131 0.05374186 0.7489903 0.9397324
#> 9      40                0         40 0.8478591 0.05580098 0.7329124 0.9276135
#> 10     40               20         60 0.8478591 0.05580098 0.7329124 0.9276135
#> 11     40               40         80 0.8478591 0.05580098 0.7329124 0.9276135
#> 12     40               60        100 0.8478591 0.05580098 0.7329124 0.9276135
#> 13     60                0         60 0.8342249 0.05771802 0.7173205 0.9153359
#> 14     60               20         80 0.8342249 0.05771802 0.7173205 0.9153359
#> 15     60               40        100 0.8342249 0.05771802 0.7173205 0.9153359
#> 16     60               60        120 0.8342249 0.05771802 0.7173205 0.9153359
#> 17     80                0         80 0.8208071 0.05950739 0.7021777 0.9029644
#> 18     80               20        100 0.8208071 0.05950739 0.7021777 0.9029644
#> 19     80               40        120 0.8208071 0.05950739 0.7021777 0.9029644
#> 20     80               60        140 0.8208071 0.05950739 0.7021777 0.9029644
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12553327 0.156885534 0.5850151
#> 2  0.30574618 0.11871037 0.126347172 0.5143330
#> 3  0.26001915 0.11199343 0.101073782 0.4526599
#> 4  0.22113100 0.10534459 0.080220687 0.3990609
#> 5  0.25589195 0.10910922 0.089716765 0.4641343
#> 6  0.21762106 0.10185713 0.070879028 0.4090229
#> 7  0.18507391 0.09502712 0.055442572 0.3612022
#> 8  0.15739448 0.08852728 0.042870913 0.3197319
#> 9  0.18213629 0.09371822 0.048572765 0.3700875
#> 10 0.15489621 0.08676952 0.037307936 0.3274378
#> 11 0.13173012 0.08033227 0.028245702 0.2904370
#> 12 0.11202872 0.07431095 0.021030694 0.2582929
#> 13 0.12963921 0.07987739 0.024280575 0.2973150
#> 14 0.11025053 0.07353876 0.017903591 0.2642728
#> 15 0.09376159 0.06770551 0.012928049 0.2355179
#> 16 0.07973872 0.06229311 0.009110984 0.2104216
#> 17 0.09227334 0.06771391 0.010810953 0.2408720
#> 18 0.07847305 0.06208553 0.007511666 0.2151009
#> 19 0.06673672 0.05691773 0.005061965 0.1925439
#> 20 0.05675565 0.05214013 0.003291847 0.1727154
```
