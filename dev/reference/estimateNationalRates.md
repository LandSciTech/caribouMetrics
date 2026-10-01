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
#> 1         0.1 0.3838513 0.03284781 0.3255800 0.4228132
#> 2         0.2 0.3832760 0.03280383 0.3250772 0.4221376
#> 3         0.3 0.3827015 0.03276005 0.3245753 0.4214631
#> 4         0.4 0.3821279 0.03271644 0.3240741 0.4207897
#> 5         0.5 0.3815551 0.03267302 0.3235737 0.4201173
#> 6         0.6 0.3809832 0.03262979 0.3230740 0.4194461
#> 7         0.7 0.3804122 0.03258674 0.3225752 0.4187759
#> 8         0.8 0.3798420 0.03254387 0.3220771 0.4181068
#> 9         0.9 0.3792726 0.03250118 0.3215797 0.4174387
#> 10        1.0 0.3787041 0.03245867 0.3210832 0.4167718

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3467345
#> 2       1        0.1        V5 0.3957949
#> 3       1        0.1        V9 0.4182675
#> 4       1        0.1        V3 0.3981073
#> 5       1        0.1        V4 0.4146378
#> 6       1        0.1        V8 0.3194384
#> 7       1        0.1        V2 0.3954128
#> 8       1        0.1        V6 0.3827003
#> 9       1        0.1        V7 0.3762855
#> 10      1        0.1       V10 0.4241329
#> 11      2        0.2        V7 0.3755909
#> 12      2        0.2        V8 0.3189556
#> 13      2        0.2        V9 0.4177001
#> 14      2        0.2       V10 0.4234259
#> 15      2        0.2        V1 0.3461629
#> 16      2        0.2        V5 0.3952378
#> 17      2        0.2        V2 0.3948370
#> 18      2        0.2        V3 0.3973916
#> 19      2        0.2        V4 0.4139971
#> 20      2        0.2        V6 0.3821447
#> 21      3        0.3        V4 0.4133573
#> 22      3        0.3        V5 0.3946816
#> 23      3        0.3        V3 0.3966772
#> 24      3        0.3        V7 0.3748977
#> 25      3        0.3        V8 0.3184735
#> 26      3        0.3        V9 0.4171334
#> 27      3        0.3        V6 0.3815899
#> 28      3        0.3       V10 0.4227201
#> 29      3        0.3        V1 0.3455924
#> 30      3        0.3        V2 0.3942621
#> 31      4        0.4        V1 0.3450228
#> 32      4        0.4        V9 0.4165676
#> 33      4        0.4        V3 0.3959641
#> 34      4        0.4        V4 0.4127185
#> 35      4        0.4        V5 0.3941261
#> 36      4        0.4        V2 0.3936880
#> 37      4        0.4        V6 0.3810359
#> 38      4        0.4        V7 0.3742057
#> 39      4        0.4        V8 0.3179922
#> 40      4        0.4       V10 0.4220154
#> 41      5        0.5        V8 0.3175116
#> 42      5        0.5        V9 0.4160025
#> 43      5        0.5       V10 0.4213120
#> 44      5        0.5        V1 0.3444541
#> 45      5        0.5        V5 0.3935714
#> 46      5        0.5        V2 0.3931147
#> 47      5        0.5        V3 0.3952523
#> 48      5        0.5        V4 0.4120807
#> 49      5        0.5        V6 0.3804827
#> 50      5        0.5        V7 0.3735150
#> 51      6        0.6        V4 0.4114439
#> 52      6        0.6        V5 0.3930175
#> 53      6        0.6        V7 0.3728256
#> 54      6        0.6        V8 0.3170318
#> 55      6        0.6        V9 0.4154381
#> 56      6        0.6        V6 0.3799303
#> 57      6        0.6       V10 0.4206097
#> 58      6        0.6        V1 0.3438863
#> 59      6        0.6        V2 0.3925423
#> 60      6        0.6        V3 0.3945418
#> 61      7        0.7        V1 0.3433195
#> 62      7        0.7        V3 0.3938325
#> 63      7        0.7        V4 0.4108081
#> 64      7        0.7        V5 0.3924643
#> 65      7        0.7        V9 0.4148746
#> 66      7        0.7        V6 0.3793787
#> 67      7        0.7        V7 0.3721374
#> 68      7        0.7        V8 0.3165526
#> 69      7        0.7        V2 0.3919707
#> 70      7        0.7       V10 0.4199085
#> 71      8        0.8        V9 0.4143117
#> 72      8        0.8       V10 0.4192086
#> 73      8        0.8        V1 0.3427536
#> 74      8        0.8        V5 0.3919120
#> 75      8        0.8        V2 0.3914000
#> 76      8        0.8        V3 0.3931246
#> 77      8        0.8        V4 0.4101732
#> 78      8        0.8        V8 0.3160742
#> 79      8        0.8        V6 0.3788279
#> 80      8        0.8        V7 0.3714506
#> 81      9        0.9        V5 0.3913604
#> 82      9        0.9        V7 0.3707650
#> 83      9        0.9        V8 0.3155965
#> 84      9        0.9        V9 0.4137497
#> 85      9        0.9        V6 0.3782779
#> 86      9        0.9       V10 0.4185098
#> 87      9        0.9        V1 0.3421887
#> 88      9        0.9        V2 0.3908301
#> 89      9        0.9        V3 0.3924179
#> 90      9        0.9        V4 0.4095394
#> 91     10        1.0        V1 0.3416247
#> 92     10        1.0        V4 0.4089065
#> 93     10        1.0        V5 0.3908096
#> 94     10        1.0        V9 0.4131884
#> 95     10        1.0        V3 0.3917124
#> 96     10        1.0        V7 0.3700806
#> 97     10        1.0        V8 0.3151195
#> 98     10        1.0        V2 0.3902610
#> 99     10        1.0        V6 0.3777286
#> 100    10        1.0       V10 0.4178121

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3467345
#> 2       1        0.1        V5 0.3957949
#> 3       1        0.1        V9 0.4182675
#> 4       1        0.1        V3 0.3981073
#> 5       1        0.1        V4 0.4146378
#> 6       1        0.1        V8 0.3194384
#> 7       1        0.1        V2 0.3954128
#> 8       1        0.1        V6 0.3827003
#> 9       1        0.1        V7 0.3762855
#> 10      1        0.1       V10 0.4241329
#> 11      2        0.2        V7 0.3755909
#> 12      2        0.2        V8 0.3189556
#> 13      2        0.2        V9 0.4177001
#> 14      2        0.2       V10 0.4234259
#> 15      2        0.2        V1 0.3461629
#> 16      2        0.2        V5 0.3952378
#> 17      2        0.2        V2 0.3948370
#> 18      2        0.2        V3 0.3973916
#> 19      2        0.2        V4 0.4139971
#> 20      2        0.2        V6 0.3821447
#> 21      3        0.3        V4 0.4133573
#> 22      3        0.3        V5 0.3946816
#> 23      3        0.3        V3 0.3966772
#> 24      3        0.3        V7 0.3748977
#> 25      3        0.3        V8 0.3184735
#> 26      3        0.3        V9 0.4171334
#> 27      3        0.3        V6 0.3815899
#> 28      3        0.3       V10 0.4227201
#> 29      3        0.3        V1 0.3455924
#> 30      3        0.3        V2 0.3942621
#> 31      4        0.4        V1 0.3450228
#> 32      4        0.4        V9 0.4165676
#> 33      4        0.4        V3 0.3959641
#> 34      4        0.4        V4 0.4127185
#> 35      4        0.4        V5 0.3941261
#> 36      4        0.4        V2 0.3936880
#> 37      4        0.4        V6 0.3810359
#> 38      4        0.4        V7 0.3742057
#> 39      4        0.4        V8 0.3179922
#> 40      4        0.4       V10 0.4220154
#> 41      5        0.5        V8 0.3175116
#> 42      5        0.5        V9 0.4160025
#> 43      5        0.5       V10 0.4213120
#> 44      5        0.5        V1 0.3444541
#> 45      5        0.5        V5 0.3935714
#> 46      5        0.5        V2 0.3931147
#> 47      5        0.5        V3 0.3952523
#> 48      5        0.5        V4 0.4120807
#> 49      5        0.5        V6 0.3804827
#> 50      5        0.5        V7 0.3735150
#> 51      6        0.6        V4 0.4114439
#> 52      6        0.6        V5 0.3930175
#> 53      6        0.6        V7 0.3728256
#> 54      6        0.6        V8 0.3170318
#> 55      6        0.6        V9 0.4154381
#> 56      6        0.6        V6 0.3799303
#> 57      6        0.6       V10 0.4206097
#> 58      6        0.6        V1 0.3438863
#> 59      6        0.6        V2 0.3925423
#> 60      6        0.6        V3 0.3945418
#> 61      7        0.7        V1 0.3433195
#> 62      7        0.7        V3 0.3938325
#> 63      7        0.7        V4 0.4108081
#> 64      7        0.7        V5 0.3924643
#> 65      7        0.7        V9 0.4148746
#> 66      7        0.7        V6 0.3793787
#> 67      7        0.7        V7 0.3721374
#> 68      7        0.7        V8 0.3165526
#> 69      7        0.7        V2 0.3919707
#> 70      7        0.7       V10 0.4199085
#> 71      8        0.8        V9 0.4143117
#> 72      8        0.8       V10 0.4192086
#> 73      8        0.8        V1 0.3427536
#> 74      8        0.8        V5 0.3919120
#> 75      8        0.8        V2 0.3914000
#> 76      8        0.8        V3 0.3931246
#> 77      8        0.8        V4 0.4101732
#> 78      8        0.8        V8 0.3160742
#> 79      8        0.8        V6 0.3788279
#> 80      8        0.8        V7 0.3714506
#> 81      9        0.9        V5 0.3913604
#> 82      9        0.9        V7 0.3707650
#> 83      9        0.9        V8 0.3155965
#> 84      9        0.9        V9 0.4137497
#> 85      9        0.9        V6 0.3782779
#> 86      9        0.9       V10 0.4185098
#> 87      9        0.9        V1 0.3421887
#> 88      9        0.9        V2 0.3908301
#> 89      9        0.9        V3 0.3924179
#> 90      9        0.9        V4 0.4095394
#> 91     10        1.0        V1 0.3416247
#> 92     10        1.0        V4 0.4089065
#> 93     10        1.0        V5 0.3908096
#> 94     10        1.0        V9 0.4131884
#> 95     10        1.0        V3 0.3917124
#> 96     10        1.0        V7 0.3700806
#> 97     10        1.0        V8 0.3151195
#> 98     10        1.0        V2 0.3902610
#> 99     10        1.0        V6 0.3777286
#> 100    10        1.0       V10 0.4178121


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
#> 1       0                0          0 0.8757906 0.04753258 0.7868981 0.9523806
#> 2       0               20         20 0.8757906 0.04753258 0.7868981 0.9523806
#> 3       0               40         40 0.8757906 0.04753258 0.7868981 0.9523806
#> 4       0               60         60 0.8757906 0.04753258 0.7868981 0.9523806
#> 5      20                0         20 0.8617131 0.05024800 0.7715274 0.9435521
#> 6      20               20         40 0.8617131 0.05024800 0.7715274 0.9435521
#> 7      20               40         60 0.8617131 0.05024800 0.7715274 0.9435521
#> 8      20               60         80 0.8617131 0.05024800 0.7715274 0.9435521
#> 9      40                0         40 0.8478591 0.05282993 0.7566063 0.9345520
#> 10     40               20         60 0.8478591 0.05282993 0.7566063 0.9345520
#> 11     40               40         80 0.8478591 0.05282993 0.7566063 0.9345520
#> 12     40               60        100 0.8478591 0.05282993 0.7566063 0.9345520
#> 13     60                0         60 0.8342249 0.05528402 0.7420977 0.9254251
#> 14     60               20         80 0.8342249 0.05528402 0.7420977 0.9254251
#> 15     60               40        100 0.8342249 0.05528402 0.7420977 0.9254251
#> 16     60               60        120 0.8342249 0.05528402 0.7420977 0.9254251
#> 17     80                0         80 0.8208071 0.05761644 0.7279715 0.9162064
#> 18     80               20        100 0.8208071 0.05761644 0.7279715 0.9162064
#> 19     80               40        120 0.8208071 0.05761644 0.7279715 0.9162064
#> 20     80               60        140 0.8208071 0.05761644 0.7279715 0.9162064
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11991352 0.169555172 0.5835928
#> 2  0.30574618 0.11583624 0.134002513 0.5172760
#> 3  0.26001915 0.11168663 0.105112068 0.4588157
#> 4  0.22113100 0.10744838 0.081711922 0.4074615
#> 5  0.25589195 0.10925252 0.097630828 0.4595658
#> 6  0.21762106 0.10434106 0.075670042 0.4081197
#> 7  0.18507391 0.09954295 0.057996302 0.3630117
#> 8  0.15739448 0.09482960 0.043868445 0.3234930
#> 9  0.18213629 0.09732478 0.053457071 0.3635895
#> 10 0.15489621 0.09219718 0.040259504 0.3239991
#> 11 0.13173012 0.08727120 0.029830575 0.2893149
#> 12 0.11202872 0.08251874 0.021683292 0.2589103
#> 13 0.12963921 0.08528892 0.027190714 0.2897593
#> 14 0.11025053 0.08029402 0.019639534 0.2593000
#> 15 0.09376159 0.07553799 0.013851861 0.2325687
#> 16 0.07973872 0.07099553 0.009498448 0.2090706
#> 17 0.09227334 0.07381123 0.012421725 0.2329116
#> 18 0.07847305 0.06913581 0.008438612 0.2093723
#> 19 0.06673672 0.06470794 0.005534695 0.1886372
#> 20 0.05675565 0.06050551 0.003481633 0.1703253
```
