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
#> 1         0.1 0.3838513 0.01996126 0.3609797 0.4129667
#> 2         0.2 0.3832760 0.01996530 0.3604268 0.4123933
#> 3         0.3 0.3827015 0.01996956 0.3598747 0.4118206
#> 4         0.4 0.3821279 0.01997405 0.3593235 0.4112488
#> 5         0.5 0.3815551 0.01997877 0.3587731 0.4107053
#> 6         0.6 0.3809832 0.01998371 0.3582235 0.4101645
#> 7         0.7 0.3804122 0.01998887 0.3576748 0.4096245
#> 8         0.8 0.3798420 0.01999424 0.3571270 0.4090851
#> 9         0.9 0.3792726 0.01999984 0.3565799 0.4085507
#> 10        1.0 0.3787041 0.02000564 0.3560338 0.4080215

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4084474
#> 2       1        0.1        V5 0.3779840
#> 3       1        0.1        V9 0.4132200
#> 4       1        0.1        V3 0.3817291
#> 5       1        0.1        V4 0.3560430
#> 6       1        0.1        V8 0.3824756
#> 7       1        0.1        V2 0.4052070
#> 8       1        0.1        V6 0.4120943
#> 9       1        0.1        V7 0.4116890
#> 10      1        0.1       V10 0.4113936
#> 11      2        0.2        V7 0.4110929
#> 12      2        0.2        V8 0.3818464
#> 13      2        0.2        V9 0.4126916
#> 14      2        0.2       V10 0.4108371
#> 15      2        0.2        V1 0.4079478
#> 16      2        0.2        V5 0.3773921
#> 17      2        0.2        V2 0.4047344
#> 18      2        0.2        V3 0.3810960
#> 19      2        0.2        V4 0.3555014
#> 20      2        0.2        V6 0.4113657
#> 21      3        0.3        V4 0.3549606
#> 22      3        0.3        V5 0.3768012
#> 23      3        0.3        V3 0.3804640
#> 24      3        0.3        V7 0.4104977
#> 25      3        0.3        V8 0.3812182
#> 26      3        0.3        V9 0.4121638
#> 27      3        0.3        V6 0.4106384
#> 28      3        0.3       V10 0.4102814
#> 29      3        0.3        V1 0.4074488
#> 30      3        0.3        V2 0.4042625
#> 31      4        0.4        V1 0.4069504
#> 32      4        0.4        V9 0.4116368
#> 33      4        0.4        V3 0.3798330
#> 34      4        0.4        V4 0.3544206
#> 35      4        0.4        V5 0.3762112
#> 36      4        0.4        V2 0.4037910
#> 37      4        0.4        V6 0.4099124
#> 38      4        0.4        V7 0.4099033
#> 39      4        0.4        V8 0.3805911
#> 40      4        0.4       V10 0.4097265
#> 41      5        0.5        V8 0.3799650
#> 42      5        0.5        V9 0.4111104
#> 43      5        0.5       V10 0.4091723
#> 44      5        0.5        V1 0.4064526
#> 45      5        0.5        V5 0.3756221
#> 46      5        0.5        V2 0.4033201
#> 47      5        0.5        V3 0.3792030
#> 48      5        0.5        V4 0.3538814
#> 49      5        0.5        V6 0.4091877
#> 50      5        0.5        V7 0.4093098
#> 51      6        0.6        V4 0.3533431
#> 52      6        0.6        V5 0.3750339
#> 53      6        0.6        V7 0.4087172
#> 54      6        0.6        V8 0.3793399
#> 55      6        0.6        V9 0.4105847
#> 56      6        0.6        V6 0.4084643
#> 57      6        0.6       V10 0.4086188
#> 58      6        0.6        V1 0.4059555
#> 59      6        0.6        V2 0.4028498
#> 60      6        0.6        V3 0.3785741
#> 61      7        0.7        V1 0.4054589
#> 62      7        0.7        V3 0.3779463
#> 63      7        0.7        V4 0.3528056
#> 64      7        0.7        V5 0.3744467
#> 65      7        0.7        V9 0.4100597
#> 66      7        0.7        V6 0.4077421
#> 67      7        0.7        V7 0.4081254
#> 68      7        0.7        V8 0.3787159
#> 69      7        0.7        V2 0.4023800
#> 70      7        0.7       V10 0.4080661
#> 71      8        0.8        V9 0.4095353
#> 72      8        0.8       V10 0.4075141
#> 73      8        0.8        V1 0.4049629
#> 74      8        0.8        V5 0.3738604
#> 75      8        0.8        V2 0.4019107
#> 76      8        0.8        V3 0.3773194
#> 77      8        0.8        V4 0.3522689
#> 78      8        0.8        V8 0.3780929
#> 79      8        0.8        V6 0.4070212
#> 80      8        0.8        V7 0.4075345
#> 81      9        0.9        V5 0.3732750
#> 82      9        0.9        V7 0.4069444
#> 83      9        0.9        V8 0.3774709
#> 84      9        0.9        V9 0.4090116
#> 85      9        0.9        V6 0.4063016
#> 86      9        0.9       V10 0.4069629
#> 87      9        0.9        V1 0.4044676
#> 88      9        0.9        V2 0.4014420
#> 89      9        0.9        V3 0.3766937
#> 90      9        0.9        V4 0.3517330
#> 91     10        1.0        V1 0.4039729
#> 92     10        1.0        V4 0.3511979
#> 93     10        1.0        V5 0.3726905
#> 94     10        1.0        V9 0.4084886
#> 95     10        1.0        V3 0.3760689
#> 96     10        1.0        V7 0.4063552
#> 97     10        1.0        V8 0.3768499
#> 98     10        1.0        V2 0.4009739
#> 99     10        1.0        V6 0.4055833
#> 100    10        1.0       V10 0.4064124

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4084474
#> 2       1        0.1        V5 0.3779840
#> 3       1        0.1        V9 0.4132200
#> 4       1        0.1        V3 0.3817291
#> 5       1        0.1        V4 0.3560430
#> 6       1        0.1        V8 0.3824756
#> 7       1        0.1        V2 0.4052070
#> 8       1        0.1        V6 0.4120943
#> 9       1        0.1        V7 0.4116890
#> 10      1        0.1       V10 0.4113936
#> 11      2        0.2        V7 0.4110929
#> 12      2        0.2        V8 0.3818464
#> 13      2        0.2        V9 0.4126916
#> 14      2        0.2       V10 0.4108371
#> 15      2        0.2        V1 0.4079478
#> 16      2        0.2        V5 0.3773921
#> 17      2        0.2        V2 0.4047344
#> 18      2        0.2        V3 0.3810960
#> 19      2        0.2        V4 0.3555014
#> 20      2        0.2        V6 0.4113657
#> 21      3        0.3        V4 0.3549606
#> 22      3        0.3        V5 0.3768012
#> 23      3        0.3        V3 0.3804640
#> 24      3        0.3        V7 0.4104977
#> 25      3        0.3        V8 0.3812182
#> 26      3        0.3        V9 0.4121638
#> 27      3        0.3        V6 0.4106384
#> 28      3        0.3       V10 0.4102814
#> 29      3        0.3        V1 0.4074488
#> 30      3        0.3        V2 0.4042625
#> 31      4        0.4        V1 0.4069504
#> 32      4        0.4        V9 0.4116368
#> 33      4        0.4        V3 0.3798330
#> 34      4        0.4        V4 0.3544206
#> 35      4        0.4        V5 0.3762112
#> 36      4        0.4        V2 0.4037910
#> 37      4        0.4        V6 0.4099124
#> 38      4        0.4        V7 0.4099033
#> 39      4        0.4        V8 0.3805911
#> 40      4        0.4       V10 0.4097265
#> 41      5        0.5        V8 0.3799650
#> 42      5        0.5        V9 0.4111104
#> 43      5        0.5       V10 0.4091723
#> 44      5        0.5        V1 0.4064526
#> 45      5        0.5        V5 0.3756221
#> 46      5        0.5        V2 0.4033201
#> 47      5        0.5        V3 0.3792030
#> 48      5        0.5        V4 0.3538814
#> 49      5        0.5        V6 0.4091877
#> 50      5        0.5        V7 0.4093098
#> 51      6        0.6        V4 0.3533431
#> 52      6        0.6        V5 0.3750339
#> 53      6        0.6        V7 0.4087172
#> 54      6        0.6        V8 0.3793399
#> 55      6        0.6        V9 0.4105847
#> 56      6        0.6        V6 0.4084643
#> 57      6        0.6       V10 0.4086188
#> 58      6        0.6        V1 0.4059555
#> 59      6        0.6        V2 0.4028498
#> 60      6        0.6        V3 0.3785741
#> 61      7        0.7        V1 0.4054589
#> 62      7        0.7        V3 0.3779463
#> 63      7        0.7        V4 0.3528056
#> 64      7        0.7        V5 0.3744467
#> 65      7        0.7        V9 0.4100597
#> 66      7        0.7        V6 0.4077421
#> 67      7        0.7        V7 0.4081254
#> 68      7        0.7        V8 0.3787159
#> 69      7        0.7        V2 0.4023800
#> 70      7        0.7       V10 0.4080661
#> 71      8        0.8        V9 0.4095353
#> 72      8        0.8       V10 0.4075141
#> 73      8        0.8        V1 0.4049629
#> 74      8        0.8        V5 0.3738604
#> 75      8        0.8        V2 0.4019107
#> 76      8        0.8        V3 0.3773194
#> 77      8        0.8        V4 0.3522689
#> 78      8        0.8        V8 0.3780929
#> 79      8        0.8        V6 0.4070212
#> 80      8        0.8        V7 0.4075345
#> 81      9        0.9        V5 0.3732750
#> 82      9        0.9        V7 0.4069444
#> 83      9        0.9        V8 0.3774709
#> 84      9        0.9        V9 0.4090116
#> 85      9        0.9        V6 0.4063016
#> 86      9        0.9       V10 0.4069629
#> 87      9        0.9        V1 0.4044676
#> 88      9        0.9        V2 0.4014420
#> 89      9        0.9        V3 0.3766937
#> 90      9        0.9        V4 0.3517330
#> 91     10        1.0        V1 0.4039729
#> 92     10        1.0        V4 0.3511979
#> 93     10        1.0        V5 0.3726905
#> 94     10        1.0        V9 0.4084886
#> 95     10        1.0        V3 0.3760689
#> 96     10        1.0        V7 0.4063552
#> 97     10        1.0        V8 0.3768499
#> 98     10        1.0        V2 0.4009739
#> 99     10        1.0        V6 0.4055833
#> 100    10        1.0       V10 0.4064124


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
#> 1       0                0          0 0.8757906 0.05048044 0.7676074 0.9389505
#> 2       0               20         20 0.8757906 0.05048044 0.7676074 0.9389505
#> 3       0               40         40 0.8757906 0.05048044 0.7676074 0.9389505
#> 4       0               60         60 0.8757906 0.05048044 0.7676074 0.9389505
#> 5      20                0         20 0.8617131 0.05348554 0.7441232 0.9286292
#> 6      20               20         40 0.8617131 0.05348554 0.7441232 0.9286292
#> 7      20               40         60 0.8617131 0.05348554 0.7441232 0.9286292
#> 8      20               60         80 0.8617131 0.05348554 0.7441232 0.9286292
#> 9      40                0         40 0.8478591 0.05630052 0.7216306 0.9182049
#> 10     40               20         60 0.8478591 0.05630052 0.7216306 0.9182049
#> 11     40               40         80 0.8478591 0.05630052 0.7216306 0.9182049
#> 12     40               60        100 0.8478591 0.05630052 0.7216306 0.9182049
#> 13     60                0         60 0.8342249 0.05894238 0.7000278 0.9077171
#> 14     60               20         80 0.8342249 0.05894238 0.7000278 0.9077171
#> 15     60               40        100 0.8342249 0.05894238 0.7000278 0.9077171
#> 16     60               60        120 0.8342249 0.05894238 0.7000278 0.9077171
#> 17     80                0         80 0.8208071 0.06142522 0.6792355 0.8971960
#> 18     80               20        100 0.8208071 0.06142522 0.6792355 0.8971960
#> 19     80               40        120 0.8208071 0.06142522 0.6792355 0.8971960
#> 20     80               60        140 0.8208071 0.06142522 0.6792355 0.8971960
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12767406 0.1509349423 0.5888872
#> 2  0.30574618 0.11870037 0.1145524708 0.5165363
#> 3  0.26001915 0.10974846 0.0859568038 0.4535310
#> 4  0.22113100 0.10103771 0.0636113470 0.3988931
#> 5  0.25589195 0.11312079 0.0778803726 0.4725003
#> 6  0.21762106 0.10391925 0.0573332431 0.4153275
#> 7  0.18507391 0.09521861 0.0414609413 0.3658219
#> 8  0.15739448 0.08705442 0.0293452785 0.3229914
#> 9  0.18213629 0.09855815 0.0370452873 0.3807082
#> 10 0.15489621 0.08993170 0.0260083431 0.3358702
#> 11 0.13173012 0.08196185 0.0177625352 0.2970718
#> 12 0.11202872 0.07461339 0.0117326949 0.2634632
#> 13 0.12963921 0.08494415 0.0155318521 0.3087405
#> 14 0.11025053 0.07721668 0.0101304907 0.2735768
#> 15 0.09376159 0.07015611 0.0063282847 0.2430785
#> 16 0.07973872 0.06370282 0.0037511247 0.2165595
#> 17 0.09227334 0.07272277 0.0053497840 0.2522612
#> 18 0.07847305 0.06597494 0.0031086919 0.2245525
#> 19 0.06673672 0.05983916 0.0016876965 0.2004064
#> 20 0.05675565 0.05425010 0.0008437663 0.1792838
```
