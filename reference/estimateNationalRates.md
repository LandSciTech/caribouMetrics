# Sample demographic rates

Apply the sampled coefficients to the disturbance covariates to
calculate expected recruitment and survival according to the beta
regression models estimated by Johnson et al.
(2020).`estimateNationalRates` is a wrapper around
`estimateNationalRate` to sample both survival and recruitment rates
based on the result of
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md)
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
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md).
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
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/simulateTrajectoriesFromPosterior.md),
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/reference/dataFromSheets.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/reference/estimateBayesianRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/reference/plotTrajectories.md),
[`popGrowthTableJohnsonECCC`](https://landscitech.github.io/caribouMetrics/reference/popGrowthTableJohnsonECCC.md),
[`simulateObservations()`](https://landscitech.github.io/caribouMetrics/reference/simulateObservations.md),
[`trajectoriesFromBayesian()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromBayesian.md),
[`trajectoriesFromNational()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromNational.md),
[`trajectoriesFromSummary()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummary.md),
[`trajectoriesFromSummaryForApp()`](https://landscitech.github.io/caribouMetrics/reference/trajectoriesFromSummaryForApp.md)

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
#> 1         0.1 0.3838513 0.01733836 0.3694991 0.4191844
#> 2         0.2 0.3832760 0.01731794 0.3689420 0.4185842
#> 3         0.3 0.3827015 0.01729773 0.3683857 0.4179848
#> 4         0.4 0.3821279 0.01727775 0.3678302 0.4173863
#> 5         0.5 0.3815551 0.01725798 0.3672756 0.4167886
#> 6         0.6 0.3809832 0.01723843 0.3667218 0.4161918
#> 7         0.7 0.3804122 0.01721909 0.3661689 0.4155959
#> 8         0.8 0.3798420 0.01719996 0.3656168 0.4150008
#> 9         0.9 0.3792726 0.01718105 0.3650655 0.4144066
#> 10        1.0 0.3787041 0.01716236 0.3645151 0.4138133

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4123812
#> 2       1        0.1        V5 0.4195674
#> 3       1        0.1        V9 0.4178653
#> 4       1        0.1        V3 0.3946122
#> 5       1        0.1        V4 0.3958369
#> 6       1        0.1        V8 0.3953822
#> 7       1        0.1        V2 0.3744726
#> 8       1        0.1        V6 0.3680552
#> 9       1        0.1        V7 0.4108720
#> 10      1        0.1       V10 0.4012677
#> 11      2        0.2        V7 0.4103809
#> 12      2        0.2        V8 0.3947965
#> 13      2        0.2        V9 0.4171442
#> 14      2        0.2       V10 0.4007205
#> 15      2        0.2        V1 0.4117892
#> 16      2        0.2        V5 0.4190022
#> 17      2        0.2        V2 0.3739523
#> 18      2        0.2        V3 0.3940482
#> 19      2        0.2        V4 0.3952325
#> 20      2        0.2        V6 0.3674874
#> 21      3        0.3        V4 0.3946290
#> 22      3        0.3        V5 0.4184378
#> 23      3        0.3        V3 0.3934849
#> 24      3        0.3        V7 0.4098903
#> 25      3        0.3        V8 0.3942118
#> 26      3        0.3        V9 0.4164244
#> 27      3        0.3        V6 0.3669204
#> 28      3        0.3       V10 0.4001740
#> 29      3        0.3        V1 0.4111980
#> 30      3        0.3        V2 0.3734327
#> 31      4        0.4        V1 0.4106077
#> 32      4        0.4        V9 0.4157058
#> 33      4        0.4        V3 0.3929225
#> 34      4        0.4        V4 0.3940265
#> 35      4        0.4        V5 0.4178741
#> 36      4        0.4        V2 0.3729139
#> 37      4        0.4        V6 0.3663543
#> 38      4        0.4        V7 0.4094004
#> 39      4        0.4        V8 0.3936279
#> 40      4        0.4       V10 0.3996283
#> 41      5        0.5        V8 0.3930448
#> 42      5        0.5        V9 0.4149885
#> 43      5        0.5       V10 0.3990833
#> 44      5        0.5        V1 0.4100183
#> 45      5        0.5        V5 0.4173112
#> 46      5        0.5        V2 0.3723957
#> 47      5        0.5        V3 0.3923609
#> 48      5        0.5        V4 0.3934248
#> 49      5        0.5        V6 0.3657891
#> 50      5        0.5        V7 0.4089110
#> 51      6        0.6        V4 0.3928241
#> 52      6        0.6        V5 0.4167491
#> 53      6        0.6        V7 0.4084222
#> 54      6        0.6        V8 0.3924627
#> 55      6        0.6        V9 0.4142724
#> 56      6        0.6        V6 0.3652248
#> 57      6        0.6       V10 0.3985391
#> 58      6        0.6        V1 0.4094297
#> 59      6        0.6        V2 0.3718783
#> 60      6        0.6        V3 0.3918001
#> 61      7        0.7        V1 0.4088419
#> 62      7        0.7        V3 0.3912401
#> 63      7        0.7        V4 0.3922243
#> 64      7        0.7        V5 0.4161877
#> 65      7        0.7        V9 0.4135575
#> 66      7        0.7        V6 0.3646613
#> 67      7        0.7        V7 0.4079340
#> 68      7        0.7        V8 0.3918813
#> 69      7        0.7        V2 0.3713616
#> 70      7        0.7       V10 0.3979956
#> 71      8        0.8        V9 0.4128439
#> 72      8        0.8       V10 0.3974528
#> 73      8        0.8        V1 0.4082550
#> 74      8        0.8        V5 0.4156270
#> 75      8        0.8        V2 0.3708456
#> 76      8        0.8        V3 0.3906809
#> 77      8        0.8        V4 0.3916254
#> 78      8        0.8        V8 0.3913009
#> 79      8        0.8        V6 0.3640987
#> 80      8        0.8        V7 0.4074464
#> 81      9        0.9        V5 0.4150671
#> 82      9        0.9        V7 0.4069593
#> 83      9        0.9        V8 0.3907213
#> 84      9        0.9        V9 0.4121315
#> 85      9        0.9        V6 0.3635370
#> 86      9        0.9       V10 0.3969108
#> 87      9        0.9        V1 0.4076689
#> 88      9        0.9        V2 0.3703303
#> 89      9        0.9        V3 0.3901225
#> 90      9        0.9        V4 0.3910274
#> 91     10        1.0        V1 0.4070837
#> 92     10        1.0        V4 0.3904303
#> 93     10        1.0        V5 0.4145080
#> 94     10        1.0        V9 0.4114203
#> 95     10        1.0        V3 0.3895648
#> 96     10        1.0        V7 0.4064729
#> 97     10        1.0        V8 0.3901426
#> 98     10        1.0        V2 0.3698157
#> 99     10        1.0        V6 0.3629761
#> 100    10        1.0       V10 0.3963695

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4123812
#> 2       1        0.1        V5 0.4195674
#> 3       1        0.1        V9 0.4178653
#> 4       1        0.1        V3 0.3946122
#> 5       1        0.1        V4 0.3958369
#> 6       1        0.1        V8 0.3953822
#> 7       1        0.1        V2 0.3744726
#> 8       1        0.1        V6 0.3680552
#> 9       1        0.1        V7 0.4108720
#> 10      1        0.1       V10 0.4012677
#> 11      2        0.2        V7 0.4103809
#> 12      2        0.2        V8 0.3947965
#> 13      2        0.2        V9 0.4171442
#> 14      2        0.2       V10 0.4007205
#> 15      2        0.2        V1 0.4117892
#> 16      2        0.2        V5 0.4190022
#> 17      2        0.2        V2 0.3739523
#> 18      2        0.2        V3 0.3940482
#> 19      2        0.2        V4 0.3952325
#> 20      2        0.2        V6 0.3674874
#> 21      3        0.3        V4 0.3946290
#> 22      3        0.3        V5 0.4184378
#> 23      3        0.3        V3 0.3934849
#> 24      3        0.3        V7 0.4098903
#> 25      3        0.3        V8 0.3942118
#> 26      3        0.3        V9 0.4164244
#> 27      3        0.3        V6 0.3669204
#> 28      3        0.3       V10 0.4001740
#> 29      3        0.3        V1 0.4111980
#> 30      3        0.3        V2 0.3734327
#> 31      4        0.4        V1 0.4106077
#> 32      4        0.4        V9 0.4157058
#> 33      4        0.4        V3 0.3929225
#> 34      4        0.4        V4 0.3940265
#> 35      4        0.4        V5 0.4178741
#> 36      4        0.4        V2 0.3729139
#> 37      4        0.4        V6 0.3663543
#> 38      4        0.4        V7 0.4094004
#> 39      4        0.4        V8 0.3936279
#> 40      4        0.4       V10 0.3996283
#> 41      5        0.5        V8 0.3930448
#> 42      5        0.5        V9 0.4149885
#> 43      5        0.5       V10 0.3990833
#> 44      5        0.5        V1 0.4100183
#> 45      5        0.5        V5 0.4173112
#> 46      5        0.5        V2 0.3723957
#> 47      5        0.5        V3 0.3923609
#> 48      5        0.5        V4 0.3934248
#> 49      5        0.5        V6 0.3657891
#> 50      5        0.5        V7 0.4089110
#> 51      6        0.6        V4 0.3928241
#> 52      6        0.6        V5 0.4167491
#> 53      6        0.6        V7 0.4084222
#> 54      6        0.6        V8 0.3924627
#> 55      6        0.6        V9 0.4142724
#> 56      6        0.6        V6 0.3652248
#> 57      6        0.6       V10 0.3985391
#> 58      6        0.6        V1 0.4094297
#> 59      6        0.6        V2 0.3718783
#> 60      6        0.6        V3 0.3918001
#> 61      7        0.7        V1 0.4088419
#> 62      7        0.7        V3 0.3912401
#> 63      7        0.7        V4 0.3922243
#> 64      7        0.7        V5 0.4161877
#> 65      7        0.7        V9 0.4135575
#> 66      7        0.7        V6 0.3646613
#> 67      7        0.7        V7 0.4079340
#> 68      7        0.7        V8 0.3918813
#> 69      7        0.7        V2 0.3713616
#> 70      7        0.7       V10 0.3979956
#> 71      8        0.8        V9 0.4128439
#> 72      8        0.8       V10 0.3974528
#> 73      8        0.8        V1 0.4082550
#> 74      8        0.8        V5 0.4156270
#> 75      8        0.8        V2 0.3708456
#> 76      8        0.8        V3 0.3906809
#> 77      8        0.8        V4 0.3916254
#> 78      8        0.8        V8 0.3913009
#> 79      8        0.8        V6 0.3640987
#> 80      8        0.8        V7 0.4074464
#> 81      9        0.9        V5 0.4150671
#> 82      9        0.9        V7 0.4069593
#> 83      9        0.9        V8 0.3907213
#> 84      9        0.9        V9 0.4121315
#> 85      9        0.9        V6 0.3635370
#> 86      9        0.9       V10 0.3969108
#> 87      9        0.9        V1 0.4076689
#> 88      9        0.9        V2 0.3703303
#> 89      9        0.9        V3 0.3901225
#> 90      9        0.9        V4 0.3910274
#> 91     10        1.0        V1 0.4070837
#> 92     10        1.0        V4 0.3904303
#> 93     10        1.0        V5 0.4145080
#> 94     10        1.0        V9 0.4114203
#> 95     10        1.0        V3 0.3895648
#> 96     10        1.0        V7 0.4064729
#> 97     10        1.0        V8 0.3901426
#> 98     10        1.0        V2 0.3698157
#> 99     10        1.0        V6 0.3629761
#> 100    10        1.0       V10 0.3963695


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
#> 1       0                0          0 0.8757906 0.05291377 0.7812288 0.9583735
#> 2       0               20         20 0.8757906 0.05291377 0.7812288 0.9583735
#> 3       0               40         40 0.8757906 0.05291377 0.7812288 0.9583735
#> 4       0               60         60 0.8757906 0.05291377 0.7812288 0.9583735
#> 5      20                0         20 0.8617131 0.05389750 0.7667851 0.9482004
#> 6      20               20         40 0.8617131 0.05389750 0.7667851 0.9482004
#> 7      20               40         60 0.8617131 0.05389750 0.7667851 0.9482004
#> 8      20               60         80 0.8617131 0.05389750 0.7667851 0.9482004
#> 9      40                0         40 0.8478591 0.05477264 0.7527429 0.9377487
#> 10     40               20         60 0.8478591 0.05477264 0.7527429 0.9377487
#> 11     40               40         80 0.8478591 0.05477264 0.7527429 0.9377487
#> 12     40               60        100 0.8478591 0.05477264 0.7527429 0.9377487
#> 13     60                0         60 0.8342249 0.05555752 0.7390709 0.9270996
#> 14     60               20         80 0.8342249 0.05555752 0.7390709 0.9270996
#> 15     60               40        100 0.8342249 0.05555752 0.7390709 0.9270996
#> 16     60               60        120 0.8342249 0.05555752 0.7390709 0.9270996
#> 17     80                0         80 0.8208071 0.05626619 0.7257429 0.9163136
#> 18     80               20        100 0.8208071 0.05626619 0.7257429 0.9163136
#> 19     80               40        120 0.8208071 0.05626619 0.7257429 0.9163136
#> 20     80               60        140 0.8208071 0.05626619 0.7257429 0.9163136
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12117747 0.158247061 0.5549855
#> 2  0.30574618 0.11569884 0.119608645 0.5119282
#> 3  0.26001915 0.11043445 0.089376631 0.4723875
#> 4  0.22113100 0.10527875 0.065855892 0.4361357
#> 5  0.25589195 0.10763944 0.094212634 0.4379584
#> 6  0.21762106 0.10169758 0.069606679 0.4046040
#> 7  0.18507391 0.09615665 0.050587590 0.3740770
#> 8  0.15739448 0.09091560 0.036039738 0.3461476
#> 9  0.18213629 0.09395964 0.053610188 0.3475509
#> 10 0.15489621 0.08816538 0.038339388 0.3218808
#> 11 0.13173012 0.08283447 0.026786525 0.2983940
#> 12 0.11202872 0.07787780 0.018193373 0.2768974
#> 13 0.12963921 0.08085081 0.028602534 0.2779779
#> 14 0.11025053 0.07553003 0.019532536 0.2582021
#> 15 0.09376159 0.07066697 0.012903516 0.2400809
#> 16 0.07973872 0.06618732 0.008186833 0.2234620
#> 17 0.09227334 0.06877563 0.013927379 0.2242983
#> 18 0.07847305 0.06406549 0.008905378 0.2089736
#> 19 0.06673672 0.05977706 0.005429683 0.1948897
#> 20 0.05675565 0.05585053 0.003123624 0.1819293
```
