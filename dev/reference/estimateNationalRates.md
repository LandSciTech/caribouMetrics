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
#> 1         0.1 0.3838513 0.02941329 0.3400259 0.4171685
#> 2         0.2 0.3832760 0.02936129 0.3395414 0.4165549
#> 3         0.3 0.3827015 0.02930941 0.3390576 0.4159422
#> 4         0.4 0.3821279 0.02925766 0.3385745 0.4153304
#> 5         0.5 0.3815551 0.02920603 0.3380920 0.4147194
#> 6         0.6 0.3809832 0.02915453 0.3376103 0.4141094
#> 7         0.7 0.3804122 0.02910315 0.3371292 0.4135003
#> 8         0.8 0.3798420 0.02905189 0.3366489 0.4128921
#> 9         0.9 0.3792726 0.02900075 0.3361692 0.4122848
#> 10        1.0 0.3787041 0.02894974 0.3356902 0.4116784

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4102632
#> 2       1        0.1        V5 0.3748773
#> 3       1        0.1        V9 0.4187528
#> 4       1        0.1        V3 0.3514791
#> 5       1        0.1        V4 0.4021971
#> 6       1        0.1        V8 0.4117115
#> 7       1        0.1        V2 0.3379513
#> 8       1        0.1        V6 0.3798400
#> 9       1        0.1        V7 0.3664042
#> 10      1        0.1       V10 0.3471716
#> 11      2        0.2        V7 0.3658279
#> 12      2        0.2        V8 0.4110430
#> 13      2        0.2        V9 0.4181551
#> 14      2        0.2       V10 0.3466291
#> 15      2        0.2        V1 0.4096009
#> 16      2        0.2        V5 0.3743088
#> 17      2        0.2        V2 0.3374837
#> 18      2        0.2        V3 0.3509628
#> 19      2        0.2        V4 0.4016218
#> 20      2        0.2        V6 0.3792865
#> 21      3        0.3        V4 0.4010472
#> 22      3        0.3        V5 0.3737412
#> 23      3        0.3        V3 0.3504473
#> 24      3        0.3        V7 0.3652525
#> 25      3        0.3        V8 0.4103757
#> 26      3        0.3        V9 0.4175582
#> 27      3        0.3        V6 0.3787338
#> 28      3        0.3       V10 0.3460875
#> 29      3        0.3        V1 0.4089395
#> 30      3        0.3        V2 0.3370166
#> 31      4        0.4        V1 0.4082793
#> 32      4        0.4        V9 0.4169622
#> 33      4        0.4        V3 0.3499325
#> 34      4        0.4        V4 0.4004735
#> 35      4        0.4        V5 0.3731745
#> 36      4        0.4        V2 0.3365503
#> 37      4        0.4        V6 0.3781819
#> 38      4        0.4        V7 0.3646781
#> 39      4        0.4        V8 0.4097094
#> 40      4        0.4       V10 0.3455467
#> 41      5        0.5        V8 0.4090442
#> 42      5        0.5        V9 0.4163671
#> 43      5        0.5       V10 0.3450067
#> 44      5        0.5        V1 0.4076201
#> 45      5        0.5        V5 0.3726086
#> 46      5        0.5        V2 0.3360845
#> 47      5        0.5        V3 0.3494185
#> 48      5        0.5        V4 0.3999006
#> 49      5        0.5        V6 0.3776308
#> 50      5        0.5        V7 0.3641045
#> 51      6        0.6        V4 0.3993285
#> 52      6        0.6        V5 0.3720435
#> 53      6        0.6        V7 0.3635318
#> 54      6        0.6        V8 0.4083801
#> 55      6        0.6        V9 0.4157728
#> 56      6        0.6        V6 0.3770805
#> 57      6        0.6       V10 0.3444676
#> 58      6        0.6        V1 0.4069620
#> 59      6        0.6        V2 0.3356195
#> 60      6        0.6        V3 0.3489052
#> 61      7        0.7        V1 0.4063050
#> 62      7        0.7        V3 0.3483927
#> 63      7        0.7        V4 0.3987572
#> 64      7        0.7        V5 0.3714793
#> 65      7        0.7        V9 0.4151793
#> 66      7        0.7        V6 0.3765310
#> 67      7        0.7        V7 0.3629600
#> 68      7        0.7        V8 0.4077171
#> 69      7        0.7        V2 0.3351550
#> 70      7        0.7       V10 0.3439293
#> 71      8        0.8        V9 0.4145867
#> 72      8        0.8       V10 0.3433919
#> 73      8        0.8        V1 0.4056490
#> 74      8        0.8        V5 0.3709160
#> 75      8        0.8        V2 0.3346912
#> 76      8        0.8        V3 0.3478809
#> 77      8        0.8        V4 0.3981868
#> 78      8        0.8        V8 0.4070552
#> 79      8        0.8        V6 0.3759823
#> 80      8        0.8        V7 0.3623891
#> 81      9        0.9        V5 0.3703535
#> 82      9        0.9        V7 0.3618191
#> 83      9        0.9        V8 0.4063943
#> 84      9        0.9        V9 0.4139949
#> 85      9        0.9        V6 0.3754344
#> 86      9        0.9       V10 0.3428553
#> 87      9        0.9        V1 0.4049941
#> 88      9        0.9        V2 0.3342280
#> 89      9        0.9        V3 0.3473699
#> 90      9        0.9        V4 0.3976171
#> 91     10        1.0        V1 0.4043402
#> 92     10        1.0        V4 0.3970483
#> 93     10        1.0        V5 0.3697919
#> 94     10        1.0        V9 0.4134040
#> 95     10        1.0        V3 0.3468597
#> 96     10        1.0        V7 0.3612500
#> 97     10        1.0        V8 0.4057345
#> 98     10        1.0        V2 0.3337655
#> 99     10        1.0        V6 0.3748873
#> 100    10        1.0       V10 0.3423195

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4102632
#> 2       1        0.1        V5 0.3748773
#> 3       1        0.1        V9 0.4187528
#> 4       1        0.1        V3 0.3514791
#> 5       1        0.1        V4 0.4021971
#> 6       1        0.1        V8 0.4117115
#> 7       1        0.1        V2 0.3379513
#> 8       1        0.1        V6 0.3798400
#> 9       1        0.1        V7 0.3664042
#> 10      1        0.1       V10 0.3471716
#> 11      2        0.2        V7 0.3658279
#> 12      2        0.2        V8 0.4110430
#> 13      2        0.2        V9 0.4181551
#> 14      2        0.2       V10 0.3466291
#> 15      2        0.2        V1 0.4096009
#> 16      2        0.2        V5 0.3743088
#> 17      2        0.2        V2 0.3374837
#> 18      2        0.2        V3 0.3509628
#> 19      2        0.2        V4 0.4016218
#> 20      2        0.2        V6 0.3792865
#> 21      3        0.3        V4 0.4010472
#> 22      3        0.3        V5 0.3737412
#> 23      3        0.3        V3 0.3504473
#> 24      3        0.3        V7 0.3652525
#> 25      3        0.3        V8 0.4103757
#> 26      3        0.3        V9 0.4175582
#> 27      3        0.3        V6 0.3787338
#> 28      3        0.3       V10 0.3460875
#> 29      3        0.3        V1 0.4089395
#> 30      3        0.3        V2 0.3370166
#> 31      4        0.4        V1 0.4082793
#> 32      4        0.4        V9 0.4169622
#> 33      4        0.4        V3 0.3499325
#> 34      4        0.4        V4 0.4004735
#> 35      4        0.4        V5 0.3731745
#> 36      4        0.4        V2 0.3365503
#> 37      4        0.4        V6 0.3781819
#> 38      4        0.4        V7 0.3646781
#> 39      4        0.4        V8 0.4097094
#> 40      4        0.4       V10 0.3455467
#> 41      5        0.5        V8 0.4090442
#> 42      5        0.5        V9 0.4163671
#> 43      5        0.5       V10 0.3450067
#> 44      5        0.5        V1 0.4076201
#> 45      5        0.5        V5 0.3726086
#> 46      5        0.5        V2 0.3360845
#> 47      5        0.5        V3 0.3494185
#> 48      5        0.5        V4 0.3999006
#> 49      5        0.5        V6 0.3776308
#> 50      5        0.5        V7 0.3641045
#> 51      6        0.6        V4 0.3993285
#> 52      6        0.6        V5 0.3720435
#> 53      6        0.6        V7 0.3635318
#> 54      6        0.6        V8 0.4083801
#> 55      6        0.6        V9 0.4157728
#> 56      6        0.6        V6 0.3770805
#> 57      6        0.6       V10 0.3444676
#> 58      6        0.6        V1 0.4069620
#> 59      6        0.6        V2 0.3356195
#> 60      6        0.6        V3 0.3489052
#> 61      7        0.7        V1 0.4063050
#> 62      7        0.7        V3 0.3483927
#> 63      7        0.7        V4 0.3987572
#> 64      7        0.7        V5 0.3714793
#> 65      7        0.7        V9 0.4151793
#> 66      7        0.7        V6 0.3765310
#> 67      7        0.7        V7 0.3629600
#> 68      7        0.7        V8 0.4077171
#> 69      7        0.7        V2 0.3351550
#> 70      7        0.7       V10 0.3439293
#> 71      8        0.8        V9 0.4145867
#> 72      8        0.8       V10 0.3433919
#> 73      8        0.8        V1 0.4056490
#> 74      8        0.8        V5 0.3709160
#> 75      8        0.8        V2 0.3346912
#> 76      8        0.8        V3 0.3478809
#> 77      8        0.8        V4 0.3981868
#> 78      8        0.8        V8 0.4070552
#> 79      8        0.8        V6 0.3759823
#> 80      8        0.8        V7 0.3623891
#> 81      9        0.9        V5 0.3703535
#> 82      9        0.9        V7 0.3618191
#> 83      9        0.9        V8 0.4063943
#> 84      9        0.9        V9 0.4139949
#> 85      9        0.9        V6 0.3754344
#> 86      9        0.9       V10 0.3428553
#> 87      9        0.9        V1 0.4049941
#> 88      9        0.9        V2 0.3342280
#> 89      9        0.9        V3 0.3473699
#> 90      9        0.9        V4 0.3976171
#> 91     10        1.0        V1 0.4043402
#> 92     10        1.0        V4 0.3970483
#> 93     10        1.0        V5 0.3697919
#> 94     10        1.0        V9 0.4134040
#> 95     10        1.0        V3 0.3468597
#> 96     10        1.0        V7 0.3612500
#> 97     10        1.0        V8 0.4057345
#> 98     10        1.0        V2 0.3337655
#> 99     10        1.0        V6 0.3748873
#> 100    10        1.0       V10 0.3423195


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
#> 1       0                0          0 0.8757906 0.05169076 0.7687993 0.9470975
#> 2       0               20         20 0.8757906 0.05169076 0.7687993 0.9470975
#> 3       0               40         40 0.8757906 0.05169076 0.7687993 0.9470975
#> 4       0               60         60 0.8757906 0.05169076 0.7687993 0.9470975
#> 5      20                0         20 0.8617131 0.05379252 0.7532037 0.9370744
#> 6      20               20         40 0.8617131 0.05379252 0.7532037 0.9370744
#> 7      20               40         60 0.8617131 0.05379252 0.7532037 0.9370744
#> 8      20               60         80 0.8617131 0.05379252 0.7532037 0.9370744
#> 9      40                0         40 0.8478591 0.05578970 0.7380611 0.9268883
#> 10     40               20         60 0.8478591 0.05578970 0.7380611 0.9268883
#> 11     40               40         80 0.8478591 0.05578970 0.7380611 0.9268883
#> 12     40               60        100 0.8478591 0.05578970 0.7380611 0.9268883
#> 13     60                0         60 0.8342249 0.05768784 0.7233368 0.9165897
#> 14     60               20         80 0.8342249 0.05768784 0.7233368 0.9165897
#> 15     60               40        100 0.8342249 0.05768784 0.7233368 0.9165897
#> 16     60               60        120 0.8342249 0.05768784 0.7233368 0.9165897
#> 17     80                0         80 0.8208071 0.05949228 0.7090020 0.9062175
#> 18     80               20        100 0.8208071 0.05949228 0.7090020 0.9062175
#> 19     80               40        120 0.8208071 0.05949228 0.7090020 0.9062175
#> 20     80               60        140 0.8208071 0.05949228 0.7090020 0.9062175
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12966341 0.157153877 0.6123377
#> 2  0.30574618 0.12578774 0.127088528 0.5568371
#> 3  0.26001915 0.12193162 0.102117711 0.5064828
#> 4  0.22113100 0.11781986 0.081437494 0.4609511
#> 5  0.25589195 0.11335565 0.093074735 0.4810473
#> 6  0.21762106 0.10882082 0.073968614 0.4379920
#> 7  0.18507391 0.10440869 0.058239525 0.3991725
#> 8  0.15739448 0.09995915 0.045363807 0.3642029
#> 9  0.18213629 0.09802151 0.052588694 0.3796255
#> 10 0.15489621 0.09342141 0.040760568 0.3465999
#> 11 0.13173012 0.08898080 0.031178268 0.3168591
#> 12 0.11202872 0.08460791 0.023487501 0.2900678
#> 13 0.12963921 0.08410137 0.027782560 0.3018853
#> 14 0.11025053 0.07971965 0.020783502 0.2765725
#> 15 0.09376159 0.07550875 0.015258919 0.2537433
#> 16 0.07973872 0.07141573 0.010962215 0.2331304
#> 17 0.09227334 0.07170547 0.013343653 0.2422296
#> 18 0.07847305 0.06767199 0.009491112 0.2227232
#> 19 0.06673672 0.06380896 0.006572494 0.2050694
#> 20 0.05675565 0.06008532 0.004412116 0.1890621
```
