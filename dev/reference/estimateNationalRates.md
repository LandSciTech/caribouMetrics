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
#> 1         0.1 0.3838513 0.02219911 0.3390991 0.4058215
#> 2         0.2 0.3832760 0.02220005 0.3386191 0.4053074
#> 3         0.3 0.3827015 0.02220112 0.3381398 0.4047940
#> 4         0.4 0.3821279 0.02220231 0.3376612 0.4042812
#> 5         0.5 0.3815551 0.02220363 0.3371833 0.4037691
#> 6         0.6 0.3809832 0.02220507 0.3367061 0.4032576
#> 7         0.7 0.3804122 0.02220664 0.3362296 0.4027468
#> 8         0.8 0.3798420 0.02220832 0.3357538 0.4022366
#> 9         0.9 0.3792726 0.02221013 0.3352786 0.4017271
#> 10        1.0 0.3787041 0.02221205 0.3348041 0.4012182

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3585267
#> 2       1        0.1        V5 0.4094526
#> 3       1        0.1        V9 0.3369756
#> 4       1        0.1        V3 0.3733557
#> 5       1        0.1        V4 0.3699310
#> 6       1        0.1        V8 0.3464132
#> 7       1        0.1        V2 0.3882812
#> 8       1        0.1        V6 0.3933144
#> 9       1        0.1        V7 0.3555791
#> 10      1        0.1       V10 0.3721737
#> 11      2        0.2        V7 0.3550518
#> 12      2        0.2        V8 0.3458169
#> 13      2        0.2        V9 0.3365294
#> 14      2        0.2       V10 0.3717430
#> 15      2        0.2        V1 0.3579471
#> 16      2        0.2        V5 0.4089375
#> 17      2        0.2        V2 0.3877439
#> 18      2        0.2        V3 0.3727586
#> 19      2        0.2        V4 0.3694909
#> 20      2        0.2        V6 0.3928037
#> 21      3        0.3        V4 0.3690513
#> 22      3        0.3        V5 0.4084231
#> 23      3        0.3        V3 0.3721624
#> 24      3        0.3        V7 0.3545252
#> 25      3        0.3        V8 0.3452216
#> 26      3        0.3        V9 0.3360838
#> 27      3        0.3        V6 0.3922937
#> 28      3        0.3       V10 0.3713128
#> 29      3        0.3        V1 0.3573684
#> 30      3        0.3        V2 0.3872073
#> 31      4        0.4        V1 0.3567907
#> 32      4        0.4        V9 0.3356388
#> 33      4        0.4        V3 0.3715671
#> 34      4        0.4        V4 0.3686123
#> 35      4        0.4        V5 0.4079093
#> 36      4        0.4        V2 0.3866715
#> 37      4        0.4        V6 0.3917844
#> 38      4        0.4        V7 0.3539995
#> 39      4        0.4        V8 0.3446273
#> 40      4        0.4       V10 0.3708832
#> 41      5        0.5        V8 0.3440340
#> 42      5        0.5        V9 0.3351944
#> 43      5        0.5       V10 0.3704540
#> 44      5        0.5        V1 0.3562139
#> 45      5        0.5        V5 0.4073962
#> 46      5        0.5        V2 0.3861364
#> 47      5        0.5        V3 0.3709728
#> 48      5        0.5        V4 0.3681738
#> 49      5        0.5        V6 0.3912757
#> 50      5        0.5        V7 0.3534745
#> 51      6        0.6        V4 0.3677358
#> 52      6        0.6        V5 0.4068837
#> 53      6        0.6        V7 0.3529503
#> 54      6        0.6        V8 0.3434418
#> 55      6        0.6        V9 0.3347506
#> 56      6        0.6        V6 0.3907677
#> 57      6        0.6       V10 0.3700253
#> 58      6        0.6        V1 0.3556380
#> 59      6        0.6        V2 0.3856021
#> 60      6        0.6        V3 0.3703795
#> 61      7        0.7        V1 0.3550631
#> 62      7        0.7        V3 0.3697871
#> 63      7        0.7        V4 0.3672983
#> 64      7        0.7        V5 0.4063719
#> 65      7        0.7        V9 0.3343074
#> 66      7        0.7        V6 0.3902604
#> 67      7        0.7        V7 0.3524268
#> 68      7        0.7        V8 0.3428506
#> 69      7        0.7        V2 0.3850685
#> 70      7        0.7       V10 0.3695971
#> 71      8        0.8        V9 0.3338647
#> 72      8        0.8       V10 0.3691694
#> 73      8        0.8        V1 0.3544891
#> 74      8        0.8        V5 0.4058607
#> 75      8        0.8        V2 0.3845356
#> 76      8        0.8        V3 0.3691956
#> 77      8        0.8        V4 0.3668613
#> 78      8        0.8        V8 0.3422604
#> 79      8        0.8        V6 0.3897537
#> 80      8        0.8        V7 0.3519042
#> 81      9        0.9        V5 0.4053501
#> 82      9        0.9        V7 0.3513823
#> 83      9        0.9        V8 0.3416712
#> 84      9        0.9        V9 0.3334227
#> 85      9        0.9        V6 0.3892477
#> 86      9        0.9       V10 0.3687422
#> 87      9        0.9        V1 0.3539160
#> 88      9        0.9        V2 0.3840035
#> 89      9        0.9        V3 0.3686051
#> 90      9        0.9        V4 0.3664249
#> 91     10        1.0        V1 0.3533438
#> 92     10        1.0        V4 0.3659890
#> 93     10        1.0        V5 0.4048402
#> 94     10        1.0        V9 0.3329812
#> 95     10        1.0        V3 0.3680156
#> 96     10        1.0        V7 0.3508612
#> 97     10        1.0        V8 0.3410830
#> 98     10        1.0        V2 0.3834721
#> 99     10        1.0        V6 0.3887423
#> 100    10        1.0       V10 0.3683155

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3585267
#> 2       1        0.1        V5 0.4094526
#> 3       1        0.1        V9 0.3369756
#> 4       1        0.1        V3 0.3733557
#> 5       1        0.1        V4 0.3699310
#> 6       1        0.1        V8 0.3464132
#> 7       1        0.1        V2 0.3882812
#> 8       1        0.1        V6 0.3933144
#> 9       1        0.1        V7 0.3555791
#> 10      1        0.1       V10 0.3721737
#> 11      2        0.2        V7 0.3550518
#> 12      2        0.2        V8 0.3458169
#> 13      2        0.2        V9 0.3365294
#> 14      2        0.2       V10 0.3717430
#> 15      2        0.2        V1 0.3579471
#> 16      2        0.2        V5 0.4089375
#> 17      2        0.2        V2 0.3877439
#> 18      2        0.2        V3 0.3727586
#> 19      2        0.2        V4 0.3694909
#> 20      2        0.2        V6 0.3928037
#> 21      3        0.3        V4 0.3690513
#> 22      3        0.3        V5 0.4084231
#> 23      3        0.3        V3 0.3721624
#> 24      3        0.3        V7 0.3545252
#> 25      3        0.3        V8 0.3452216
#> 26      3        0.3        V9 0.3360838
#> 27      3        0.3        V6 0.3922937
#> 28      3        0.3       V10 0.3713128
#> 29      3        0.3        V1 0.3573684
#> 30      3        0.3        V2 0.3872073
#> 31      4        0.4        V1 0.3567907
#> 32      4        0.4        V9 0.3356388
#> 33      4        0.4        V3 0.3715671
#> 34      4        0.4        V4 0.3686123
#> 35      4        0.4        V5 0.4079093
#> 36      4        0.4        V2 0.3866715
#> 37      4        0.4        V6 0.3917844
#> 38      4        0.4        V7 0.3539995
#> 39      4        0.4        V8 0.3446273
#> 40      4        0.4       V10 0.3708832
#> 41      5        0.5        V8 0.3440340
#> 42      5        0.5        V9 0.3351944
#> 43      5        0.5       V10 0.3704540
#> 44      5        0.5        V1 0.3562139
#> 45      5        0.5        V5 0.4073962
#> 46      5        0.5        V2 0.3861364
#> 47      5        0.5        V3 0.3709728
#> 48      5        0.5        V4 0.3681738
#> 49      5        0.5        V6 0.3912757
#> 50      5        0.5        V7 0.3534745
#> 51      6        0.6        V4 0.3677358
#> 52      6        0.6        V5 0.4068837
#> 53      6        0.6        V7 0.3529503
#> 54      6        0.6        V8 0.3434418
#> 55      6        0.6        V9 0.3347506
#> 56      6        0.6        V6 0.3907677
#> 57      6        0.6       V10 0.3700253
#> 58      6        0.6        V1 0.3556380
#> 59      6        0.6        V2 0.3856021
#> 60      6        0.6        V3 0.3703795
#> 61      7        0.7        V1 0.3550631
#> 62      7        0.7        V3 0.3697871
#> 63      7        0.7        V4 0.3672983
#> 64      7        0.7        V5 0.4063719
#> 65      7        0.7        V9 0.3343074
#> 66      7        0.7        V6 0.3902604
#> 67      7        0.7        V7 0.3524268
#> 68      7        0.7        V8 0.3428506
#> 69      7        0.7        V2 0.3850685
#> 70      7        0.7       V10 0.3695971
#> 71      8        0.8        V9 0.3338647
#> 72      8        0.8       V10 0.3691694
#> 73      8        0.8        V1 0.3544891
#> 74      8        0.8        V5 0.4058607
#> 75      8        0.8        V2 0.3845356
#> 76      8        0.8        V3 0.3691956
#> 77      8        0.8        V4 0.3668613
#> 78      8        0.8        V8 0.3422604
#> 79      8        0.8        V6 0.3897537
#> 80      8        0.8        V7 0.3519042
#> 81      9        0.9        V5 0.4053501
#> 82      9        0.9        V7 0.3513823
#> 83      9        0.9        V8 0.3416712
#> 84      9        0.9        V9 0.3334227
#> 85      9        0.9        V6 0.3892477
#> 86      9        0.9       V10 0.3687422
#> 87      9        0.9        V1 0.3539160
#> 88      9        0.9        V2 0.3840035
#> 89      9        0.9        V3 0.3686051
#> 90      9        0.9        V4 0.3664249
#> 91     10        1.0        V1 0.3533438
#> 92     10        1.0        V4 0.3659890
#> 93     10        1.0        V5 0.4048402
#> 94     10        1.0        V9 0.3329812
#> 95     10        1.0        V3 0.3680156
#> 96     10        1.0        V7 0.3508612
#> 97     10        1.0        V8 0.3410830
#> 98     10        1.0        V2 0.3834721
#> 99     10        1.0        V6 0.3887423
#> 100    10        1.0       V10 0.3683155


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
#> 1       0                0          0 0.8757906 0.04221558 0.8068574 0.9487313
#> 2       0               20         20 0.8757906 0.04221558 0.8068574 0.9487313
#> 3       0               40         40 0.8757906 0.04221558 0.8068574 0.9487313
#> 4       0               60         60 0.8757906 0.04221558 0.8068574 0.9487313
#> 5      20                0         20 0.8617131 0.04543158 0.7848729 0.9402485
#> 6      20               20         40 0.8617131 0.04543158 0.7848729 0.9402485
#> 7      20               40         60 0.8617131 0.04543158 0.7848729 0.9402485
#> 8      20               60         80 0.8617131 0.04543158 0.7848729 0.9402485
#> 9      40                0         40 0.8478591 0.04848639 0.7638002 0.9316363
#> 10     40               20         60 0.8478591 0.04848639 0.7638002 0.9316363
#> 11     40               40         80 0.8478591 0.04848639 0.7638002 0.9316363
#> 12     40               60        100 0.8478591 0.04848639 0.7638002 0.9316363
#> 13     60                0         60 0.8342249 0.05138459 0.7435324 0.9229291
#> 14     60               20         80 0.8342249 0.05138459 0.7435324 0.9229291
#> 15     60               40        100 0.8342249 0.05138459 0.7435324 0.9229291
#> 16     60               60        120 0.8342249 0.05138459 0.7435324 0.9229291
#> 17     80                0         80 0.8208071 0.05413301 0.7239894 0.9141537
#> 18     80               20        100 0.8208071 0.05413301 0.7239894 0.9141537
#> 19     80               40        120 0.8208071 0.05413301 0.7239894 0.9141537
#> 20     80               60        140 0.8208071 0.05413301 0.7239894 0.9141537
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11067001 0.180001216 0.5658347
#> 2  0.30574618 0.10247519 0.141987265 0.4912172
#> 3  0.26001915 0.09504251 0.111149845 0.4271100
#> 4  0.22113100 0.08818215 0.086209450 0.3722325
#> 5  0.25589195 0.10460060 0.100873896 0.4674801
#> 6  0.21762106 0.09516827 0.077924350 0.4067756
#> 7  0.18507391 0.08682520 0.059492994 0.3548451
#> 8  0.15739448 0.07934380 0.044796275 0.3104463
#> 9  0.18213629 0.09548146 0.053408951 0.3875312
#> 10 0.15489621 0.08609078 0.039975133 0.3383928
#> 11 0.13173012 0.07782870 0.029405659 0.2963739
#> 12 0.11202872 0.07048150 0.021193894 0.2603933
#> 13 0.12963921 0.08530912 0.025977711 0.3228243
#> 14 0.11025053 0.07654133 0.018558846 0.2830508
#> 15 0.09376159 0.06882831 0.012923802 0.2489666
#> 16 0.07973872 0.06198394 0.008732390 0.2196721
#> 17 0.09227334 0.07517709 0.011149813 0.2704336
#> 18 0.07847305 0.06726969 0.007436032 0.2381341
#> 19 0.06673672 0.06030250 0.004774717 0.2103377
#> 20 0.05675565 0.05411863 0.002931782 0.1863126
```
