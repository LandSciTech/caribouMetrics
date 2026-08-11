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
#> 1         0.1 0.3838513 0.02838848 0.3375059 0.4260076
#> 2         0.2 0.3832760 0.02831648 0.3370853 0.4253734
#> 3         0.3 0.3827015 0.02824468 0.3366652 0.4247401
#> 4         0.4 0.3821279 0.02817307 0.3362457 0.4241078
#> 5         0.5 0.3815551 0.02810167 0.3358267 0.4234764
#> 6         0.6 0.3809832 0.02803047 0.3354082 0.4228460
#> 7         0.7 0.3804122 0.02795947 0.3349902 0.4222165
#> 8         0.8 0.3798420 0.02788867 0.3345727 0.4215879
#> 9         0.9 0.3792726 0.02781807 0.3341558 0.4209603
#> 10        1.0 0.3787041 0.02774767 0.3337394 0.4203336

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4023259
#> 2       1        0.1        V5 0.3467962
#> 3       1        0.1        V9 0.3348087
#> 4       1        0.1        V3 0.3763585
#> 5       1        0.1        V4 0.3946340
#> 6       1        0.1        V8 0.3716749
#> 7       1        0.1        V2 0.3816536
#> 8       1        0.1        V6 0.4065694
#> 9       1        0.1        V7 0.4316510
#> 10      1        0.1       V10 0.3869708
#> 11      2        0.2        V7 0.4309979
#> 12      2        0.2        V8 0.3711981
#> 13      2        0.2        V9 0.3343970
#> 14      2        0.2       V10 0.3863633
#> 15      2        0.2        V1 0.4016983
#> 16      2        0.2        V5 0.3463450
#> 17      2        0.2        V2 0.3810258
#> 18      2        0.2        V3 0.3757998
#> 19      2        0.2        V4 0.3940566
#> 20      2        0.2        V6 0.4060000
#> 21      3        0.3        V4 0.3934800
#> 22      3        0.3        V5 0.3458943
#> 23      3        0.3        V3 0.3752420
#> 24      3        0.3        V7 0.4303459
#> 25      3        0.3        V8 0.3707219
#> 26      3        0.3        V9 0.3339858
#> 27      3        0.3        V6 0.4054314
#> 28      3        0.3       V10 0.3857567
#> 29      3        0.3        V1 0.4010717
#> 30      3        0.3        V2 0.3803991
#> 31      4        0.4        V1 0.4004460
#> 32      4        0.4        V9 0.3335751
#> 33      4        0.4        V3 0.3746850
#> 34      4        0.4        V4 0.3929042
#> 35      4        0.4        V5 0.3454443
#> 36      4        0.4        V2 0.3797733
#> 37      4        0.4        V6 0.4048635
#> 38      4        0.4        V7 0.4296949
#> 39      4        0.4        V8 0.3702464
#> 40      4        0.4       V10 0.3851511
#> 41      5        0.5        V8 0.3697714
#> 42      5        0.5        V9 0.3331649
#> 43      5        0.5       V10 0.3845464
#> 44      5        0.5        V1 0.3998213
#> 45      5        0.5        V5 0.3449948
#> 46      5        0.5        V2 0.3791487
#> 47      5        0.5        V3 0.3741289
#> 48      5        0.5        V4 0.3923293
#> 49      5        0.5        V6 0.4042965
#> 50      5        0.5        V7 0.4290448
#> 51      6        0.6        V4 0.3917553
#> 52      6        0.6        V5 0.3445459
#> 53      6        0.6        V7 0.4283957
#> 54      6        0.6        V8 0.3692971
#> 55      6        0.6        V9 0.3327553
#> 56      6        0.6        V6 0.4037303
#> 57      6        0.6       V10 0.3839427
#> 58      6        0.6        V1 0.3991976
#> 59      6        0.6        V2 0.3785250
#> 60      6        0.6        V3 0.3735735
#> 61      7        0.7        V1 0.3985748
#> 62      7        0.7        V3 0.3730190
#> 63      7        0.7        V4 0.3911820
#> 64      7        0.7        V5 0.3440976
#> 65      7        0.7        V9 0.3323461
#> 66      7        0.7        V6 0.4031649
#> 67      7        0.7        V7 0.4277476
#> 68      7        0.7        V8 0.3688233
#> 69      7        0.7        V2 0.3779024
#> 70      7        0.7       V10 0.3833400
#> 71      8        0.8        V9 0.3319374
#> 72      8        0.8       V10 0.3827381
#> 73      8        0.8        V1 0.3979530
#> 74      8        0.8        V5 0.3436499
#> 75      8        0.8        V2 0.3772807
#> 76      8        0.8        V3 0.3724653
#> 77      8        0.8        V4 0.3906096
#> 78      8        0.8        V8 0.3683502
#> 79      8        0.8        V6 0.4026002
#> 80      8        0.8        V7 0.4271005
#> 81      9        0.9        V5 0.3432028
#> 82      9        0.9        V7 0.4264544
#> 83      9        0.9        V8 0.3678777
#> 84      9        0.9        V9 0.3315293
#> 85      9        0.9        V6 0.4020364
#> 86      9        0.9       V10 0.3821373
#> 87      9        0.9        V1 0.3973322
#> 88      9        0.9        V2 0.3766602
#> 89      9        0.9        V3 0.3719124
#> 90      9        0.9        V4 0.3900381
#> 91     10        1.0        V1 0.3967124
#> 92     10        1.0        V4 0.3894674
#> 93     10        1.0        V5 0.3427562
#> 94     10        1.0        V9 0.3311216
#> 95     10        1.0        V3 0.3713604
#> 96     10        1.0        V7 0.4258092
#> 97     10        1.0        V8 0.3674057
#> 98     10        1.0        V2 0.3760406
#> 99     10        1.0        V6 0.4014733
#> 100    10        1.0       V10 0.3815373

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4023259
#> 2       1        0.1        V5 0.3467962
#> 3       1        0.1        V9 0.3348087
#> 4       1        0.1        V3 0.3763585
#> 5       1        0.1        V4 0.3946340
#> 6       1        0.1        V8 0.3716749
#> 7       1        0.1        V2 0.3816536
#> 8       1        0.1        V6 0.4065694
#> 9       1        0.1        V7 0.4316510
#> 10      1        0.1       V10 0.3869708
#> 11      2        0.2        V7 0.4309979
#> 12      2        0.2        V8 0.3711981
#> 13      2        0.2        V9 0.3343970
#> 14      2        0.2       V10 0.3863633
#> 15      2        0.2        V1 0.4016983
#> 16      2        0.2        V5 0.3463450
#> 17      2        0.2        V2 0.3810258
#> 18      2        0.2        V3 0.3757998
#> 19      2        0.2        V4 0.3940566
#> 20      2        0.2        V6 0.4060000
#> 21      3        0.3        V4 0.3934800
#> 22      3        0.3        V5 0.3458943
#> 23      3        0.3        V3 0.3752420
#> 24      3        0.3        V7 0.4303459
#> 25      3        0.3        V8 0.3707219
#> 26      3        0.3        V9 0.3339858
#> 27      3        0.3        V6 0.4054314
#> 28      3        0.3       V10 0.3857567
#> 29      3        0.3        V1 0.4010717
#> 30      3        0.3        V2 0.3803991
#> 31      4        0.4        V1 0.4004460
#> 32      4        0.4        V9 0.3335751
#> 33      4        0.4        V3 0.3746850
#> 34      4        0.4        V4 0.3929042
#> 35      4        0.4        V5 0.3454443
#> 36      4        0.4        V2 0.3797733
#> 37      4        0.4        V6 0.4048635
#> 38      4        0.4        V7 0.4296949
#> 39      4        0.4        V8 0.3702464
#> 40      4        0.4       V10 0.3851511
#> 41      5        0.5        V8 0.3697714
#> 42      5        0.5        V9 0.3331649
#> 43      5        0.5       V10 0.3845464
#> 44      5        0.5        V1 0.3998213
#> 45      5        0.5        V5 0.3449948
#> 46      5        0.5        V2 0.3791487
#> 47      5        0.5        V3 0.3741289
#> 48      5        0.5        V4 0.3923293
#> 49      5        0.5        V6 0.4042965
#> 50      5        0.5        V7 0.4290448
#> 51      6        0.6        V4 0.3917553
#> 52      6        0.6        V5 0.3445459
#> 53      6        0.6        V7 0.4283957
#> 54      6        0.6        V8 0.3692971
#> 55      6        0.6        V9 0.3327553
#> 56      6        0.6        V6 0.4037303
#> 57      6        0.6       V10 0.3839427
#> 58      6        0.6        V1 0.3991976
#> 59      6        0.6        V2 0.3785250
#> 60      6        0.6        V3 0.3735735
#> 61      7        0.7        V1 0.3985748
#> 62      7        0.7        V3 0.3730190
#> 63      7        0.7        V4 0.3911820
#> 64      7        0.7        V5 0.3440976
#> 65      7        0.7        V9 0.3323461
#> 66      7        0.7        V6 0.4031649
#> 67      7        0.7        V7 0.4277476
#> 68      7        0.7        V8 0.3688233
#> 69      7        0.7        V2 0.3779024
#> 70      7        0.7       V10 0.3833400
#> 71      8        0.8        V9 0.3319374
#> 72      8        0.8       V10 0.3827381
#> 73      8        0.8        V1 0.3979530
#> 74      8        0.8        V5 0.3436499
#> 75      8        0.8        V2 0.3772807
#> 76      8        0.8        V3 0.3724653
#> 77      8        0.8        V4 0.3906096
#> 78      8        0.8        V8 0.3683502
#> 79      8        0.8        V6 0.4026002
#> 80      8        0.8        V7 0.4271005
#> 81      9        0.9        V5 0.3432028
#> 82      9        0.9        V7 0.4264544
#> 83      9        0.9        V8 0.3678777
#> 84      9        0.9        V9 0.3315293
#> 85      9        0.9        V6 0.4020364
#> 86      9        0.9       V10 0.3821373
#> 87      9        0.9        V1 0.3973322
#> 88      9        0.9        V2 0.3766602
#> 89      9        0.9        V3 0.3719124
#> 90      9        0.9        V4 0.3900381
#> 91     10        1.0        V1 0.3967124
#> 92     10        1.0        V4 0.3894674
#> 93     10        1.0        V5 0.3427562
#> 94     10        1.0        V9 0.3311216
#> 95     10        1.0        V3 0.3713604
#> 96     10        1.0        V7 0.4258092
#> 97     10        1.0        V8 0.3674057
#> 98     10        1.0        V2 0.3760406
#> 99     10        1.0        V6 0.4014733
#> 100    10        1.0       V10 0.3815373


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
#> 1       0                0          0 0.8757906 0.04696881 0.7968480 0.9453410
#> 2       0               20         20 0.8757906 0.04696881 0.7968480 0.9453410
#> 3       0               40         40 0.8757906 0.04696881 0.7968480 0.9453410
#> 4       0               60         60 0.8757906 0.04696881 0.7968480 0.9453410
#> 5      20                0         20 0.8617131 0.04882297 0.7795709 0.9337497
#> 6      20               20         40 0.8617131 0.04882297 0.7795709 0.9337497
#> 7      20               40         60 0.8617131 0.04882297 0.7795709 0.9337497
#> 8      20               60         80 0.8617131 0.04882297 0.7795709 0.9337497
#> 9      40                0         40 0.8478591 0.05055101 0.7628507 0.9219877
#> 10     40               20         60 0.8478591 0.05055101 0.7628507 0.9219877
#> 11     40               40         80 0.8478591 0.05055101 0.7628507 0.9219877
#> 12     40               60        100 0.8478591 0.05055101 0.7628507 0.9219877
#> 13     60                0         60 0.8342249 0.05216538 0.7466368 0.9101207
#> 14     60               20         80 0.8342249 0.05216538 0.7466368 0.9101207
#> 15     60               40        100 0.8342249 0.05216538 0.7466368 0.9101207
#> 16     60               60        120 0.8342249 0.05216538 0.7466368 0.9101207
#> 17     80                0         80 0.8208071 0.05367620 0.7308890 0.8981976
#> 18     80               20        100 0.8208071 0.05367620 0.7308890 0.8981976
#> 19     80               40        120 0.8208071 0.05367620 0.7308890 0.8981976
#> 20     80               60        140 0.8208071 0.05367620 0.7308890 0.8981976
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12209410 0.148956991 0.5698200
#> 2  0.30574618 0.10949294 0.116464504 0.4797944
#> 3  0.26001915 0.09907536 0.090242220 0.4049630
#> 4  0.22113100 0.09024440 0.069174263 0.3519595
#> 5  0.25589195 0.10990466 0.078890085 0.4533865
#> 6  0.21762106 0.09778642 0.060094798 0.3830921
#> 7  0.18507391 0.08770583 0.045144863 0.3249551
#> 8  0.15739448 0.07913988 0.033361642 0.2768694
#> 9  0.18213629 0.09589668 0.038763735 0.3625586
#> 10 0.15489621 0.08501440 0.028375012 0.3079781
#> 11 0.13173012 0.07587590 0.020336650 0.2628097
#> 12 0.11202872 0.06806017 0.014216126 0.2253210
#> 13 0.12963921 0.08208176 0.016992735 0.2920353
#> 14 0.11025053 0.07263630 0.011707741 0.2495935
#> 15 0.09376159 0.06462700 0.007810352 0.2143195
#> 16 0.07973872 0.05772510 0.005014624 0.1848477
#> 17 0.09227334 0.06942078 0.006259842 0.2371646
#> 18 0.07847305 0.06137280 0.003930629 0.2039551
#> 19 0.06673672 0.05448357 0.002348620 0.1761494
#> 20 0.05675565 0.04849865 0.001323179 0.1526930
```
