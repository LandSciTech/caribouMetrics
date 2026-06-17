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
#> 1         0.1 0.3838513 0.03804201 0.3468195 0.4573875
#> 2         0.2 0.3832760 0.03801147 0.3462180 0.4567306
#> 3         0.3 0.3827015 0.03798100 0.3456176 0.4560746
#> 4         0.4 0.3821279 0.03795062 0.3450183 0.4554197
#> 5         0.5 0.3815551 0.03792032 0.3444200 0.4547656
#> 6         0.6 0.3809832 0.03789009 0.3438227 0.4541125
#> 7         0.7 0.3804122 0.03785994 0.3432265 0.4534603
#> 8         0.8 0.3798420 0.03782987 0.3426313 0.4528091
#> 9         0.9 0.3792726 0.03779988 0.3420372 0.4521588
#> 10        1.0 0.3787041 0.03776996 0.3414441 0.4515095

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4357627
#> 2       1        0.1        V5 0.3678916
#> 3       1        0.1        V9 0.4065422
#> 4       1        0.1        V3 0.4171930
#> 5       1        0.1        V4 0.3677322
#> 6       1        0.1        V8 0.3811218
#> 7       1        0.1        V2 0.4077981
#> 8       1        0.1        V6 0.3446304
#> 9       1        0.1        V7 0.4636656
#> 10      1        0.1       V10 0.3543596
#> 11      2        0.2        V7 0.4630097
#> 12      2        0.2        V8 0.3805937
#> 13      2        0.2        V9 0.4059970
#> 14      2        0.2       V10 0.3538606
#> 15      2        0.2        V1 0.4351025
#> 16      2        0.2        V5 0.3673795
#> 17      2        0.2        V2 0.4071769
#> 18      2        0.2        V3 0.4165882
#> 19      2        0.2        V4 0.3670790
#> 20      2        0.2        V6 0.3439992
#> 21      3        0.3        V4 0.3664269
#> 22      3        0.3        V5 0.3668682
#> 23      3        0.3        V3 0.4159843
#> 24      3        0.3        V7 0.4623547
#> 25      3        0.3        V8 0.3800663
#> 26      3        0.3        V9 0.4054526
#> 27      3        0.3        V6 0.3433692
#> 28      3        0.3       V10 0.3533623
#> 29      3        0.3        V1 0.4344432
#> 30      3        0.3        V2 0.4065566
#> 31      4        0.4        V1 0.4337850
#> 32      4        0.4        V9 0.4049088
#> 33      4        0.4        V3 0.4153812
#> 34      4        0.4        V4 0.3657760
#> 35      4        0.4        V5 0.3663575
#> 36      4        0.4        V2 0.4059374
#> 37      4        0.4        V6 0.3427403
#> 38      4        0.4        V7 0.4617007
#> 39      4        0.4        V8 0.3795396
#> 40      4        0.4       V10 0.3528646
#> 41      5        0.5        V8 0.3790137
#> 42      5        0.5        V9 0.4043659
#> 43      5        0.5       V10 0.3523677
#> 44      5        0.5        V1 0.4331277
#> 45      5        0.5        V5 0.3658476
#> 46      5        0.5        V2 0.4053190
#> 47      5        0.5        V3 0.4147790
#> 48      5        0.5        V4 0.3651262
#> 49      5        0.5        V6 0.3421126
#> 50      5        0.5        V7 0.4610476
#> 51      6        0.6        V4 0.3644776
#> 52      6        0.6        V5 0.3653384
#> 53      6        0.6        V7 0.4603954
#> 54      6        0.6        V8 0.3784884
#> 55      6        0.6        V9 0.4038236
#> 56      6        0.6        V6 0.3414860
#> 57      6        0.6       V10 0.3518715
#> 58      6        0.6        V1 0.4324715
#> 59      6        0.6        V2 0.4047016
#> 60      6        0.6        V3 0.4141777
#> 61      7        0.7        V1 0.4318162
#> 62      7        0.7        V3 0.4135773
#> 63      7        0.7        V4 0.3638302
#> 64      7        0.7        V5 0.3648299
#> 65      7        0.7        V9 0.4032821
#> 66      7        0.7        V6 0.3408605
#> 67      7        0.7        V7 0.4597441
#> 68      7        0.7        V8 0.3779640
#> 69      7        0.7        V2 0.4040851
#> 70      7        0.7       V10 0.3513760
#> 71      8        0.8        V9 0.4027412
#> 72      8        0.8       V10 0.3508811
#> 73      8        0.8        V1 0.4311619
#> 74      8        0.8        V5 0.3643221
#> 75      8        0.8        V2 0.4034696
#> 76      8        0.8        V3 0.4129777
#> 77      8        0.8        V4 0.3631839
#> 78      8        0.8        V8 0.3774402
#> 79      8        0.8        V6 0.3402362
#> 80      8        0.8        V7 0.4590938
#> 81      9        0.9        V5 0.3638150
#> 82      9        0.9        V7 0.4584444
#> 83      9        0.9        V8 0.3769172
#> 84      9        0.9        V9 0.4022012
#> 85      9        0.9        V6 0.3396131
#> 86      9        0.9       V10 0.3503870
#> 87      9        0.9        V1 0.4305086
#> 88      9        0.9        V2 0.4028550
#> 89      9        0.9        V3 0.4123790
#> 90      9        0.9        V4 0.3625387
#> 91     10        1.0        V1 0.4298563
#> 92     10        1.0        V4 0.3618947
#> 93     10        1.0        V5 0.3633086
#> 94     10        1.0        V9 0.4016618
#> 95     10        1.0        V3 0.4117812
#> 96     10        1.0        V7 0.4577958
#> 97     10        1.0        V8 0.3763949
#> 98     10        1.0        V2 0.4022414
#> 99     10        1.0        V6 0.3389911
#> 100    10        1.0       V10 0.3498936

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4357627
#> 2       1        0.1        V5 0.3678916
#> 3       1        0.1        V9 0.4065422
#> 4       1        0.1        V3 0.4171930
#> 5       1        0.1        V4 0.3677322
#> 6       1        0.1        V8 0.3811218
#> 7       1        0.1        V2 0.4077981
#> 8       1        0.1        V6 0.3446304
#> 9       1        0.1        V7 0.4636656
#> 10      1        0.1       V10 0.3543596
#> 11      2        0.2        V7 0.4630097
#> 12      2        0.2        V8 0.3805937
#> 13      2        0.2        V9 0.4059970
#> 14      2        0.2       V10 0.3538606
#> 15      2        0.2        V1 0.4351025
#> 16      2        0.2        V5 0.3673795
#> 17      2        0.2        V2 0.4071769
#> 18      2        0.2        V3 0.4165882
#> 19      2        0.2        V4 0.3670790
#> 20      2        0.2        V6 0.3439992
#> 21      3        0.3        V4 0.3664269
#> 22      3        0.3        V5 0.3668682
#> 23      3        0.3        V3 0.4159843
#> 24      3        0.3        V7 0.4623547
#> 25      3        0.3        V8 0.3800663
#> 26      3        0.3        V9 0.4054526
#> 27      3        0.3        V6 0.3433692
#> 28      3        0.3       V10 0.3533623
#> 29      3        0.3        V1 0.4344432
#> 30      3        0.3        V2 0.4065566
#> 31      4        0.4        V1 0.4337850
#> 32      4        0.4        V9 0.4049088
#> 33      4        0.4        V3 0.4153812
#> 34      4        0.4        V4 0.3657760
#> 35      4        0.4        V5 0.3663575
#> 36      4        0.4        V2 0.4059374
#> 37      4        0.4        V6 0.3427403
#> 38      4        0.4        V7 0.4617007
#> 39      4        0.4        V8 0.3795396
#> 40      4        0.4       V10 0.3528646
#> 41      5        0.5        V8 0.3790137
#> 42      5        0.5        V9 0.4043659
#> 43      5        0.5       V10 0.3523677
#> 44      5        0.5        V1 0.4331277
#> 45      5        0.5        V5 0.3658476
#> 46      5        0.5        V2 0.4053190
#> 47      5        0.5        V3 0.4147790
#> 48      5        0.5        V4 0.3651262
#> 49      5        0.5        V6 0.3421126
#> 50      5        0.5        V7 0.4610476
#> 51      6        0.6        V4 0.3644776
#> 52      6        0.6        V5 0.3653384
#> 53      6        0.6        V7 0.4603954
#> 54      6        0.6        V8 0.3784884
#> 55      6        0.6        V9 0.4038236
#> 56      6        0.6        V6 0.3414860
#> 57      6        0.6       V10 0.3518715
#> 58      6        0.6        V1 0.4324715
#> 59      6        0.6        V2 0.4047016
#> 60      6        0.6        V3 0.4141777
#> 61      7        0.7        V1 0.4318162
#> 62      7        0.7        V3 0.4135773
#> 63      7        0.7        V4 0.3638302
#> 64      7        0.7        V5 0.3648299
#> 65      7        0.7        V9 0.4032821
#> 66      7        0.7        V6 0.3408605
#> 67      7        0.7        V7 0.4597441
#> 68      7        0.7        V8 0.3779640
#> 69      7        0.7        V2 0.4040851
#> 70      7        0.7       V10 0.3513760
#> 71      8        0.8        V9 0.4027412
#> 72      8        0.8       V10 0.3508811
#> 73      8        0.8        V1 0.4311619
#> 74      8        0.8        V5 0.3643221
#> 75      8        0.8        V2 0.4034696
#> 76      8        0.8        V3 0.4129777
#> 77      8        0.8        V4 0.3631839
#> 78      8        0.8        V8 0.3774402
#> 79      8        0.8        V6 0.3402362
#> 80      8        0.8        V7 0.4590938
#> 81      9        0.9        V5 0.3638150
#> 82      9        0.9        V7 0.4584444
#> 83      9        0.9        V8 0.3769172
#> 84      9        0.9        V9 0.4022012
#> 85      9        0.9        V6 0.3396131
#> 86      9        0.9       V10 0.3503870
#> 87      9        0.9        V1 0.4305086
#> 88      9        0.9        V2 0.4028550
#> 89      9        0.9        V3 0.4123790
#> 90      9        0.9        V4 0.3625387
#> 91     10        1.0        V1 0.4298563
#> 92     10        1.0        V4 0.3618947
#> 93     10        1.0        V5 0.3633086
#> 94     10        1.0        V9 0.4016618
#> 95     10        1.0        V3 0.4117812
#> 96     10        1.0        V7 0.4577958
#> 97     10        1.0        V8 0.3763949
#> 98     10        1.0        V2 0.4022414
#> 99     10        1.0        V6 0.3389911
#> 100    10        1.0       V10 0.3498936


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
#> 1       0                0          0 0.8757906 0.04937462 0.7804384 0.9573003
#> 2       0               20         20 0.8757906 0.04937462 0.7804384 0.9573003
#> 3       0               40         40 0.8757906 0.04937462 0.7804384 0.9573003
#> 4       0               60         60 0.8757906 0.04937462 0.7804384 0.9573003
#> 5      20                0         20 0.8617131 0.05139258 0.7640129 0.9465428
#> 6      20               20         40 0.8617131 0.05139258 0.7640129 0.9465428
#> 7      20               40         60 0.8617131 0.05139258 0.7640129 0.9465428
#> 8      20               60         80 0.8617131 0.05139258 0.7640129 0.9465428
#> 9      40                0         40 0.8478591 0.05322694 0.7480822 0.9355227
#> 10     40               20         60 0.8478591 0.05322694 0.7480822 0.9355227
#> 11     40               40         80 0.8478591 0.05322694 0.7480822 0.9355227
#> 12     40               60        100 0.8478591 0.05322694 0.7480822 0.9355227
#> 13     60                0         60 0.8342249 0.05490272 0.7326063 0.9243247
#> 14     60               20         80 0.8342249 0.05490272 0.7326063 0.9243247
#> 15     60               40        100 0.8342249 0.05490272 0.7326063 0.9243247
#> 16     60               60        120 0.8342249 0.05490272 0.7326063 0.9243247
#> 17     80                0         80 0.8208071 0.05643934 0.7175528 0.9130109
#> 18     80               20        100 0.8208071 0.05643934 0.7175528 0.9130109
#> 19     80               40        120 0.8208071 0.05643934 0.7175528 0.9130109
#> 20     80               60        140 0.8208071 0.05643934 0.7175528 0.9130109
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13651958 0.139307863 0.5675289
#> 2  0.30574618 0.12831944 0.111035190 0.5068030
#> 3  0.26001915 0.12051868 0.087811674 0.4529661
#> 4  0.22113100 0.11297960 0.068810622 0.4053760
#> 5  0.25589195 0.11602258 0.081766772 0.4460869
#> 6  0.21762106 0.10814084 0.063881027 0.3993014
#> 7  0.18507391 0.10088906 0.049350704 0.3580107
#> 8  0.15739448 0.09408410 0.037630714 0.3215815
#> 9  0.18213629 0.09791230 0.045602329 0.3527424
#> 10 0.15489621 0.09083219 0.034624318 0.3169329
#> 11 0.13173012 0.08438978 0.025874495 0.2853226
#> 12 0.11202872 0.07841685 0.018981141 0.2573869
#> 13 0.12963921 0.08229095 0.023650773 0.2812864
#> 14 0.11025053 0.07612162 0.017244997 0.2538166
#> 15 0.09376159 0.07051925 0.012293188 0.2294917
#> 16 0.07973872 0.06534570 0.008534909 0.2079014
#> 17 0.09227334 0.06894298 0.011064327 0.2263785
#> 18 0.07847305 0.06364316 0.007615489 0.2051339
#> 19 0.06673672 0.05881876 0.005074856 0.1862161
#> 20 0.05675565 0.05436216 0.003256250 0.1693127
```
