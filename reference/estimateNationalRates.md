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
#> 1         0.1 0.3838513 0.02176614 0.3503156 0.4185517
#> 2         0.2 0.3832760 0.02174222 0.3497442 0.4179032
#> 3         0.3 0.3827015 0.02171839 0.3491738 0.4172557
#> 4         0.4 0.3821279 0.02169465 0.3486044 0.4166093
#> 5         0.5 0.3815551 0.02167100 0.3480358 0.4159638
#> 6         0.6 0.3809832 0.02164744 0.3474682 0.4153193
#> 7         0.7 0.3804122 0.02162397 0.3469015 0.4146759
#> 8         0.8 0.3798420 0.02160058 0.3463358 0.4140334
#> 9         0.9 0.3792726 0.02157727 0.3457709 0.4133919
#> 10        1.0 0.3787041 0.02155406 0.3452070 0.4127514

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4221954
#> 2       1        0.1        V5 0.3989312
#> 3       1        0.1        V9 0.3479295
#> 4       1        0.1        V3 0.3857562
#> 5       1        0.1        V4 0.4060012
#> 6       1        0.1        V8 0.3739189
#> 7       1        0.1        V2 0.3585344
#> 8       1        0.1        V6 0.3820225
#> 9       1        0.1        V7 0.3872709
#> 10      1        0.1       V10 0.3785931
#> 11      2        0.2        V7 0.3866518
#> 12      2        0.2        V8 0.3732744
#> 13      2        0.2        V9 0.3473628
#> 14      2        0.2       V10 0.3779948
#> 15      2        0.2        V1 0.4215420
#> 16      2        0.2        V5 0.3982793
#> 17      2        0.2        V2 0.3579470
#> 18      2        0.2        V3 0.3852135
#> 19      2        0.2        V4 0.4053696
#> 20      2        0.2        V6 0.3814944
#> 21      3        0.3        V4 0.4047390
#> 22      3        0.3        V5 0.3976285
#> 23      3        0.3        V3 0.3846716
#> 24      3        0.3        V7 0.3860338
#> 25      3        0.3        V8 0.3726310
#> 26      3        0.3        V9 0.3467971
#> 27      3        0.3        V6 0.3809671
#> 28      3        0.3       V10 0.3773975
#> 29      3        0.3        V1 0.4208896
#> 30      3        0.3        V2 0.3573604
#> 31      4        0.4        V1 0.4202382
#> 32      4        0.4        V9 0.3462323
#> 33      4        0.4        V3 0.3841304
#> 34      4        0.4        V4 0.4041094
#> 35      4        0.4        V5 0.3969788
#> 36      4        0.4        V2 0.3567749
#> 37      4        0.4        V6 0.3804404
#> 38      4        0.4        V7 0.3854167
#> 39      4        0.4        V8 0.3719887
#> 40      4        0.4       V10 0.3768010
#> 41      5        0.5        V8 0.3713476
#> 42      5        0.5        V9 0.3456684
#> 43      5        0.5       V10 0.3762056
#> 44      5        0.5        V1 0.4195879
#> 45      5        0.5        V5 0.3963301
#> 46      5        0.5        V2 0.3561903
#> 47      5        0.5        V3 0.3835900
#> 48      5        0.5        V4 0.4034808
#> 49      5        0.5        V6 0.3799145
#> 50      5        0.5        V7 0.3848006
#> 51      6        0.6        V4 0.4028532
#> 52      6        0.6        V5 0.3956825
#> 53      6        0.6        V7 0.3841855
#> 54      6        0.6        V8 0.3707075
#> 55      6        0.6        V9 0.3451055
#> 56      6        0.6        V6 0.3793893
#> 57      6        0.6       V10 0.3756110
#> 58      6        0.6        V1 0.4189385
#> 59      6        0.6        V2 0.3556066
#> 60      6        0.6        V3 0.3830504
#> 61      7        0.7        V1 0.4182902
#> 62      7        0.7        V3 0.3825115
#> 63      7        0.7        V4 0.4022265
#> 64      7        0.7        V5 0.3950359
#> 65      7        0.7        V9 0.3445434
#> 66      7        0.7        V6 0.3788649
#> 67      7        0.7        V7 0.3835714
#> 68      7        0.7        V8 0.3700686
#> 69      7        0.7        V2 0.3550239
#> 70      7        0.7       V10 0.3750174
#> 71      8        0.8        V9 0.3439823
#> 72      8        0.8       V10 0.3744248
#> 73      8        0.8        V1 0.4176429
#> 74      8        0.8        V5 0.3943904
#> 75      8        0.8        V2 0.3544422
#> 76      8        0.8        V3 0.3819733
#> 77      8        0.8        V4 0.4016008
#> 78      8        0.8        V8 0.3694307
#> 79      8        0.8        V6 0.3783412
#> 80      8        0.8        V7 0.3829583
#> 81      9        0.9        V5 0.3937460
#> 82      9        0.9        V7 0.3823461
#> 83      9        0.9        V8 0.3687940
#> 84      9        0.9        V9 0.3434221
#> 85      9        0.9        V6 0.3778182
#> 86      9        0.9       V10 0.3738331
#> 87      9        0.9        V1 0.4169965
#> 88      9        0.9        V2 0.3538614
#> 89      9        0.9        V3 0.3814360
#> 90      9        0.9        V4 0.4009761
#> 91     10        1.0        V1 0.4163512
#> 92     10        1.0        V4 0.4003524
#> 93     10        1.0        V5 0.3931026
#> 94     10        1.0        V9 0.3428628
#> 95     10        1.0        V3 0.3808993
#> 96     10        1.0        V7 0.3817349
#> 97     10        1.0        V8 0.3681583
#> 98     10        1.0        V2 0.3532816
#> 99     10        1.0        V6 0.3772959
#> 100    10        1.0       V10 0.3732423

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4221954
#> 2       1        0.1        V5 0.3989312
#> 3       1        0.1        V9 0.3479295
#> 4       1        0.1        V3 0.3857562
#> 5       1        0.1        V4 0.4060012
#> 6       1        0.1        V8 0.3739189
#> 7       1        0.1        V2 0.3585344
#> 8       1        0.1        V6 0.3820225
#> 9       1        0.1        V7 0.3872709
#> 10      1        0.1       V10 0.3785931
#> 11      2        0.2        V7 0.3866518
#> 12      2        0.2        V8 0.3732744
#> 13      2        0.2        V9 0.3473628
#> 14      2        0.2       V10 0.3779948
#> 15      2        0.2        V1 0.4215420
#> 16      2        0.2        V5 0.3982793
#> 17      2        0.2        V2 0.3579470
#> 18      2        0.2        V3 0.3852135
#> 19      2        0.2        V4 0.4053696
#> 20      2        0.2        V6 0.3814944
#> 21      3        0.3        V4 0.4047390
#> 22      3        0.3        V5 0.3976285
#> 23      3        0.3        V3 0.3846716
#> 24      3        0.3        V7 0.3860338
#> 25      3        0.3        V8 0.3726310
#> 26      3        0.3        V9 0.3467971
#> 27      3        0.3        V6 0.3809671
#> 28      3        0.3       V10 0.3773975
#> 29      3        0.3        V1 0.4208896
#> 30      3        0.3        V2 0.3573604
#> 31      4        0.4        V1 0.4202382
#> 32      4        0.4        V9 0.3462323
#> 33      4        0.4        V3 0.3841304
#> 34      4        0.4        V4 0.4041094
#> 35      4        0.4        V5 0.3969788
#> 36      4        0.4        V2 0.3567749
#> 37      4        0.4        V6 0.3804404
#> 38      4        0.4        V7 0.3854167
#> 39      4        0.4        V8 0.3719887
#> 40      4        0.4       V10 0.3768010
#> 41      5        0.5        V8 0.3713476
#> 42      5        0.5        V9 0.3456684
#> 43      5        0.5       V10 0.3762056
#> 44      5        0.5        V1 0.4195879
#> 45      5        0.5        V5 0.3963301
#> 46      5        0.5        V2 0.3561903
#> 47      5        0.5        V3 0.3835900
#> 48      5        0.5        V4 0.4034808
#> 49      5        0.5        V6 0.3799145
#> 50      5        0.5        V7 0.3848006
#> 51      6        0.6        V4 0.4028532
#> 52      6        0.6        V5 0.3956825
#> 53      6        0.6        V7 0.3841855
#> 54      6        0.6        V8 0.3707075
#> 55      6        0.6        V9 0.3451055
#> 56      6        0.6        V6 0.3793893
#> 57      6        0.6       V10 0.3756110
#> 58      6        0.6        V1 0.4189385
#> 59      6        0.6        V2 0.3556066
#> 60      6        0.6        V3 0.3830504
#> 61      7        0.7        V1 0.4182902
#> 62      7        0.7        V3 0.3825115
#> 63      7        0.7        V4 0.4022265
#> 64      7        0.7        V5 0.3950359
#> 65      7        0.7        V9 0.3445434
#> 66      7        0.7        V6 0.3788649
#> 67      7        0.7        V7 0.3835714
#> 68      7        0.7        V8 0.3700686
#> 69      7        0.7        V2 0.3550239
#> 70      7        0.7       V10 0.3750174
#> 71      8        0.8        V9 0.3439823
#> 72      8        0.8       V10 0.3744248
#> 73      8        0.8        V1 0.4176429
#> 74      8        0.8        V5 0.3943904
#> 75      8        0.8        V2 0.3544422
#> 76      8        0.8        V3 0.3819733
#> 77      8        0.8        V4 0.4016008
#> 78      8        0.8        V8 0.3694307
#> 79      8        0.8        V6 0.3783412
#> 80      8        0.8        V7 0.3829583
#> 81      9        0.9        V5 0.3937460
#> 82      9        0.9        V7 0.3823461
#> 83      9        0.9        V8 0.3687940
#> 84      9        0.9        V9 0.3434221
#> 85      9        0.9        V6 0.3778182
#> 86      9        0.9       V10 0.3738331
#> 87      9        0.9        V1 0.4169965
#> 88      9        0.9        V2 0.3538614
#> 89      9        0.9        V3 0.3814360
#> 90      9        0.9        V4 0.4009761
#> 91     10        1.0        V1 0.4163512
#> 92     10        1.0        V4 0.4003524
#> 93     10        1.0        V5 0.3931026
#> 94     10        1.0        V9 0.3428628
#> 95     10        1.0        V3 0.3808993
#> 96     10        1.0        V7 0.3817349
#> 97     10        1.0        V8 0.3681583
#> 98     10        1.0        V2 0.3532816
#> 99     10        1.0        V6 0.3772959
#> 100    10        1.0       V10 0.3732423


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
#> 1       0                0          0 0.8757906 0.04418423 0.7955674 0.9418946
#> 2       0               20         20 0.8757906 0.04418423 0.7955674 0.9418946
#> 3       0               40         40 0.8757906 0.04418423 0.7955674 0.9418946
#> 4       0               60         60 0.8757906 0.04418423 0.7955674 0.9418946
#> 5      20                0         20 0.8617131 0.04603990 0.7782885 0.9308895
#> 6      20               20         40 0.8617131 0.04603990 0.7782885 0.9308895
#> 7      20               40         60 0.8617131 0.04603990 0.7782885 0.9308895
#> 8      20               60         80 0.8617131 0.04603990 0.7782885 0.9308895
#> 9      40                0         40 0.8478591 0.04770651 0.7615637 0.9197407
#> 10     40               20         60 0.8478591 0.04770651 0.7615637 0.9197407
#> 11     40               40         80 0.8478591 0.04770651 0.7615637 0.9197407
#> 12     40               60        100 0.8478591 0.04770651 0.7615637 0.9197407
#> 13     60                0         60 0.8342249 0.04921081 0.7453435 0.9085011
#> 14     60               20         80 0.8342249 0.04921081 0.7453435 0.9085011
#> 15     60               40        100 0.8342249 0.04921081 0.7453435 0.9085011
#> 16     60               60        120 0.8342249 0.04921081 0.7453435 0.9085011
#> 17     80                0         80 0.8208071 0.05057382 0.7295880 0.8972112
#> 18     80               20        100 0.8208071 0.05057382 0.7295880 0.8972112
#> 19     80               40        120 0.8208071 0.05057382 0.7295880 0.8972112
#> 20     80               60        140 0.8208071 0.05057382 0.7295880 0.8972112
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12072951 0.196466193 0.5637835
#> 2  0.30574618 0.11796633 0.140625441 0.4991194
#> 3  0.26001915 0.11356126 0.099187739 0.4423002
#> 4  0.22113100 0.10792177 0.068641406 0.3925263
#> 5  0.25589195 0.10771695 0.120146726 0.4483047
#> 6  0.21762106 0.10278203 0.084058459 0.3977815
#> 7  0.18507391 0.09716002 0.057572664 0.3535844
#> 8  0.15739448 0.09106798 0.038384049 0.3149418
#> 9  0.18213629 0.09424452 0.070922508 0.3582493
#> 10 0.15489621 0.08853367 0.048019089 0.3190205
#> 11 0.13173012 0.08264121 0.031551283 0.2847137
#> 12 0.11202872 0.07668681 0.019957881 0.2546836
#> 13 0.12963921 0.08141768 0.039801886 0.2883358
#> 14 0.11025053 0.07565326 0.025730800 0.2578561
#> 15 0.09376159 0.06997670 0.015943627 0.2311408
#> 16 0.07973872 0.06445464 0.009358954 0.2076774
#> 17 0.09227334 0.06970579 0.020799986 0.2339652
#> 18 0.07847305 0.06425369 0.012594445 0.2101607
#> 19 0.06673672 0.05903080 0.007177809 0.1892067
#> 20 0.05675565 0.05407482 0.003784608 0.1707056
```
