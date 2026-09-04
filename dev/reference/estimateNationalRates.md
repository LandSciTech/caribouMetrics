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
#> 1         0.1 0.3838513 0.02486065 0.3380595 0.4061899
#> 2         0.2 0.3832760 0.02478857 0.3375432 0.4054337
#> 3         0.3 0.3827015 0.02471688 0.3370277 0.4046788
#> 4         0.4 0.3821279 0.02464558 0.3365130 0.4039254
#> 5         0.5 0.3815551 0.02457467 0.3359991 0.4031734
#> 6         0.6 0.3809832 0.02450415 0.3354860 0.4024228
#> 7         0.7 0.3804122 0.02443402 0.3349737 0.4016736
#> 8         0.8 0.3798420 0.02436427 0.3344622 0.4009258
#> 9         0.9 0.3792726 0.02429492 0.3339514 0.4001794
#> 10        1.0 0.3787041 0.02422595 0.3334414 0.3994344

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3838964
#> 2       1        0.1        V5 0.3986645
#> 3       1        0.1        V9 0.4083748
#> 4       1        0.1        V3 0.3822683
#> 5       1        0.1        V4 0.3401447
#> 6       1        0.1        V8 0.3798209
#> 7       1        0.1        V2 0.3828968
#> 8       1        0.1        V6 0.3374541
#> 9       1        0.1        V7 0.3491250
#> 10      1        0.1       V10 0.3945469
#> 11      2        0.2        V7 0.3486790
#> 12      2        0.2        V8 0.3793289
#> 13      2        0.2        V9 0.4075804
#> 14      2        0.2       V10 0.3938913
#> 15      2        0.2        V1 0.3834065
#> 16      2        0.2        V5 0.3980395
#> 17      2        0.2        V2 0.3823158
#> 18      2        0.2        V3 0.3817411
#> 19      2        0.2        V4 0.3396228
#> 20      2        0.2        V6 0.3369394
#> 21      3        0.3        V4 0.3391017
#> 22      3        0.3        V5 0.3974155
#> 23      3        0.3        V3 0.3812146
#> 24      3        0.3        V7 0.3482336
#> 25      3        0.3        V8 0.3788376
#> 26      3        0.3        V9 0.4067875
#> 27      3        0.3        V6 0.3364256
#> 28      3        0.3       V10 0.3932367
#> 29      3        0.3        V1 0.3829173
#> 30      3        0.3        V2 0.3817356
#> 31      4        0.4        V1 0.3824287
#> 32      4        0.4        V9 0.4059963
#> 33      4        0.4        V3 0.3806888
#> 34      4        0.4        V4 0.3385814
#> 35      4        0.4        V5 0.3967924
#> 36      4        0.4        V2 0.3811564
#> 37      4        0.4        V6 0.3359126
#> 38      4        0.4        V7 0.3477887
#> 39      4        0.4        V8 0.3783470
#> 40      4        0.4       V10 0.3925832
#> 41      5        0.5        V8 0.3778569
#> 42      5        0.5        V9 0.4052065
#> 43      5        0.5       V10 0.3919308
#> 44      5        0.5        V1 0.3819407
#> 45      5        0.5        V5 0.3961704
#> 46      5        0.5        V2 0.3805780
#> 47      5        0.5        V3 0.3801637
#> 48      5        0.5        V4 0.3380619
#> 49      5        0.5        V6 0.3354003
#> 50      5        0.5        V7 0.3473444
#> 51      6        0.6        V4 0.3375432
#> 52      6        0.6        V5 0.3955493
#> 53      6        0.6        V7 0.3469006
#> 54      6        0.6        V8 0.3773676
#> 55      6        0.6        V9 0.4044183
#> 56      6        0.6        V6 0.3348888
#> 57      6        0.6       V10 0.3912795
#> 58      6        0.6        V1 0.3814534
#> 59      6        0.6        V2 0.3800005
#> 60      6        0.6        V3 0.3796394
#> 61      7        0.7        V1 0.3809667
#> 62      7        0.7        V3 0.3791158
#> 63      7        0.7        V4 0.3370253
#> 64      7        0.7        V5 0.3949292
#> 65      7        0.7        V9 0.4036316
#> 66      7        0.7        V6 0.3343781
#> 67      7        0.7        V7 0.3464575
#> 68      7        0.7        V8 0.3768788
#> 69      7        0.7        V2 0.3794238
#> 70      7        0.7       V10 0.3906293
#> 71      8        0.8        V9 0.4028465
#> 72      8        0.8       V10 0.3899802
#> 73      8        0.8        V1 0.3804806
#> 74      8        0.8        V5 0.3943100
#> 75      8        0.8        V2 0.3788481
#> 76      8        0.8        V3 0.3785929
#> 77      8        0.8        V4 0.3365081
#> 78      8        0.8        V8 0.3763907
#> 79      8        0.8        V6 0.3338682
#> 80      8        0.8        V7 0.3460148
#> 81      9        0.9        V5 0.3936919
#> 82      9        0.9        V7 0.3455728
#> 83      9        0.9        V8 0.3759032
#> 84      9        0.9        V9 0.4020629
#> 85      9        0.9        V6 0.3333590
#> 86      9        0.9       V10 0.3893321
#> 87      9        0.9        V1 0.3799951
#> 88      9        0.9        V2 0.3782732
#> 89      9        0.9        V3 0.3780707
#> 90      9        0.9        V4 0.3359918
#> 91     10        1.0        V1 0.3795102
#> 92     10        1.0        V4 0.3354763
#> 93     10        1.0        V5 0.3930747
#> 94     10        1.0        V9 0.4012808
#> 95     10        1.0        V3 0.3775493
#> 96     10        1.0        V7 0.3451313
#> 97     10        1.0        V8 0.3754163
#> 98     10        1.0        V2 0.3776992
#> 99     10        1.0        V6 0.3328506
#> 100    10        1.0       V10 0.3886851

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3838964
#> 2       1        0.1        V5 0.3986645
#> 3       1        0.1        V9 0.4083748
#> 4       1        0.1        V3 0.3822683
#> 5       1        0.1        V4 0.3401447
#> 6       1        0.1        V8 0.3798209
#> 7       1        0.1        V2 0.3828968
#> 8       1        0.1        V6 0.3374541
#> 9       1        0.1        V7 0.3491250
#> 10      1        0.1       V10 0.3945469
#> 11      2        0.2        V7 0.3486790
#> 12      2        0.2        V8 0.3793289
#> 13      2        0.2        V9 0.4075804
#> 14      2        0.2       V10 0.3938913
#> 15      2        0.2        V1 0.3834065
#> 16      2        0.2        V5 0.3980395
#> 17      2        0.2        V2 0.3823158
#> 18      2        0.2        V3 0.3817411
#> 19      2        0.2        V4 0.3396228
#> 20      2        0.2        V6 0.3369394
#> 21      3        0.3        V4 0.3391017
#> 22      3        0.3        V5 0.3974155
#> 23      3        0.3        V3 0.3812146
#> 24      3        0.3        V7 0.3482336
#> 25      3        0.3        V8 0.3788376
#> 26      3        0.3        V9 0.4067875
#> 27      3        0.3        V6 0.3364256
#> 28      3        0.3       V10 0.3932367
#> 29      3        0.3        V1 0.3829173
#> 30      3        0.3        V2 0.3817356
#> 31      4        0.4        V1 0.3824287
#> 32      4        0.4        V9 0.4059963
#> 33      4        0.4        V3 0.3806888
#> 34      4        0.4        V4 0.3385814
#> 35      4        0.4        V5 0.3967924
#> 36      4        0.4        V2 0.3811564
#> 37      4        0.4        V6 0.3359126
#> 38      4        0.4        V7 0.3477887
#> 39      4        0.4        V8 0.3783470
#> 40      4        0.4       V10 0.3925832
#> 41      5        0.5        V8 0.3778569
#> 42      5        0.5        V9 0.4052065
#> 43      5        0.5       V10 0.3919308
#> 44      5        0.5        V1 0.3819407
#> 45      5        0.5        V5 0.3961704
#> 46      5        0.5        V2 0.3805780
#> 47      5        0.5        V3 0.3801637
#> 48      5        0.5        V4 0.3380619
#> 49      5        0.5        V6 0.3354003
#> 50      5        0.5        V7 0.3473444
#> 51      6        0.6        V4 0.3375432
#> 52      6        0.6        V5 0.3955493
#> 53      6        0.6        V7 0.3469006
#> 54      6        0.6        V8 0.3773676
#> 55      6        0.6        V9 0.4044183
#> 56      6        0.6        V6 0.3348888
#> 57      6        0.6       V10 0.3912795
#> 58      6        0.6        V1 0.3814534
#> 59      6        0.6        V2 0.3800005
#> 60      6        0.6        V3 0.3796394
#> 61      7        0.7        V1 0.3809667
#> 62      7        0.7        V3 0.3791158
#> 63      7        0.7        V4 0.3370253
#> 64      7        0.7        V5 0.3949292
#> 65      7        0.7        V9 0.4036316
#> 66      7        0.7        V6 0.3343781
#> 67      7        0.7        V7 0.3464575
#> 68      7        0.7        V8 0.3768788
#> 69      7        0.7        V2 0.3794238
#> 70      7        0.7       V10 0.3906293
#> 71      8        0.8        V9 0.4028465
#> 72      8        0.8       V10 0.3899802
#> 73      8        0.8        V1 0.3804806
#> 74      8        0.8        V5 0.3943100
#> 75      8        0.8        V2 0.3788481
#> 76      8        0.8        V3 0.3785929
#> 77      8        0.8        V4 0.3365081
#> 78      8        0.8        V8 0.3763907
#> 79      8        0.8        V6 0.3338682
#> 80      8        0.8        V7 0.3460148
#> 81      9        0.9        V5 0.3936919
#> 82      9        0.9        V7 0.3455728
#> 83      9        0.9        V8 0.3759032
#> 84      9        0.9        V9 0.4020629
#> 85      9        0.9        V6 0.3333590
#> 86      9        0.9       V10 0.3893321
#> 87      9        0.9        V1 0.3799951
#> 88      9        0.9        V2 0.3782732
#> 89      9        0.9        V3 0.3780707
#> 90      9        0.9        V4 0.3359918
#> 91     10        1.0        V1 0.3795102
#> 92     10        1.0        V4 0.3354763
#> 93     10        1.0        V5 0.3930747
#> 94     10        1.0        V9 0.4012808
#> 95     10        1.0        V3 0.3775493
#> 96     10        1.0        V7 0.3451313
#> 97     10        1.0        V8 0.3754163
#> 98     10        1.0        V2 0.3776992
#> 99     10        1.0        V6 0.3328506
#> 100    10        1.0       V10 0.3886851


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
#> 1       0                0          0 0.8757906 0.05177242 0.7744392 0.9430759
#> 2       0               20         20 0.8757906 0.05177242 0.7744392 0.9430759
#> 3       0               40         40 0.8757906 0.05177242 0.7744392 0.9430759
#> 4       0               60         60 0.8757906 0.05177242 0.7744392 0.9430759
#> 5      20                0         20 0.8617131 0.05345036 0.7643739 0.9365390
#> 6      20               20         40 0.8617131 0.05345036 0.7643739 0.9365390
#> 7      20               40         60 0.8617131 0.05345036 0.7643739 0.9365390
#> 8      20               60         80 0.8617131 0.05345036 0.7643739 0.9365390
#> 9      40                0         40 0.8478591 0.05522261 0.7545059 0.9299413
#> 10     40               20         60 0.8478591 0.05522261 0.7545059 0.9299413
#> 11     40               40         80 0.8478591 0.05522261 0.7545059 0.9299413
#> 12     40               60        100 0.8478591 0.05522261 0.7545059 0.9299413
#> 13     60                0         60 0.8342249 0.05706383 0.7448242 0.9232963
#> 14     60               20         80 0.8342249 0.05706383 0.7448242 0.9232963
#> 15     60               40        100 0.8342249 0.05706383 0.7448242 0.9232963
#> 16     60               60        120 0.8342249 0.05706383 0.7448242 0.9232963
#> 17     80                0         80 0.8208071 0.05895265 0.7353189 0.9166155
#> 18     80               20        100 0.8208071 0.05895265 0.7353189 0.9166155
#> 19     80               40        120 0.8208071 0.05895265 0.7353189 0.9166155
#> 20     80               60        140 0.8208071 0.05895265 0.7353189 0.9166155
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12366095 0.173606006 0.5623918
#> 2  0.30574618 0.11427514 0.128708031 0.4934027
#> 3  0.26001915 0.10582206 0.094199323 0.4333893
#> 4  0.22113100 0.09797626 0.067836995 0.3813537
#> 5  0.25589195 0.11306828 0.102257123 0.4502102
#> 6  0.21762106 0.10319778 0.073973230 0.3959275
#> 7  0.18507391 0.09437147 0.052509539 0.3489160
#> 8  0.15739448 0.08634316 0.036417410 0.3082232
#> 9  0.18213629 0.10041214 0.057489509 0.3620799
#> 10 0.15489621 0.09105864 0.040129885 0.3196182
#> 11 0.13173012 0.08270500 0.027265278 0.2828544
#> 12 0.11202872 0.07516278 0.017920846 0.2509894
#> 13 0.12963921 0.08749769 0.030217130 0.2931515
#> 14 0.11025053 0.07906002 0.020045117 0.2599194
#> 15 0.09376159 0.07152285 0.012793863 0.2310803
#> 16 0.07973872 0.06473929 0.007787083 0.2059950
#> 17 0.09227334 0.07523645 0.014427955 0.2391667
#> 18 0.07847305 0.06783088 0.008898660 0.2130358
#> 19 0.06673672 0.06121137 0.005192528 0.1902596
#> 20 0.05675565 0.05526189 0.002829990 0.1703381
```
