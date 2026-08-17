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
#> 1         0.1 0.3838513 0.01443161 0.3627471 0.4040784
#> 2         0.2 0.3832760 0.01442847 0.3621762 0.4035010
#> 3         0.3 0.3827015 0.01442538 0.3616063 0.4029244
#> 4         0.4 0.3821279 0.01442235 0.3610372 0.4023486
#> 5         0.5 0.3815551 0.01441938 0.3604691 0.4017736
#> 6         0.6 0.3809832 0.01441647 0.3599018 0.4011995
#> 7         0.7 0.3804122 0.01441361 0.3593354 0.4006262
#> 8         0.8 0.3798420 0.01441080 0.3587699 0.4000537
#> 9         0.9 0.3792726 0.01440805 0.3582053 0.3994820
#> 10        1.0 0.3787041 0.01440536 0.3576416 0.3989111

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4024413
#> 2       1        0.1        V5 0.3883481
#> 3       1        0.1        V9 0.3775318
#> 4       1        0.1        V3 0.4015041
#> 5       1        0.1        V4 0.3799598
#> 6       1        0.1        V8 0.3791045
#> 7       1        0.1        V2 0.4045537
#> 8       1        0.1        V6 0.3596899
#> 9       1        0.1        V7 0.3732773
#> 10      1        0.1       V10 0.3817307
#> 11      2        0.2        V7 0.3726746
#> 12      2        0.2        V8 0.3785900
#> 13      2        0.2        V9 0.3769518
#> 14      2        0.2       V10 0.3810949
#> 15      2        0.2        V1 0.4018815
#> 16      2        0.2        V5 0.3877666
#> 17      2        0.2        V2 0.4039711
#> 18      2        0.2        V3 0.4008948
#> 19      2        0.2        V4 0.3793432
#> 20      2        0.2        V6 0.3591283
#> 21      3        0.3        V4 0.3787275
#> 22      3        0.3        V5 0.3871860
#> 23      3        0.3        V3 0.4002865
#> 24      3        0.3        V7 0.3720729
#> 25      3        0.3        V8 0.3780762
#> 26      3        0.3        V9 0.3763727
#> 27      3        0.3        V6 0.3585676
#> 28      3        0.3       V10 0.3804602
#> 29      3        0.3        V1 0.4013225
#> 30      3        0.3        V2 0.4033894
#> 31      4        0.4        V1 0.4007642
#> 32      4        0.4        V9 0.3757945
#> 33      4        0.4        V3 0.3996790
#> 34      4        0.4        V4 0.3781128
#> 35      4        0.4        V5 0.3866062
#> 36      4        0.4        V2 0.4028086
#> 37      4        0.4        V6 0.3580077
#> 38      4        0.4        V7 0.3714721
#> 39      4        0.4        V8 0.3775632
#> 40      4        0.4       V10 0.3798265
#> 41      5        0.5        V8 0.3770508
#> 42      5        0.5        V9 0.3752172
#> 43      5        0.5       V10 0.3791938
#> 44      5        0.5        V1 0.4002067
#> 45      5        0.5        V5 0.3860273
#> 46      5        0.5        V2 0.4022285
#> 47      5        0.5        V3 0.3990725
#> 48      5        0.5        V4 0.3774991
#> 49      5        0.5        V6 0.3574488
#> 50      5        0.5        V7 0.3708723
#> 51      6        0.6        V4 0.3768865
#> 52      6        0.6        V5 0.3854493
#> 53      6        0.6        V7 0.3702735
#> 54      6        0.6        V8 0.3765391
#> 55      6        0.6        V9 0.3746408
#> 56      6        0.6        V6 0.3568906
#> 57      6        0.6       V10 0.3785623
#> 58      6        0.6        V1 0.3996500
#> 59      6        0.6        V2 0.4016493
#> 60      6        0.6        V3 0.3984669
#> 61      7        0.7        V1 0.3990940
#> 62      7        0.7        V3 0.3978622
#> 63      7        0.7        V4 0.3762748
#> 64      7        0.7        V5 0.3848721
#> 65      7        0.7        V9 0.3740652
#> 66      7        0.7        V6 0.3563334
#> 67      7        0.7        V7 0.3696757
#> 68      7        0.7        V8 0.3760282
#> 69      7        0.7        V2 0.4010710
#> 70      7        0.7       V10 0.3779318
#> 71      8        0.8        V9 0.3734906
#> 72      8        0.8       V10 0.3773023
#> 73      8        0.8        V1 0.3985389
#> 74      8        0.8        V5 0.3842958
#> 75      8        0.8        V2 0.4004934
#> 76      8        0.8        V3 0.3972584
#> 77      8        0.8        V4 0.3756641
#> 78      8        0.8        V8 0.3755179
#> 79      8        0.8        V6 0.3557770
#> 80      8        0.8        V7 0.3690788
#> 81      9        0.9        V5 0.3837204
#> 82      9        0.9        V7 0.3684829
#> 83      9        0.9        V8 0.3750083
#> 84      9        0.9        V9 0.3729168
#> 85      9        0.9        V6 0.3552215
#> 86      9        0.9       V10 0.3766739
#> 87      9        0.9        V1 0.3979845
#> 88      9        0.9        V2 0.3999168
#> 89      9        0.9        V3 0.3966556
#> 90      9        0.9        V4 0.3750544
#> 91     10        1.0        V1 0.3974308
#> 92     10        1.0        V4 0.3744457
#> 93     10        1.0        V5 0.3831458
#> 94     10        1.0        V9 0.3723439
#> 95     10        1.0        V3 0.3960536
#> 96     10        1.0        V7 0.3678879
#> 97     10        1.0        V8 0.3744994
#> 98     10        1.0        V2 0.3993409
#> 99     10        1.0        V6 0.3546669
#> 100    10        1.0       V10 0.3760465

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4024413
#> 2       1        0.1        V5 0.3883481
#> 3       1        0.1        V9 0.3775318
#> 4       1        0.1        V3 0.4015041
#> 5       1        0.1        V4 0.3799598
#> 6       1        0.1        V8 0.3791045
#> 7       1        0.1        V2 0.4045537
#> 8       1        0.1        V6 0.3596899
#> 9       1        0.1        V7 0.3732773
#> 10      1        0.1       V10 0.3817307
#> 11      2        0.2        V7 0.3726746
#> 12      2        0.2        V8 0.3785900
#> 13      2        0.2        V9 0.3769518
#> 14      2        0.2       V10 0.3810949
#> 15      2        0.2        V1 0.4018815
#> 16      2        0.2        V5 0.3877666
#> 17      2        0.2        V2 0.4039711
#> 18      2        0.2        V3 0.4008948
#> 19      2        0.2        V4 0.3793432
#> 20      2        0.2        V6 0.3591283
#> 21      3        0.3        V4 0.3787275
#> 22      3        0.3        V5 0.3871860
#> 23      3        0.3        V3 0.4002865
#> 24      3        0.3        V7 0.3720729
#> 25      3        0.3        V8 0.3780762
#> 26      3        0.3        V9 0.3763727
#> 27      3        0.3        V6 0.3585676
#> 28      3        0.3       V10 0.3804602
#> 29      3        0.3        V1 0.4013225
#> 30      3        0.3        V2 0.4033894
#> 31      4        0.4        V1 0.4007642
#> 32      4        0.4        V9 0.3757945
#> 33      4        0.4        V3 0.3996790
#> 34      4        0.4        V4 0.3781128
#> 35      4        0.4        V5 0.3866062
#> 36      4        0.4        V2 0.4028086
#> 37      4        0.4        V6 0.3580077
#> 38      4        0.4        V7 0.3714721
#> 39      4        0.4        V8 0.3775632
#> 40      4        0.4       V10 0.3798265
#> 41      5        0.5        V8 0.3770508
#> 42      5        0.5        V9 0.3752172
#> 43      5        0.5       V10 0.3791938
#> 44      5        0.5        V1 0.4002067
#> 45      5        0.5        V5 0.3860273
#> 46      5        0.5        V2 0.4022285
#> 47      5        0.5        V3 0.3990725
#> 48      5        0.5        V4 0.3774991
#> 49      5        0.5        V6 0.3574488
#> 50      5        0.5        V7 0.3708723
#> 51      6        0.6        V4 0.3768865
#> 52      6        0.6        V5 0.3854493
#> 53      6        0.6        V7 0.3702735
#> 54      6        0.6        V8 0.3765391
#> 55      6        0.6        V9 0.3746408
#> 56      6        0.6        V6 0.3568906
#> 57      6        0.6       V10 0.3785623
#> 58      6        0.6        V1 0.3996500
#> 59      6        0.6        V2 0.4016493
#> 60      6        0.6        V3 0.3984669
#> 61      7        0.7        V1 0.3990940
#> 62      7        0.7        V3 0.3978622
#> 63      7        0.7        V4 0.3762748
#> 64      7        0.7        V5 0.3848721
#> 65      7        0.7        V9 0.3740652
#> 66      7        0.7        V6 0.3563334
#> 67      7        0.7        V7 0.3696757
#> 68      7        0.7        V8 0.3760282
#> 69      7        0.7        V2 0.4010710
#> 70      7        0.7       V10 0.3779318
#> 71      8        0.8        V9 0.3734906
#> 72      8        0.8       V10 0.3773023
#> 73      8        0.8        V1 0.3985389
#> 74      8        0.8        V5 0.3842958
#> 75      8        0.8        V2 0.4004934
#> 76      8        0.8        V3 0.3972584
#> 77      8        0.8        V4 0.3756641
#> 78      8        0.8        V8 0.3755179
#> 79      8        0.8        V6 0.3557770
#> 80      8        0.8        V7 0.3690788
#> 81      9        0.9        V5 0.3837204
#> 82      9        0.9        V7 0.3684829
#> 83      9        0.9        V8 0.3750083
#> 84      9        0.9        V9 0.3729168
#> 85      9        0.9        V6 0.3552215
#> 86      9        0.9       V10 0.3766739
#> 87      9        0.9        V1 0.3979845
#> 88      9        0.9        V2 0.3999168
#> 89      9        0.9        V3 0.3966556
#> 90      9        0.9        V4 0.3750544
#> 91     10        1.0        V1 0.3974308
#> 92     10        1.0        V4 0.3744457
#> 93     10        1.0        V5 0.3831458
#> 94     10        1.0        V9 0.3723439
#> 95     10        1.0        V3 0.3960536
#> 96     10        1.0        V7 0.3678879
#> 97     10        1.0        V8 0.3744994
#> 98     10        1.0        V2 0.3993409
#> 99     10        1.0        V6 0.3546669
#> 100    10        1.0       V10 0.3760465


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
#> 1       0                0          0 0.8757906 0.04364040 0.7865739 0.9458194
#> 2       0               20         20 0.8757906 0.04364040 0.7865739 0.9458194
#> 3       0               40         40 0.8757906 0.04364040 0.7865739 0.9458194
#> 4       0               60         60 0.8757906 0.04364040 0.7865739 0.9458194
#> 5      20                0         20 0.8617131 0.04629427 0.7675205 0.9358998
#> 6      20               20         40 0.8617131 0.04629427 0.7675205 0.9358998
#> 7      20               40         60 0.8617131 0.04629427 0.7675205 0.9358998
#> 8      20               60         80 0.8617131 0.04629427 0.7675205 0.9358998
#> 9      40                0         40 0.8478591 0.04872846 0.7491540 0.9258222
#> 10     40               20         60 0.8478591 0.04872846 0.7491540 0.9258222
#> 11     40               40         80 0.8478591 0.04872846 0.7491540 0.9258222
#> 12     40               60        100 0.8478591 0.04872846 0.7491540 0.9258222
#> 13     60                0         60 0.8342249 0.05097052 0.7314066 0.9156345
#> 14     60               20         80 0.8342249 0.05097052 0.7314066 0.9156345
#> 15     60               40        100 0.8342249 0.05097052 0.7314066 0.9156345
#> 16     60               60        120 0.8342249 0.05097052 0.7314066 0.9156345
#> 17     80                0         80 0.8208071 0.05304237 0.7142250 0.9053740
#> 18     80               20        100 0.8208071 0.05304237 0.7142250 0.9053740
#> 19     80               40        120 0.8208071 0.05304237 0.7142250 0.9053740
#> 20     80               60        140 0.8208071 0.05304237 0.7142250 0.9053740
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12956798 0.1494724414 0.5767944
#> 2  0.30574618 0.11883715 0.1080558177 0.5060157
#> 3  0.26001915 0.10897311 0.0768385165 0.4444023
#> 4  0.22113100 0.09985918 0.0535199150 0.3909671
#> 5  0.25589195 0.11480440 0.0794807578 0.4505910
#> 6  0.21762106 0.10431102 0.0554827226 0.3963287
#> 7  0.18507391 0.09486573 0.0377650072 0.3493440
#> 8  0.15739448 0.08631360 0.0249116372 0.3086797
#> 9  0.18213629 0.09910318 0.0392461562 0.3540567
#> 10 0.15489621 0.08964350 0.0259752663 0.3127590
#> 11 0.13173012 0.08119113 0.0165479140 0.2769991
#> 12 0.11202872 0.07360182 0.0100523353 0.2459918
#> 13 0.12963921 0.08420421 0.0173186253 0.2805880
#> 14 0.11025053 0.07600770 0.0105740468 0.2491065
#> 15 0.09376159 0.06870546 0.0060914041 0.2217572
#> 16 0.07973872 0.06217350 0.0032630377 0.1979311
#> 17 0.09227334 0.07080369 0.0064439524 0.2245074
#> 18 0.07847305 0.06386145 0.0034787875 0.2003305
#> 19 0.06673672 0.05768105 0.0017182771 0.1792023
#> 20 0.05675565 0.05215571 0.0007596943 0.1606615
```
