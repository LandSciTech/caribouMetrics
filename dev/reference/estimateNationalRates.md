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
#> 1         0.1 0.3838513 0.02700043 0.3392418 0.4110038
#> 2         0.2 0.3832760 0.02694019 0.3387087 0.4103296
#> 3         0.3 0.3827015 0.02688019 0.3381764 0.4096564
#> 4         0.4 0.3821279 0.02682041 0.3376449 0.4089843
#> 5         0.5 0.3815551 0.02676086 0.3371142 0.4083134
#> 6         0.6 0.3809832 0.02670154 0.3365844 0.4076435
#> 7         0.7 0.3804122 0.02664245 0.3360554 0.4069748
#> 8         0.8 0.3798420 0.02658358 0.3355273 0.4063071
#> 9         0.9 0.3792726 0.02652495 0.3350000 0.4056406
#> 10        1.0 0.3787041 0.02646653 0.3344735 0.4049751

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3928838
#> 2       1        0.1        V5 0.4115061
#> 3       1        0.1        V9 0.3373241
#> 4       1        0.1        V3 0.3975685
#> 5       1        0.1        V4 0.3955646
#> 6       1        0.1        V8 0.3591620
#> 7       1        0.1        V2 0.4092737
#> 8       1        0.1        V6 0.3458473
#> 9       1        0.1        V7 0.3802398
#> 10      1        0.1       V10 0.4064830
#> 11      2        0.2        V7 0.3797007
#> 12      2        0.2        V8 0.3586391
#> 13      2        0.2        V9 0.3367853
#> 14      2        0.2       V10 0.4057308
#> 15      2        0.2        V1 0.3923128
#> 16      2        0.2        V5 0.4108237
#> 17      2        0.2        V2 0.4086276
#> 18      2        0.2        V3 0.3970195
#> 19      2        0.2        V4 0.3949241
#> 20      2        0.2        V6 0.3453337
#> 21      3        0.3        V4 0.3942846
#> 22      3        0.3        V5 0.4101423
#> 23      3        0.3        V3 0.3964713
#> 24      3        0.3        V7 0.3791624
#> 25      3        0.3        V8 0.3581170
#> 26      3        0.3        V9 0.3362473
#> 27      3        0.3        V6 0.3448209
#> 28      3        0.3       V10 0.4049800
#> 29      3        0.3        V1 0.3917427
#> 30      3        0.3        V2 0.4079825
#> 31      4        0.4        V1 0.3911734
#> 32      4        0.4        V9 0.3357102
#> 33      4        0.4        V3 0.3959239
#> 34      4        0.4        V4 0.3936461
#> 35      4        0.4        V5 0.4094622
#> 36      4        0.4        V2 0.4073385
#> 37      4        0.4        V6 0.3443088
#> 38      4        0.4        V7 0.3786248
#> 39      4        0.4        V8 0.3575956
#> 40      4        0.4       V10 0.4042306
#> 41      5        0.5        V8 0.3570750
#> 42      5        0.5        V9 0.3351739
#> 43      5        0.5       V10 0.4034826
#> 44      5        0.5        V1 0.3906049
#> 45      5        0.5        V5 0.4087831
#> 46      5        0.5        V2 0.4066954
#> 47      5        0.5        V3 0.3953772
#> 48      5        0.5        V4 0.3930087
#> 49      5        0.5        V6 0.3437975
#> 50      5        0.5        V7 0.3780880
#> 51      6        0.6        V4 0.3923723
#> 52      6        0.6        V5 0.4081052
#> 53      6        0.6        V7 0.3775520
#> 54      6        0.6        V8 0.3565551
#> 55      6        0.6        V9 0.3346385
#> 56      6        0.6        V6 0.3432870
#> 57      6        0.6       V10 0.4027359
#> 58      6        0.6        V1 0.3900373
#> 59      6        0.6        V2 0.4060534
#> 60      6        0.6        V3 0.3948312
#> 61      7        0.7        V1 0.3894705
#> 62      7        0.7        V3 0.3942861
#> 63      7        0.7        V4 0.3917369
#> 64      7        0.7        V5 0.4074284
#> 65      7        0.7        V9 0.3341040
#> 66      7        0.7        V6 0.3427772
#> 67      7        0.7        V7 0.3770167
#> 68      7        0.7        V8 0.3560360
#> 69      7        0.7        V2 0.4054124
#> 70      7        0.7       V10 0.4019907
#> 71      8        0.8        V9 0.3335703
#> 72      8        0.8       V10 0.4012468
#> 73      8        0.8        V1 0.3889045
#> 74      8        0.8        V5 0.4067527
#> 75      8        0.8        V2 0.4047724
#> 76      8        0.8        V3 0.3937416
#> 77      8        0.8        V4 0.3911026
#> 78      8        0.8        V8 0.3555176
#> 79      8        0.8        V6 0.3422682
#> 80      8        0.8        V7 0.3764822
#> 81      9        0.9        V5 0.4060781
#> 82      9        0.9        V7 0.3759485
#> 83      9        0.9        V8 0.3550000
#> 84      9        0.9        V9 0.3330374
#> 85      9        0.9        V6 0.3417599
#> 86      9        0.9       V10 0.4005043
#> 87      9        0.9        V1 0.3883393
#> 88      9        0.9        V2 0.4041334
#> 89      9        0.9        V3 0.3931980
#> 90      9        0.9        V4 0.3904693
#> 91     10        1.0        V1 0.3877750
#> 92     10        1.0        V4 0.3898370
#> 93     10        1.0        V5 0.4054047
#> 94     10        1.0        V9 0.3325054
#> 95     10        1.0        V3 0.3926550
#> 96     10        1.0        V7 0.3754155
#> 97     10        1.0        V8 0.3544832
#> 98     10        1.0        V2 0.4034954
#> 99     10        1.0        V6 0.3412524
#> 100    10        1.0       V10 0.3997632

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3928838
#> 2       1        0.1        V5 0.4115061
#> 3       1        0.1        V9 0.3373241
#> 4       1        0.1        V3 0.3975685
#> 5       1        0.1        V4 0.3955646
#> 6       1        0.1        V8 0.3591620
#> 7       1        0.1        V2 0.4092737
#> 8       1        0.1        V6 0.3458473
#> 9       1        0.1        V7 0.3802398
#> 10      1        0.1       V10 0.4064830
#> 11      2        0.2        V7 0.3797007
#> 12      2        0.2        V8 0.3586391
#> 13      2        0.2        V9 0.3367853
#> 14      2        0.2       V10 0.4057308
#> 15      2        0.2        V1 0.3923128
#> 16      2        0.2        V5 0.4108237
#> 17      2        0.2        V2 0.4086276
#> 18      2        0.2        V3 0.3970195
#> 19      2        0.2        V4 0.3949241
#> 20      2        0.2        V6 0.3453337
#> 21      3        0.3        V4 0.3942846
#> 22      3        0.3        V5 0.4101423
#> 23      3        0.3        V3 0.3964713
#> 24      3        0.3        V7 0.3791624
#> 25      3        0.3        V8 0.3581170
#> 26      3        0.3        V9 0.3362473
#> 27      3        0.3        V6 0.3448209
#> 28      3        0.3       V10 0.4049800
#> 29      3        0.3        V1 0.3917427
#> 30      3        0.3        V2 0.4079825
#> 31      4        0.4        V1 0.3911734
#> 32      4        0.4        V9 0.3357102
#> 33      4        0.4        V3 0.3959239
#> 34      4        0.4        V4 0.3936461
#> 35      4        0.4        V5 0.4094622
#> 36      4        0.4        V2 0.4073385
#> 37      4        0.4        V6 0.3443088
#> 38      4        0.4        V7 0.3786248
#> 39      4        0.4        V8 0.3575956
#> 40      4        0.4       V10 0.4042306
#> 41      5        0.5        V8 0.3570750
#> 42      5        0.5        V9 0.3351739
#> 43      5        0.5       V10 0.4034826
#> 44      5        0.5        V1 0.3906049
#> 45      5        0.5        V5 0.4087831
#> 46      5        0.5        V2 0.4066954
#> 47      5        0.5        V3 0.3953772
#> 48      5        0.5        V4 0.3930087
#> 49      5        0.5        V6 0.3437975
#> 50      5        0.5        V7 0.3780880
#> 51      6        0.6        V4 0.3923723
#> 52      6        0.6        V5 0.4081052
#> 53      6        0.6        V7 0.3775520
#> 54      6        0.6        V8 0.3565551
#> 55      6        0.6        V9 0.3346385
#> 56      6        0.6        V6 0.3432870
#> 57      6        0.6       V10 0.4027359
#> 58      6        0.6        V1 0.3900373
#> 59      6        0.6        V2 0.4060534
#> 60      6        0.6        V3 0.3948312
#> 61      7        0.7        V1 0.3894705
#> 62      7        0.7        V3 0.3942861
#> 63      7        0.7        V4 0.3917369
#> 64      7        0.7        V5 0.4074284
#> 65      7        0.7        V9 0.3341040
#> 66      7        0.7        V6 0.3427772
#> 67      7        0.7        V7 0.3770167
#> 68      7        0.7        V8 0.3560360
#> 69      7        0.7        V2 0.4054124
#> 70      7        0.7       V10 0.4019907
#> 71      8        0.8        V9 0.3335703
#> 72      8        0.8       V10 0.4012468
#> 73      8        0.8        V1 0.3889045
#> 74      8        0.8        V5 0.4067527
#> 75      8        0.8        V2 0.4047724
#> 76      8        0.8        V3 0.3937416
#> 77      8        0.8        V4 0.3911026
#> 78      8        0.8        V8 0.3555176
#> 79      8        0.8        V6 0.3422682
#> 80      8        0.8        V7 0.3764822
#> 81      9        0.9        V5 0.4060781
#> 82      9        0.9        V7 0.3759485
#> 83      9        0.9        V8 0.3550000
#> 84      9        0.9        V9 0.3330374
#> 85      9        0.9        V6 0.3417599
#> 86      9        0.9       V10 0.4005043
#> 87      9        0.9        V1 0.3883393
#> 88      9        0.9        V2 0.4041334
#> 89      9        0.9        V3 0.3931980
#> 90      9        0.9        V4 0.3904693
#> 91     10        1.0        V1 0.3877750
#> 92     10        1.0        V4 0.3898370
#> 93     10        1.0        V5 0.4054047
#> 94     10        1.0        V9 0.3325054
#> 95     10        1.0        V3 0.3926550
#> 96     10        1.0        V7 0.3754155
#> 97     10        1.0        V8 0.3544832
#> 98     10        1.0        V2 0.4034954
#> 99     10        1.0        V6 0.3412524
#> 100    10        1.0       V10 0.3997632


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
#> 1       0                0          0 0.8757906 0.05033323 0.7682514 0.9470793
#> 2       0               20         20 0.8757906 0.05033323 0.7682514 0.9470793
#> 3       0               40         40 0.8757906 0.05033323 0.7682514 0.9470793
#> 4       0               60         60 0.8757906 0.05033323 0.7682514 0.9470793
#> 5      20                0         20 0.8617131 0.05169933 0.7529647 0.9378029
#> 6      20               20         40 0.8617131 0.05169933 0.7529647 0.9378029
#> 7      20               40         60 0.8617131 0.05169933 0.7529647 0.9378029
#> 8      20               60         80 0.8617131 0.05169933 0.7529647 0.9378029
#> 9      40                0         40 0.8478591 0.05295522 0.7381065 0.9283857
#> 10     40               20         60 0.8478591 0.05295522 0.7381065 0.9283857
#> 11     40               40         80 0.8478591 0.05295522 0.7381065 0.9283857
#> 12     40               60        100 0.8478591 0.05295522 0.7381065 0.9283857
#> 13     60                0         60 0.8342249 0.05411553 0.7236455 0.9188686
#> 14     60               20         80 0.8342249 0.05411553 0.7236455 0.9188686
#> 15     60               40        100 0.8342249 0.05411553 0.7236455 0.9188686
#> 16     60               60        120 0.8342249 0.05411553 0.7236455 0.9188686
#> 17     80                0         80 0.8208071 0.05519170 0.7095557 0.9092835
#> 18     80               20        100 0.8208071 0.05519170 0.7095557 0.9092835
#> 19     80               40        120 0.8208071 0.05519170 0.7095557 0.9092835
#> 20     80               60        140 0.8208071 0.05519170 0.7095557 0.9092835
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11758045 0.167855197 0.5703711
#> 2  0.30574618 0.11993429 0.122880782 0.5311427
#> 3  0.26001915 0.12058855 0.088693845 0.4947243
#> 4  0.22113100 0.11981670 0.062889122 0.4609673
#> 5  0.25589195 0.10823272 0.099388301 0.4606061
#> 6  0.21762106 0.10744394 0.070936075 0.4293775
#> 7  0.18507391 0.10579429 0.049599947 0.4004850
#> 8  0.15739448 0.10342594 0.033810283 0.3737664
#> 9  0.18213629 0.09672312 0.056234952 0.3734808
#> 10 0.15489621 0.09435935 0.038693507 0.3488003
#> 11 0.13173012 0.09160591 0.025856088 0.3259840
#> 12 0.11202872 0.08853861 0.016662284 0.3048887
#> 13 0.12963921 0.08473356 0.029808001 0.3046631
#> 14 0.11025053 0.08170381 0.019467320 0.2851711
#> 15 0.09376159 0.07854911 0.012193045 0.2671380
#> 16 0.07973872 0.07531021 0.007249129 0.2504460
#> 17 0.09227334 0.07320931 0.014395992 0.2502674
#> 18 0.07847305 0.07002741 0.008724953 0.2348203
#> 19 0.06673672 0.06686853 0.004978590 0.2205024
#> 20 0.05675565 0.06375369 0.002634778 0.2072198
```
