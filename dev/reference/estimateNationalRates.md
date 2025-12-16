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
[`bbouNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/bbouNationalPriors.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`getNationalCoefficients()`](https://landscitech.github.io/caribouMetrics/dev/reference/getNationalCoefficients.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
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
#> 1         0.1 0.3838513 0.01722816 0.3670972 0.4188928
#> 2         0.2 0.3832760 0.01719940 0.3665884 0.4183568
#> 3         0.3 0.3827015 0.01717100 0.3660803 0.4178215
#> 4         0.4 0.3821279 0.01714294 0.3655729 0.4172869
#> 5         0.5 0.3815551 0.01711524 0.3650662 0.4167530
#> 6         0.6 0.3809832 0.01708788 0.3645602 0.4162198
#> 7         0.7 0.3804122 0.01706087 0.3640549 0.4156873
#> 8         0.8 0.3798420 0.01703420 0.3635504 0.4151556
#> 9         0.9 0.3792726 0.01700787 0.3630465 0.4146245
#> 10        1.0 0.3787041 0.01698189 0.3625433 0.4140941

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3816591
#> 2       1        0.1        V5 0.3637750
#> 3       1        0.1        V9 0.4101887
#> 4       1        0.1        V3 0.4004175
#> 5       1        0.1        V4 0.3785400
#> 6       1        0.1        V8 0.3875027
#> 7       1        0.1        V2 0.4214197
#> 8       1        0.1        V6 0.4042018
#> 9       1        0.1        V7 0.4071746
#> 10      1        0.1       V10 0.3944883
#> 11      2        0.2        V7 0.4064747
#> 12      2        0.2        V8 0.3869299
#> 13      2        0.2        V9 0.4095115
#> 14      2        0.2       V10 0.3938556
#> 15      2        0.2        V1 0.3810750
#> 16      2        0.2        V5 0.3632891
#> 17      2        0.2        V2 0.4209247
#> 18      2        0.2        V3 0.3998544
#> 19      2        0.2        V4 0.3779526
#> 20      2        0.2        V6 0.4035041
#> 21      3        0.3        V4 0.3773661
#> 22      3        0.3        V5 0.3628037
#> 23      3        0.3        V3 0.3992920
#> 24      3        0.3        V7 0.4057759
#> 25      3        0.3        V8 0.3863580
#> 26      3        0.3        V9 0.4088355
#> 27      3        0.3        V6 0.4028075
#> 28      3        0.3       V10 0.3932240
#> 29      3        0.3        V1 0.3804918
#> 30      3        0.3        V2 0.4204303
#> 31      4        0.4        V1 0.3799095
#> 32      4        0.4        V9 0.4081606
#> 33      4        0.4        V3 0.3987304
#> 34      4        0.4        V4 0.3767805
#> 35      4        0.4        V5 0.3623190
#> 36      4        0.4        V2 0.4199365
#> 37      4        0.4        V6 0.4021122
#> 38      4        0.4        V7 0.4050783
#> 39      4        0.4        V8 0.3857869
#> 40      4        0.4       V10 0.3925933
#> 41      5        0.5        V8 0.3852166
#> 42      5        0.5        V9 0.4074869
#> 43      5        0.5       V10 0.3919637
#> 44      5        0.5        V1 0.3793281
#> 45      5        0.5        V5 0.3618350
#> 46      5        0.5        V2 0.4194432
#> 47      5        0.5        V3 0.3981696
#> 48      5        0.5        V4 0.3761958
#> 49      5        0.5        V6 0.4014181
#> 50      5        0.5        V7 0.4043820
#> 51      6        0.6        V4 0.3756120
#> 52      6        0.6        V5 0.3613516
#> 53      6        0.6        V7 0.4036868
#> 54      6        0.6        V8 0.3846472
#> 55      6        0.6        V9 0.4068142
#> 56      6        0.6        V6 0.4007251
#> 57      6        0.6       V10 0.3913350
#> 58      6        0.6        V1 0.3787476
#> 59      6        0.6        V2 0.4189505
#> 60      6        0.6        V3 0.3976096
#> 61      7        0.7        V1 0.3781680
#> 62      7        0.7        V3 0.3970503
#> 63      7        0.7        V4 0.3750291
#> 64      7        0.7        V5 0.3608689
#> 65      7        0.7        V9 0.4061426
#> 66      7        0.7        V6 0.4000334
#> 67      7        0.7        V7 0.4029928
#> 68      7        0.7        V8 0.3840786
#> 69      7        0.7        V2 0.4184584
#> 70      7        0.7       V10 0.3907074
#> 71      8        0.8        V9 0.4054722
#> 72      8        0.8       V10 0.3900808
#> 73      8        0.8        V1 0.3775892
#> 74      8        0.8        V5 0.3603868
#> 75      8        0.8        V2 0.4179669
#> 76      8        0.8        V3 0.3964919
#> 77      8        0.8        V4 0.3744472
#> 78      8        0.8        V8 0.3835109
#> 79      8        0.8        V6 0.3993428
#> 80      8        0.8        V7 0.4023001
#> 81      9        0.9        V5 0.3599053
#> 82      9        0.9        V7 0.4016085
#> 83      9        0.9        V8 0.3829440
#> 84      9        0.9        V9 0.4048029
#> 85      9        0.9        V6 0.3986535
#> 86      9        0.9       V10 0.3894551
#> 87      9        0.9        V1 0.3770114
#> 88      9        0.9        V2 0.4174759
#> 89      9        0.9        V3 0.3959342
#> 90      9        0.9        V4 0.3738661
#> 91     10        1.0        V1 0.3764344
#> 92     10        1.0        V4 0.3732859
#> 93     10        1.0        V5 0.3594245
#> 94     10        1.0        V9 0.4041346
#> 95     10        1.0        V3 0.3953773
#> 96     10        1.0        V7 0.4009181
#> 97     10        1.0        V8 0.3823780
#> 98     10        1.0        V2 0.4169855
#> 99     10        1.0        V6 0.3979653
#> 100    10        1.0       V10 0.3888305

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3816591
#> 2       1        0.1        V5 0.3637750
#> 3       1        0.1        V9 0.4101887
#> 4       1        0.1        V3 0.4004175
#> 5       1        0.1        V4 0.3785400
#> 6       1        0.1        V8 0.3875027
#> 7       1        0.1        V2 0.4214197
#> 8       1        0.1        V6 0.4042018
#> 9       1        0.1        V7 0.4071746
#> 10      1        0.1       V10 0.3944883
#> 11      2        0.2        V7 0.4064747
#> 12      2        0.2        V8 0.3869299
#> 13      2        0.2        V9 0.4095115
#> 14      2        0.2       V10 0.3938556
#> 15      2        0.2        V1 0.3810750
#> 16      2        0.2        V5 0.3632891
#> 17      2        0.2        V2 0.4209247
#> 18      2        0.2        V3 0.3998544
#> 19      2        0.2        V4 0.3779526
#> 20      2        0.2        V6 0.4035041
#> 21      3        0.3        V4 0.3773661
#> 22      3        0.3        V5 0.3628037
#> 23      3        0.3        V3 0.3992920
#> 24      3        0.3        V7 0.4057759
#> 25      3        0.3        V8 0.3863580
#> 26      3        0.3        V9 0.4088355
#> 27      3        0.3        V6 0.4028075
#> 28      3        0.3       V10 0.3932240
#> 29      3        0.3        V1 0.3804918
#> 30      3        0.3        V2 0.4204303
#> 31      4        0.4        V1 0.3799095
#> 32      4        0.4        V9 0.4081606
#> 33      4        0.4        V3 0.3987304
#> 34      4        0.4        V4 0.3767805
#> 35      4        0.4        V5 0.3623190
#> 36      4        0.4        V2 0.4199365
#> 37      4        0.4        V6 0.4021122
#> 38      4        0.4        V7 0.4050783
#> 39      4        0.4        V8 0.3857869
#> 40      4        0.4       V10 0.3925933
#> 41      5        0.5        V8 0.3852166
#> 42      5        0.5        V9 0.4074869
#> 43      5        0.5       V10 0.3919637
#> 44      5        0.5        V1 0.3793281
#> 45      5        0.5        V5 0.3618350
#> 46      5        0.5        V2 0.4194432
#> 47      5        0.5        V3 0.3981696
#> 48      5        0.5        V4 0.3761958
#> 49      5        0.5        V6 0.4014181
#> 50      5        0.5        V7 0.4043820
#> 51      6        0.6        V4 0.3756120
#> 52      6        0.6        V5 0.3613516
#> 53      6        0.6        V7 0.4036868
#> 54      6        0.6        V8 0.3846472
#> 55      6        0.6        V9 0.4068142
#> 56      6        0.6        V6 0.4007251
#> 57      6        0.6       V10 0.3913350
#> 58      6        0.6        V1 0.3787476
#> 59      6        0.6        V2 0.4189505
#> 60      6        0.6        V3 0.3976096
#> 61      7        0.7        V1 0.3781680
#> 62      7        0.7        V3 0.3970503
#> 63      7        0.7        V4 0.3750291
#> 64      7        0.7        V5 0.3608689
#> 65      7        0.7        V9 0.4061426
#> 66      7        0.7        V6 0.4000334
#> 67      7        0.7        V7 0.4029928
#> 68      7        0.7        V8 0.3840786
#> 69      7        0.7        V2 0.4184584
#> 70      7        0.7       V10 0.3907074
#> 71      8        0.8        V9 0.4054722
#> 72      8        0.8       V10 0.3900808
#> 73      8        0.8        V1 0.3775892
#> 74      8        0.8        V5 0.3603868
#> 75      8        0.8        V2 0.4179669
#> 76      8        0.8        V3 0.3964919
#> 77      8        0.8        V4 0.3744472
#> 78      8        0.8        V8 0.3835109
#> 79      8        0.8        V6 0.3993428
#> 80      8        0.8        V7 0.4023001
#> 81      9        0.9        V5 0.3599053
#> 82      9        0.9        V7 0.4016085
#> 83      9        0.9        V8 0.3829440
#> 84      9        0.9        V9 0.4048029
#> 85      9        0.9        V6 0.3986535
#> 86      9        0.9       V10 0.3894551
#> 87      9        0.9        V1 0.3770114
#> 88      9        0.9        V2 0.4174759
#> 89      9        0.9        V3 0.3959342
#> 90      9        0.9        V4 0.3738661
#> 91     10        1.0        V1 0.3764344
#> 92     10        1.0        V4 0.3732859
#> 93     10        1.0        V5 0.3594245
#> 94     10        1.0        V9 0.4041346
#> 95     10        1.0        V3 0.3953773
#> 96     10        1.0        V7 0.4009181
#> 97     10        1.0        V8 0.3823780
#> 98     10        1.0        V2 0.4169855
#> 99     10        1.0        V6 0.3979653
#> 100    10        1.0       V10 0.3888305


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
#> 1       0                0          0 0.8757906 0.04360871 0.8026298 0.9379030
#> 2       0               20         20 0.8757906 0.04360871 0.8026298 0.9379030
#> 3       0               40         40 0.8757906 0.04360871 0.8026298 0.9379030
#> 4       0               60         60 0.8757906 0.04360871 0.8026298 0.9379030
#> 5      20                0         20 0.8617131 0.04537344 0.7852669 0.9272629
#> 6      20               20         40 0.8617131 0.04537344 0.7852669 0.9272629
#> 7      20               40         60 0.8617131 0.04537344 0.7852669 0.9272629
#> 8      20               60         80 0.8617131 0.04537344 0.7852669 0.9272629
#> 9      40                0         40 0.8478591 0.04698492 0.7684805 0.9165138
#> 10     40               20         60 0.8478591 0.04698492 0.7684805 0.9165138
#> 11     40               40         80 0.8478591 0.04698492 0.7684805 0.9165138
#> 12     40               60        100 0.8478591 0.04698492 0.7684805 0.9165138
#> 13     60                0         60 0.8342249 0.04846361 0.7522148 0.9056978
#> 14     60               20         80 0.8342249 0.04846361 0.7522148 0.9056978
#> 15     60               40        100 0.8342249 0.04846361 0.7522148 0.9056978
#> 16     60               60        120 0.8342249 0.04846361 0.7522148 0.9056978
#> 17     80                0         80 0.8208071 0.04982562 0.7364259 0.8948477
#> 18     80               20        100 0.8208071 0.04982562 0.7364259 0.8948477
#> 19     80               40        120 0.8208071 0.04982562 0.7364259 0.8948477
#> 20     80               60        140 0.8208071 0.04982562 0.7364259 0.8948477
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11450260 0.185631272 0.5374343
#> 2  0.30574618 0.10830360 0.136999918 0.4639916
#> 3  0.26001915 0.10217391 0.099821229 0.4013743
#> 4  0.22113100 0.09607683 0.071559917 0.3481353
#> 5  0.25589195 0.10569271 0.107634867 0.4445307
#> 6  0.21762106 0.09790492 0.077480674 0.3848186
#> 7  0.18507391 0.09071332 0.054713098 0.3340688
#> 8  0.15739448 0.08399716 0.037729878 0.2909424
#> 9  0.18213629 0.09526496 0.059466636 0.3690363
#> 10 0.15489621 0.08715128 0.041254701 0.3206598
#> 11 0.13173012 0.07982681 0.027834494 0.2795400
#> 12 0.11202872 0.07315174 0.018146436 0.2445310
#> 13 0.12963921 0.08451634 0.030603407 0.3078764
#> 14 0.11025053 0.07671564 0.020125572 0.2686645
#> 15 0.09376159 0.06973970 0.012712007 0.2352561
#> 16 0.07973872 0.06345311 0.007638694 0.2067039
#> 17 0.09227334 0.07417228 0.014211845 0.2582892
#> 18 0.07847305 0.06698594 0.008648551 0.2263997
#> 19 0.06673672 0.06058535 0.004962681 0.1991153
#> 20 0.05675565 0.05484697 0.002647718 0.1756665
```
