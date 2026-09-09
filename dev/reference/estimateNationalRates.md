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
#> 1         0.1 0.3838513 0.03101878 0.3386941 0.4250180
#> 2         0.2 0.3832760 0.03098793 0.3381806 0.4243993
#> 3         0.3 0.3827015 0.03095719 0.3376679 0.4237815
#> 4         0.4 0.3821279 0.03092657 0.3371559 0.4231646
#> 5         0.5 0.3815551 0.03089606 0.3366448 0.4225486
#> 6         0.6 0.3809832 0.03086567 0.3361344 0.4219335
#> 7         0.7 0.3804122 0.03083539 0.3356248 0.4213193
#> 8         0.8 0.3798420 0.03080522 0.3351159 0.4207060
#> 9         0.9 0.3792726 0.03077517 0.3346079 0.4200936
#> 10        1.0 0.3787041 0.03074523 0.3341006 0.4194820

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3964821
#> 2       1        0.1        V5 0.4005309
#> 3       1        0.1        V9 0.3535829
#> 4       1        0.1        V3 0.4062681
#> 5       1        0.1        V4 0.4081948
#> 6       1        0.1        V8 0.3521568
#> 7       1        0.1        V2 0.3347856
#> 8       1        0.1        V6 0.4299021
#> 9       1        0.1        V7 0.3618433
#> 10      1        0.1       V10 0.3977919
#> 11      2        0.2        V7 0.3613235
#> 12      2        0.2        V8 0.3516683
#> 13      2        0.2        V9 0.3530737
#> 14      2        0.2       V10 0.3973803
#> 15      2        0.2        V1 0.3959154
#> 16      2        0.2        V5 0.3999890
#> 17      2        0.2        V2 0.3342648
#> 18      2        0.2        V3 0.4056951
#> 19      2        0.2        V4 0.4075943
#> 20      2        0.2        V6 0.4292781
#> 21      3        0.3        V4 0.4069947
#> 22      3        0.3        V5 0.3994478
#> 23      3        0.3        V3 0.4051230
#> 24      3        0.3        V7 0.3608044
#> 25      3        0.3        V8 0.3511804
#> 26      3        0.3        V9 0.3525653
#> 27      3        0.3        V6 0.4286550
#> 28      3        0.3       V10 0.3969692
#> 29      3        0.3        V1 0.3953496
#> 30      3        0.3        V2 0.3337449
#> 31      4        0.4        V1 0.3947845
#> 32      4        0.4        V9 0.3520577
#> 33      4        0.4        V3 0.4045516
#> 34      4        0.4        V4 0.4063960
#> 35      4        0.4        V5 0.3989073
#> 36      4        0.4        V2 0.3332258
#> 37      4        0.4        V6 0.4280329
#> 38      4        0.4        V7 0.3602861
#> 39      4        0.4        V8 0.3506933
#> 40      4        0.4       V10 0.3965585
#> 41      5        0.5        V8 0.3502068
#> 42      5        0.5        V9 0.3515507
#> 43      5        0.5       V10 0.3961482
#> 44      5        0.5        V1 0.3942202
#> 45      5        0.5        V5 0.3983676
#> 46      5        0.5        V2 0.3327074
#> 47      5        0.5        V3 0.4039810
#> 48      5        0.5        V4 0.4057981
#> 49      5        0.5        V6 0.4274116
#> 50      5        0.5        V7 0.3597685
#> 51      6        0.6        V4 0.4052012
#> 52      6        0.6        V5 0.3978286
#> 53      6        0.6        V7 0.3592516
#> 54      6        0.6        V8 0.3497209
#> 55      6        0.6        V9 0.3510445
#> 56      6        0.6        V6 0.4267913
#> 57      6        0.6       V10 0.3957384
#> 58      6        0.6        V1 0.3936568
#> 59      6        0.6        V2 0.3321899
#> 60      6        0.6        V3 0.4034113
#> 61      7        0.7        V1 0.3930941
#> 62      7        0.7        V3 0.4028423
#> 63      7        0.7        V4 0.4046051
#> 64      7        0.7        V5 0.3972903
#> 65      7        0.7        V9 0.3505390
#> 66      7        0.7        V6 0.4261718
#> 67      7        0.7        V7 0.3587355
#> 68      7        0.7        V8 0.3492358
#> 69      7        0.7        V2 0.3316732
#> 70      7        0.7       V10 0.3953289
#> 71      8        0.8        V9 0.3500343
#> 72      8        0.8       V10 0.3949199
#> 73      8        0.8        V1 0.3925323
#> 74      8        0.8        V5 0.3967528
#> 75      8        0.8        V2 0.3311573
#> 76      8        0.8        V3 0.4022741
#> 77      8        0.8        V4 0.4040099
#> 78      8        0.8        V8 0.3487513
#> 79      8        0.8        V6 0.4255532
#> 80      8        0.8        V7 0.3582201
#> 81      9        0.9        V5 0.3962160
#> 82      9        0.9        V7 0.3577055
#> 83      9        0.9        V8 0.3482675
#> 84      9        0.9        V9 0.3495303
#> 85      9        0.9        V6 0.4249356
#> 86      9        0.9       V10 0.3945114
#> 87      9        0.9        V1 0.3919713
#> 88      9        0.9        V2 0.3306422
#> 89      9        0.9        V3 0.4017068
#> 90      9        0.9        V4 0.4034155
#> 91     10        1.0        V1 0.3914110
#> 92     10        1.0        V4 0.4028221
#> 93     10        1.0        V5 0.3956799
#> 94     10        1.0        V9 0.3490270
#> 95     10        1.0        V3 0.4011402
#> 96     10        1.0        V7 0.3571916
#> 97     10        1.0        V8 0.3477844
#> 98     10        1.0        V2 0.3301278
#> 99     10        1.0        V6 0.4243188
#> 100    10        1.0       V10 0.3941032

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3964821
#> 2       1        0.1        V5 0.4005309
#> 3       1        0.1        V9 0.3535829
#> 4       1        0.1        V3 0.4062681
#> 5       1        0.1        V4 0.4081948
#> 6       1        0.1        V8 0.3521568
#> 7       1        0.1        V2 0.3347856
#> 8       1        0.1        V6 0.4299021
#> 9       1        0.1        V7 0.3618433
#> 10      1        0.1       V10 0.3977919
#> 11      2        0.2        V7 0.3613235
#> 12      2        0.2        V8 0.3516683
#> 13      2        0.2        V9 0.3530737
#> 14      2        0.2       V10 0.3973803
#> 15      2        0.2        V1 0.3959154
#> 16      2        0.2        V5 0.3999890
#> 17      2        0.2        V2 0.3342648
#> 18      2        0.2        V3 0.4056951
#> 19      2        0.2        V4 0.4075943
#> 20      2        0.2        V6 0.4292781
#> 21      3        0.3        V4 0.4069947
#> 22      3        0.3        V5 0.3994478
#> 23      3        0.3        V3 0.4051230
#> 24      3        0.3        V7 0.3608044
#> 25      3        0.3        V8 0.3511804
#> 26      3        0.3        V9 0.3525653
#> 27      3        0.3        V6 0.4286550
#> 28      3        0.3       V10 0.3969692
#> 29      3        0.3        V1 0.3953496
#> 30      3        0.3        V2 0.3337449
#> 31      4        0.4        V1 0.3947845
#> 32      4        0.4        V9 0.3520577
#> 33      4        0.4        V3 0.4045516
#> 34      4        0.4        V4 0.4063960
#> 35      4        0.4        V5 0.3989073
#> 36      4        0.4        V2 0.3332258
#> 37      4        0.4        V6 0.4280329
#> 38      4        0.4        V7 0.3602861
#> 39      4        0.4        V8 0.3506933
#> 40      4        0.4       V10 0.3965585
#> 41      5        0.5        V8 0.3502068
#> 42      5        0.5        V9 0.3515507
#> 43      5        0.5       V10 0.3961482
#> 44      5        0.5        V1 0.3942202
#> 45      5        0.5        V5 0.3983676
#> 46      5        0.5        V2 0.3327074
#> 47      5        0.5        V3 0.4039810
#> 48      5        0.5        V4 0.4057981
#> 49      5        0.5        V6 0.4274116
#> 50      5        0.5        V7 0.3597685
#> 51      6        0.6        V4 0.4052012
#> 52      6        0.6        V5 0.3978286
#> 53      6        0.6        V7 0.3592516
#> 54      6        0.6        V8 0.3497209
#> 55      6        0.6        V9 0.3510445
#> 56      6        0.6        V6 0.4267913
#> 57      6        0.6       V10 0.3957384
#> 58      6        0.6        V1 0.3936568
#> 59      6        0.6        V2 0.3321899
#> 60      6        0.6        V3 0.4034113
#> 61      7        0.7        V1 0.3930941
#> 62      7        0.7        V3 0.4028423
#> 63      7        0.7        V4 0.4046051
#> 64      7        0.7        V5 0.3972903
#> 65      7        0.7        V9 0.3505390
#> 66      7        0.7        V6 0.4261718
#> 67      7        0.7        V7 0.3587355
#> 68      7        0.7        V8 0.3492358
#> 69      7        0.7        V2 0.3316732
#> 70      7        0.7       V10 0.3953289
#> 71      8        0.8        V9 0.3500343
#> 72      8        0.8       V10 0.3949199
#> 73      8        0.8        V1 0.3925323
#> 74      8        0.8        V5 0.3967528
#> 75      8        0.8        V2 0.3311573
#> 76      8        0.8        V3 0.4022741
#> 77      8        0.8        V4 0.4040099
#> 78      8        0.8        V8 0.3487513
#> 79      8        0.8        V6 0.4255532
#> 80      8        0.8        V7 0.3582201
#> 81      9        0.9        V5 0.3962160
#> 82      9        0.9        V7 0.3577055
#> 83      9        0.9        V8 0.3482675
#> 84      9        0.9        V9 0.3495303
#> 85      9        0.9        V6 0.4249356
#> 86      9        0.9       V10 0.3945114
#> 87      9        0.9        V1 0.3919713
#> 88      9        0.9        V2 0.3306422
#> 89      9        0.9        V3 0.4017068
#> 90      9        0.9        V4 0.4034155
#> 91     10        1.0        V1 0.3914110
#> 92     10        1.0        V4 0.4028221
#> 93     10        1.0        V5 0.3956799
#> 94     10        1.0        V9 0.3490270
#> 95     10        1.0        V3 0.4011402
#> 96     10        1.0        V7 0.3571916
#> 97     10        1.0        V8 0.3477844
#> 98     10        1.0        V2 0.3301278
#> 99     10        1.0        V6 0.4243188
#> 100    10        1.0       V10 0.3941032


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
#> 1       0                0          0 0.8757906 0.04759047 0.8029766 0.9539037
#> 2       0               20         20 0.8757906 0.04759047 0.8029766 0.9539037
#> 3       0               40         40 0.8757906 0.04759047 0.8029766 0.9539037
#> 4       0               60         60 0.8757906 0.04759047 0.8029766 0.9539037
#> 5      20                0         20 0.8617131 0.04945080 0.7856744 0.9455893
#> 6      20               20         40 0.8617131 0.04945080 0.7856744 0.9455893
#> 7      20               40         60 0.8617131 0.04945080 0.7856744 0.9455893
#> 8      20               60         80 0.8617131 0.04945080 0.7856744 0.9455893
#> 9      40                0         40 0.8478591 0.05126806 0.7689456 0.9371094
#> 10     40               20         60 0.8478591 0.05126806 0.7689456 0.9371094
#> 11     40               40         80 0.8478591 0.05126806 0.7689456 0.9371094
#> 12     40               60        100 0.8478591 0.05126806 0.7689456 0.9371094
#> 13     60                0         60 0.8342249 0.05304108 0.7527346 0.9285045
#> 14     60               20         80 0.8342249 0.05304108 0.7527346 0.9285045
#> 15     60               40        100 0.8342249 0.05304108 0.7527346 0.9285045
#> 16     60               60        120 0.8342249 0.05304108 0.7527346 0.9285045
#> 17     80                0         80 0.8208071 0.05476866 0.7369978 0.9198068
#> 18     80               20        100 0.8208071 0.05476866 0.7369978 0.9198068
#> 19     80               40        120 0.8208071 0.05476866 0.7369978 0.9198068
#> 20     80               60        140 0.8208071 0.05476866 0.7369978 0.9198068
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11147100 0.203121684 0.5902968
#> 2  0.30574618 0.11322052 0.165435563 0.5461081
#> 3  0.26001915 0.11343068 0.134079609 0.5053322
#> 4  0.22113100 0.11230186 0.108029655 0.4677852
#> 5  0.25589195 0.09632765 0.121678082 0.4511566
#> 6  0.21762106 0.09602223 0.097744453 0.4179834
#> 7  0.18507391 0.09480386 0.077939518 0.3875201
#> 8  0.15739448 0.09280234 0.061615686 0.3595572
#> 9  0.18213629 0.08234899 0.070149721 0.3471920
#> 10 0.15489621 0.08096952 0.055217976 0.3225447
#> 11 0.13173012 0.07905478 0.043005098 0.2999196
#> 12 0.11202872 0.07669022 0.033084165 0.2791432
#> 13 0.12963921 0.06988100 0.038250920 0.2699491
#> 14 0.11025053 0.06797768 0.029245181 0.2516022
#> 15 0.09376159 0.06577228 0.022021923 0.2347266
#> 16 0.07973872 0.06332157 0.016292216 0.2191897
#> 17 0.09227334 0.05889805 0.019257555 0.2122986
#> 18 0.07847305 0.05679074 0.014120483 0.1985124
#> 19 0.06673672 0.05452687 0.010128105 0.1857820
#> 20 0.05675565 0.05214510 0.007080158 0.1740091
```
