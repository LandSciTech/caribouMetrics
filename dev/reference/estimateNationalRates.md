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
#> 1         0.1 0.3838513 0.02665852 0.3458537 0.4211757
#> 2         0.2 0.3832760 0.02659709 0.3454144 0.4205164
#> 3         0.3 0.3827015 0.02653583 0.3449757 0.4198581
#> 4         0.4 0.3821279 0.02647474 0.3445375 0.4192009
#> 5         0.5 0.3815551 0.02641382 0.3441000 0.4185447
#> 6         0.6 0.3809832 0.02635307 0.3436629 0.4178895
#> 7         0.7 0.3804122 0.02629249 0.3432265 0.4172353
#> 8         0.8 0.3798420 0.02623207 0.3427906 0.4165822
#> 9         0.9 0.3792726 0.02617182 0.3423553 0.4159301
#> 10        1.0 0.3787041 0.02611174 0.3419205 0.4152790

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3587570
#> 2       1        0.1        V5 0.4253353
#> 3       1        0.1        V9 0.3489233
#> 4       1        0.1        V3 0.3999275
#> 5       1        0.1        V4 0.3953103
#> 6       1        0.1        V8 0.3696357
#> 7       1        0.1        V2 0.3449625
#> 8       1        0.1        V6 0.3970299
#> 9       1        0.1        V7 0.4068481
#> 10      1        0.1       V10 0.3883201
#> 11      2        0.2        V7 0.4062314
#> 12      2        0.2        V8 0.3690957
#> 13      2        0.2        V9 0.3483805
#> 14      2        0.2       V10 0.3877326
#> 15      2        0.2        V1 0.3581930
#> 16      2        0.2        V5 0.4246637
#> 17      2        0.2        V2 0.3445533
#> 18      2        0.2        V3 0.3992929
#> 19      2        0.2        V4 0.3946977
#> 20      2        0.2        V6 0.3964561
#> 21      3        0.3        V4 0.3940861
#> 22      3        0.3        V5 0.4239931
#> 23      3        0.3        V3 0.3986593
#> 24      3        0.3        V7 0.4056155
#> 25      3        0.3        V8 0.3685564
#> 26      3        0.3        V9 0.3478386
#> 27      3        0.3        V6 0.3958832
#> 28      3        0.3       V10 0.3871460
#> 29      3        0.3        V1 0.3576299
#> 30      3        0.3        V2 0.3441445
#> 31      4        0.4        V1 0.3570676
#> 32      4        0.4        V9 0.3472975
#> 33      4        0.4        V3 0.3980267
#> 34      4        0.4        V4 0.3934754
#> 35      4        0.4        V5 0.4233235
#> 36      4        0.4        V2 0.3437363
#> 37      4        0.4        V6 0.3953111
#> 38      4        0.4        V7 0.4050006
#> 39      4        0.4        V8 0.3680179
#> 40      4        0.4       V10 0.3865602
#> 41      5        0.5        V8 0.3674803
#> 42      5        0.5        V9 0.3467573
#> 43      5        0.5       V10 0.3859754
#> 44      5        0.5        V1 0.3565062
#> 45      5        0.5        V5 0.4226551
#> 46      5        0.5        V2 0.3433285
#> 47      5        0.5        V3 0.3973951
#> 48      5        0.5        V4 0.3928657
#> 49      5        0.5        V6 0.3947398
#> 50      5        0.5        V7 0.4043867
#> 51      6        0.6        V4 0.3922569
#> 52      6        0.6        V5 0.4219876
#> 53      6        0.6        V7 0.4037736
#> 54      6        0.6        V8 0.3669434
#> 55      6        0.6        V9 0.3462179
#> 56      6        0.6        V6 0.3941694
#> 57      6        0.6       V10 0.3853914
#> 58      6        0.6        V1 0.3559457
#> 59      6        0.6        V2 0.3429212
#> 60      6        0.6        V3 0.3967644
#> 61      7        0.7        V1 0.3553861
#> 62      7        0.7        V3 0.3961348
#> 63      7        0.7        V4 0.3916490
#> 64      7        0.7        V5 0.4213213
#> 65      7        0.7        V9 0.3456793
#> 66      7        0.7        V6 0.3935998
#> 67      7        0.7        V7 0.4031615
#> 68      7        0.7        V8 0.3664073
#> 69      7        0.7        V2 0.3425144
#> 70      7        0.7       V10 0.3848083
#> 71      8        0.8        V9 0.3451416
#> 72      8        0.8       V10 0.3842261
#> 73      8        0.8        V1 0.3548274
#> 74      8        0.8        V5 0.4206560
#> 75      8        0.8        V2 0.3421080
#> 76      8        0.8        V3 0.3955062
#> 77      8        0.8        V4 0.3910421
#> 78      8        0.8        V8 0.3658719
#> 79      8        0.8        V6 0.3930310
#> 80      8        0.8        V7 0.4025503
#> 81      9        0.9        V5 0.4199917
#> 82      9        0.9        V7 0.4019401
#> 83      9        0.9        V8 0.3653374
#> 84      9        0.9        V9 0.3446048
#> 85      9        0.9        V6 0.3924630
#> 86      9        0.9       V10 0.3836448
#> 87      9        0.9        V1 0.3542695
#> 88      9        0.9        V2 0.3417022
#> 89      9        0.9        V3 0.3948786
#> 90      9        0.9        V4 0.3904361
#> 91     10        1.0        V1 0.3537125
#> 92     10        1.0        V4 0.3898311
#> 93     10        1.0        V5 0.4193285
#> 94     10        1.0        V9 0.3440687
#> 95     10        1.0        V3 0.3942520
#> 96     10        1.0        V7 0.4013308
#> 97     10        1.0        V8 0.3648036
#> 98     10        1.0        V2 0.3412968
#> 99     10        1.0        V6 0.3918958
#> 100    10        1.0       V10 0.3830643

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3587570
#> 2       1        0.1        V5 0.4253353
#> 3       1        0.1        V9 0.3489233
#> 4       1        0.1        V3 0.3999275
#> 5       1        0.1        V4 0.3953103
#> 6       1        0.1        V8 0.3696357
#> 7       1        0.1        V2 0.3449625
#> 8       1        0.1        V6 0.3970299
#> 9       1        0.1        V7 0.4068481
#> 10      1        0.1       V10 0.3883201
#> 11      2        0.2        V7 0.4062314
#> 12      2        0.2        V8 0.3690957
#> 13      2        0.2        V9 0.3483805
#> 14      2        0.2       V10 0.3877326
#> 15      2        0.2        V1 0.3581930
#> 16      2        0.2        V5 0.4246637
#> 17      2        0.2        V2 0.3445533
#> 18      2        0.2        V3 0.3992929
#> 19      2        0.2        V4 0.3946977
#> 20      2        0.2        V6 0.3964561
#> 21      3        0.3        V4 0.3940861
#> 22      3        0.3        V5 0.4239931
#> 23      3        0.3        V3 0.3986593
#> 24      3        0.3        V7 0.4056155
#> 25      3        0.3        V8 0.3685564
#> 26      3        0.3        V9 0.3478386
#> 27      3        0.3        V6 0.3958832
#> 28      3        0.3       V10 0.3871460
#> 29      3        0.3        V1 0.3576299
#> 30      3        0.3        V2 0.3441445
#> 31      4        0.4        V1 0.3570676
#> 32      4        0.4        V9 0.3472975
#> 33      4        0.4        V3 0.3980267
#> 34      4        0.4        V4 0.3934754
#> 35      4        0.4        V5 0.4233235
#> 36      4        0.4        V2 0.3437363
#> 37      4        0.4        V6 0.3953111
#> 38      4        0.4        V7 0.4050006
#> 39      4        0.4        V8 0.3680179
#> 40      4        0.4       V10 0.3865602
#> 41      5        0.5        V8 0.3674803
#> 42      5        0.5        V9 0.3467573
#> 43      5        0.5       V10 0.3859754
#> 44      5        0.5        V1 0.3565062
#> 45      5        0.5        V5 0.4226551
#> 46      5        0.5        V2 0.3433285
#> 47      5        0.5        V3 0.3973951
#> 48      5        0.5        V4 0.3928657
#> 49      5        0.5        V6 0.3947398
#> 50      5        0.5        V7 0.4043867
#> 51      6        0.6        V4 0.3922569
#> 52      6        0.6        V5 0.4219876
#> 53      6        0.6        V7 0.4037736
#> 54      6        0.6        V8 0.3669434
#> 55      6        0.6        V9 0.3462179
#> 56      6        0.6        V6 0.3941694
#> 57      6        0.6       V10 0.3853914
#> 58      6        0.6        V1 0.3559457
#> 59      6        0.6        V2 0.3429212
#> 60      6        0.6        V3 0.3967644
#> 61      7        0.7        V1 0.3553861
#> 62      7        0.7        V3 0.3961348
#> 63      7        0.7        V4 0.3916490
#> 64      7        0.7        V5 0.4213213
#> 65      7        0.7        V9 0.3456793
#> 66      7        0.7        V6 0.3935998
#> 67      7        0.7        V7 0.4031615
#> 68      7        0.7        V8 0.3664073
#> 69      7        0.7        V2 0.3425144
#> 70      7        0.7       V10 0.3848083
#> 71      8        0.8        V9 0.3451416
#> 72      8        0.8       V10 0.3842261
#> 73      8        0.8        V1 0.3548274
#> 74      8        0.8        V5 0.4206560
#> 75      8        0.8        V2 0.3421080
#> 76      8        0.8        V3 0.3955062
#> 77      8        0.8        V4 0.3910421
#> 78      8        0.8        V8 0.3658719
#> 79      8        0.8        V6 0.3930310
#> 80      8        0.8        V7 0.4025503
#> 81      9        0.9        V5 0.4199917
#> 82      9        0.9        V7 0.4019401
#> 83      9        0.9        V8 0.3653374
#> 84      9        0.9        V9 0.3446048
#> 85      9        0.9        V6 0.3924630
#> 86      9        0.9       V10 0.3836448
#> 87      9        0.9        V1 0.3542695
#> 88      9        0.9        V2 0.3417022
#> 89      9        0.9        V3 0.3948786
#> 90      9        0.9        V4 0.3904361
#> 91     10        1.0        V1 0.3537125
#> 92     10        1.0        V4 0.3898311
#> 93     10        1.0        V5 0.4193285
#> 94     10        1.0        V9 0.3440687
#> 95     10        1.0        V3 0.3942520
#> 96     10        1.0        V7 0.4013308
#> 97     10        1.0        V8 0.3648036
#> 98     10        1.0        V2 0.3412968
#> 99     10        1.0        V6 0.3918958
#> 100    10        1.0       V10 0.3830643


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
#> 1       0                0          0 0.8757906 0.04605283 0.7872648 0.9482641
#> 2       0               20         20 0.8757906 0.04605283 0.7872648 0.9482641
#> 3       0               40         40 0.8757906 0.04605283 0.7872648 0.9482641
#> 4       0               60         60 0.8757906 0.04605283 0.7872648 0.9482641
#> 5      20                0         20 0.8617131 0.04684459 0.7745061 0.9386587
#> 6      20               20         40 0.8617131 0.04684459 0.7745061 0.9386587
#> 7      20               40         60 0.8617131 0.04684459 0.7745061 0.9386587
#> 8      20               60         80 0.8617131 0.04684459 0.7745061 0.9386587
#> 9      40                0         40 0.8478591 0.04760881 0.7620828 0.9289056
#> 10     40               20         60 0.8478591 0.04760881 0.7620828 0.9289056
#> 11     40               40         80 0.8478591 0.04760881 0.7620828 0.9289056
#> 12     40               60        100 0.8478591 0.04760881 0.7620828 0.9289056
#> 13     60                0         60 0.8342249 0.04834979 0.7499681 0.9190500
#> 14     60               20         80 0.8342249 0.04834979 0.7499681 0.9190500
#> 15     60               40        100 0.8342249 0.04834979 0.7499681 0.9190500
#> 16     60               60        120 0.8342249 0.04834979 0.7499681 0.9190500
#> 17     80                0         80 0.8208071 0.04907052 0.7381400 0.9091268
#> 18     80               20        100 0.8208071 0.04907052 0.7381400 0.9091268
#> 19     80               40        120 0.8208071 0.04907052 0.7381400 0.9091268
#> 20     80               60        140 0.8208071 0.04907052 0.7381400 0.9091268
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11017164 0.185407054 0.5720626
#> 2  0.30574618 0.10738287 0.164024862 0.5258196
#> 3  0.26001915 0.10525977 0.144854324 0.4834949
#> 4  0.22113100 0.10306470 0.112947288 0.4448363
#> 5  0.25589195 0.10396763 0.101661554 0.4598584
#> 6  0.21762106 0.10007328 0.089024959 0.4232717
#> 7  0.18507391 0.09667064 0.077742378 0.3899164
#> 8  0.15739448 0.09334906 0.067683611 0.3595232
#> 9  0.18213629 0.09441768 0.052585139 0.3713273
#> 10 0.15489621 0.09022373 0.045330029 0.3425876
#> 11 0.13173012 0.08641174 0.038911673 0.3164028
#> 12 0.11202872 0.08274243 0.033249032 0.2925384
#> 13 0.12963921 0.08368143 0.024888822 0.3018086
#> 14 0.11025053 0.07960303 0.020950442 0.2792317
#> 15 0.09376159 0.07582689 0.017522808 0.2586345
#> 16 0.07973872 0.07220937 0.014553761 0.2398262
#> 17 0.09227334 0.07299197 0.010294546 0.2471374
#> 18 0.07847305 0.06922011 0.008356007 0.2293182
#> 19 0.06673672 0.06568872 0.006715942 0.2130152
#> 20 0.05675565 0.06231030 0.005339701 0.1980777
```
