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
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
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
#> 1         0.1 0.3838513 0.02109705 0.3603521 0.4147716
#> 2         0.2 0.3832760 0.02105503 0.3597916 0.4141219
#> 3         0.3 0.3827015 0.02101325 0.3592320 0.4134732
#> 4         0.4 0.3821279 0.02097171 0.3586732 0.4128255
#> 5         0.5 0.3815551 0.02093042 0.3581154 0.4121789
#> 6         0.6 0.3809832 0.02088936 0.3575583 0.4115332
#> 7         0.7 0.3804122 0.02084854 0.3570022 0.4108886
#> 8         0.8 0.3798420 0.02080796 0.3564469 0.4102450
#> 9         0.9 0.3792726 0.02076761 0.3558925 0.4096024
#> 10        1.0 0.3787041 0.02072751 0.3553390 0.4089608

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3814906
#> 2       1        0.1        V5 0.3779617
#> 3       1        0.1        V9 0.4147863
#> 4       1        0.1        V3 0.3783502
#> 5       1        0.1        V4 0.3616693
#> 6       1        0.1        V8 0.4107025
#> 7       1        0.1        V2 0.3814234
#> 8       1        0.1        V6 0.3599698
#> 9       1        0.1        V7 0.4147212
#> 10      1        0.1       V10 0.4054525
#> 11      2        0.2        V7 0.4140305
#> 12      2        0.2        V8 0.4101182
#> 13      2        0.2        V9 0.4141484
#> 14      2        0.2       V10 0.4048534
#> 15      2        0.2        V1 0.3810373
#> 16      2        0.2        V5 0.3773203
#> 17      2        0.2        V2 0.3808746
#> 18      2        0.2        V3 0.3777778
#> 19      2        0.2        V4 0.3611914
#> 20      2        0.2        V6 0.3593852
#> 21      3        0.3        V4 0.3607142
#> 22      3        0.3        V5 0.3766801
#> 23      3        0.3        V3 0.3772063
#> 24      3        0.3        V7 0.4133410
#> 25      3        0.3        V8 0.4095347
#> 26      3        0.3        V9 0.4135116
#> 27      3        0.3        V6 0.3588017
#> 28      3        0.3       V10 0.4042551
#> 29      3        0.3        V1 0.3805844
#> 30      3        0.3        V2 0.3803265
#> 31      4        0.4        V1 0.3801322
#> 32      4        0.4        V9 0.4128757
#> 33      4        0.4        V3 0.3766356
#> 34      4        0.4        V4 0.3602377
#> 35      4        0.4        V5 0.3760409
#> 36      4        0.4        V2 0.3797792
#> 37      4        0.4        V6 0.3582191
#> 38      4        0.4        V7 0.4126527
#> 39      4        0.4        V8 0.4089521
#> 40      4        0.4       V10 0.4036578
#> 41      5        0.5        V8 0.4083703
#> 42      5        0.5        V9 0.4122408
#> 43      5        0.5       V10 0.4030613
#> 44      5        0.5        V1 0.3796804
#> 45      5        0.5        V5 0.3754028
#> 46      5        0.5        V2 0.3792327
#> 47      5        0.5        V3 0.3760658
#> 48      5        0.5        V4 0.3597617
#> 49      5        0.5        V6 0.3576374
#> 50      5        0.5        V7 0.4119654
#> 51      6        0.6        V4 0.3592864
#> 52      6        0.6        V5 0.3747658
#> 53      6        0.6        V7 0.4112794
#> 54      6        0.6        V8 0.4077893
#> 55      6        0.6        V9 0.4116069
#> 56      6        0.6        V6 0.3570567
#> 57      6        0.6       V10 0.4024657
#> 58      6        0.6        V1 0.3792292
#> 59      6        0.6        V2 0.3786870
#> 60      6        0.6        V3 0.3754969
#> 61      7        0.7        V1 0.3787785
#> 62      7        0.7        V3 0.3749288
#> 63      7        0.7        V4 0.3588117
#> 64      7        0.7        V5 0.3741299
#> 65      7        0.7        V9 0.4109740
#> 66      7        0.7        V6 0.3564769
#> 67      7        0.7        V7 0.4105945
#> 68      7        0.7        V8 0.4072091
#> 69      7        0.7        V2 0.3781421
#> 70      7        0.7       V10 0.4018710
#> 71      8        0.8        V9 0.4103420
#> 72      8        0.8       V10 0.4012772
#> 73      8        0.8        V1 0.3783284
#> 74      8        0.8        V5 0.3734951
#> 75      8        0.8        V2 0.3775980
#> 76      8        0.8        V3 0.3743616
#> 77      8        0.8        V4 0.3583377
#> 78      8        0.8        V8 0.4066298
#> 79      8        0.8        V6 0.3558980
#> 80      8        0.8        V7 0.4099107
#> 81      9        0.9        V5 0.3728613
#> 82      9        0.9        V7 0.4092280
#> 83      9        0.9        V8 0.4060513
#> 84      9        0.9        V9 0.4097111
#> 85      9        0.9        V6 0.3553201
#> 86      9        0.9       V10 0.4006843
#> 87      9        0.9        V1 0.3778788
#> 88      9        0.9        V2 0.3770546
#> 89      9        0.9        V3 0.3737952
#> 90      9        0.9        V4 0.3578642
#> 91     10        1.0        V1 0.3774297
#> 92     10        1.0        V4 0.3573914
#> 93     10        1.0        V5 0.3722286
#> 94     10        1.0        V9 0.4090810
#> 95     10        1.0        V3 0.3732297
#> 96     10        1.0        V7 0.4085465
#> 97     10        1.0        V8 0.4054736
#> 98     10        1.0        V2 0.3765120
#> 99     10        1.0        V6 0.3547432
#> 100    10        1.0       V10 0.4000922

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3814906
#> 2       1        0.1        V5 0.3779617
#> 3       1        0.1        V9 0.4147863
#> 4       1        0.1        V3 0.3783502
#> 5       1        0.1        V4 0.3616693
#> 6       1        0.1        V8 0.4107025
#> 7       1        0.1        V2 0.3814234
#> 8       1        0.1        V6 0.3599698
#> 9       1        0.1        V7 0.4147212
#> 10      1        0.1       V10 0.4054525
#> 11      2        0.2        V7 0.4140305
#> 12      2        0.2        V8 0.4101182
#> 13      2        0.2        V9 0.4141484
#> 14      2        0.2       V10 0.4048534
#> 15      2        0.2        V1 0.3810373
#> 16      2        0.2        V5 0.3773203
#> 17      2        0.2        V2 0.3808746
#> 18      2        0.2        V3 0.3777778
#> 19      2        0.2        V4 0.3611914
#> 20      2        0.2        V6 0.3593852
#> 21      3        0.3        V4 0.3607142
#> 22      3        0.3        V5 0.3766801
#> 23      3        0.3        V3 0.3772063
#> 24      3        0.3        V7 0.4133410
#> 25      3        0.3        V8 0.4095347
#> 26      3        0.3        V9 0.4135116
#> 27      3        0.3        V6 0.3588017
#> 28      3        0.3       V10 0.4042551
#> 29      3        0.3        V1 0.3805844
#> 30      3        0.3        V2 0.3803265
#> 31      4        0.4        V1 0.3801322
#> 32      4        0.4        V9 0.4128757
#> 33      4        0.4        V3 0.3766356
#> 34      4        0.4        V4 0.3602377
#> 35      4        0.4        V5 0.3760409
#> 36      4        0.4        V2 0.3797792
#> 37      4        0.4        V6 0.3582191
#> 38      4        0.4        V7 0.4126527
#> 39      4        0.4        V8 0.4089521
#> 40      4        0.4       V10 0.4036578
#> 41      5        0.5        V8 0.4083703
#> 42      5        0.5        V9 0.4122408
#> 43      5        0.5       V10 0.4030613
#> 44      5        0.5        V1 0.3796804
#> 45      5        0.5        V5 0.3754028
#> 46      5        0.5        V2 0.3792327
#> 47      5        0.5        V3 0.3760658
#> 48      5        0.5        V4 0.3597617
#> 49      5        0.5        V6 0.3576374
#> 50      5        0.5        V7 0.4119654
#> 51      6        0.6        V4 0.3592864
#> 52      6        0.6        V5 0.3747658
#> 53      6        0.6        V7 0.4112794
#> 54      6        0.6        V8 0.4077893
#> 55      6        0.6        V9 0.4116069
#> 56      6        0.6        V6 0.3570567
#> 57      6        0.6       V10 0.4024657
#> 58      6        0.6        V1 0.3792292
#> 59      6        0.6        V2 0.3786870
#> 60      6        0.6        V3 0.3754969
#> 61      7        0.7        V1 0.3787785
#> 62      7        0.7        V3 0.3749288
#> 63      7        0.7        V4 0.3588117
#> 64      7        0.7        V5 0.3741299
#> 65      7        0.7        V9 0.4109740
#> 66      7        0.7        V6 0.3564769
#> 67      7        0.7        V7 0.4105945
#> 68      7        0.7        V8 0.4072091
#> 69      7        0.7        V2 0.3781421
#> 70      7        0.7       V10 0.4018710
#> 71      8        0.8        V9 0.4103420
#> 72      8        0.8       V10 0.4012772
#> 73      8        0.8        V1 0.3783284
#> 74      8        0.8        V5 0.3734951
#> 75      8        0.8        V2 0.3775980
#> 76      8        0.8        V3 0.3743616
#> 77      8        0.8        V4 0.3583377
#> 78      8        0.8        V8 0.4066298
#> 79      8        0.8        V6 0.3558980
#> 80      8        0.8        V7 0.4099107
#> 81      9        0.9        V5 0.3728613
#> 82      9        0.9        V7 0.4092280
#> 83      9        0.9        V8 0.4060513
#> 84      9        0.9        V9 0.4097111
#> 85      9        0.9        V6 0.3553201
#> 86      9        0.9       V10 0.4006843
#> 87      9        0.9        V1 0.3778788
#> 88      9        0.9        V2 0.3770546
#> 89      9        0.9        V3 0.3737952
#> 90      9        0.9        V4 0.3578642
#> 91     10        1.0        V1 0.3774297
#> 92     10        1.0        V4 0.3573914
#> 93     10        1.0        V5 0.3722286
#> 94     10        1.0        V9 0.4090810
#> 95     10        1.0        V3 0.3732297
#> 96     10        1.0        V7 0.4085465
#> 97     10        1.0        V8 0.4054736
#> 98     10        1.0        V2 0.3765120
#> 99     10        1.0        V6 0.3547432
#> 100    10        1.0       V10 0.4000922


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
#> 1       0                0          0 0.8757906 0.04501074 0.7856236 0.9370466
#> 2       0               20         20 0.8757906 0.04501074 0.7856236 0.9370466
#> 3       0               40         40 0.8757906 0.04501074 0.7856236 0.9370466
#> 4       0               60         60 0.8757906 0.04501074 0.7856236 0.9370466
#> 5      20                0         20 0.8617131 0.04753264 0.7673869 0.9286851
#> 6      20               20         40 0.8617131 0.04753264 0.7673869 0.9286851
#> 7      20               40         60 0.8617131 0.04753264 0.7673869 0.9286851
#> 8      20               60         80 0.8617131 0.04753264 0.7673869 0.9286851
#> 9      40                0         40 0.8478591 0.04989644 0.7497816 0.9202505
#> 10     40               20         60 0.8478591 0.04989644 0.7497816 0.9202505
#> 11     40               40         80 0.8478591 0.04989644 0.7497816 0.9202505
#> 12     40               60        100 0.8478591 0.04989644 0.7497816 0.9202505
#> 13     60                0         60 0.8342249 0.05211881 0.7327475 0.9117649
#> 14     60               20         80 0.8342249 0.05211881 0.7327475 0.9117649
#> 15     60               40        100 0.8342249 0.05211881 0.7327475 0.9117649
#> 16     60               60        120 0.8342249 0.05211881 0.7327475 0.9117649
#> 17     80                0         80 0.8208071 0.05421316 0.7162369 0.9032462
#> 18     80               20        100 0.8208071 0.05421316 0.7162369 0.9032462
#> 19     80               40        120 0.8208071 0.05421316 0.7162369 0.9032462
#> 20     80               60        140 0.8208071 0.05421316 0.7162369 0.9032462
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.10714007 0.159313967 0.5161944
#> 2  0.30574618 0.10607971 0.120251680 0.4577656
#> 3  0.26001915 0.10484895 0.089720623 0.4247636
#> 4  0.22113100 0.10312676 0.065994120 0.3999485
#> 5  0.25589195 0.09827849 0.089017899 0.4269194
#> 6  0.21762106 0.09524890 0.065450194 0.3795023
#> 7  0.18507391 0.09232702 0.047291556 0.3380065
#> 8  0.15739448 0.08930244 0.033457919 0.3079403
#> 9  0.18213629 0.08884926 0.046877619 0.3545208
#> 10 0.15489621 0.08477429 0.033144867 0.3161475
#> 11 0.13173012 0.08097500 0.022841277 0.2825560
#> 12 0.11202872 0.07730331 0.015258092 0.2531164
#> 13 0.12963921 0.07924352 0.022610376 0.2959286
#> 14 0.11025053 0.07478000 0.015090239 0.2648415
#> 15 0.09376159 0.07067018 0.009695070 0.2375701
#> 16 0.07973872 0.06680557 0.005946388 0.2135946
#> 17 0.09227334 0.06993747 0.009577591 0.2484355
#> 18 0.07847305 0.06549088 0.005866432 0.2231539
#> 19 0.06673672 0.06141875 0.003394446 0.2008937
#> 20 0.05675565 0.05763826 0.001832310 0.1812334
```
