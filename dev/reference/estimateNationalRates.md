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
#> 1         0.1 0.3838513 0.02001853 0.3371173 0.3926348
#> 2         0.2 0.3832760 0.02001762 0.3365501 0.3920938
#> 3         0.3 0.3827015 0.02001693 0.3359839 0.3915535
#> 4         0.4 0.3821279 0.02001646 0.3354186 0.3910140
#> 5         0.5 0.3815551 0.02001618 0.3348543 0.3904752
#> 6         0.6 0.3809832 0.02001612 0.3342910 0.3899372
#> 7         0.7 0.3804122 0.02001626 0.3337286 0.3894000
#> 8         0.8 0.3798420 0.02001660 0.3331671 0.3888635
#> 9         0.9 0.3792726 0.02001714 0.3326066 0.3883277
#> 10        1.0 0.3787041 0.02001789 0.3320471 0.3877927

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3601960
#> 2       1        0.1        V5 0.3817508
#> 3       1        0.1        V9 0.3881706
#> 4       1        0.1        V3 0.3482474
#> 5       1        0.1        V4 0.3338859
#> 6       1        0.1        V8 0.3582465
#> 7       1        0.1        V2 0.3605270
#> 8       1        0.1        V6 0.3841487
#> 9       1        0.1        V7 0.3863168
#> 10      1        0.1       V10 0.3939309
#> 11      2        0.2        V7 0.3857215
#> 12      2        0.2        V8 0.3577104
#> 13      2        0.2        V9 0.3877667
#> 14      2        0.2       V10 0.3933500
#> 15      2        0.2        V1 0.3597027
#> 16      2        0.2        V5 0.3810915
#> 17      2        0.2        V2 0.3599483
#> 18      2        0.2        V3 0.3477527
#> 19      2        0.2        V4 0.3332977
#> 20      2        0.2        V6 0.3836114
#> 21      3        0.3        V4 0.3327105
#> 22      3        0.3        V5 0.3804332
#> 23      3        0.3        V3 0.3472587
#> 24      3        0.3        V7 0.3851271
#> 25      3        0.3        V8 0.3571751
#> 26      3        0.3        V9 0.3873633
#> 27      3        0.3        V6 0.3830748
#> 28      3        0.3       V10 0.3927700
#> 29      3        0.3        V1 0.3592102
#> 30      3        0.3        V2 0.3593706
#> 31      4        0.4        V1 0.3587183
#> 32      4        0.4        V9 0.3869603
#> 33      4        0.4        V3 0.3467655
#> 34      4        0.4        V4 0.3321243
#> 35      4        0.4        V5 0.3797762
#> 36      4        0.4        V2 0.3587938
#> 37      4        0.4        V6 0.3825390
#> 38      4        0.4        V7 0.3845336
#> 39      4        0.4        V8 0.3566406
#> 40      4        0.4       V10 0.3921909
#> 41      5        0.5        V8 0.3561069
#> 42      5        0.5        V9 0.3865577
#> 43      5        0.5       V10 0.3916126
#> 44      5        0.5        V1 0.3582271
#> 45      5        0.5        V5 0.3791202
#> 46      5        0.5        V2 0.3582179
#> 47      5        0.5        V3 0.3462729
#> 48      5        0.5        V4 0.3315392
#> 49      5        0.5        V6 0.3820039
#> 50      5        0.5        V7 0.3839410
#> 51      6        0.6        V4 0.3309551
#> 52      6        0.6        V5 0.3784654
#> 53      6        0.6        V7 0.3833493
#> 54      6        0.6        V8 0.3555740
#> 55      6        0.6        V9 0.3861555
#> 56      6        0.6        V6 0.3814696
#> 57      6        0.6       V10 0.3910351
#> 58      6        0.6        V1 0.3577366
#> 59      6        0.6        V2 0.3576429
#> 60      6        0.6        V3 0.3457811
#> 61      7        0.7        V1 0.3572467
#> 62      7        0.7        V3 0.3452899
#> 63      7        0.7        V4 0.3303720
#> 64      7        0.7        V5 0.3778117
#> 65      7        0.7        V9 0.3857537
#> 66      7        0.7        V6 0.3809360
#> 67      7        0.7        V7 0.3827585
#> 68      7        0.7        V8 0.3550419
#> 69      7        0.7        V2 0.3570689
#> 70      7        0.7       V10 0.3904586
#> 71      8        0.8        V9 0.3853524
#> 72      8        0.8       V10 0.3898828
#> 73      8        0.8        V1 0.3567575
#> 74      8        0.8        V5 0.3771591
#> 75      8        0.8        V2 0.3564957
#> 76      8        0.8        V3 0.3447994
#> 77      8        0.8        V4 0.3297900
#> 78      8        0.8        V8 0.3545106
#> 79      8        0.8        V6 0.3804032
#> 80      8        0.8        V7 0.3821687
#> 81      9        0.9        V5 0.3765077
#> 82      9        0.9        V7 0.3815797
#> 83      9        0.9        V8 0.3539801
#> 84      9        0.9        V9 0.3849515
#> 85      9        0.9        V6 0.3798711
#> 86      9        0.9       V10 0.3893079
#> 87      9        0.9        V1 0.3562690
#> 88      9        0.9        V2 0.3559235
#> 89      9        0.9        V3 0.3443097
#> 90      9        0.9        V4 0.3292090
#> 91     10        1.0        V1 0.3557811
#> 92     10        1.0        V4 0.3286290
#> 93     10        1.0        V5 0.3758574
#> 94     10        1.0        V9 0.3845510
#> 95     10        1.0        V3 0.3438206
#> 96     10        1.0        V7 0.3809917
#> 97     10        1.0        V8 0.3534503
#> 98     10        1.0        V2 0.3553522
#> 99     10        1.0        V6 0.3793398
#> 100    10        1.0       V10 0.3887339

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3601960
#> 2       1        0.1        V5 0.3817508
#> 3       1        0.1        V9 0.3881706
#> 4       1        0.1        V3 0.3482474
#> 5       1        0.1        V4 0.3338859
#> 6       1        0.1        V8 0.3582465
#> 7       1        0.1        V2 0.3605270
#> 8       1        0.1        V6 0.3841487
#> 9       1        0.1        V7 0.3863168
#> 10      1        0.1       V10 0.3939309
#> 11      2        0.2        V7 0.3857215
#> 12      2        0.2        V8 0.3577104
#> 13      2        0.2        V9 0.3877667
#> 14      2        0.2       V10 0.3933500
#> 15      2        0.2        V1 0.3597027
#> 16      2        0.2        V5 0.3810915
#> 17      2        0.2        V2 0.3599483
#> 18      2        0.2        V3 0.3477527
#> 19      2        0.2        V4 0.3332977
#> 20      2        0.2        V6 0.3836114
#> 21      3        0.3        V4 0.3327105
#> 22      3        0.3        V5 0.3804332
#> 23      3        0.3        V3 0.3472587
#> 24      3        0.3        V7 0.3851271
#> 25      3        0.3        V8 0.3571751
#> 26      3        0.3        V9 0.3873633
#> 27      3        0.3        V6 0.3830748
#> 28      3        0.3       V10 0.3927700
#> 29      3        0.3        V1 0.3592102
#> 30      3        0.3        V2 0.3593706
#> 31      4        0.4        V1 0.3587183
#> 32      4        0.4        V9 0.3869603
#> 33      4        0.4        V3 0.3467655
#> 34      4        0.4        V4 0.3321243
#> 35      4        0.4        V5 0.3797762
#> 36      4        0.4        V2 0.3587938
#> 37      4        0.4        V6 0.3825390
#> 38      4        0.4        V7 0.3845336
#> 39      4        0.4        V8 0.3566406
#> 40      4        0.4       V10 0.3921909
#> 41      5        0.5        V8 0.3561069
#> 42      5        0.5        V9 0.3865577
#> 43      5        0.5       V10 0.3916126
#> 44      5        0.5        V1 0.3582271
#> 45      5        0.5        V5 0.3791202
#> 46      5        0.5        V2 0.3582179
#> 47      5        0.5        V3 0.3462729
#> 48      5        0.5        V4 0.3315392
#> 49      5        0.5        V6 0.3820039
#> 50      5        0.5        V7 0.3839410
#> 51      6        0.6        V4 0.3309551
#> 52      6        0.6        V5 0.3784654
#> 53      6        0.6        V7 0.3833493
#> 54      6        0.6        V8 0.3555740
#> 55      6        0.6        V9 0.3861555
#> 56      6        0.6        V6 0.3814696
#> 57      6        0.6       V10 0.3910351
#> 58      6        0.6        V1 0.3577366
#> 59      6        0.6        V2 0.3576429
#> 60      6        0.6        V3 0.3457811
#> 61      7        0.7        V1 0.3572467
#> 62      7        0.7        V3 0.3452899
#> 63      7        0.7        V4 0.3303720
#> 64      7        0.7        V5 0.3778117
#> 65      7        0.7        V9 0.3857537
#> 66      7        0.7        V6 0.3809360
#> 67      7        0.7        V7 0.3827585
#> 68      7        0.7        V8 0.3550419
#> 69      7        0.7        V2 0.3570689
#> 70      7        0.7       V10 0.3904586
#> 71      8        0.8        V9 0.3853524
#> 72      8        0.8       V10 0.3898828
#> 73      8        0.8        V1 0.3567575
#> 74      8        0.8        V5 0.3771591
#> 75      8        0.8        V2 0.3564957
#> 76      8        0.8        V3 0.3447994
#> 77      8        0.8        V4 0.3297900
#> 78      8        0.8        V8 0.3545106
#> 79      8        0.8        V6 0.3804032
#> 80      8        0.8        V7 0.3821687
#> 81      9        0.9        V5 0.3765077
#> 82      9        0.9        V7 0.3815797
#> 83      9        0.9        V8 0.3539801
#> 84      9        0.9        V9 0.3849515
#> 85      9        0.9        V6 0.3798711
#> 86      9        0.9       V10 0.3893079
#> 87      9        0.9        V1 0.3562690
#> 88      9        0.9        V2 0.3559235
#> 89      9        0.9        V3 0.3443097
#> 90      9        0.9        V4 0.3292090
#> 91     10        1.0        V1 0.3557811
#> 92     10        1.0        V4 0.3286290
#> 93     10        1.0        V5 0.3758574
#> 94     10        1.0        V9 0.3845510
#> 95     10        1.0        V3 0.3438206
#> 96     10        1.0        V7 0.3809917
#> 97     10        1.0        V8 0.3534503
#> 98     10        1.0        V2 0.3553522
#> 99     10        1.0        V6 0.3793398
#> 100    10        1.0       V10 0.3887339


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
#> 1       0                0          0 0.8757906 0.04475125 0.7905220 0.9466890
#> 2       0               20         20 0.8757906 0.04475125 0.7905220 0.9466890
#> 3       0               40         40 0.8757906 0.04475125 0.7905220 0.9466890
#> 4       0               60         60 0.8757906 0.04475125 0.7905220 0.9466890
#> 5      20                0         20 0.8617131 0.04747926 0.7707350 0.9387013
#> 6      20               20         40 0.8617131 0.04747926 0.7707350 0.9387013
#> 7      20               40         60 0.8617131 0.04747926 0.7707350 0.9387013
#> 8      20               60         80 0.8617131 0.04747926 0.7707350 0.9387013
#> 9      40                0         40 0.8478591 0.05003841 0.7516781 0.9306103
#> 10     40               20         60 0.8478591 0.05003841 0.7516781 0.9306103
#> 11     40               40         80 0.8478591 0.05003841 0.7516781 0.9306103
#> 12     40               60        100 0.8478591 0.05003841 0.7516781 0.9306103
#> 13     60                0         60 0.8342249 0.05244757 0.7332779 0.9224427
#> 14     60               20         80 0.8342249 0.05244757 0.7332779 0.9224427
#> 15     60               40        100 0.8342249 0.05244757 0.7332779 0.9224427
#> 16     60               60        120 0.8342249 0.05244757 0.7332779 0.9224427
#> 17     80                0         80 0.8208071 0.05472148 0.7154771 0.9142198
#> 18     80               20        100 0.8208071 0.05472148 0.7154771 0.9142198
#> 19     80               40        120 0.8208071 0.05472148 0.7154771 0.9142198
#> 20     80               60        140 0.8208071 0.05472148 0.7154771 0.9142198
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12973474 0.139948434 0.5643049
#> 2  0.30574618 0.12910393 0.109094124 0.5074612
#> 3  0.26001915 0.12838002 0.084233990 0.4567363
#> 4  0.22113100 0.12710942 0.064302363 0.4298461
#> 5  0.25589195 0.11742111 0.074670871 0.4622868
#> 6  0.21762106 0.11417583 0.056671401 0.4165264
#> 7  0.18507391 0.11115298 0.042389699 0.3758466
#> 8  0.15739448 0.10805787 0.031167414 0.3396966
#> 9  0.18213629 0.10345691 0.036973601 0.3802933
#> 10 0.15489621 0.09911587 0.026948475 0.3436485
#> 11 0.13173012 0.09511516 0.019216906 0.3110730
#> 12 0.11202872 0.09126108 0.013353742 0.2820847
#> 13 0.12963921 0.08969216 0.016358349 0.3146354
#> 14 0.11025053 0.08506376 0.011218249 0.2852569
#> 15 0.09376159 0.08080312 0.007442880 0.2590762
#> 16 0.07973872 0.07678020 0.004747959 0.2356949
#> 17 0.09227334 0.07699649 0.006106430 0.2619433
#> 18 0.07847305 0.07250251 0.003817993 0.2382584
#> 19 0.06673672 0.06836231 0.002270008 0.2170563
#> 20 0.05675565 0.06448587 0.001271529 0.1980170
```
