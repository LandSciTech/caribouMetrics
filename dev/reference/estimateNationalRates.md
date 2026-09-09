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
#> 1         0.1 0.3838513 0.02706093 0.3392684 0.4247380
#> 2         0.2 0.3832760 0.02702685 0.3387788 0.4241534
#> 3         0.3 0.3827015 0.02699293 0.3382900 0.4235695
#> 4         0.4 0.3821279 0.02695916 0.3378018 0.4229865
#> 5         0.5 0.3815551 0.02692555 0.3373144 0.4224042
#> 6         0.6 0.3809832 0.02689210 0.3368277 0.4218228
#> 7         0.7 0.3804122 0.02685881 0.3363417 0.4212422
#> 8         0.8 0.3798420 0.02682567 0.3358563 0.4206623
#> 9         0.9 0.3792726 0.02679269 0.3353717 0.4200833
#> 10        1.0 0.3787041 0.02675986 0.3348878 0.4195051

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4259400
#> 2       1        0.1        V5 0.3322579
#> 3       1        0.1        V9 0.3889414
#> 4       1        0.1        V3 0.3734788
#> 5       1        0.1        V4 0.3859903
#> 6       1        0.1        V8 0.4205979
#> 7       1        0.1        V2 0.3898619
#> 8       1        0.1        V6 0.3796976
#> 9       1        0.1        V7 0.3634157
#> 10      1        0.1       V10 0.3994205
#> 11      2        0.2        V7 0.3628603
#> 12      2        0.2        V8 0.4200064
#> 13      2        0.2        V9 0.3883811
#> 14      2        0.2       V10 0.3989226
#> 15      2        0.2        V1 0.4253573
#> 16      2        0.2        V5 0.3317874
#> 17      2        0.2        V2 0.3892310
#> 18      2        0.2        V3 0.3730411
#> 19      2        0.2        V4 0.3853702
#> 20      2        0.2        V6 0.3792067
#> 21      3        0.3        V4 0.3847512
#> 22      3        0.3        V5 0.3313176
#> 23      3        0.3        V3 0.3726040
#> 24      3        0.3        V7 0.3623058
#> 25      3        0.3        V8 0.4194157
#> 26      3        0.3        V9 0.3878216
#> 27      3        0.3        V6 0.3787165
#> 28      3        0.3       V10 0.3984253
#> 29      3        0.3        V1 0.4247755
#> 30      3        0.3        V2 0.3886012
#> 31      4        0.4        V1 0.4241944
#> 32      4        0.4        V9 0.3872629
#> 33      4        0.4        V3 0.3721674
#> 34      4        0.4        V4 0.3841331
#> 35      4        0.4        V5 0.3308485
#> 36      4        0.4        V2 0.3879724
#> 37      4        0.4        V6 0.3782269
#> 38      4        0.4        V7 0.3617522
#> 39      4        0.4        V8 0.4188258
#> 40      4        0.4       V10 0.3979287
#> 41      5        0.5        V8 0.4182368
#> 42      5        0.5        V9 0.3867050
#> 43      5        0.5       V10 0.3974326
#> 44      5        0.5        V1 0.4236141
#> 45      5        0.5        V5 0.3303801
#> 46      5        0.5        V2 0.3873446
#> 47      5        0.5        V3 0.3717313
#> 48      5        0.5        V4 0.3835160
#> 49      5        0.5        V6 0.3777379
#> 50      5        0.5        V7 0.3611994
#> 51      6        0.6        V4 0.3828999
#> 52      6        0.6        V5 0.3299123
#> 53      6        0.6        V7 0.3606475
#> 54      6        0.6        V8 0.4176486
#> 55      6        0.6        V9 0.3861479
#> 56      6        0.6        V6 0.3772496
#> 57      6        0.6       V10 0.3969372
#> 58      6        0.6        V1 0.4230347
#> 59      6        0.6        V2 0.3867179
#> 60      6        0.6        V3 0.3712957
#> 61      7        0.7        V1 0.4224560
#> 62      7        0.7        V3 0.3708607
#> 63      7        0.7        V4 0.3822848
#> 64      7        0.7        V5 0.3294451
#> 65      7        0.7        V9 0.3855916
#> 66      7        0.7        V6 0.3767619
#> 67      7        0.7        V7 0.3600964
#> 68      7        0.7        V8 0.4170613
#> 69      7        0.7        V2 0.3860921
#> 70      7        0.7       V10 0.3964423
#> 71      8        0.8        V9 0.3850361
#> 72      8        0.8       V10 0.3959481
#> 73      8        0.8        V1 0.4218781
#> 74      8        0.8        V5 0.3289787
#> 75      8        0.8        V2 0.3854674
#> 76      8        0.8        V3 0.3704261
#> 77      8        0.8        V4 0.3816707
#> 78      8        0.8        V8 0.4164747
#> 79      8        0.8        V6 0.3762749
#> 80      8        0.8        V7 0.3595461
#> 81      9        0.9        V5 0.3285129
#> 82      9        0.9        V7 0.3589967
#> 83      9        0.9        V8 0.4158890
#> 84      9        0.9        V9 0.3844814
#> 85      9        0.9        V6 0.3757884
#> 86      9        0.9       V10 0.3954546
#> 87      9        0.9        V1 0.4213010
#> 88      9        0.9        V2 0.3848437
#> 89      9        0.9        V3 0.3699920
#> 90      9        0.9        V4 0.3810575
#> 91     10        1.0        V1 0.4207247
#> 92     10        1.0        V4 0.3804454
#> 93     10        1.0        V5 0.3280477
#> 94     10        1.0        V9 0.3839275
#> 95     10        1.0        V3 0.3695585
#> 96     10        1.0        V7 0.3584481
#> 97     10        1.0        V8 0.4153041
#> 98     10        1.0        V2 0.3842209
#> 99     10        1.0        V6 0.3753026
#> 100    10        1.0       V10 0.3949616

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4259400
#> 2       1        0.1        V5 0.3322579
#> 3       1        0.1        V9 0.3889414
#> 4       1        0.1        V3 0.3734788
#> 5       1        0.1        V4 0.3859903
#> 6       1        0.1        V8 0.4205979
#> 7       1        0.1        V2 0.3898619
#> 8       1        0.1        V6 0.3796976
#> 9       1        0.1        V7 0.3634157
#> 10      1        0.1       V10 0.3994205
#> 11      2        0.2        V7 0.3628603
#> 12      2        0.2        V8 0.4200064
#> 13      2        0.2        V9 0.3883811
#> 14      2        0.2       V10 0.3989226
#> 15      2        0.2        V1 0.4253573
#> 16      2        0.2        V5 0.3317874
#> 17      2        0.2        V2 0.3892310
#> 18      2        0.2        V3 0.3730411
#> 19      2        0.2        V4 0.3853702
#> 20      2        0.2        V6 0.3792067
#> 21      3        0.3        V4 0.3847512
#> 22      3        0.3        V5 0.3313176
#> 23      3        0.3        V3 0.3726040
#> 24      3        0.3        V7 0.3623058
#> 25      3        0.3        V8 0.4194157
#> 26      3        0.3        V9 0.3878216
#> 27      3        0.3        V6 0.3787165
#> 28      3        0.3       V10 0.3984253
#> 29      3        0.3        V1 0.4247755
#> 30      3        0.3        V2 0.3886012
#> 31      4        0.4        V1 0.4241944
#> 32      4        0.4        V9 0.3872629
#> 33      4        0.4        V3 0.3721674
#> 34      4        0.4        V4 0.3841331
#> 35      4        0.4        V5 0.3308485
#> 36      4        0.4        V2 0.3879724
#> 37      4        0.4        V6 0.3782269
#> 38      4        0.4        V7 0.3617522
#> 39      4        0.4        V8 0.4188258
#> 40      4        0.4       V10 0.3979287
#> 41      5        0.5        V8 0.4182368
#> 42      5        0.5        V9 0.3867050
#> 43      5        0.5       V10 0.3974326
#> 44      5        0.5        V1 0.4236141
#> 45      5        0.5        V5 0.3303801
#> 46      5        0.5        V2 0.3873446
#> 47      5        0.5        V3 0.3717313
#> 48      5        0.5        V4 0.3835160
#> 49      5        0.5        V6 0.3777379
#> 50      5        0.5        V7 0.3611994
#> 51      6        0.6        V4 0.3828999
#> 52      6        0.6        V5 0.3299123
#> 53      6        0.6        V7 0.3606475
#> 54      6        0.6        V8 0.4176486
#> 55      6        0.6        V9 0.3861479
#> 56      6        0.6        V6 0.3772496
#> 57      6        0.6       V10 0.3969372
#> 58      6        0.6        V1 0.4230347
#> 59      6        0.6        V2 0.3867179
#> 60      6        0.6        V3 0.3712957
#> 61      7        0.7        V1 0.4224560
#> 62      7        0.7        V3 0.3708607
#> 63      7        0.7        V4 0.3822848
#> 64      7        0.7        V5 0.3294451
#> 65      7        0.7        V9 0.3855916
#> 66      7        0.7        V6 0.3767619
#> 67      7        0.7        V7 0.3600964
#> 68      7        0.7        V8 0.4170613
#> 69      7        0.7        V2 0.3860921
#> 70      7        0.7       V10 0.3964423
#> 71      8        0.8        V9 0.3850361
#> 72      8        0.8       V10 0.3959481
#> 73      8        0.8        V1 0.4218781
#> 74      8        0.8        V5 0.3289787
#> 75      8        0.8        V2 0.3854674
#> 76      8        0.8        V3 0.3704261
#> 77      8        0.8        V4 0.3816707
#> 78      8        0.8        V8 0.4164747
#> 79      8        0.8        V6 0.3762749
#> 80      8        0.8        V7 0.3595461
#> 81      9        0.9        V5 0.3285129
#> 82      9        0.9        V7 0.3589967
#> 83      9        0.9        V8 0.4158890
#> 84      9        0.9        V9 0.3844814
#> 85      9        0.9        V6 0.3757884
#> 86      9        0.9       V10 0.3954546
#> 87      9        0.9        V1 0.4213010
#> 88      9        0.9        V2 0.3848437
#> 89      9        0.9        V3 0.3699920
#> 90      9        0.9        V4 0.3810575
#> 91     10        1.0        V1 0.4207247
#> 92     10        1.0        V4 0.3804454
#> 93     10        1.0        V5 0.3280477
#> 94     10        1.0        V9 0.3839275
#> 95     10        1.0        V3 0.3695585
#> 96     10        1.0        V7 0.3584481
#> 97     10        1.0        V8 0.4153041
#> 98     10        1.0        V2 0.3842209
#> 99     10        1.0        V6 0.3753026
#> 100    10        1.0       V10 0.3949616


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
#> 1       0                0          0 0.8757906 0.04706552 0.7899735 0.9481530
#> 2       0               20         20 0.8757906 0.04706552 0.7899735 0.9481530
#> 3       0               40         40 0.8757906 0.04706552 0.7899735 0.9481530
#> 4       0               60         60 0.8757906 0.04706552 0.7899735 0.9481530
#> 5      20                0         20 0.8617131 0.04836202 0.7741422 0.9371927
#> 6      20               20         40 0.8617131 0.04836202 0.7741422 0.9371927
#> 7      20               40         60 0.8617131 0.04836202 0.7741422 0.9371927
#> 8      20               60         80 0.8617131 0.04836202 0.7741422 0.9371927
#> 9      40                0         40 0.8478591 0.04953320 0.7587866 0.9260375
#> 10     40               20         60 0.8478591 0.04953320 0.7587866 0.9260375
#> 11     40               40         80 0.8478591 0.04953320 0.7587866 0.9260375
#> 12     40               60        100 0.8478591 0.04953320 0.7587866 0.9260375
#> 13     60                0         60 0.8342249 0.05059791 0.7438664 0.9147528
#> 14     60               20         80 0.8342249 0.05059791 0.7438664 0.9147528
#> 15     60               40        100 0.8342249 0.05059791 0.7438664 0.9147528
#> 16     60               60        120 0.8342249 0.05059791 0.7438664 0.9147528
#> 17     80                0         80 0.8208071 0.05157070 0.7293486 0.9033880
#> 18     80               20        100 0.8208071 0.05157070 0.7293486 0.9033880
#> 19     80               40        120 0.8208071 0.05157070 0.7293486 0.9033880
#> 20     80               60        140 0.8208071 0.05157070 0.7293486 0.9033880
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11942096 0.181505429 0.5843110
#> 2  0.30574618 0.11301517 0.143267274 0.5094531
#> 3  0.26001915 0.10748773 0.112250650 0.4447264
#> 4  0.22113100 0.10234833 0.087165984 0.3889952
#> 5  0.25589195 0.10815095 0.103056569 0.4648546
#> 6  0.21762106 0.10031658 0.079751667 0.4063102
#> 7  0.18507391 0.09351863 0.061024839 0.3559723
#> 8  0.15739448 0.08739084 0.046077885 0.3127159
#> 9  0.18213629 0.09513516 0.055520703 0.3716081
#> 10 0.15489621 0.08725930 0.041709446 0.3261532
#> 11 0.13173012 0.08035405 0.030818273 0.2870778
#> 12 0.11202872 0.07415491 0.022329870 0.2534353
#> 13 0.12963921 0.08224135 0.027666720 0.2992199
#> 14 0.11025053 0.07491316 0.019897083 0.2638969
#> 15 0.09376159 0.06842675 0.013963564 0.2334363
#> 16 0.07973872 0.06259030 0.009519740 0.2070852
#> 17 0.09227334 0.07029954 0.012291134 0.2429149
#> 18 0.07847305 0.06373583 0.008287083 0.2152955
#> 19 0.06673672 0.05787746 0.005387638 0.1913389
#> 20 0.05675565 0.05258372 0.003354480 0.1704617
```
