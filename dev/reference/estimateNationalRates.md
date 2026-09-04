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
#> 1         0.1 0.3838513 0.01824662 0.3463318 0.3917782
#> 2         0.2 0.3832760 0.01823970 0.3457956 0.3912539
#> 3         0.3 0.3827015 0.01823287 0.3452602 0.3907302
#> 4         0.4 0.3821279 0.01822613 0.3447256 0.3902072
#> 5         0.5 0.3815551 0.01821947 0.3441919 0.3896850
#> 6         0.6 0.3809832 0.01821289 0.3436590 0.3891634
#> 7         0.7 0.3804122 0.01820640 0.3431270 0.3886425
#> 8         0.8 0.3798420 0.01820000 0.3425957 0.3881224
#> 9         0.9 0.3792726 0.01819367 0.3420653 0.3876029
#> 10        1.0 0.3787041 0.01818743 0.3415357 0.3870841

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3492871
#> 2       1        0.1        V5 0.3454738
#> 3       1        0.1        V9 0.3901837
#> 4       1        0.1        V3 0.3682380
#> 5       1        0.1        V4 0.3539687
#> 6       1        0.1        V8 0.3509185
#> 7       1        0.1        V2 0.3890108
#> 8       1        0.1        V6 0.3558933
#> 9       1        0.1        V7 0.3622389
#> 10      1        0.1       V10 0.3922411
#> 11      2        0.2        V7 0.3617147
#> 12      2        0.2        V8 0.3503614
#> 13      2        0.2        V9 0.3896019
#> 14      2        0.2       V10 0.3917334
#> 15      2        0.2        V1 0.3487873
#> 16      2        0.2        V5 0.3449270
#> 17      2        0.2        V2 0.3884669
#> 18      2        0.2        V3 0.3677662
#> 19      2        0.2        V4 0.3534055
#> 20      2        0.2        V6 0.3554553
#> 21      3        0.3        V4 0.3528432
#> 22      3        0.3        V5 0.3443811
#> 23      3        0.3        V3 0.3672950
#> 24      3        0.3        V7 0.3611913
#> 25      3        0.3        V8 0.3498052
#> 26      3        0.3        V9 0.3890210
#> 27      3        0.3        V6 0.3550178
#> 28      3        0.3       V10 0.3912264
#> 29      3        0.3        V1 0.3482882
#> 30      3        0.3        V2 0.3879237
#> 31      4        0.4        V1 0.3477898
#> 32      4        0.4        V9 0.3884409
#> 33      4        0.4        V3 0.3668244
#> 34      4        0.4        V4 0.3522818
#> 35      4        0.4        V5 0.3438360
#> 36      4        0.4        V2 0.3873813
#> 37      4        0.4        V6 0.3545809
#> 38      4        0.4        V7 0.3606686
#> 39      4        0.4        V8 0.3492498
#> 40      4        0.4       V10 0.3907200
#> 41      5        0.5        V8 0.3486953
#> 42      5        0.5        V9 0.3878617
#> 43      5        0.5       V10 0.3902143
#> 44      5        0.5        V1 0.3472922
#> 45      5        0.5        V5 0.3432919
#> 46      5        0.5        V2 0.3868397
#> 47      5        0.5        V3 0.3663543
#> 48      5        0.5        V4 0.3517213
#> 49      5        0.5        V6 0.3541445
#> 50      5        0.5        V7 0.3601467
#> 51      6        0.6        V4 0.3511616
#> 52      6        0.6        V5 0.3427485
#> 53      6        0.6        V7 0.3596256
#> 54      6        0.6        V8 0.3481417
#> 55      6        0.6        V9 0.3872834
#> 56      6        0.6        V6 0.3537086
#> 57      6        0.6       V10 0.3897092
#> 58      6        0.6        V1 0.3467952
#> 59      6        0.6        V2 0.3862988
#> 60      6        0.6        V3 0.3658849
#> 61      7        0.7        V1 0.3462990
#> 62      7        0.7        V3 0.3654161
#> 63      7        0.7        V4 0.3506029
#> 64      7        0.7        V5 0.3422060
#> 65      7        0.7        V9 0.3867059
#> 66      7        0.7        V6 0.3532733
#> 67      7        0.7        V7 0.3591052
#> 68      7        0.7        V8 0.3475890
#> 69      7        0.7        V2 0.3857587
#> 70      7        0.7       V10 0.3892048
#> 71      8        0.8        V9 0.3861293
#> 72      8        0.8       V10 0.3887010
#> 73      8        0.8        V1 0.3458035
#> 74      8        0.8        V5 0.3416644
#> 75      8        0.8        V2 0.3852193
#> 76      8        0.8        V3 0.3649479
#> 77      8        0.8        V4 0.3500451
#> 78      8        0.8        V8 0.3470372
#> 79      8        0.8        V6 0.3528385
#> 80      8        0.8        V7 0.3585855
#> 81      9        0.9        V5 0.3411237
#> 82      9        0.9        V7 0.3580666
#> 83      9        0.9        V8 0.3464862
#> 84      9        0.9        V9 0.3855535
#> 85      9        0.9        V6 0.3524042
#> 86      9        0.9       V10 0.3881979
#> 87      9        0.9        V1 0.3453087
#> 88      9        0.9        V2 0.3846807
#> 89      9        0.9        V3 0.3644803
#> 90      9        0.9        V4 0.3494881
#> 91     10        1.0        V1 0.3448146
#> 92     10        1.0        V4 0.3489321
#> 93     10        1.0        V5 0.3405838
#> 94     10        1.0        V9 0.3849786
#> 95     10        1.0        V3 0.3640133
#> 96     10        1.0        V7 0.3575485
#> 97     10        1.0        V8 0.3459361
#> 98     10        1.0        V2 0.3841428
#> 99     10        1.0        V6 0.3519705
#> 100    10        1.0       V10 0.3876954

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3492871
#> 2       1        0.1        V5 0.3454738
#> 3       1        0.1        V9 0.3901837
#> 4       1        0.1        V3 0.3682380
#> 5       1        0.1        V4 0.3539687
#> 6       1        0.1        V8 0.3509185
#> 7       1        0.1        V2 0.3890108
#> 8       1        0.1        V6 0.3558933
#> 9       1        0.1        V7 0.3622389
#> 10      1        0.1       V10 0.3922411
#> 11      2        0.2        V7 0.3617147
#> 12      2        0.2        V8 0.3503614
#> 13      2        0.2        V9 0.3896019
#> 14      2        0.2       V10 0.3917334
#> 15      2        0.2        V1 0.3487873
#> 16      2        0.2        V5 0.3449270
#> 17      2        0.2        V2 0.3884669
#> 18      2        0.2        V3 0.3677662
#> 19      2        0.2        V4 0.3534055
#> 20      2        0.2        V6 0.3554553
#> 21      3        0.3        V4 0.3528432
#> 22      3        0.3        V5 0.3443811
#> 23      3        0.3        V3 0.3672950
#> 24      3        0.3        V7 0.3611913
#> 25      3        0.3        V8 0.3498052
#> 26      3        0.3        V9 0.3890210
#> 27      3        0.3        V6 0.3550178
#> 28      3        0.3       V10 0.3912264
#> 29      3        0.3        V1 0.3482882
#> 30      3        0.3        V2 0.3879237
#> 31      4        0.4        V1 0.3477898
#> 32      4        0.4        V9 0.3884409
#> 33      4        0.4        V3 0.3668244
#> 34      4        0.4        V4 0.3522818
#> 35      4        0.4        V5 0.3438360
#> 36      4        0.4        V2 0.3873813
#> 37      4        0.4        V6 0.3545809
#> 38      4        0.4        V7 0.3606686
#> 39      4        0.4        V8 0.3492498
#> 40      4        0.4       V10 0.3907200
#> 41      5        0.5        V8 0.3486953
#> 42      5        0.5        V9 0.3878617
#> 43      5        0.5       V10 0.3902143
#> 44      5        0.5        V1 0.3472922
#> 45      5        0.5        V5 0.3432919
#> 46      5        0.5        V2 0.3868397
#> 47      5        0.5        V3 0.3663543
#> 48      5        0.5        V4 0.3517213
#> 49      5        0.5        V6 0.3541445
#> 50      5        0.5        V7 0.3601467
#> 51      6        0.6        V4 0.3511616
#> 52      6        0.6        V5 0.3427485
#> 53      6        0.6        V7 0.3596256
#> 54      6        0.6        V8 0.3481417
#> 55      6        0.6        V9 0.3872834
#> 56      6        0.6        V6 0.3537086
#> 57      6        0.6       V10 0.3897092
#> 58      6        0.6        V1 0.3467952
#> 59      6        0.6        V2 0.3862988
#> 60      6        0.6        V3 0.3658849
#> 61      7        0.7        V1 0.3462990
#> 62      7        0.7        V3 0.3654161
#> 63      7        0.7        V4 0.3506029
#> 64      7        0.7        V5 0.3422060
#> 65      7        0.7        V9 0.3867059
#> 66      7        0.7        V6 0.3532733
#> 67      7        0.7        V7 0.3591052
#> 68      7        0.7        V8 0.3475890
#> 69      7        0.7        V2 0.3857587
#> 70      7        0.7       V10 0.3892048
#> 71      8        0.8        V9 0.3861293
#> 72      8        0.8       V10 0.3887010
#> 73      8        0.8        V1 0.3458035
#> 74      8        0.8        V5 0.3416644
#> 75      8        0.8        V2 0.3852193
#> 76      8        0.8        V3 0.3649479
#> 77      8        0.8        V4 0.3500451
#> 78      8        0.8        V8 0.3470372
#> 79      8        0.8        V6 0.3528385
#> 80      8        0.8        V7 0.3585855
#> 81      9        0.9        V5 0.3411237
#> 82      9        0.9        V7 0.3580666
#> 83      9        0.9        V8 0.3464862
#> 84      9        0.9        V9 0.3855535
#> 85      9        0.9        V6 0.3524042
#> 86      9        0.9       V10 0.3881979
#> 87      9        0.9        V1 0.3453087
#> 88      9        0.9        V2 0.3846807
#> 89      9        0.9        V3 0.3644803
#> 90      9        0.9        V4 0.3494881
#> 91     10        1.0        V1 0.3448146
#> 92     10        1.0        V4 0.3489321
#> 93     10        1.0        V5 0.3405838
#> 94     10        1.0        V9 0.3849786
#> 95     10        1.0        V3 0.3640133
#> 96     10        1.0        V7 0.3575485
#> 97     10        1.0        V8 0.3459361
#> 98     10        1.0        V2 0.3841428
#> 99     10        1.0        V6 0.3519705
#> 100    10        1.0       V10 0.3876954


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
#> 1       0                0          0 0.8757906 0.04695874 0.7915018 0.9475659
#> 2       0               20         20 0.8757906 0.04695874 0.7915018 0.9475659
#> 3       0               40         40 0.8757906 0.04695874 0.7915018 0.9475659
#> 4       0               60         60 0.8757906 0.04695874 0.7915018 0.9475659
#> 5      20                0         20 0.8617131 0.04951946 0.7743969 0.9401207
#> 6      20               20         40 0.8617131 0.04951946 0.7743969 0.9401207
#> 7      20               40         60 0.8617131 0.04951946 0.7743969 0.9401207
#> 8      20               60         80 0.8617131 0.04951946 0.7743969 0.9401207
#> 9      40                0         40 0.8478591 0.05198585 0.7578379 0.9325709
#> 10     40               20         60 0.8478591 0.05198585 0.7578379 0.9325709
#> 11     40               40         80 0.8478591 0.05198585 0.7578379 0.9325709
#> 12     40               60        100 0.8478591 0.05198585 0.7578379 0.9325709
#> 13     60                0         60 0.8342249 0.05435966 0.7417765 0.9249404
#> 14     60               20         80 0.8342249 0.05435966 0.7417765 0.9249404
#> 15     60               40        100 0.8342249 0.05435966 0.7417765 0.9249404
#> 16     60               60        120 0.8342249 0.05435966 0.7417765 0.9249404
#> 17     80                0         80 0.8208071 0.05664312 0.7261737 0.9172487
#> 18     80               20        100 0.8208071 0.05664312 0.7261737 0.9172487
#> 19     80               40        120 0.8208071 0.05664312 0.7261737 0.9172487
#> 20     80               60        140 0.8208071 0.05664312 0.7261737 0.9172487
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11815304 0.184329287 0.5638835
#> 2  0.30574618 0.11433382 0.140630963 0.5115142
#> 3  0.26001915 0.11089331 0.106235145 0.4642971
#> 4  0.22113100 0.10741594 0.079274115 0.4218256
#> 5  0.25589195 0.10874363 0.103645434 0.4531401
#> 6  0.21762106 0.10356709 0.077250877 0.4117998
#> 7  0.18507391 0.09881763 0.056711246 0.3746724
#> 8  0.15739448 0.09426212 0.040876920 0.3413437
#> 9  0.18213629 0.09719905 0.055177859 0.3659129
#> 10 0.15489621 0.09162307 0.039702407 0.3334805
#> 11 0.13173012 0.08649179 0.027933165 0.3043611
#> 12 0.11202872 0.08166145 0.019127613 0.2781992
#> 13 0.12963921 0.08518781 0.027068192 0.2974884
#> 14 0.11025053 0.07974800 0.018487563 0.2720207
#> 15 0.09376159 0.07474047 0.012213632 0.2491086
#> 16 0.07973872 0.07006919 0.007749898 0.2284647
#> 17 0.09227334 0.07365532 0.011764640 0.2436918
#> 18 0.07847305 0.06862273 0.007436323 0.2235787
#> 19 0.06673672 0.06399050 0.004472275 0.2054142
#> 20 0.05675565 0.05968948 0.002532023 0.1889725
```
