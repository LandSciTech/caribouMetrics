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
#> 1         0.1 0.3838513 0.01877244 0.3506951 0.4002366
#> 2         0.2 0.3832760 0.01875502 0.3501962 0.3996818
#> 3         0.3 0.3827015 0.01873767 0.3496980 0.3991277
#> 4         0.4 0.3821279 0.01872038 0.3492005 0.3985744
#> 5         0.5 0.3815551 0.01870315 0.3487037 0.3980218
#> 6         0.6 0.3809832 0.01868599 0.3482077 0.3974700
#> 7         0.7 0.3804122 0.01866890 0.3477123 0.3969190
#> 8         0.8 0.3798420 0.01865187 0.3472176 0.3963688
#> 9         0.9 0.3792726 0.01863491 0.3467237 0.3958193
#> 10        1.0 0.3787041 0.01861800 0.3462304 0.3952705

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3926426
#> 2       1        0.1        V5 0.3638412
#> 3       1        0.1        V9 0.3920499
#> 4       1        0.1        V3 0.3924543
#> 5       1        0.1        V4 0.3901944
#> 6       1        0.1        V8 0.3559097
#> 7       1        0.1        V2 0.3491812
#> 8       1        0.1        V6 0.3828409
#> 9       1        0.1        V7 0.3617668
#> 10      1        0.1       V10 0.4024413
#> 11      2        0.2        V7 0.3612133
#> 12      2        0.2        V8 0.3553617
#> 13      2        0.2        V9 0.3914571
#> 14      2        0.2       V10 0.4018815
#> 15      2        0.2        V1 0.3921049
#> 16      2        0.2        V5 0.3632376
#> 17      2        0.2        V2 0.3486965
#> 18      2        0.2        V3 0.3918563
#> 19      2        0.2        V4 0.3896208
#> 20      2        0.2        V6 0.3822408
#> 21      3        0.3        V4 0.3890479
#> 22      3        0.3        V5 0.3626350
#> 23      3        0.3        V3 0.3912593
#> 24      3        0.3        V7 0.3606607
#> 25      3        0.3        V8 0.3548146
#> 26      3        0.3        V9 0.3908652
#> 27      3        0.3        V6 0.3816417
#> 28      3        0.3       V10 0.4013225
#> 29      3        0.3        V1 0.3915679
#> 30      3        0.3        V2 0.3482125
#> 31      4        0.4        V1 0.3910317
#> 32      4        0.4        V9 0.3902741
#> 33      4        0.4        V3 0.3906631
#> 34      4        0.4        V4 0.3884760
#> 35      4        0.4        V5 0.3620333
#> 36      4        0.4        V2 0.3477292
#> 37      4        0.4        V6 0.3810435
#> 38      4        0.4        V7 0.3601089
#> 39      4        0.4        V8 0.3542683
#> 40      4        0.4       V10 0.4007642
#> 41      5        0.5        V8 0.3537228
#> 42      5        0.5        V9 0.3896840
#> 43      5        0.5       V10 0.4002067
#> 44      5        0.5        V1 0.3904961
#> 45      5        0.5        V5 0.3614327
#> 46      5        0.5        V2 0.3472466
#> 47      5        0.5        V3 0.3900679
#> 48      5        0.5        V4 0.3879049
#> 49      5        0.5        V6 0.3804462
#> 50      5        0.5        V7 0.3595580
#> 51      6        0.6        V4 0.3873346
#> 52      6        0.6        V5 0.3608331
#> 53      6        0.6        V7 0.3590079
#> 54      6        0.6        V8 0.3531782
#> 55      6        0.6        V9 0.3890947
#> 56      6        0.6        V6 0.3798499
#> 57      6        0.6       V10 0.3996500
#> 58      6        0.6        V1 0.3899614
#> 59      6        0.6        V2 0.3467646
#> 60      6        0.6        V3 0.3894735
#> 61      7        0.7        V1 0.3894273
#> 62      7        0.7        V3 0.3888801
#> 63      7        0.7        V4 0.3867651
#> 64      7        0.7        V5 0.3602345
#> 65      7        0.7        V9 0.3885064
#> 66      7        0.7        V6 0.3792545
#> 67      7        0.7        V7 0.3584587
#> 68      7        0.7        V8 0.3526344
#> 69      7        0.7        V2 0.3462833
#> 70      7        0.7       V10 0.3990940
#> 71      8        0.8        V9 0.3879189
#> 72      8        0.8       V10 0.3985389
#> 73      8        0.8        V1 0.3888940
#> 74      8        0.8        V5 0.3596369
#> 75      8        0.8        V2 0.3458027
#> 76      8        0.8        V3 0.3882875
#> 77      8        0.8        V4 0.3861965
#> 78      8        0.8        V8 0.3520914
#> 79      8        0.8        V6 0.3786600
#> 80      8        0.8        V7 0.3579103
#> 81      9        0.9        V5 0.3590402
#> 82      9        0.9        V7 0.3573627
#> 83      9        0.9        V8 0.3515493
#> 84      9        0.9        V9 0.3873323
#> 85      9        0.9        V6 0.3780665
#> 86      9        0.9       V10 0.3979845
#> 87      9        0.9        V1 0.3883614
#> 88      9        0.9        V2 0.3453227
#> 89      9        0.9        V3 0.3876959
#> 90      9        0.9        V4 0.3856287
#> 91     10        1.0        V1 0.3878295
#> 92     10        1.0        V4 0.3850618
#> 93     10        1.0        V5 0.3584446
#> 94     10        1.0        V9 0.3867466
#> 95     10        1.0        V3 0.3871052
#> 96     10        1.0        V7 0.3568160
#> 97     10        1.0        V8 0.3510080
#> 98     10        1.0        V2 0.3448434
#> 99     10        1.0        V6 0.3774739
#> 100    10        1.0       V10 0.3974308

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3926426
#> 2       1        0.1        V5 0.3638412
#> 3       1        0.1        V9 0.3920499
#> 4       1        0.1        V3 0.3924543
#> 5       1        0.1        V4 0.3901944
#> 6       1        0.1        V8 0.3559097
#> 7       1        0.1        V2 0.3491812
#> 8       1        0.1        V6 0.3828409
#> 9       1        0.1        V7 0.3617668
#> 10      1        0.1       V10 0.4024413
#> 11      2        0.2        V7 0.3612133
#> 12      2        0.2        V8 0.3553617
#> 13      2        0.2        V9 0.3914571
#> 14      2        0.2       V10 0.4018815
#> 15      2        0.2        V1 0.3921049
#> 16      2        0.2        V5 0.3632376
#> 17      2        0.2        V2 0.3486965
#> 18      2        0.2        V3 0.3918563
#> 19      2        0.2        V4 0.3896208
#> 20      2        0.2        V6 0.3822408
#> 21      3        0.3        V4 0.3890479
#> 22      3        0.3        V5 0.3626350
#> 23      3        0.3        V3 0.3912593
#> 24      3        0.3        V7 0.3606607
#> 25      3        0.3        V8 0.3548146
#> 26      3        0.3        V9 0.3908652
#> 27      3        0.3        V6 0.3816417
#> 28      3        0.3       V10 0.4013225
#> 29      3        0.3        V1 0.3915679
#> 30      3        0.3        V2 0.3482125
#> 31      4        0.4        V1 0.3910317
#> 32      4        0.4        V9 0.3902741
#> 33      4        0.4        V3 0.3906631
#> 34      4        0.4        V4 0.3884760
#> 35      4        0.4        V5 0.3620333
#> 36      4        0.4        V2 0.3477292
#> 37      4        0.4        V6 0.3810435
#> 38      4        0.4        V7 0.3601089
#> 39      4        0.4        V8 0.3542683
#> 40      4        0.4       V10 0.4007642
#> 41      5        0.5        V8 0.3537228
#> 42      5        0.5        V9 0.3896840
#> 43      5        0.5       V10 0.4002067
#> 44      5        0.5        V1 0.3904961
#> 45      5        0.5        V5 0.3614327
#> 46      5        0.5        V2 0.3472466
#> 47      5        0.5        V3 0.3900679
#> 48      5        0.5        V4 0.3879049
#> 49      5        0.5        V6 0.3804462
#> 50      5        0.5        V7 0.3595580
#> 51      6        0.6        V4 0.3873346
#> 52      6        0.6        V5 0.3608331
#> 53      6        0.6        V7 0.3590079
#> 54      6        0.6        V8 0.3531782
#> 55      6        0.6        V9 0.3890947
#> 56      6        0.6        V6 0.3798499
#> 57      6        0.6       V10 0.3996500
#> 58      6        0.6        V1 0.3899614
#> 59      6        0.6        V2 0.3467646
#> 60      6        0.6        V3 0.3894735
#> 61      7        0.7        V1 0.3894273
#> 62      7        0.7        V3 0.3888801
#> 63      7        0.7        V4 0.3867651
#> 64      7        0.7        V5 0.3602345
#> 65      7        0.7        V9 0.3885064
#> 66      7        0.7        V6 0.3792545
#> 67      7        0.7        V7 0.3584587
#> 68      7        0.7        V8 0.3526344
#> 69      7        0.7        V2 0.3462833
#> 70      7        0.7       V10 0.3990940
#> 71      8        0.8        V9 0.3879189
#> 72      8        0.8       V10 0.3985389
#> 73      8        0.8        V1 0.3888940
#> 74      8        0.8        V5 0.3596369
#> 75      8        0.8        V2 0.3458027
#> 76      8        0.8        V3 0.3882875
#> 77      8        0.8        V4 0.3861965
#> 78      8        0.8        V8 0.3520914
#> 79      8        0.8        V6 0.3786600
#> 80      8        0.8        V7 0.3579103
#> 81      9        0.9        V5 0.3590402
#> 82      9        0.9        V7 0.3573627
#> 83      9        0.9        V8 0.3515493
#> 84      9        0.9        V9 0.3873323
#> 85      9        0.9        V6 0.3780665
#> 86      9        0.9       V10 0.3979845
#> 87      9        0.9        V1 0.3883614
#> 88      9        0.9        V2 0.3453227
#> 89      9        0.9        V3 0.3876959
#> 90      9        0.9        V4 0.3856287
#> 91     10        1.0        V1 0.3878295
#> 92     10        1.0        V4 0.3850618
#> 93     10        1.0        V5 0.3584446
#> 94     10        1.0        V9 0.3867466
#> 95     10        1.0        V3 0.3871052
#> 96     10        1.0        V7 0.3568160
#> 97     10        1.0        V8 0.3510080
#> 98     10        1.0        V2 0.3448434
#> 99     10        1.0        V6 0.3774739
#> 100    10        1.0       V10 0.3974308


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
#> 1       0                0          0 0.8757906 0.04945237 0.7794039 0.9488849
#> 2       0               20         20 0.8757906 0.04945237 0.7794039 0.9488849
#> 3       0               40         40 0.8757906 0.04945237 0.7794039 0.9488849
#> 4       0               60         60 0.8757906 0.04945237 0.7794039 0.9488849
#> 5      20                0         20 0.8617131 0.05202644 0.7600923 0.9389068
#> 6      20               20         40 0.8617131 0.05202644 0.7600923 0.9389068
#> 7      20               40         60 0.8617131 0.05202644 0.7600923 0.9389068
#> 8      20               60         80 0.8617131 0.05202644 0.7600923 0.9389068
#> 9      40                0         40 0.8478591 0.05440967 0.7414893 0.9287542
#> 10     40               20         60 0.8478591 0.05440967 0.7414893 0.9287542
#> 11     40               40         80 0.8478591 0.05440967 0.7414893 0.9287542
#> 12     40               60        100 0.8478591 0.05440967 0.7414893 0.9287542
#> 13     60                0         60 0.8342249 0.05662432 0.7235250 0.9184801
#> 14     60               20         80 0.8342249 0.05662432 0.7235250 0.9184801
#> 15     60               40        100 0.8342249 0.05662432 0.7235250 0.9184801
#> 16     60               60        120 0.8342249 0.05662432 0.7235250 0.9184801
#> 17     80                0         80 0.8208071 0.05868821 0.7061443 0.9081253
#> 18     80               20        100 0.8208071 0.05868821 0.7061443 0.9081253
#> 19     80               40        120 0.8208071 0.05868821 0.7061443 0.9081253
#> 20     80               60        140 0.8208071 0.05868821 0.7061443 0.9081253
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13201655 0.165976862 0.6306567
#> 2  0.30574618 0.12221708 0.133080933 0.5557119
#> 3  0.26001915 0.11330065 0.105980572 0.4898732
#> 4  0.22113100 0.10512312 0.083719559 0.4323577
#> 5  0.25589195 0.11893239 0.088064382 0.4981518
#> 6  0.21762106 0.10926789 0.069057155 0.4395778
#> 7  0.18507391 0.10060718 0.053575767 0.3885450
#> 8  0.15739448 0.09276787 0.041050597 0.3441471
#> 9  0.18213629 0.10405834 0.043477559 0.3949467
#> 10 0.15489621 0.09527344 0.032941169 0.3497149
#> 11 0.13173012 0.08741575 0.024554679 0.3103737
#> 12 0.11202872 0.08032266 0.017959524 0.2761348
#> 13 0.12963921 0.08946422 0.019220633 0.3153081
#> 14 0.11025053 0.08177812 0.013819418 0.2804315
#> 15 0.09376159 0.07488952 0.009691428 0.2500420
#> 16 0.07973872 0.06866555 0.006601708 0.2235080
#> 17 0.09227334 0.07606638 0.007178991 0.2538583
#> 18 0.07847305 0.06947909 0.004763208 0.2268438
#> 19 0.06673672 0.06355439 0.003041458 0.2032014
#> 20 0.05675565 0.05818661 0.001857122 0.1824415
```
