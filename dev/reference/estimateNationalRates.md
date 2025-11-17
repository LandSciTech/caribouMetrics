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
[`simTrajectory()`](https://landscitech.github.io/caribouMetrics/dev/reference/simTrajectory.md),
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
#> 1         0.1 0.3838513 0.01269751 0.3563401 0.3923086
#> 2         0.2 0.3832760 0.01269390 0.3557687 0.3916549
#> 3         0.3 0.3827015 0.01269048 0.3551982 0.3910022
#> 4         0.4 0.3821279 0.01268724 0.3546286 0.3903655
#> 5         0.5 0.3815551 0.01268418 0.3540585 0.3897365
#> 6         0.6 0.3809832 0.01268130 0.3534684 0.3891085
#> 7         0.7 0.3804122 0.01267860 0.3528794 0.3884816
#> 8         0.8 0.3798420 0.01267608 0.3522913 0.3878557
#> 9         0.9 0.3792726 0.01267373 0.3517042 0.3872308
#> 10        1.0 0.3787041 0.01267155 0.3511180 0.3866069

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3765566
#> 2       1        0.1        V5 0.3907568
#> 3       1        0.1        V9 0.3709793
#> 4       1        0.1        V3 0.3909790
#> 5       1        0.1        V4 0.3781865
#> 6       1        0.1        V8 0.3891790
#> 7       1        0.1        V2 0.3713723
#> 8       1        0.1        V6 0.3926946
#> 9       1        0.1        V7 0.3520901
#> 10      1        0.1       V10 0.3739041
#> 11      2        0.2        V7 0.3515116
#> 12      2        0.2        V8 0.3886442
#> 13      2        0.2        V9 0.3704319
#> 14      2        0.2       V10 0.3733324
#> 15      2        0.2        V1 0.3760272
#> 16      2        0.2        V5 0.3902353
#> 17      2        0.2        V2 0.3707246
#> 18      2        0.2        V3 0.3903613
#> 19      2        0.2        V4 0.3776402
#> 20      2        0.2        V6 0.3920304
#> 21      3        0.3        V4 0.3770947
#> 22      3        0.3        V5 0.3897146
#> 23      3        0.3        V3 0.3897445
#> 24      3        0.3        V7 0.3509342
#> 25      3        0.3        V8 0.3881102
#> 26      3        0.3        V9 0.3698854
#> 27      3        0.3        V6 0.3913674
#> 28      3        0.3       V10 0.3727616
#> 29      3        0.3        V1 0.3754985
#> 30      3        0.3        V2 0.3700781
#> 31      4        0.4        V1 0.3749706
#> 32      4        0.4        V9 0.3693396
#> 33      4        0.4        V3 0.3891287
#> 34      4        0.4        V4 0.3765499
#> 35      4        0.4        V5 0.3891945
#> 36      4        0.4        V2 0.3694326
#> 37      4        0.4        V6 0.3907054
#> 38      4        0.4        V7 0.3503577
#> 39      4        0.4        V8 0.3875769
#> 40      4        0.4       V10 0.3721916
#> 41      5        0.5        V8 0.3870443
#> 42      5        0.5        V9 0.3687946
#> 43      5        0.5       V10 0.3716225
#> 44      5        0.5        V1 0.3744435
#> 45      5        0.5        V5 0.3886752
#> 46      5        0.5        V2 0.3687883
#> 47      5        0.5        V3 0.3885139
#> 48      5        0.5        V4 0.3760059
#> 49      5        0.5        V6 0.3900446
#> 50      5        0.5        V7 0.3497821
#> 51      6        0.6        V4 0.3754628
#> 52      6        0.6        V5 0.3881565
#> 53      6        0.6        V7 0.3492075
#> 54      6        0.6        V8 0.3865124
#> 55      6        0.6        V9 0.3682505
#> 56      6        0.6        V6 0.3893849
#> 57      6        0.6       V10 0.3710543
#> 58      6        0.6        V1 0.3739171
#> 59      6        0.6        V2 0.3681451
#> 60      6        0.6        V3 0.3879001
#> 61      7        0.7        V1 0.3733914
#> 62      7        0.7        V3 0.3872872
#> 63      7        0.7        V4 0.3749204
#> 64      7        0.7        V5 0.3876386
#> 65      7        0.7        V9 0.3677071
#> 66      7        0.7        V6 0.3887264
#> 67      7        0.7        V7 0.3486338
#> 68      7        0.7        V8 0.3859813
#> 69      7        0.7        V2 0.3675030
#> 70      7        0.7       V10 0.3704870
#> 71      8        0.8        V9 0.3671646
#> 72      8        0.8       V10 0.3699205
#> 73      8        0.8        V1 0.3728664
#> 74      8        0.8        V5 0.3871213
#> 75      8        0.8        V2 0.3668621
#> 76      8        0.8        V3 0.3866753
#> 77      8        0.8        V4 0.3743787
#> 78      8        0.8        V8 0.3854509
#> 79      8        0.8        V6 0.3880689
#> 80      8        0.8        V7 0.3480610
#> 81      9        0.9        V5 0.3866047
#> 82      9        0.9        V7 0.3474892
#> 83      9        0.9        V8 0.3849213
#> 84      9        0.9        V9 0.3666228
#> 85      9        0.9        V6 0.3874125
#> 86      9        0.9       V10 0.3693548
#> 87      9        0.9        V1 0.3723423
#> 88      9        0.9        V2 0.3662223
#> 89      9        0.9        V3 0.3860644
#> 90      9        0.9        V4 0.3738379
#> 91     10        1.0        V1 0.3718188
#> 92     10        1.0        V4 0.3732979
#> 93     10        1.0        V5 0.3860888
#> 94     10        1.0        V9 0.3660819
#> 95     10        1.0        V3 0.3854544
#> 96     10        1.0        V7 0.3469184
#> 97     10        1.0        V8 0.3843923
#> 98     10        1.0        V2 0.3655835
#> 99     10        1.0        V6 0.3867573
#> 100    10        1.0       V10 0.3687901

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3765566
#> 2       1        0.1        V5 0.3907568
#> 3       1        0.1        V9 0.3709793
#> 4       1        0.1        V3 0.3909790
#> 5       1        0.1        V4 0.3781865
#> 6       1        0.1        V8 0.3891790
#> 7       1        0.1        V2 0.3713723
#> 8       1        0.1        V6 0.3926946
#> 9       1        0.1        V7 0.3520901
#> 10      1        0.1       V10 0.3739041
#> 11      2        0.2        V7 0.3515116
#> 12      2        0.2        V8 0.3886442
#> 13      2        0.2        V9 0.3704319
#> 14      2        0.2       V10 0.3733324
#> 15      2        0.2        V1 0.3760272
#> 16      2        0.2        V5 0.3902353
#> 17      2        0.2        V2 0.3707246
#> 18      2        0.2        V3 0.3903613
#> 19      2        0.2        V4 0.3776402
#> 20      2        0.2        V6 0.3920304
#> 21      3        0.3        V4 0.3770947
#> 22      3        0.3        V5 0.3897146
#> 23      3        0.3        V3 0.3897445
#> 24      3        0.3        V7 0.3509342
#> 25      3        0.3        V8 0.3881102
#> 26      3        0.3        V9 0.3698854
#> 27      3        0.3        V6 0.3913674
#> 28      3        0.3       V10 0.3727616
#> 29      3        0.3        V1 0.3754985
#> 30      3        0.3        V2 0.3700781
#> 31      4        0.4        V1 0.3749706
#> 32      4        0.4        V9 0.3693396
#> 33      4        0.4        V3 0.3891287
#> 34      4        0.4        V4 0.3765499
#> 35      4        0.4        V5 0.3891945
#> 36      4        0.4        V2 0.3694326
#> 37      4        0.4        V6 0.3907054
#> 38      4        0.4        V7 0.3503577
#> 39      4        0.4        V8 0.3875769
#> 40      4        0.4       V10 0.3721916
#> 41      5        0.5        V8 0.3870443
#> 42      5        0.5        V9 0.3687946
#> 43      5        0.5       V10 0.3716225
#> 44      5        0.5        V1 0.3744435
#> 45      5        0.5        V5 0.3886752
#> 46      5        0.5        V2 0.3687883
#> 47      5        0.5        V3 0.3885139
#> 48      5        0.5        V4 0.3760059
#> 49      5        0.5        V6 0.3900446
#> 50      5        0.5        V7 0.3497821
#> 51      6        0.6        V4 0.3754628
#> 52      6        0.6        V5 0.3881565
#> 53      6        0.6        V7 0.3492075
#> 54      6        0.6        V8 0.3865124
#> 55      6        0.6        V9 0.3682505
#> 56      6        0.6        V6 0.3893849
#> 57      6        0.6       V10 0.3710543
#> 58      6        0.6        V1 0.3739171
#> 59      6        0.6        V2 0.3681451
#> 60      6        0.6        V3 0.3879001
#> 61      7        0.7        V1 0.3733914
#> 62      7        0.7        V3 0.3872872
#> 63      7        0.7        V4 0.3749204
#> 64      7        0.7        V5 0.3876386
#> 65      7        0.7        V9 0.3677071
#> 66      7        0.7        V6 0.3887264
#> 67      7        0.7        V7 0.3486338
#> 68      7        0.7        V8 0.3859813
#> 69      7        0.7        V2 0.3675030
#> 70      7        0.7       V10 0.3704870
#> 71      8        0.8        V9 0.3671646
#> 72      8        0.8       V10 0.3699205
#> 73      8        0.8        V1 0.3728664
#> 74      8        0.8        V5 0.3871213
#> 75      8        0.8        V2 0.3668621
#> 76      8        0.8        V3 0.3866753
#> 77      8        0.8        V4 0.3743787
#> 78      8        0.8        V8 0.3854509
#> 79      8        0.8        V6 0.3880689
#> 80      8        0.8        V7 0.3480610
#> 81      9        0.9        V5 0.3866047
#> 82      9        0.9        V7 0.3474892
#> 83      9        0.9        V8 0.3849213
#> 84      9        0.9        V9 0.3666228
#> 85      9        0.9        V6 0.3874125
#> 86      9        0.9       V10 0.3693548
#> 87      9        0.9        V1 0.3723423
#> 88      9        0.9        V2 0.3662223
#> 89      9        0.9        V3 0.3860644
#> 90      9        0.9        V4 0.3738379
#> 91     10        1.0        V1 0.3718188
#> 92     10        1.0        V4 0.3732979
#> 93     10        1.0        V5 0.3860888
#> 94     10        1.0        V9 0.3660819
#> 95     10        1.0        V3 0.3854544
#> 96     10        1.0        V7 0.3469184
#> 97     10        1.0        V8 0.3843923
#> 98     10        1.0        V2 0.3655835
#> 99     10        1.0        V6 0.3867573
#> 100    10        1.0       V10 0.3687901


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
#> 1       0                0          0 0.8757906 0.04408131 0.7948010 0.9409652
#> 2       0               20         20 0.8757906 0.04408131 0.7948010 0.9409652
#> 3       0               40         40 0.8757906 0.04408131 0.7948010 0.9409652
#> 4       0               60         60 0.8757906 0.04408131 0.7948010 0.9409652
#> 5      20                0         20 0.8617131 0.04656881 0.7751327 0.9305502
#> 6      20               20         40 0.8617131 0.04656881 0.7751327 0.9305502
#> 7      20               40         60 0.8617131 0.04656881 0.7751327 0.9305502
#> 8      20               60         80 0.8617131 0.04656881 0.7751327 0.9305502
#> 9      40                0         40 0.8478591 0.04893804 0.7561703 0.9200039
#> 10     40               20         60 0.8478591 0.04893804 0.7561703 0.9200039
#> 11     40               40         80 0.8478591 0.04893804 0.7561703 0.9200039
#> 12     40               60        100 0.8478591 0.04893804 0.7561703 0.9200039
#> 13     60                0         60 0.8342249 0.05119994 0.7378450 0.9093715
#> 14     60               20         80 0.8342249 0.05119994 0.7378450 0.9093715
#> 15     60               40        100 0.8342249 0.05119994 0.7378450 0.9093715
#> 16     60               60        120 0.8342249 0.05119994 0.7378450 0.9093715
#> 17     80                0         80 0.8208071 0.05336292 0.7201032 0.8986883
#> 18     80               20        100 0.8208071 0.05336292 0.7201032 0.8986883
#> 19     80               40        120 0.8208071 0.05336292 0.7201032 0.8986883
#> 20     80               60        140 0.8208071 0.05336292 0.7201032 0.8986883
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11326886 0.172688381 0.5646447
#> 2  0.30574618 0.10504952 0.132716458 0.5042689
#> 3  0.26001915 0.09826318 0.101033423 0.4507460
#> 4  0.22113100 0.09234200 0.076027277 0.4034331
#> 5  0.25589195 0.10568712 0.094815677 0.4686054
#> 6  0.21762106 0.09720830 0.071137825 0.4192100
#> 7  0.18507391 0.08994241 0.052599792 0.3755922
#> 8  0.15739448 0.08353065 0.038217393 0.3371047
#> 9  0.18213629 0.09589652 0.048998165 0.3901338
#> 10 0.15489621 0.08785650 0.035442691 0.3499349
#> 11 0.13173012 0.08084123 0.025081758 0.3144650
#> 12 0.11202872 0.07461021 0.017288037 0.2831503
#> 13 0.12963921 0.08552265 0.023105556 0.3262902
#> 14 0.11025053 0.07822533 0.015819499 0.2935935
#> 15 0.09376159 0.07179278 0.010477699 0.2647068
#> 16 0.07973872 0.06605553 0.006666802 0.2391469
#> 17 0.09227334 0.07544369 0.009490737 0.2743432
#> 18 0.07847305 0.06899045 0.005977300 0.2476788
#> 19 0.06673672 0.06326492 0.003578117 0.2240546
#> 20 0.05675565 0.05814191 0.002014015 0.2030735
```
