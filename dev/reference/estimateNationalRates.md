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
#> 1         0.1 0.3838513 0.02260511 0.3556453 0.4193930
#> 2         0.2 0.3832760 0.02260739 0.3550414 0.4187706
#> 3         0.3 0.3827015 0.02260972 0.3544386 0.4181491
#> 4         0.4 0.3821279 0.02261210 0.3538367 0.4175286
#> 5         0.5 0.3815551 0.02261452 0.3532359 0.4169089
#> 6         0.6 0.3809832 0.02261699 0.3526361 0.4162902
#> 7         0.7 0.3804122 0.02261949 0.3520374 0.4156724
#> 8         0.8 0.3798420 0.02262204 0.3514396 0.4150556
#> 9         0.9 0.3792726 0.02262463 0.3508429 0.4144396
#> 10        1.0 0.3787041 0.02262727 0.3502471 0.4138246

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3546157
#> 2       1        0.1        V5 0.3832138
#> 3       1        0.1        V9 0.3635982
#> 4       1        0.1        V3 0.4221020
#> 5       1        0.1        V4 0.4100619
#> 6       1        0.1        V8 0.3920838
#> 7       1        0.1        V2 0.3591919
#> 8       1        0.1        V6 0.3992131
#> 9       1        0.1        V7 0.4035323
#> 10      1        0.1       V10 0.3851688
#> 11      2        0.2        V7 0.4030061
#> 12      2        0.2        V8 0.3915224
#> 13      2        0.2        V9 0.3630568
#> 14      2        0.2       V10 0.3844928
#> 15      2        0.2        V1 0.3540236
#> 16      2        0.2        V5 0.3826725
#> 17      2        0.2        V2 0.3585472
#> 18      2        0.2        V3 0.4214684
#> 19      2        0.2        V4 0.4094782
#> 20      2        0.2        V6 0.3986417
#> 21      3        0.3        V4 0.4088953
#> 22      3        0.3        V5 0.3821320
#> 23      3        0.3        V3 0.4208357
#> 24      3        0.3        V7 0.4024806
#> 25      3        0.3        V8 0.3909618
#> 26      3        0.3        V9 0.3625163
#> 27      3        0.3        V6 0.3980710
#> 28      3        0.3       V10 0.3838181
#> 29      3        0.3        V1 0.3534326
#> 30      3        0.3        V2 0.3579037
#> 31      4        0.4        V1 0.3528425
#> 32      4        0.4        V9 0.3619765
#> 33      4        0.4        V3 0.4202040
#> 34      4        0.4        V4 0.4083132
#> 35      4        0.4        V5 0.3815922
#> 36      4        0.4        V2 0.3572614
#> 37      4        0.4        V6 0.3975012
#> 38      4        0.4        V7 0.4019558
#> 39      4        0.4        V8 0.3904020
#> 40      4        0.4       V10 0.3831446
#> 41      5        0.5        V8 0.3898431
#> 42      5        0.5        V9 0.3614376
#> 43      5        0.5       V10 0.3824722
#> 44      5        0.5        V1 0.3522534
#> 45      5        0.5        V5 0.3810532
#> 46      5        0.5        V2 0.3566202
#> 47      5        0.5        V3 0.4195732
#> 48      5        0.5        V4 0.4077320
#> 49      5        0.5        V6 0.3969322
#> 50      5        0.5        V7 0.4014317
#> 51      6        0.6        V4 0.4071516
#> 52      6        0.6        V5 0.3805149
#> 53      6        0.6        V7 0.4009083
#> 54      6        0.6        V8 0.3892849
#> 55      6        0.6        V9 0.3608995
#> 56      6        0.6        V6 0.3963640
#> 57      6        0.6       V10 0.3818011
#> 58      6        0.6        V1 0.3516653
#> 59      6        0.6        V2 0.3559801
#> 60      6        0.6        V3 0.4189434
#> 61      7        0.7        V1 0.3510782
#> 62      7        0.7        V3 0.4183145
#> 63      7        0.7        V4 0.4065720
#> 64      7        0.7        V5 0.3799775
#> 65      7        0.7        V9 0.3603621
#> 66      7        0.7        V6 0.3957966
#> 67      7        0.7        V7 0.4003855
#> 68      7        0.7        V8 0.3887275
#> 69      7        0.7        V2 0.3553412
#> 70      7        0.7       V10 0.3811311
#> 71      8        0.8        V9 0.3598256
#> 72      8        0.8       V10 0.3804623
#> 73      8        0.8        V1 0.3504920
#> 74      8        0.8        V5 0.3794407
#> 75      8        0.8        V2 0.3547034
#> 76      8        0.8        V3 0.4176866
#> 77      8        0.8        V4 0.4059933
#> 78      8        0.8        V8 0.3881710
#> 79      8        0.8        V6 0.3952301
#> 80      8        0.8        V7 0.3998635
#> 81      9        0.9        V5 0.3789048
#> 82      9        0.9        V7 0.3993421
#> 83      9        0.9        V8 0.3876152
#> 84      9        0.9        V9 0.3592899
#> 85      9        0.9        V6 0.3946643
#> 86      9        0.9       V10 0.3797946
#> 87      9        0.9        V1 0.3499069
#> 88      9        0.9        V2 0.3540668
#> 89      9        0.9        V3 0.4170596
#> 90      9        0.9        V4 0.4054153
#> 91     10        1.0        V1 0.3493227
#> 92     10        1.0        V4 0.4048382
#> 93     10        1.0        V5 0.3783695
#> 94     10        1.0        V9 0.3587550
#> 95     10        1.0        V3 0.4164335
#> 96     10        1.0        V7 0.3988214
#> 97     10        1.0        V8 0.3870602
#> 98     10        1.0        V2 0.3534313
#> 99     10        1.0        V6 0.3940994
#> 100    10        1.0       V10 0.3791281

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3546157
#> 2       1        0.1        V5 0.3832138
#> 3       1        0.1        V9 0.3635982
#> 4       1        0.1        V3 0.4221020
#> 5       1        0.1        V4 0.4100619
#> 6       1        0.1        V8 0.3920838
#> 7       1        0.1        V2 0.3591919
#> 8       1        0.1        V6 0.3992131
#> 9       1        0.1        V7 0.4035323
#> 10      1        0.1       V10 0.3851688
#> 11      2        0.2        V7 0.4030061
#> 12      2        0.2        V8 0.3915224
#> 13      2        0.2        V9 0.3630568
#> 14      2        0.2       V10 0.3844928
#> 15      2        0.2        V1 0.3540236
#> 16      2        0.2        V5 0.3826725
#> 17      2        0.2        V2 0.3585472
#> 18      2        0.2        V3 0.4214684
#> 19      2        0.2        V4 0.4094782
#> 20      2        0.2        V6 0.3986417
#> 21      3        0.3        V4 0.4088953
#> 22      3        0.3        V5 0.3821320
#> 23      3        0.3        V3 0.4208357
#> 24      3        0.3        V7 0.4024806
#> 25      3        0.3        V8 0.3909618
#> 26      3        0.3        V9 0.3625163
#> 27      3        0.3        V6 0.3980710
#> 28      3        0.3       V10 0.3838181
#> 29      3        0.3        V1 0.3534326
#> 30      3        0.3        V2 0.3579037
#> 31      4        0.4        V1 0.3528425
#> 32      4        0.4        V9 0.3619765
#> 33      4        0.4        V3 0.4202040
#> 34      4        0.4        V4 0.4083132
#> 35      4        0.4        V5 0.3815922
#> 36      4        0.4        V2 0.3572614
#> 37      4        0.4        V6 0.3975012
#> 38      4        0.4        V7 0.4019558
#> 39      4        0.4        V8 0.3904020
#> 40      4        0.4       V10 0.3831446
#> 41      5        0.5        V8 0.3898431
#> 42      5        0.5        V9 0.3614376
#> 43      5        0.5       V10 0.3824722
#> 44      5        0.5        V1 0.3522534
#> 45      5        0.5        V5 0.3810532
#> 46      5        0.5        V2 0.3566202
#> 47      5        0.5        V3 0.4195732
#> 48      5        0.5        V4 0.4077320
#> 49      5        0.5        V6 0.3969322
#> 50      5        0.5        V7 0.4014317
#> 51      6        0.6        V4 0.4071516
#> 52      6        0.6        V5 0.3805149
#> 53      6        0.6        V7 0.4009083
#> 54      6        0.6        V8 0.3892849
#> 55      6        0.6        V9 0.3608995
#> 56      6        0.6        V6 0.3963640
#> 57      6        0.6       V10 0.3818011
#> 58      6        0.6        V1 0.3516653
#> 59      6        0.6        V2 0.3559801
#> 60      6        0.6        V3 0.4189434
#> 61      7        0.7        V1 0.3510782
#> 62      7        0.7        V3 0.4183145
#> 63      7        0.7        V4 0.4065720
#> 64      7        0.7        V5 0.3799775
#> 65      7        0.7        V9 0.3603621
#> 66      7        0.7        V6 0.3957966
#> 67      7        0.7        V7 0.4003855
#> 68      7        0.7        V8 0.3887275
#> 69      7        0.7        V2 0.3553412
#> 70      7        0.7       V10 0.3811311
#> 71      8        0.8        V9 0.3598256
#> 72      8        0.8       V10 0.3804623
#> 73      8        0.8        V1 0.3504920
#> 74      8        0.8        V5 0.3794407
#> 75      8        0.8        V2 0.3547034
#> 76      8        0.8        V3 0.4176866
#> 77      8        0.8        V4 0.4059933
#> 78      8        0.8        V8 0.3881710
#> 79      8        0.8        V6 0.3952301
#> 80      8        0.8        V7 0.3998635
#> 81      9        0.9        V5 0.3789048
#> 82      9        0.9        V7 0.3993421
#> 83      9        0.9        V8 0.3876152
#> 84      9        0.9        V9 0.3592899
#> 85      9        0.9        V6 0.3946643
#> 86      9        0.9       V10 0.3797946
#> 87      9        0.9        V1 0.3499069
#> 88      9        0.9        V2 0.3540668
#> 89      9        0.9        V3 0.4170596
#> 90      9        0.9        V4 0.4054153
#> 91     10        1.0        V1 0.3493227
#> 92     10        1.0        V4 0.4048382
#> 93     10        1.0        V5 0.3783695
#> 94     10        1.0        V9 0.3587550
#> 95     10        1.0        V3 0.4164335
#> 96     10        1.0        V7 0.3988214
#> 97     10        1.0        V8 0.3870602
#> 98     10        1.0        V2 0.3534313
#> 99     10        1.0        V6 0.3940994
#> 100    10        1.0       V10 0.3791281


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
#> 1       0                0          0 0.8757906 0.05426859 0.7519684 0.9378056
#> 2       0               20         20 0.8757906 0.05426859 0.7519684 0.9378056
#> 3       0               40         40 0.8757906 0.05426859 0.7519684 0.9378056
#> 4       0               60         60 0.8757906 0.05426859 0.7519684 0.9378056
#> 5      20                0         20 0.8617131 0.05583377 0.7351913 0.9260388
#> 6      20               20         40 0.8617131 0.05583377 0.7351913 0.9260388
#> 7      20               40         60 0.8617131 0.05583377 0.7351913 0.9260388
#> 8      20               60         80 0.8617131 0.05583377 0.7351913 0.9260388
#> 9      40                0         40 0.8478591 0.05730703 0.7189543 0.9141195
#> 10     40               20         60 0.8478591 0.05730703 0.7189543 0.9141195
#> 11     40               40         80 0.8478591 0.05730703 0.7189543 0.9141195
#> 12     40               60        100 0.8478591 0.05730703 0.7189543 0.9141195
#> 13     60                0         60 0.8342249 0.05869990 0.7032127 0.9021065
#> 14     60               20         80 0.8342249 0.05869990 0.7032127 0.9021065
#> 15     60               40        100 0.8342249 0.05869990 0.7032127 0.9021065
#> 16     60               60        120 0.8342249 0.05869990 0.7032127 0.9021065
#> 17     80                0         80 0.8208071 0.06002103 0.6879298 0.8900449
#> 18     80               20        100 0.8208071 0.06002103 0.6879298 0.8900449
#> 19     80               40        120 0.8208071 0.06002103 0.6879298 0.8900449
#> 20     80               60        140 0.8208071 0.06002103 0.6879298 0.8900449
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.10556748 0.178840189 0.5465219
#> 2  0.30574618 0.09942157 0.133128720 0.4887196
#> 3  0.26001915 0.09386875 0.097912304 0.4374563
#> 4  0.22113100 0.08862287 0.070931239 0.3920992
#> 5  0.25589195 0.09735836 0.101276819 0.4443288
#> 6  0.21762106 0.09012759 0.073499977 0.3981758
#> 7  0.18507391 0.08373755 0.052377354 0.3573827
#> 8  0.15739448 0.07794055 0.036498712 0.3213405
#> 9  0.18213629 0.08799821 0.054379569 0.3628467
#> 10 0.15489621 0.08079731 0.037994039 0.3261683
#> 11 0.13173012 0.07444145 0.025842269 0.2937532
#> 12 0.11202872 0.06873769 0.017006542 0.2650811
#> 13 0.12963921 0.07828709 0.026977936 0.2980963
#> 14 0.11025053 0.07161152 0.017823188 0.2689247
#> 15 0.09376159 0.06570377 0.011316258 0.2430915
#> 16 0.07973872 0.06041604 0.006840727 0.2201728
#> 17 0.09227334 0.06882880 0.011910003 0.2465568
#> 18 0.07847305 0.06286650 0.007241575 0.2232500
#> 19 0.06673672 0.05757199 0.004149250 0.2025323
#> 20 0.05675565 0.05283067 0.002208896 0.1840660
```
