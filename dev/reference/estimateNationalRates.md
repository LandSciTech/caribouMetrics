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
#> 1         0.1 0.3838513 0.02066132 0.3630189 0.4163998
#> 2         0.2 0.3832760 0.02064758 0.3624795 0.4158327
#> 3         0.3 0.3827015 0.02063399 0.3619410 0.4152664
#> 4         0.4 0.3821279 0.02062053 0.3614032 0.4147009
#> 5         0.5 0.3815551 0.02060722 0.3608663 0.4141362
#> 6         0.6 0.3809832 0.02059405 0.3603301 0.4135722
#> 7         0.7 0.3804122 0.02058102 0.3597948 0.4130090
#> 8         0.8 0.3798420 0.02056812 0.3592602 0.4124465
#> 9         0.9 0.3792726 0.02055537 0.3587264 0.4118848
#> 10        1.0 0.3787041 0.02054275 0.3581934 0.4113239

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4055959
#> 2       1        0.1        V5 0.4107628
#> 3       1        0.1        V9 0.3940517
#> 4       1        0.1        V3 0.3812394
#> 5       1        0.1        V4 0.3672101
#> 6       1        0.1        V8 0.3626400
#> 7       1        0.1        V2 0.4180364
#> 8       1        0.1        V6 0.3727637
#> 9       1        0.1        V7 0.3643239
#> 10      1        0.1       V10 0.3718402
#> 11      2        0.2        V7 0.3638155
#> 12      2        0.2        V8 0.3620917
#> 13      2        0.2        V9 0.3933747
#> 14      2        0.2       V10 0.3712703
#> 15      2        0.2        V1 0.4049995
#> 16      2        0.2        V5 0.4102080
#> 17      2        0.2        V2 0.4174657
#> 18      2        0.2        V3 0.3806846
#> 19      2        0.2        V4 0.3665625
#> 20      2        0.2        V6 0.3722745
#> 21      3        0.3        V4 0.3659159
#> 22      3        0.3        V5 0.4096540
#> 23      3        0.3        V3 0.3801306
#> 24      3        0.3        V7 0.3633079
#> 25      3        0.3        V8 0.3615442
#> 26      3        0.3        V9 0.3926988
#> 27      3        0.3        V6 0.3717859
#> 28      3        0.3       V10 0.3707012
#> 29      3        0.3        V1 0.4044039
#> 30      3        0.3        V2 0.4168959
#> 31      4        0.4        V1 0.4038092
#> 32      4        0.4        V9 0.3920241
#> 33      4        0.4        V3 0.3795774
#> 34      4        0.4        V4 0.3652706
#> 35      4        0.4        V5 0.4091007
#> 36      4        0.4        V2 0.4163268
#> 37      4        0.4        V6 0.3712980
#> 38      4        0.4        V7 0.3628009
#> 39      4        0.4        V8 0.3609975
#> 40      4        0.4       V10 0.3701331
#> 41      5        0.5        V8 0.3604516
#> 42      5        0.5        V9 0.3913505
#> 43      5        0.5       V10 0.3695658
#> 44      5        0.5        V1 0.4032154
#> 45      5        0.5        V5 0.4085482
#> 46      5        0.5        V2 0.4157585
#> 47      5        0.5        V3 0.3790251
#> 48      5        0.5        V4 0.3646263
#> 49      5        0.5        V6 0.3708107
#> 50      5        0.5        V7 0.3622946
#> 51      6        0.6        V4 0.3639832
#> 52      6        0.6        V5 0.4079964
#> 53      6        0.6        V7 0.3617891
#> 54      6        0.6        V8 0.3599065
#> 55      6        0.6        V9 0.3906781
#> 56      6        0.6        V6 0.3703240
#> 57      6        0.6       V10 0.3689993
#> 58      6        0.6        V1 0.4026224
#> 59      6        0.6        V2 0.4151910
#> 60      6        0.6        V3 0.3784735
#> 61      7        0.7        V1 0.4020303
#> 62      7        0.7        V3 0.3779228
#> 63      7        0.7        V4 0.3633413
#> 64      7        0.7        V5 0.4074454
#> 65      7        0.7        V9 0.3900069
#> 66      7        0.7        V6 0.3698380
#> 67      7        0.7        V7 0.3612842
#> 68      7        0.7        V8 0.3593623
#> 69      7        0.7        V2 0.4146242
#> 70      7        0.7       V10 0.3684338
#> 71      8        0.8        V9 0.3893368
#> 72      8        0.8       V10 0.3678691
#> 73      8        0.8        V1 0.4014391
#> 74      8        0.8        V5 0.4068951
#> 75      8        0.8        V2 0.4140582
#> 76      8        0.8        V3 0.3773728
#> 77      8        0.8        V4 0.3627004
#> 78      8        0.8        V8 0.3588189
#> 79      8        0.8        V6 0.3693527
#> 80      8        0.8        V7 0.3607801
#> 81      9        0.9        V5 0.4063456
#> 82      9        0.9        V7 0.3602767
#> 83      9        0.9        V8 0.3582764
#> 84      9        0.9        V9 0.3886678
#> 85      9        0.9        V6 0.3688679
#> 86      9        0.9       V10 0.3673053
#> 87      9        0.9        V1 0.4008488
#> 88      9        0.9        V2 0.4134930
#> 89      9        0.9        V3 0.3768237
#> 90      9        0.9        V4 0.3620607
#> 91     10        1.0        V1 0.4002593
#> 92     10        1.0        V4 0.3614221
#> 93     10        1.0        V5 0.4057968
#> 94     10        1.0        V9 0.3880000
#> 95     10        1.0        V3 0.3762753
#> 96     10        1.0        V7 0.3597739
#> 97     10        1.0        V8 0.3577346
#> 98     10        1.0        V2 0.4129286
#> 99     10        1.0        V6 0.3683838
#> 100    10        1.0       V10 0.3667423

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4055959
#> 2       1        0.1        V5 0.4107628
#> 3       1        0.1        V9 0.3940517
#> 4       1        0.1        V3 0.3812394
#> 5       1        0.1        V4 0.3672101
#> 6       1        0.1        V8 0.3626400
#> 7       1        0.1        V2 0.4180364
#> 8       1        0.1        V6 0.3727637
#> 9       1        0.1        V7 0.3643239
#> 10      1        0.1       V10 0.3718402
#> 11      2        0.2        V7 0.3638155
#> 12      2        0.2        V8 0.3620917
#> 13      2        0.2        V9 0.3933747
#> 14      2        0.2       V10 0.3712703
#> 15      2        0.2        V1 0.4049995
#> 16      2        0.2        V5 0.4102080
#> 17      2        0.2        V2 0.4174657
#> 18      2        0.2        V3 0.3806846
#> 19      2        0.2        V4 0.3665625
#> 20      2        0.2        V6 0.3722745
#> 21      3        0.3        V4 0.3659159
#> 22      3        0.3        V5 0.4096540
#> 23      3        0.3        V3 0.3801306
#> 24      3        0.3        V7 0.3633079
#> 25      3        0.3        V8 0.3615442
#> 26      3        0.3        V9 0.3926988
#> 27      3        0.3        V6 0.3717859
#> 28      3        0.3       V10 0.3707012
#> 29      3        0.3        V1 0.4044039
#> 30      3        0.3        V2 0.4168959
#> 31      4        0.4        V1 0.4038092
#> 32      4        0.4        V9 0.3920241
#> 33      4        0.4        V3 0.3795774
#> 34      4        0.4        V4 0.3652706
#> 35      4        0.4        V5 0.4091007
#> 36      4        0.4        V2 0.4163268
#> 37      4        0.4        V6 0.3712980
#> 38      4        0.4        V7 0.3628009
#> 39      4        0.4        V8 0.3609975
#> 40      4        0.4       V10 0.3701331
#> 41      5        0.5        V8 0.3604516
#> 42      5        0.5        V9 0.3913505
#> 43      5        0.5       V10 0.3695658
#> 44      5        0.5        V1 0.4032154
#> 45      5        0.5        V5 0.4085482
#> 46      5        0.5        V2 0.4157585
#> 47      5        0.5        V3 0.3790251
#> 48      5        0.5        V4 0.3646263
#> 49      5        0.5        V6 0.3708107
#> 50      5        0.5        V7 0.3622946
#> 51      6        0.6        V4 0.3639832
#> 52      6        0.6        V5 0.4079964
#> 53      6        0.6        V7 0.3617891
#> 54      6        0.6        V8 0.3599065
#> 55      6        0.6        V9 0.3906781
#> 56      6        0.6        V6 0.3703240
#> 57      6        0.6       V10 0.3689993
#> 58      6        0.6        V1 0.4026224
#> 59      6        0.6        V2 0.4151910
#> 60      6        0.6        V3 0.3784735
#> 61      7        0.7        V1 0.4020303
#> 62      7        0.7        V3 0.3779228
#> 63      7        0.7        V4 0.3633413
#> 64      7        0.7        V5 0.4074454
#> 65      7        0.7        V9 0.3900069
#> 66      7        0.7        V6 0.3698380
#> 67      7        0.7        V7 0.3612842
#> 68      7        0.7        V8 0.3593623
#> 69      7        0.7        V2 0.4146242
#> 70      7        0.7       V10 0.3684338
#> 71      8        0.8        V9 0.3893368
#> 72      8        0.8       V10 0.3678691
#> 73      8        0.8        V1 0.4014391
#> 74      8        0.8        V5 0.4068951
#> 75      8        0.8        V2 0.4140582
#> 76      8        0.8        V3 0.3773728
#> 77      8        0.8        V4 0.3627004
#> 78      8        0.8        V8 0.3588189
#> 79      8        0.8        V6 0.3693527
#> 80      8        0.8        V7 0.3607801
#> 81      9        0.9        V5 0.4063456
#> 82      9        0.9        V7 0.3602767
#> 83      9        0.9        V8 0.3582764
#> 84      9        0.9        V9 0.3886678
#> 85      9        0.9        V6 0.3688679
#> 86      9        0.9       V10 0.3673053
#> 87      9        0.9        V1 0.4008488
#> 88      9        0.9        V2 0.4134930
#> 89      9        0.9        V3 0.3768237
#> 90      9        0.9        V4 0.3620607
#> 91     10        1.0        V1 0.4002593
#> 92     10        1.0        V4 0.3614221
#> 93     10        1.0        V5 0.4057968
#> 94     10        1.0        V9 0.3880000
#> 95     10        1.0        V3 0.3762753
#> 96     10        1.0        V7 0.3597739
#> 97     10        1.0        V8 0.3577346
#> 98     10        1.0        V2 0.4129286
#> 99     10        1.0        V6 0.3683838
#> 100    10        1.0       V10 0.3667423


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
#> 1       0                0          0 0.8757906 0.04105968 0.8022922 0.9378847
#> 2       0               20         20 0.8757906 0.04105968 0.8022922 0.9378847
#> 3       0               40         40 0.8757906 0.04105968 0.8022922 0.9378847
#> 4       0               60         60 0.8757906 0.04105968 0.8022922 0.9378847
#> 5      20                0         20 0.8617131 0.04260535 0.7864968 0.9288917
#> 6      20               20         40 0.8617131 0.04260535 0.7864968 0.9288917
#> 7      20               40         60 0.8617131 0.04260535 0.7864968 0.9288917
#> 8      20               60         80 0.8617131 0.04260535 0.7864968 0.9288917
#> 9      40                0         40 0.8478591 0.04402573 0.7711843 0.9198106
#> 10     40               20         60 0.8478591 0.04402573 0.7711843 0.9198106
#> 11     40               40         80 0.8478591 0.04402573 0.7711843 0.9198106
#> 12     40               60        100 0.8478591 0.04402573 0.7711843 0.9198106
#> 13     60                0         60 0.8342249 0.04534002 0.7563108 0.9106690
#> 14     60               20         80 0.8342249 0.04534002 0.7563108 0.9106690
#> 15     60               40        100 0.8342249 0.04534002 0.7563108 0.9106690
#> 16     60               60        120 0.8342249 0.04534002 0.7563108 0.9106690
#> 17     80                0         80 0.8208071 0.04656263 0.7418415 0.9014894
#> 18     80               20        100 0.8208071 0.04656263 0.7418415 0.9014894
#> 19     80               40        120 0.8208071 0.04656263 0.7418415 0.9014894
#> 20     80               60        140 0.8208071 0.04656263 0.7418415 0.9014894
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12288387 0.188226179 0.5911356
#> 2  0.30574618 0.11635698 0.139300338 0.5208247
#> 3  0.26001915 0.10965162 0.101860639 0.4593274
#> 4  0.22113100 0.10283637 0.073362490 0.4057611
#> 5  0.25589195 0.10846091 0.109816051 0.4648663
#> 6  0.21762106 0.10126970 0.079400498 0.4105800
#> 7  0.18507391 0.09432958 0.056390304 0.3633783
#> 8  0.15739448 0.08761855 0.039172047 0.3223557
#> 9  0.18213629 0.09422499 0.061250604 0.3676232
#> 10 0.15489621 0.08726056 0.042789598 0.3260454
#> 11 0.13173012 0.08068798 0.029125901 0.2898882
#> 12 0.11202872 0.07446578 0.019198274 0.2583932
#> 13 0.12963921 0.08099771 0.031981524 0.2931421
#> 14 0.11025053 0.07461996 0.021254682 0.2612305
#> 15 0.09376159 0.06865828 0.013600470 0.2333714
#> 16 0.07973872 0.06307179 0.008301699 0.2089702
#> 17 0.09227334 0.06909965 0.015172184 0.2358843
#> 18 0.07847305 0.06341993 0.009373847 0.2111751
#> 19 0.06673672 0.05813313 0.005476623 0.1894538
#> 20 0.05675565 0.05320389 0.002983540 0.1702680
```
