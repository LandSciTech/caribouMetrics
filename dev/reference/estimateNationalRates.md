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
#> 1         0.1 0.3838513 0.02543683 0.3607070 0.4323489
#> 2         0.2 0.3832760 0.02543606 0.3601190 0.4317187
#> 3         0.3 0.3827015 0.02543533 0.3595319 0.4310894
#> 4         0.4 0.3821279 0.02543464 0.3589458 0.4304610
#> 5         0.5 0.3815551 0.02543400 0.3583606 0.4298336
#> 6         0.6 0.3809832 0.02543339 0.3577764 0.4292070
#> 7         0.7 0.3804122 0.02543283 0.3571931 0.4285814
#> 8         0.8 0.3798420 0.02543231 0.3566108 0.4279567
#> 9         0.9 0.3792726 0.02543183 0.3560295 0.4273329
#> 10        1.0 0.3787041 0.02543138 0.3554491 0.4267101

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3943932
#> 2       1        0.1        V5 0.3614875
#> 3       1        0.1        V9 0.3916625
#> 4       1        0.1        V3 0.3755326
#> 5       1        0.1        V4 0.3604804
#> 6       1        0.1        V8 0.3746962
#> 7       1        0.1        V2 0.3705368
#> 8       1        0.1        V6 0.4367345
#> 9       1        0.1        V7 0.4172429
#> 10      1        0.1       V10 0.3645971
#> 11      2        0.2        V7 0.4166840
#> 12      2        0.2        V8 0.3740474
#> 13      2        0.2        V9 0.3911834
#> 14      2        0.2       V10 0.3639912
#> 15      2        0.2        V1 0.3937931
#> 16      2        0.2        V5 0.3608575
#> 17      2        0.2        V2 0.3699458
#> 18      2        0.2        V3 0.3749662
#> 19      2        0.2        V4 0.3599046
#> 20      2        0.2        V6 0.4360835
#> 21      3        0.3        V4 0.3593297
#> 22      3        0.3        V5 0.3602285
#> 23      3        0.3        V3 0.3744007
#> 24      3        0.3        V7 0.4161259
#> 25      3        0.3        V8 0.3733996
#> 26      3        0.3        V9 0.3907049
#> 27      3        0.3        V6 0.4354336
#> 28      3        0.3       V10 0.3633862
#> 29      3        0.3        V1 0.3931938
#> 30      3        0.3        V2 0.3693559
#> 31      4        0.4        V1 0.3925954
#> 32      4        0.4        V9 0.3902270
#> 33      4        0.4        V3 0.3738360
#> 34      4        0.4        V4 0.3587557
#> 35      4        0.4        V5 0.3596006
#> 36      4        0.4        V2 0.3687668
#> 37      4        0.4        V6 0.4347846
#> 38      4        0.4        V7 0.4155686
#> 39      4        0.4        V8 0.3727530
#> 40      4        0.4       V10 0.3627823
#> 41      5        0.5        V8 0.3721075
#> 42      5        0.5        V9 0.3897496
#> 43      5        0.5       V10 0.3621793
#> 44      5        0.5        V1 0.3919980
#> 45      5        0.5        V5 0.3589738
#> 46      5        0.5        V2 0.3681787
#> 47      5        0.5        V3 0.3732721
#> 48      5        0.5        V4 0.3581826
#> 49      5        0.5        V6 0.4341366
#> 50      5        0.5        V7 0.4150120
#> 51      6        0.6        V4 0.3576104
#> 52      6        0.6        V5 0.3583481
#> 53      6        0.6        V7 0.4144561
#> 54      6        0.6        V8 0.3714632
#> 55      6        0.6        V9 0.3892728
#> 56      6        0.6        V6 0.4334896
#> 57      6        0.6       V10 0.3615774
#> 58      6        0.6        V1 0.3914015
#> 59      6        0.6        V2 0.3675916
#> 60      6        0.6        V3 0.3727091
#> 61      7        0.7        V1 0.3908058
#> 62      7        0.7        V3 0.3721470
#> 63      7        0.7        V4 0.3570391
#> 64      7        0.7        V5 0.3577235
#> 65      7        0.7        V9 0.3887967
#> 66      7        0.7        V6 0.4328435
#> 67      7        0.7        V7 0.4139010
#> 68      7        0.7        V8 0.3708199
#> 69      7        0.7        V2 0.3670054
#> 70      7        0.7       V10 0.3609764
#> 71      8        0.8        V9 0.3883211
#> 72      8        0.8       V10 0.3603765
#> 73      8        0.8        V1 0.3902111
#> 74      8        0.8        V5 0.3571000
#> 75      8        0.8        V2 0.3664201
#> 76      8        0.8        V3 0.3715857
#> 77      8        0.8        V4 0.3564688
#> 78      8        0.8        V8 0.3701778
#> 79      8        0.8        V6 0.4321984
#> 80      8        0.8        V7 0.4133466
#> 81      9        0.9        V5 0.3564776
#> 82      9        0.9        V7 0.4127930
#> 83      9        0.9        V8 0.3695367
#> 84      9        0.9        V9 0.3878460
#> 85      9        0.9        V6 0.4315542
#> 86      9        0.9       V10 0.3597776
#> 87      9        0.9        V1 0.3896173
#> 88      9        0.9        V2 0.3658357
#> 89      9        0.9        V3 0.3710252
#> 90      9        0.9        V4 0.3558994
#> 91     10        1.0        V1 0.3890244
#> 92     10        1.0        V4 0.3553308
#> 93     10        1.0        V5 0.3558563
#> 94     10        1.0        V9 0.3873716
#> 95     10        1.0        V3 0.3704656
#> 96     10        1.0        V7 0.4122401
#> 97     10        1.0        V8 0.3688968
#> 98     10        1.0        V2 0.3652523
#> 99     10        1.0        V6 0.4309110
#> 100    10        1.0       V10 0.3591796

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3943932
#> 2       1        0.1        V5 0.3614875
#> 3       1        0.1        V9 0.3916625
#> 4       1        0.1        V3 0.3755326
#> 5       1        0.1        V4 0.3604804
#> 6       1        0.1        V8 0.3746962
#> 7       1        0.1        V2 0.3705368
#> 8       1        0.1        V6 0.4367345
#> 9       1        0.1        V7 0.4172429
#> 10      1        0.1       V10 0.3645971
#> 11      2        0.2        V7 0.4166840
#> 12      2        0.2        V8 0.3740474
#> 13      2        0.2        V9 0.3911834
#> 14      2        0.2       V10 0.3639912
#> 15      2        0.2        V1 0.3937931
#> 16      2        0.2        V5 0.3608575
#> 17      2        0.2        V2 0.3699458
#> 18      2        0.2        V3 0.3749662
#> 19      2        0.2        V4 0.3599046
#> 20      2        0.2        V6 0.4360835
#> 21      3        0.3        V4 0.3593297
#> 22      3        0.3        V5 0.3602285
#> 23      3        0.3        V3 0.3744007
#> 24      3        0.3        V7 0.4161259
#> 25      3        0.3        V8 0.3733996
#> 26      3        0.3        V9 0.3907049
#> 27      3        0.3        V6 0.4354336
#> 28      3        0.3       V10 0.3633862
#> 29      3        0.3        V1 0.3931938
#> 30      3        0.3        V2 0.3693559
#> 31      4        0.4        V1 0.3925954
#> 32      4        0.4        V9 0.3902270
#> 33      4        0.4        V3 0.3738360
#> 34      4        0.4        V4 0.3587557
#> 35      4        0.4        V5 0.3596006
#> 36      4        0.4        V2 0.3687668
#> 37      4        0.4        V6 0.4347846
#> 38      4        0.4        V7 0.4155686
#> 39      4        0.4        V8 0.3727530
#> 40      4        0.4       V10 0.3627823
#> 41      5        0.5        V8 0.3721075
#> 42      5        0.5        V9 0.3897496
#> 43      5        0.5       V10 0.3621793
#> 44      5        0.5        V1 0.3919980
#> 45      5        0.5        V5 0.3589738
#> 46      5        0.5        V2 0.3681787
#> 47      5        0.5        V3 0.3732721
#> 48      5        0.5        V4 0.3581826
#> 49      5        0.5        V6 0.4341366
#> 50      5        0.5        V7 0.4150120
#> 51      6        0.6        V4 0.3576104
#> 52      6        0.6        V5 0.3583481
#> 53      6        0.6        V7 0.4144561
#> 54      6        0.6        V8 0.3714632
#> 55      6        0.6        V9 0.3892728
#> 56      6        0.6        V6 0.4334896
#> 57      6        0.6       V10 0.3615774
#> 58      6        0.6        V1 0.3914015
#> 59      6        0.6        V2 0.3675916
#> 60      6        0.6        V3 0.3727091
#> 61      7        0.7        V1 0.3908058
#> 62      7        0.7        V3 0.3721470
#> 63      7        0.7        V4 0.3570391
#> 64      7        0.7        V5 0.3577235
#> 65      7        0.7        V9 0.3887967
#> 66      7        0.7        V6 0.4328435
#> 67      7        0.7        V7 0.4139010
#> 68      7        0.7        V8 0.3708199
#> 69      7        0.7        V2 0.3670054
#> 70      7        0.7       V10 0.3609764
#> 71      8        0.8        V9 0.3883211
#> 72      8        0.8       V10 0.3603765
#> 73      8        0.8        V1 0.3902111
#> 74      8        0.8        V5 0.3571000
#> 75      8        0.8        V2 0.3664201
#> 76      8        0.8        V3 0.3715857
#> 77      8        0.8        V4 0.3564688
#> 78      8        0.8        V8 0.3701778
#> 79      8        0.8        V6 0.4321984
#> 80      8        0.8        V7 0.4133466
#> 81      9        0.9        V5 0.3564776
#> 82      9        0.9        V7 0.4127930
#> 83      9        0.9        V8 0.3695367
#> 84      9        0.9        V9 0.3878460
#> 85      9        0.9        V6 0.4315542
#> 86      9        0.9       V10 0.3597776
#> 87      9        0.9        V1 0.3896173
#> 88      9        0.9        V2 0.3658357
#> 89      9        0.9        V3 0.3710252
#> 90      9        0.9        V4 0.3558994
#> 91     10        1.0        V1 0.3890244
#> 92     10        1.0        V4 0.3553308
#> 93     10        1.0        V5 0.3558563
#> 94     10        1.0        V9 0.3873716
#> 95     10        1.0        V3 0.3704656
#> 96     10        1.0        V7 0.4122401
#> 97     10        1.0        V8 0.3688968
#> 98     10        1.0        V2 0.3652523
#> 99     10        1.0        V6 0.4309110
#> 100    10        1.0       V10 0.3591796


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
#> 1       0                0          0 0.8757906 0.04611602 0.7855611 0.9435791
#> 2       0               20         20 0.8757906 0.04611602 0.7855611 0.9435791
#> 3       0               40         40 0.8757906 0.04611602 0.7855611 0.9435791
#> 4       0               60         60 0.8757906 0.04611602 0.7855611 0.9435791
#> 5      20                0         20 0.8617131 0.04928462 0.7643885 0.9353165
#> 6      20               20         40 0.8617131 0.04928462 0.7643885 0.9353165
#> 7      20               40         60 0.8617131 0.04928462 0.7643885 0.9353165
#> 8      20               60         80 0.8617131 0.04928462 0.7643885 0.9353165
#> 9      40                0         40 0.8478591 0.05223555 0.7440467 0.9269537
#> 10     40               20         60 0.8478591 0.05223555 0.7440467 0.9269537
#> 11     40               40         80 0.8478591 0.05223555 0.7440467 0.9269537
#> 12     40               60        100 0.8478591 0.05223555 0.7440467 0.9269537
#> 13     60                0         60 0.8342249 0.05499461 0.7244491 0.9185175
#> 14     60               20         80 0.8342249 0.05499461 0.7244491 0.9185175
#> 15     60               40        100 0.8342249 0.05499461 0.7244491 0.9185175
#> 16     60               60        120 0.8342249 0.05499461 0.7244491 0.9185175
#> 17     80                0         80 0.8208071 0.05758195 0.7055289 0.9100297
#> 18     80               20        100 0.8208071 0.05758195 0.7055289 0.9100297
#> 19     80               40        120 0.8208071 0.05758195 0.7055289 0.9100297
#> 20     80               60        140 0.8208071 0.05758195 0.7055289 0.9100297
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.10413368 0.189686050 0.5371633
#> 2  0.30574618 0.10396977 0.147640859 0.4957900
#> 3  0.26001915 0.10343040 0.113998269 0.4577938
#> 4  0.22113100 0.10214444 0.087163824 0.4229467
#> 5  0.25589195 0.09576146 0.109514364 0.4260157
#> 6  0.21762106 0.09373177 0.083596768 0.3938279
#> 7  0.18507391 0.09156881 0.063044565 0.3643497
#> 8  0.15739448 0.08906836 0.046862333 0.3373605
#> 9  0.18213629 0.08539400 0.060324589 0.3397362
#> 10 0.15489621 0.08252186 0.044731666 0.3148265
#> 11 0.13173012 0.07966133 0.032586918 0.2920181
#> 12 0.11202872 0.07669394 0.023242938 0.2711270
#> 13 0.12963921 0.07463685 0.031000179 0.2729665
#> 14 0.11025053 0.07150401 0.022032674 0.2536688
#> 15 0.09376159 0.06846461 0.015257509 0.2359745
#> 16 0.07973872 0.06544579 0.010241365 0.2197378
#> 17 0.09227334 0.06433921 0.014391209 0.2211690
#> 18 0.07847305 0.06126680 0.009609160 0.2061401
#> 19 0.06673672 0.05832972 0.006173284 0.1923234
#> 20 0.05675565 0.05548065 0.003785714 0.1796060
```
