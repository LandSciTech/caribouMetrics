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
#> 1         0.1 0.3838513 0.02142434 0.3489456 0.4099852
#> 2         0.2 0.3832760 0.02137357 0.3484536 0.4093213
#> 3         0.3 0.3827015 0.02132296 0.3479622 0.4086585
#> 4         0.4 0.3821279 0.02127250 0.3474715 0.4079968
#> 5         0.5 0.3815551 0.02122219 0.3469815 0.4073361
#> 6         0.6 0.3809832 0.02117203 0.3464922 0.4066765
#> 7         0.7 0.3804122 0.02112202 0.3460036 0.4060180
#> 8         0.8 0.3798420 0.02107216 0.3455157 0.4053606
#> 9         0.9 0.3792726 0.02102245 0.3450284 0.4047042
#> 10        1.0 0.3787041 0.02097290 0.3445419 0.4040489

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3869836
#> 2       1        0.1        V5 0.3825259
#> 3       1        0.1        V9 0.4103797
#> 4       1        0.1        V3 0.3806388
#> 5       1        0.1        V4 0.3602686
#> 6       1        0.1        V8 0.3934771
#> 7       1        0.1        V2 0.4086264
#> 8       1        0.1        V6 0.3456583
#> 9       1        0.1        V7 0.3667118
#> 10      1        0.1       V10 0.4036079
#> 11      2        0.2        V7 0.3661474
#> 12      2        0.2        V8 0.3928815
#> 13      2        0.2        V9 0.4097011
#> 14      2        0.2       V10 0.4029876
#> 15      2        0.2        V1 0.3863598
#> 16      2        0.2        V5 0.3818690
#> 17      2        0.2        V2 0.4080131
#> 18      2        0.2        V3 0.3800785
#> 19      2        0.2        V4 0.3597457
#> 20      2        0.2        V6 0.3451752
#> 21      3        0.3        V4 0.3592236
#> 22      3        0.3        V5 0.3812132
#> 23      3        0.3        V3 0.3795189
#> 24      3        0.3        V7 0.3655838
#> 25      3        0.3        V8 0.3922868
#> 26      3        0.3        V9 0.4090237
#> 27      3        0.3        V6 0.3446927
#> 28      3        0.3       V10 0.4023683
#> 29      3        0.3        V1 0.3857370
#> 30      3        0.3        V2 0.4074006
#> 31      4        0.4        V1 0.3851153
#> 32      4        0.4        V9 0.4083474
#> 33      4        0.4        V3 0.3789602
#> 34      4        0.4        V4 0.3587023
#> 35      4        0.4        V5 0.3805586
#> 36      4        0.4        V2 0.4067891
#> 37      4        0.4        V6 0.3442109
#> 38      4        0.4        V7 0.3650211
#> 39      4        0.4        V8 0.3916930
#> 40      4        0.4       V10 0.4017500
#> 41      5        0.5        V8 0.3911002
#> 42      5        0.5        V9 0.4076722
#> 43      5        0.5       V10 0.4011326
#> 44      5        0.5        V1 0.3844945
#> 45      5        0.5        V5 0.3799051
#> 46      5        0.5        V2 0.4061785
#> 47      5        0.5        V3 0.3784023
#> 48      5        0.5        V4 0.3581817
#> 49      5        0.5        V6 0.3437298
#> 50      5        0.5        V7 0.3644592
#> 51      6        0.6        V4 0.3576619
#> 52      6        0.6        V5 0.3792526
#> 53      6        0.6        V7 0.3638983
#> 54      6        0.6        V8 0.3905082
#> 55      6        0.6        V9 0.4069981
#> 56      6        0.6        V6 0.3432494
#> 57      6        0.6       V10 0.4005162
#> 58      6        0.6        V1 0.3838748
#> 59      6        0.6        V2 0.4055688
#> 60      6        0.6        V3 0.3778453
#> 61      7        0.7        V1 0.3832560
#> 62      7        0.7        V3 0.3772890
#> 63      7        0.7        V4 0.3571428
#> 64      7        0.7        V5 0.3786014
#> 65      7        0.7        V9 0.4063252
#> 66      7        0.7        V6 0.3427696
#> 67      7        0.7        V7 0.3633381
#> 68      7        0.7        V8 0.3899171
#> 69      7        0.7        V2 0.4049600
#> 70      7        0.7       V10 0.3999007
#> 71      8        0.8        V9 0.4056533
#> 72      8        0.8       V10 0.3992862
#> 73      8        0.8        V1 0.3826383
#> 74      8        0.8        V5 0.3779512
#> 75      8        0.8        V2 0.4043521
#> 76      8        0.8        V3 0.3767336
#> 77      8        0.8        V4 0.3566245
#> 78      8        0.8        V8 0.3893269
#> 79      8        0.8        V6 0.3422905
#> 80      8        0.8        V7 0.3627789
#> 81      9        0.9        V5 0.3773021
#> 82      9        0.9        V7 0.3622205
#> 83      9        0.9        V8 0.3887377
#> 84      9        0.9        V9 0.4049826
#> 85      9        0.9        V6 0.3418121
#> 86      9        0.9       V10 0.3986726
#> 87      9        0.9        V1 0.3820215
#> 88      9        0.9        V2 0.4037452
#> 89      9        0.9        V3 0.3761790
#> 90      9        0.9        V4 0.3561069
#> 91     10        1.0        V1 0.3814058
#> 92     10        1.0        V4 0.3555901
#> 93     10        1.0        V5 0.3766542
#> 94     10        1.0        V9 0.4043130
#> 95     10        1.0        V3 0.3756252
#> 96     10        1.0        V7 0.3616630
#> 97     10        1.0        V8 0.3881493
#> 98     10        1.0        V2 0.4031392
#> 99     10        1.0        V6 0.3413343
#> 100    10        1.0       V10 0.3980599

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3869836
#> 2       1        0.1        V5 0.3825259
#> 3       1        0.1        V9 0.4103797
#> 4       1        0.1        V3 0.3806388
#> 5       1        0.1        V4 0.3602686
#> 6       1        0.1        V8 0.3934771
#> 7       1        0.1        V2 0.4086264
#> 8       1        0.1        V6 0.3456583
#> 9       1        0.1        V7 0.3667118
#> 10      1        0.1       V10 0.4036079
#> 11      2        0.2        V7 0.3661474
#> 12      2        0.2        V8 0.3928815
#> 13      2        0.2        V9 0.4097011
#> 14      2        0.2       V10 0.4029876
#> 15      2        0.2        V1 0.3863598
#> 16      2        0.2        V5 0.3818690
#> 17      2        0.2        V2 0.4080131
#> 18      2        0.2        V3 0.3800785
#> 19      2        0.2        V4 0.3597457
#> 20      2        0.2        V6 0.3451752
#> 21      3        0.3        V4 0.3592236
#> 22      3        0.3        V5 0.3812132
#> 23      3        0.3        V3 0.3795189
#> 24      3        0.3        V7 0.3655838
#> 25      3        0.3        V8 0.3922868
#> 26      3        0.3        V9 0.4090237
#> 27      3        0.3        V6 0.3446927
#> 28      3        0.3       V10 0.4023683
#> 29      3        0.3        V1 0.3857370
#> 30      3        0.3        V2 0.4074006
#> 31      4        0.4        V1 0.3851153
#> 32      4        0.4        V9 0.4083474
#> 33      4        0.4        V3 0.3789602
#> 34      4        0.4        V4 0.3587023
#> 35      4        0.4        V5 0.3805586
#> 36      4        0.4        V2 0.4067891
#> 37      4        0.4        V6 0.3442109
#> 38      4        0.4        V7 0.3650211
#> 39      4        0.4        V8 0.3916930
#> 40      4        0.4       V10 0.4017500
#> 41      5        0.5        V8 0.3911002
#> 42      5        0.5        V9 0.4076722
#> 43      5        0.5       V10 0.4011326
#> 44      5        0.5        V1 0.3844945
#> 45      5        0.5        V5 0.3799051
#> 46      5        0.5        V2 0.4061785
#> 47      5        0.5        V3 0.3784023
#> 48      5        0.5        V4 0.3581817
#> 49      5        0.5        V6 0.3437298
#> 50      5        0.5        V7 0.3644592
#> 51      6        0.6        V4 0.3576619
#> 52      6        0.6        V5 0.3792526
#> 53      6        0.6        V7 0.3638983
#> 54      6        0.6        V8 0.3905082
#> 55      6        0.6        V9 0.4069981
#> 56      6        0.6        V6 0.3432494
#> 57      6        0.6       V10 0.4005162
#> 58      6        0.6        V1 0.3838748
#> 59      6        0.6        V2 0.4055688
#> 60      6        0.6        V3 0.3778453
#> 61      7        0.7        V1 0.3832560
#> 62      7        0.7        V3 0.3772890
#> 63      7        0.7        V4 0.3571428
#> 64      7        0.7        V5 0.3786014
#> 65      7        0.7        V9 0.4063252
#> 66      7        0.7        V6 0.3427696
#> 67      7        0.7        V7 0.3633381
#> 68      7        0.7        V8 0.3899171
#> 69      7        0.7        V2 0.4049600
#> 70      7        0.7       V10 0.3999007
#> 71      8        0.8        V9 0.4056533
#> 72      8        0.8       V10 0.3992862
#> 73      8        0.8        V1 0.3826383
#> 74      8        0.8        V5 0.3779512
#> 75      8        0.8        V2 0.4043521
#> 76      8        0.8        V3 0.3767336
#> 77      8        0.8        V4 0.3566245
#> 78      8        0.8        V8 0.3893269
#> 79      8        0.8        V6 0.3422905
#> 80      8        0.8        V7 0.3627789
#> 81      9        0.9        V5 0.3773021
#> 82      9        0.9        V7 0.3622205
#> 83      9        0.9        V8 0.3887377
#> 84      9        0.9        V9 0.4049826
#> 85      9        0.9        V6 0.3418121
#> 86      9        0.9       V10 0.3986726
#> 87      9        0.9        V1 0.3820215
#> 88      9        0.9        V2 0.4037452
#> 89      9        0.9        V3 0.3761790
#> 90      9        0.9        V4 0.3561069
#> 91     10        1.0        V1 0.3814058
#> 92     10        1.0        V4 0.3555901
#> 93     10        1.0        V5 0.3766542
#> 94     10        1.0        V9 0.4043130
#> 95     10        1.0        V3 0.3756252
#> 96     10        1.0        V7 0.3616630
#> 97     10        1.0        V8 0.3881493
#> 98     10        1.0        V2 0.4031392
#> 99     10        1.0        V6 0.3413343
#> 100    10        1.0       V10 0.3980599


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
#> 1       0                0          0 0.8757906 0.04701223 0.7772624 0.9393645
#> 2       0               20         20 0.8757906 0.04701223 0.7772624 0.9393645
#> 3       0               40         40 0.8757906 0.04701223 0.7772624 0.9393645
#> 4       0               60         60 0.8757906 0.04701223 0.7772624 0.9393645
#> 5      20                0         20 0.8617131 0.04865952 0.7610652 0.9317019
#> 6      20               20         40 0.8617131 0.04865952 0.7610652 0.9317019
#> 7      20               40         60 0.8617131 0.04865952 0.7610652 0.9317019
#> 8      20               60         80 0.8617131 0.04865952 0.7610652 0.9317019
#> 9      40                0         40 0.8478591 0.05029089 0.7453708 0.9239751
#> 10     40               20         60 0.8478591 0.05029089 0.7453708 0.9239751
#> 11     40               40         80 0.8478591 0.05029089 0.7453708 0.9239751
#> 12     40               60        100 0.8478591 0.05029089 0.7453708 0.9239751
#> 13     60                0         60 0.8342249 0.05190562 0.7301362 0.9162021
#> 14     60               20         80 0.8342249 0.05190562 0.7301362 0.9162021
#> 15     60               40        100 0.8342249 0.05190562 0.7301362 0.9162021
#> 16     60               60        120 0.8342249 0.05190562 0.7301362 0.9162021
#> 17     80                0         80 0.8208071 0.05350188 0.7153267 0.9083977
#> 18     80               20        100 0.8208071 0.05350188 0.7153267 0.9083977
#> 19     80               40        120 0.8208071 0.05350188 0.7153267 0.9083977
#> 20     80               60        140 0.8208071 0.05350188 0.7153267 0.9083977
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11323634 0.182707448 0.5684821
#> 2  0.30574618 0.10999872 0.139286176 0.5025934
#> 3  0.26001915 0.10583907 0.105137514 0.4448416
#> 4  0.22113100 0.10083442 0.078395159 0.3943859
#> 5  0.25589195 0.10247388 0.106804550 0.4585653
#> 6  0.21762106 0.09776786 0.079696922 0.4063664
#> 7  0.18507391 0.09265879 0.058598991 0.3608163
#> 8  0.15739448 0.08721447 0.042324833 0.3210835
#> 9  0.18213629 0.09079832 0.059622104 0.3716298
#> 10 0.15489621 0.08560631 0.043109827 0.3305170
#> 11 0.13173012 0.08030105 0.030513583 0.2946397
#> 12 0.11202872 0.07493274 0.021048980 0.2632885
#> 13 0.12963921 0.07941266 0.031117216 0.3031606
#> 14 0.11025053 0.07423697 0.021498536 0.2707398
#> 15 0.09376159 0.06910845 0.014401984 0.2423648
#> 16 0.07973872 0.06405960 0.009293430 0.2174626
#> 17 0.09227334 0.06884342 0.014735470 0.2491137
#> 18 0.07847305 0.06393425 0.009530043 0.2233928
#> 19 0.06673672 0.05915626 0.005891911 0.2007623
#> 20 0.05675565 0.05452756 0.003448084 0.1807711
```
