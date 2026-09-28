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
[`addMissingYears()`](https://landscitech.github.io/caribouMetrics/dev/reference/addMissingYears.md),
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
[`getCaribouYear()`](https://landscitech.github.io/caribouMetrics/dev/reference/getCaribouYear.md),
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
#> 1         0.1 0.3838513 0.01840283 0.3606373 0.4125699
#> 2         0.2 0.3832760 0.01834177 0.3601943 0.4119290
#> 3         0.3 0.3827015 0.01828106 0.3597518 0.4112890
#> 4         0.4 0.3821279 0.01822070 0.3593099 0.4106501
#> 5         0.5 0.3815551 0.01816069 0.3588686 0.4100121
#> 6         0.6 0.3809832 0.01810102 0.3584279 0.4093751
#> 7         0.7 0.3804122 0.01804170 0.3579877 0.4087392
#> 8         0.8 0.3798420 0.01798273 0.3575481 0.4081042
#> 9         0.9 0.3792726 0.01792411 0.3571091 0.4074702
#> 10        1.0 0.3787041 0.01786583 0.3566706 0.4068371

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3594988
#> 2       1        0.1        V5 0.3920200
#> 3       1        0.1        V9 0.4032897
#> 4       1        0.1        V3 0.3726516
#> 5       1        0.1        V4 0.3756742
#> 6       1        0.1        V8 0.3676891
#> 7       1        0.1        V2 0.3959049
#> 8       1        0.1        V6 0.4152642
#> 9       1        0.1        V7 0.3747931
#> 10      1        0.1       V10 0.3645588
#> 11      2        0.2        V7 0.3741903
#> 12      2        0.2        V8 0.3672205
#> 13      2        0.2        V9 0.4026581
#> 14      2        0.2       V10 0.3639307
#> 15      2        0.2        V1 0.3591095
#> 16      2        0.2        V5 0.3913761
#> 17      2        0.2        V2 0.3952859
#> 18      2        0.2        V3 0.3721105
#> 19      2        0.2        V4 0.3751684
#> 20      2        0.2        V6 0.4146205
#> 21      3        0.3        V4 0.3746633
#> 22      3        0.3        V5 0.3907333
#> 23      3        0.3        V3 0.3715701
#> 24      3        0.3        V7 0.3735885
#> 25      3        0.3        V8 0.3667524
#> 26      3        0.3        V9 0.4020276
#> 27      3        0.3        V6 0.4139778
#> 28      3        0.3       V10 0.3633036
#> 29      3        0.3        V1 0.3587206
#> 30      3        0.3        V2 0.3946679
#> 31      4        0.4        V1 0.3583322
#> 32      4        0.4        V9 0.4013980
#> 33      4        0.4        V3 0.3710306
#> 34      4        0.4        V4 0.3741588
#> 35      4        0.4        V5 0.3900915
#> 36      4        0.4        V2 0.3940508
#> 37      4        0.4        V6 0.4133362
#> 38      4        0.4        V7 0.3729877
#> 39      4        0.4        V8 0.3662850
#> 40      4        0.4       V10 0.3626777
#> 41      5        0.5        V8 0.3658181
#> 42      5        0.5        V9 0.4007694
#> 43      5        0.5       V10 0.3620528
#> 44      5        0.5        V1 0.3579442
#> 45      5        0.5        V5 0.3894508
#> 46      5        0.5        V2 0.3934347
#> 47      5        0.5        V3 0.3704919
#> 48      5        0.5        V4 0.3736551
#> 49      5        0.5        V6 0.4126955
#> 50      5        0.5        V7 0.3723878
#> 51      6        0.6        V4 0.3731520
#> 52      6        0.6        V5 0.3888112
#> 53      6        0.6        V7 0.3717888
#> 54      6        0.6        V8 0.3653519
#> 55      6        0.6        V9 0.4001417
#> 56      6        0.6        V6 0.4120558
#> 57      6        0.6       V10 0.3614290
#> 58      6        0.6        V1 0.3575566
#> 59      6        0.6        V2 0.3928196
#> 60      6        0.6        V3 0.3699539
#> 61      7        0.7        V1 0.3571695
#> 62      7        0.7        V3 0.3694167
#> 63      7        0.7        V4 0.3726495
#> 64      7        0.7        V5 0.3881725
#> 65      7        0.7        V9 0.3995151
#> 66      7        0.7        V6 0.4114171
#> 67      7        0.7        V7 0.3711909
#> 68      7        0.7        V8 0.3648863
#> 69      7        0.7        V2 0.3922054
#> 70      7        0.7       V10 0.3608063
#> 71      8        0.8        V9 0.3988894
#> 72      8        0.8       V10 0.3601846
#> 73      8        0.8        V1 0.3567827
#> 74      8        0.8        V5 0.3875350
#> 75      8        0.8        V2 0.3915921
#> 76      8        0.8        V3 0.3688803
#> 77      8        0.8        V4 0.3721478
#> 78      8        0.8        V8 0.3644212
#> 79      8        0.8        V6 0.4107794
#> 80      8        0.8        V7 0.3705939
#> 81      9        0.9        V5 0.3868985
#> 82      9        0.9        V7 0.3699978
#> 83      9        0.9        V8 0.3639567
#> 84      9        0.9        V9 0.3982648
#> 85      9        0.9        V6 0.4101427
#> 86      9        0.9       V10 0.3595640
#> 87      9        0.9        V1 0.3563964
#> 88      9        0.9        V2 0.3909799
#> 89      9        0.9        V3 0.3683447
#> 90      9        0.9        V4 0.3716467
#> 91     10        1.0        V1 0.3560105
#> 92     10        1.0        V4 0.3711464
#> 93     10        1.0        V5 0.3862630
#> 94     10        1.0        V9 0.3976410
#> 95     10        1.0        V3 0.3678099
#> 96     10        1.0        V7 0.3694028
#> 97     10        1.0        V8 0.3634929
#> 98     10        1.0        V2 0.3903686
#> 99     10        1.0        V6 0.4095070
#> 100    10        1.0       V10 0.3589445

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3594988
#> 2       1        0.1        V5 0.3920200
#> 3       1        0.1        V9 0.4032897
#> 4       1        0.1        V3 0.3726516
#> 5       1        0.1        V4 0.3756742
#> 6       1        0.1        V8 0.3676891
#> 7       1        0.1        V2 0.3959049
#> 8       1        0.1        V6 0.4152642
#> 9       1        0.1        V7 0.3747931
#> 10      1        0.1       V10 0.3645588
#> 11      2        0.2        V7 0.3741903
#> 12      2        0.2        V8 0.3672205
#> 13      2        0.2        V9 0.4026581
#> 14      2        0.2       V10 0.3639307
#> 15      2        0.2        V1 0.3591095
#> 16      2        0.2        V5 0.3913761
#> 17      2        0.2        V2 0.3952859
#> 18      2        0.2        V3 0.3721105
#> 19      2        0.2        V4 0.3751684
#> 20      2        0.2        V6 0.4146205
#> 21      3        0.3        V4 0.3746633
#> 22      3        0.3        V5 0.3907333
#> 23      3        0.3        V3 0.3715701
#> 24      3        0.3        V7 0.3735885
#> 25      3        0.3        V8 0.3667524
#> 26      3        0.3        V9 0.4020276
#> 27      3        0.3        V6 0.4139778
#> 28      3        0.3       V10 0.3633036
#> 29      3        0.3        V1 0.3587206
#> 30      3        0.3        V2 0.3946679
#> 31      4        0.4        V1 0.3583322
#> 32      4        0.4        V9 0.4013980
#> 33      4        0.4        V3 0.3710306
#> 34      4        0.4        V4 0.3741588
#> 35      4        0.4        V5 0.3900915
#> 36      4        0.4        V2 0.3940508
#> 37      4        0.4        V6 0.4133362
#> 38      4        0.4        V7 0.3729877
#> 39      4        0.4        V8 0.3662850
#> 40      4        0.4       V10 0.3626777
#> 41      5        0.5        V8 0.3658181
#> 42      5        0.5        V9 0.4007694
#> 43      5        0.5       V10 0.3620528
#> 44      5        0.5        V1 0.3579442
#> 45      5        0.5        V5 0.3894508
#> 46      5        0.5        V2 0.3934347
#> 47      5        0.5        V3 0.3704919
#> 48      5        0.5        V4 0.3736551
#> 49      5        0.5        V6 0.4126955
#> 50      5        0.5        V7 0.3723878
#> 51      6        0.6        V4 0.3731520
#> 52      6        0.6        V5 0.3888112
#> 53      6        0.6        V7 0.3717888
#> 54      6        0.6        V8 0.3653519
#> 55      6        0.6        V9 0.4001417
#> 56      6        0.6        V6 0.4120558
#> 57      6        0.6       V10 0.3614290
#> 58      6        0.6        V1 0.3575566
#> 59      6        0.6        V2 0.3928196
#> 60      6        0.6        V3 0.3699539
#> 61      7        0.7        V1 0.3571695
#> 62      7        0.7        V3 0.3694167
#> 63      7        0.7        V4 0.3726495
#> 64      7        0.7        V5 0.3881725
#> 65      7        0.7        V9 0.3995151
#> 66      7        0.7        V6 0.4114171
#> 67      7        0.7        V7 0.3711909
#> 68      7        0.7        V8 0.3648863
#> 69      7        0.7        V2 0.3922054
#> 70      7        0.7       V10 0.3608063
#> 71      8        0.8        V9 0.3988894
#> 72      8        0.8       V10 0.3601846
#> 73      8        0.8        V1 0.3567827
#> 74      8        0.8        V5 0.3875350
#> 75      8        0.8        V2 0.3915921
#> 76      8        0.8        V3 0.3688803
#> 77      8        0.8        V4 0.3721478
#> 78      8        0.8        V8 0.3644212
#> 79      8        0.8        V6 0.4107794
#> 80      8        0.8        V7 0.3705939
#> 81      9        0.9        V5 0.3868985
#> 82      9        0.9        V7 0.3699978
#> 83      9        0.9        V8 0.3639567
#> 84      9        0.9        V9 0.3982648
#> 85      9        0.9        V6 0.4101427
#> 86      9        0.9       V10 0.3595640
#> 87      9        0.9        V1 0.3563964
#> 88      9        0.9        V2 0.3909799
#> 89      9        0.9        V3 0.3683447
#> 90      9        0.9        V4 0.3716467
#> 91     10        1.0        V1 0.3560105
#> 92     10        1.0        V4 0.3711464
#> 93     10        1.0        V5 0.3862630
#> 94     10        1.0        V9 0.3976410
#> 95     10        1.0        V3 0.3678099
#> 96     10        1.0        V7 0.3694028
#> 97     10        1.0        V8 0.3634929
#> 98     10        1.0        V2 0.3903686
#> 99     10        1.0        V6 0.4095070
#> 100    10        1.0       V10 0.3589445


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
#> 1       0                0          0 0.8757906 0.05333972 0.7733898 0.9474673
#> 2       0               20         20 0.8757906 0.05333972 0.7733898 0.9474673
#> 3       0               40         40 0.8757906 0.05333972 0.7733898 0.9474673
#> 4       0               60         60 0.8757906 0.05333972 0.7733898 0.9474673
#> 5      20                0         20 0.8617131 0.05541037 0.7547268 0.9386083
#> 6      20               20         40 0.8617131 0.05541037 0.7547268 0.9386083
#> 7      20               40         60 0.8617131 0.05541037 0.7547268 0.9386083
#> 8      20               60         80 0.8617131 0.05541037 0.7547268 0.9386083
#> 9      40                0         40 0.8478591 0.05732475 0.7367157 0.9296266
#> 10     40               20         60 0.8478591 0.05732475 0.7367157 0.9296266
#> 11     40               40         80 0.8478591 0.05732475 0.7367157 0.9296266
#> 12     40               60        100 0.8478591 0.05732475 0.7367157 0.9296266
#> 13     60                0         60 0.8342249 0.05910173 0.7192966 0.9205575
#> 14     60               20         80 0.8342249 0.05910173 0.7192966 0.9205575
#> 15     60               40        100 0.8342249 0.05910173 0.7192966 0.9205575
#> 16     60               60        120 0.8342249 0.05910173 0.7192966 0.9205575
#> 17     80                0         80 0.8208071 0.06075627 0.7024216 0.9114289
#> 18     80               20        100 0.8208071 0.06075627 0.7024216 0.9114289
#> 19     80               40        120 0.8208071 0.06075627 0.7024216 0.9114289
#> 20     80               60        140 0.8208071 0.06075627 0.7024216 0.9114289
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12313236 0.149906493 0.5417560
#> 2  0.30574618 0.11968610 0.110065421 0.4933593
#> 3  0.26001915 0.11539799 0.079611767 0.4495915
#> 4  0.22113100 0.11048593 0.056520056 0.4100815
#> 5  0.25589195 0.11036514 0.083501523 0.4360910
#> 6  0.21762106 0.10546467 0.059455805 0.3979038
#> 7  0.18507391 0.10029207 0.041400743 0.3634734
#> 8  0.15739448 0.09495835 0.028050092 0.3324379
#> 9  0.18213629 0.09674472 0.043683774 0.3528655
#> 10 0.15489621 0.09142581 0.029724274 0.3228753
#> 11 0.13173012 0.08613467 0.019576919 0.2958329
#> 12 0.11202872 0.08092490 0.012387866 0.2714316
#> 13 0.12963921 0.08355398 0.020837726 0.2874965
#> 14 0.11025053 0.07835666 0.013268792 0.2639045
#> 15 0.09376159 0.07334065 0.008055149 0.2425863
#> 16 0.07973872 0.06852614 0.004611100 0.2232949
#> 17 0.09227334 0.07139561 0.008684350 0.2360025
#> 18 0.07847305 0.06659055 0.005017364 0.2173304
#> 19 0.06673672 0.06203599 0.002702475 0.2003935
#> 20 0.05675565 0.05773124 0.001335563 0.1849980
```
