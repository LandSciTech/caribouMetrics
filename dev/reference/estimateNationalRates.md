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
#> 1         0.1 0.3838513 0.01501732 0.3758418 0.4204279
#> 2         0.2 0.3832760 0.01501084 0.3752250 0.4198051
#> 3         0.3 0.3827015 0.01500464 0.3746093 0.4191832
#> 4         0.4 0.3821279 0.01499872 0.3739946 0.4185622
#> 5         0.5 0.3815551 0.01499308 0.3733809 0.4179422
#> 6         0.6 0.3809832 0.01498771 0.3727682 0.4173231
#> 7         0.7 0.3804122 0.01498261 0.3721565 0.4167049
#> 8         0.8 0.3798420 0.01497779 0.3715458 0.4160877
#> 9         0.9 0.3792726 0.01497323 0.3709361 0.4154713
#> 10        1.0 0.3787041 0.01496893 0.3703275 0.4148559

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3879088
#> 2       1        0.1        V5 0.3801421
#> 3       1        0.1        V9 0.3813512
#> 4       1        0.1        V3 0.3991321
#> 5       1        0.1        V4 0.4266105
#> 6       1        0.1        V8 0.3814825
#> 7       1        0.1        V2 0.3746566
#> 8       1        0.1        V6 0.3927174
#> 9       1        0.1        V7 0.3936485
#> 10      1        0.1       V10 0.3799239
#> 11      2        0.2        V7 0.3931379
#> 12      2        0.2        V8 0.3808474
#> 13      2        0.2        V9 0.3806964
#> 14      2        0.2       V10 0.3793021
#> 15      2        0.2        V1 0.3873474
#> 16      2        0.2        V5 0.3796382
#> 17      2        0.2        V2 0.3740414
#> 18      2        0.2        V3 0.3986173
#> 19      2        0.2        V4 0.4259564
#> 20      2        0.2        V6 0.3920337
#> 21      3        0.3        V4 0.4253032
#> 22      3        0.3        V5 0.3791350
#> 23      3        0.3        V3 0.3981032
#> 24      3        0.3        V7 0.3926280
#> 25      3        0.3        V8 0.3802133
#> 26      3        0.3        V9 0.3800428
#> 27      3        0.3        V6 0.3913513
#> 28      3        0.3       V10 0.3786812
#> 29      3        0.3        V1 0.3867867
#> 30      3        0.3        V2 0.3734272
#> 31      4        0.4        V1 0.3862269
#> 32      4        0.4        V9 0.3793903
#> 33      4        0.4        V3 0.3975898
#> 34      4        0.4        V4 0.4246510
#> 35      4        0.4        V5 0.3786324
#> 36      4        0.4        V2 0.3728139
#> 37      4        0.4        V6 0.3906700
#> 38      4        0.4        V7 0.3921188
#> 39      4        0.4        V8 0.3795803
#> 40      4        0.4       V10 0.3780614
#> 41      5        0.5        V8 0.3789483
#> 42      5        0.5        V9 0.3787389
#> 43      5        0.5       V10 0.3774426
#> 44      5        0.5        V1 0.3856679
#> 45      5        0.5        V5 0.3781305
#> 46      5        0.5        V2 0.3722017
#> 47      5        0.5        V3 0.3970770
#> 48      5        0.5        V4 0.4239998
#> 49      5        0.5        V6 0.3899899
#> 50      5        0.5        V7 0.3916102
#> 51      6        0.6        V4 0.4233497
#> 52      6        0.6        V5 0.3776292
#> 53      6        0.6        V7 0.3911023
#> 54      6        0.6        V8 0.3783174
#> 55      6        0.6        V9 0.3780886
#> 56      6        0.6        V6 0.3893110
#> 57      6        0.6       V10 0.3768247
#> 58      6        0.6        V1 0.3851097
#> 59      6        0.6        V2 0.3715905
#> 60      6        0.6        V3 0.3965649
#> 61      7        0.7        V1 0.3845523
#> 62      7        0.7        V3 0.3960535
#> 63      7        0.7        V4 0.4227005
#> 64      7        0.7        V5 0.3771286
#> 65      7        0.7        V9 0.3774395
#> 66      7        0.7        V6 0.3886333
#> 67      7        0.7        V7 0.3905951
#> 68      7        0.7        V8 0.3776876
#> 69      7        0.7        V2 0.3709803
#> 70      7        0.7       V10 0.3762080
#> 71      8        0.8        V9 0.3767915
#> 72      8        0.8       V10 0.3755922
#> 73      8        0.8        V1 0.3839957
#> 74      8        0.8        V5 0.3766287
#> 75      8        0.8        V2 0.3703711
#> 76      8        0.8        V3 0.3955427
#> 77      8        0.8        V4 0.4220523
#> 78      8        0.8        V8 0.3770588
#> 79      8        0.8        V6 0.3879568
#> 80      8        0.8        V7 0.3900885
#> 81      9        0.9        V5 0.3761294
#> 82      9        0.9        V7 0.3895826
#> 83      9        0.9        V8 0.3764310
#> 84      9        0.9        V9 0.3761445
#> 85      9        0.9        V6 0.3872815
#> 86      9        0.9       V10 0.3749774
#> 87      9        0.9        V1 0.3834399
#> 88      9        0.9        V2 0.3697629
#> 89      9        0.9        V3 0.3950325
#> 90      9        0.9        V4 0.4214051
#> 91     10        1.0        V1 0.3828849
#> 92     10        1.0        V4 0.4207589
#> 93     10        1.0        V5 0.3756308
#> 94     10        1.0        V9 0.3754987
#> 95     10        1.0        V3 0.3945231
#> 96     10        1.0        V7 0.3890773
#> 97     10        1.0        V8 0.3758043
#> 98     10        1.0        V2 0.3691557
#> 99     10        1.0        V6 0.3866073
#> 100    10        1.0       V10 0.3743636

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3879088
#> 2       1        0.1        V5 0.3801421
#> 3       1        0.1        V9 0.3813512
#> 4       1        0.1        V3 0.3991321
#> 5       1        0.1        V4 0.4266105
#> 6       1        0.1        V8 0.3814825
#> 7       1        0.1        V2 0.3746566
#> 8       1        0.1        V6 0.3927174
#> 9       1        0.1        V7 0.3936485
#> 10      1        0.1       V10 0.3799239
#> 11      2        0.2        V7 0.3931379
#> 12      2        0.2        V8 0.3808474
#> 13      2        0.2        V9 0.3806964
#> 14      2        0.2       V10 0.3793021
#> 15      2        0.2        V1 0.3873474
#> 16      2        0.2        V5 0.3796382
#> 17      2        0.2        V2 0.3740414
#> 18      2        0.2        V3 0.3986173
#> 19      2        0.2        V4 0.4259564
#> 20      2        0.2        V6 0.3920337
#> 21      3        0.3        V4 0.4253032
#> 22      3        0.3        V5 0.3791350
#> 23      3        0.3        V3 0.3981032
#> 24      3        0.3        V7 0.3926280
#> 25      3        0.3        V8 0.3802133
#> 26      3        0.3        V9 0.3800428
#> 27      3        0.3        V6 0.3913513
#> 28      3        0.3       V10 0.3786812
#> 29      3        0.3        V1 0.3867867
#> 30      3        0.3        V2 0.3734272
#> 31      4        0.4        V1 0.3862269
#> 32      4        0.4        V9 0.3793903
#> 33      4        0.4        V3 0.3975898
#> 34      4        0.4        V4 0.4246510
#> 35      4        0.4        V5 0.3786324
#> 36      4        0.4        V2 0.3728139
#> 37      4        0.4        V6 0.3906700
#> 38      4        0.4        V7 0.3921188
#> 39      4        0.4        V8 0.3795803
#> 40      4        0.4       V10 0.3780614
#> 41      5        0.5        V8 0.3789483
#> 42      5        0.5        V9 0.3787389
#> 43      5        0.5       V10 0.3774426
#> 44      5        0.5        V1 0.3856679
#> 45      5        0.5        V5 0.3781305
#> 46      5        0.5        V2 0.3722017
#> 47      5        0.5        V3 0.3970770
#> 48      5        0.5        V4 0.4239998
#> 49      5        0.5        V6 0.3899899
#> 50      5        0.5        V7 0.3916102
#> 51      6        0.6        V4 0.4233497
#> 52      6        0.6        V5 0.3776292
#> 53      6        0.6        V7 0.3911023
#> 54      6        0.6        V8 0.3783174
#> 55      6        0.6        V9 0.3780886
#> 56      6        0.6        V6 0.3893110
#> 57      6        0.6       V10 0.3768247
#> 58      6        0.6        V1 0.3851097
#> 59      6        0.6        V2 0.3715905
#> 60      6        0.6        V3 0.3965649
#> 61      7        0.7        V1 0.3845523
#> 62      7        0.7        V3 0.3960535
#> 63      7        0.7        V4 0.4227005
#> 64      7        0.7        V5 0.3771286
#> 65      7        0.7        V9 0.3774395
#> 66      7        0.7        V6 0.3886333
#> 67      7        0.7        V7 0.3905951
#> 68      7        0.7        V8 0.3776876
#> 69      7        0.7        V2 0.3709803
#> 70      7        0.7       V10 0.3762080
#> 71      8        0.8        V9 0.3767915
#> 72      8        0.8       V10 0.3755922
#> 73      8        0.8        V1 0.3839957
#> 74      8        0.8        V5 0.3766287
#> 75      8        0.8        V2 0.3703711
#> 76      8        0.8        V3 0.3955427
#> 77      8        0.8        V4 0.4220523
#> 78      8        0.8        V8 0.3770588
#> 79      8        0.8        V6 0.3879568
#> 80      8        0.8        V7 0.3900885
#> 81      9        0.9        V5 0.3761294
#> 82      9        0.9        V7 0.3895826
#> 83      9        0.9        V8 0.3764310
#> 84      9        0.9        V9 0.3761445
#> 85      9        0.9        V6 0.3872815
#> 86      9        0.9       V10 0.3749774
#> 87      9        0.9        V1 0.3834399
#> 88      9        0.9        V2 0.3697629
#> 89      9        0.9        V3 0.3950325
#> 90      9        0.9        V4 0.4214051
#> 91     10        1.0        V1 0.3828849
#> 92     10        1.0        V4 0.4207589
#> 93     10        1.0        V5 0.3756308
#> 94     10        1.0        V9 0.3754987
#> 95     10        1.0        V3 0.3945231
#> 96     10        1.0        V7 0.3890773
#> 97     10        1.0        V8 0.3758043
#> 98     10        1.0        V2 0.3691557
#> 99     10        1.0        V6 0.3866073
#> 100    10        1.0       V10 0.3743636


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
#> 1       0                0          0 0.8757906 0.05065070 0.7730320 0.9393025
#> 2       0               20         20 0.8757906 0.05065070 0.7730320 0.9393025
#> 3       0               40         40 0.8757906 0.05065070 0.7730320 0.9393025
#> 4       0               60         60 0.8757906 0.05065070 0.7730320 0.9393025
#> 5      20                0         20 0.8617131 0.05231800 0.7556914 0.9276194
#> 6      20               20         40 0.8617131 0.05231800 0.7556914 0.9276194
#> 7      20               40         60 0.8617131 0.05231800 0.7556914 0.9276194
#> 8      20               60         80 0.8617131 0.05231800 0.7556914 0.9276194
#> 9      40                0         40 0.8478591 0.05381202 0.7389281 0.9157914
#> 10     40               20         60 0.8478591 0.05381202 0.7389281 0.9157914
#> 11     40               40         80 0.8478591 0.05381202 0.7389281 0.9157914
#> 12     40               60        100 0.8478591 0.05381202 0.7389281 0.9157914
#> 13     60                0         60 0.8342249 0.05515853 0.7226903 0.9038763
#> 14     60               20         80 0.8342249 0.05515853 0.7226903 0.9038763
#> 15     60               40        100 0.8342249 0.05515853 0.7226903 0.9038763
#> 16     60               60        120 0.8342249 0.05515853 0.7226903 0.9038763
#> 17     80                0         80 0.8208071 0.05637758 0.7069361 0.8919181
#> 18     80               20        100 0.8208071 0.05637758 0.7069361 0.8919181
#> 19     80               40        120 0.8208071 0.05637758 0.7069361 0.8919181
#> 20     80               60        140 0.8208071 0.05637758 0.7069361 0.8919181
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11538114 0.150670892 0.5460294
#> 2  0.30574618 0.11528275 0.114797266 0.5149970
#> 3  0.26001915 0.11517993 0.086500622 0.4858313
#> 4  0.22113100 0.11439555 0.064306008 0.4584459
#> 5  0.25589195 0.10777168 0.079171664 0.4502386
#> 6  0.21762106 0.10568368 0.058586010 0.4250522
#> 7  0.18507391 0.10378512 0.042614127 0.4014345
#> 8  0.15739448 0.10162694 0.030363075 0.3792943
#> 9  0.18213629 0.09760703 0.038535057 0.3726642
#> 10 0.15489621 0.09464412 0.027263473 0.3523282
#> 11 0.13173012 0.09193999 0.018789003 0.3332679
#> 12 0.11202872 0.08919535 0.012545458 0.3154016
#> 13 0.12963921 0.08671573 0.016679544 0.3100507
#> 14 0.11025053 0.08349670 0.011016775 0.2936338
#> 15 0.09376159 0.08054872 0.006987944 0.2782367
#> 16 0.07973872 0.07766979 0.004220667 0.2637906
#> 17 0.09227334 0.07610677 0.006029252 0.2594607
#> 18 0.07847305 0.07296621 0.003580878 0.2461650
#> 19 0.06673672 0.07008078 0.001997352 0.2336761
#> 20 0.05675565 0.06731141 0.001032743 0.2219379
```
