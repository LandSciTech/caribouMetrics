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
#> 1         0.1 0.3838513 0.02010149 0.3658234 0.4223580
#> 2         0.2 0.3832760 0.02006614 0.3652684 0.4216806
#> 3         0.3 0.3827015 0.02003091 0.3647143 0.4210043
#> 4         0.4 0.3821279 0.01999581 0.3641610 0.4203291
#> 5         0.5 0.3815551 0.01996083 0.3636085 0.4196550
#> 6         0.6 0.3809832 0.01992597 0.3630569 0.4189819
#> 7         0.7 0.3804122 0.01989123 0.3625061 0.4183099
#> 8         0.8 0.3798420 0.01985662 0.3619562 0.4176390
#> 9         0.9 0.3792726 0.01982212 0.3614071 0.4169692
#> 10        1.0 0.3787041 0.01978775 0.3608588 0.4163005

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3655479
#> 2       1        0.1        V5 0.4177900
#> 3       1        0.1        V9 0.3941529
#> 4       1        0.1        V3 0.3937916
#> 5       1        0.1        V4 0.3667722
#> 6       1        0.1        V8 0.4068546
#> 7       1        0.1        V2 0.4071452
#> 8       1        0.1        V6 0.3971648
#> 9       1        0.1        V7 0.4236842
#> 10      1        0.1       V10 0.3766876
#> 11      2        0.2        V7 0.4229918
#> 12      2        0.2        V8 0.4062172
#> 13      2        0.2        V9 0.3935735
#> 14      2        0.2       V10 0.3761104
#> 15      2        0.2        V1 0.3649894
#> 16      2        0.2        V5 0.4171644
#> 17      2        0.2        V2 0.4065754
#> 18      2        0.2        V3 0.3931527
#> 19      2        0.2        V4 0.3662292
#> 20      2        0.2        V6 0.3966245
#> 21      3        0.3        V4 0.3656870
#> 22      3        0.3        V5 0.4165396
#> 23      3        0.3        V3 0.3925149
#> 24      3        0.3        V7 0.4223005
#> 25      3        0.3        V8 0.4055808
#> 26      3        0.3        V9 0.3929949
#> 27      3        0.3        V6 0.3960850
#> 28      3        0.3       V10 0.3755341
#> 29      3        0.3        V1 0.3644318
#> 30      3        0.3        V2 0.4060065
#> 31      4        0.4        V1 0.3638751
#> 32      4        0.4        V9 0.3924172
#> 33      4        0.4        V3 0.3918781
#> 34      4        0.4        V4 0.3651456
#> 35      4        0.4        V5 0.4159158
#> 36      4        0.4        V2 0.4054383
#> 37      4        0.4        V6 0.3955462
#> 38      4        0.4        V7 0.4216104
#> 39      4        0.4        V8 0.4049455
#> 40      4        0.4       V10 0.3749587
#> 41      5        0.5        V8 0.4043111
#> 42      5        0.5        V9 0.3918403
#> 43      5        0.5       V10 0.3743841
#> 44      5        0.5        V1 0.3633192
#> 45      5        0.5        V5 0.4152929
#> 46      5        0.5        V2 0.4048710
#> 47      5        0.5        V3 0.3912423
#> 48      5        0.5        V4 0.3646050
#> 49      5        0.5        V6 0.3950081
#> 50      5        0.5        V7 0.4209213
#> 51      6        0.6        V4 0.3640653
#> 52      6        0.6        V5 0.4146710
#> 53      6        0.6        V7 0.4202335
#> 54      6        0.6        V8 0.4036777
#> 55      6        0.6        V9 0.3912642
#> 56      6        0.6        V6 0.3944708
#> 57      6        0.6       V10 0.3738104
#> 58      6        0.6        V1 0.3627641
#> 59      6        0.6        V2 0.4043044
#> 60      6        0.6        V3 0.3906076
#> 61      7        0.7        V1 0.3622099
#> 62      7        0.7        V3 0.3899739
#> 63      7        0.7        V4 0.3635263
#> 64      7        0.7        V5 0.4140500
#> 65      7        0.7        V9 0.3906891
#> 66      7        0.7        V6 0.3939342
#> 67      7        0.7        V7 0.4195467
#> 68      7        0.7        V8 0.4030453
#> 69      7        0.7        V2 0.4037386
#> 70      7        0.7       V10 0.3732376
#> 71      8        0.8        V9 0.3901147
#> 72      8        0.8       V10 0.3726657
#> 73      8        0.8        V1 0.3616566
#> 74      8        0.8        V5 0.4134299
#> 75      8        0.8        V2 0.4031736
#> 76      8        0.8        V3 0.3893412
#> 77      8        0.8        V4 0.3629881
#> 78      8        0.8        V8 0.4024138
#> 79      8        0.8        V6 0.3933983
#> 80      8        0.8        V7 0.4188610
#> 81      9        0.9        V5 0.4128108
#> 82      9        0.9        V7 0.4181765
#> 83      9        0.9        V8 0.4017834
#> 84      9        0.9        V9 0.3895412
#> 85      9        0.9        V6 0.3928632
#> 86      9        0.9       V10 0.3720947
#> 87      9        0.9        V1 0.3611041
#> 88      9        0.9        V2 0.4026095
#> 89      9        0.9        V3 0.3887095
#> 90      9        0.9        V4 0.3624507
#> 91     10        1.0        V1 0.3605524
#> 92     10        1.0        V4 0.3619141
#> 93     10        1.0        V5 0.4121925
#> 94     10        1.0        V9 0.3889685
#> 95     10        1.0        V3 0.3880789
#> 96     10        1.0        V7 0.4174931
#> 97     10        1.0        V8 0.4011540
#> 98     10        1.0        V2 0.4020461
#> 99     10        1.0        V6 0.3923288
#> 100    10        1.0       V10 0.3715245

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3655479
#> 2       1        0.1        V5 0.4177900
#> 3       1        0.1        V9 0.3941529
#> 4       1        0.1        V3 0.3937916
#> 5       1        0.1        V4 0.3667722
#> 6       1        0.1        V8 0.4068546
#> 7       1        0.1        V2 0.4071452
#> 8       1        0.1        V6 0.3971648
#> 9       1        0.1        V7 0.4236842
#> 10      1        0.1       V10 0.3766876
#> 11      2        0.2        V7 0.4229918
#> 12      2        0.2        V8 0.4062172
#> 13      2        0.2        V9 0.3935735
#> 14      2        0.2       V10 0.3761104
#> 15      2        0.2        V1 0.3649894
#> 16      2        0.2        V5 0.4171644
#> 17      2        0.2        V2 0.4065754
#> 18      2        0.2        V3 0.3931527
#> 19      2        0.2        V4 0.3662292
#> 20      2        0.2        V6 0.3966245
#> 21      3        0.3        V4 0.3656870
#> 22      3        0.3        V5 0.4165396
#> 23      3        0.3        V3 0.3925149
#> 24      3        0.3        V7 0.4223005
#> 25      3        0.3        V8 0.4055808
#> 26      3        0.3        V9 0.3929949
#> 27      3        0.3        V6 0.3960850
#> 28      3        0.3       V10 0.3755341
#> 29      3        0.3        V1 0.3644318
#> 30      3        0.3        V2 0.4060065
#> 31      4        0.4        V1 0.3638751
#> 32      4        0.4        V9 0.3924172
#> 33      4        0.4        V3 0.3918781
#> 34      4        0.4        V4 0.3651456
#> 35      4        0.4        V5 0.4159158
#> 36      4        0.4        V2 0.4054383
#> 37      4        0.4        V6 0.3955462
#> 38      4        0.4        V7 0.4216104
#> 39      4        0.4        V8 0.4049455
#> 40      4        0.4       V10 0.3749587
#> 41      5        0.5        V8 0.4043111
#> 42      5        0.5        V9 0.3918403
#> 43      5        0.5       V10 0.3743841
#> 44      5        0.5        V1 0.3633192
#> 45      5        0.5        V5 0.4152929
#> 46      5        0.5        V2 0.4048710
#> 47      5        0.5        V3 0.3912423
#> 48      5        0.5        V4 0.3646050
#> 49      5        0.5        V6 0.3950081
#> 50      5        0.5        V7 0.4209213
#> 51      6        0.6        V4 0.3640653
#> 52      6        0.6        V5 0.4146710
#> 53      6        0.6        V7 0.4202335
#> 54      6        0.6        V8 0.4036777
#> 55      6        0.6        V9 0.3912642
#> 56      6        0.6        V6 0.3944708
#> 57      6        0.6       V10 0.3738104
#> 58      6        0.6        V1 0.3627641
#> 59      6        0.6        V2 0.4043044
#> 60      6        0.6        V3 0.3906076
#> 61      7        0.7        V1 0.3622099
#> 62      7        0.7        V3 0.3899739
#> 63      7        0.7        V4 0.3635263
#> 64      7        0.7        V5 0.4140500
#> 65      7        0.7        V9 0.3906891
#> 66      7        0.7        V6 0.3939342
#> 67      7        0.7        V7 0.4195467
#> 68      7        0.7        V8 0.4030453
#> 69      7        0.7        V2 0.4037386
#> 70      7        0.7       V10 0.3732376
#> 71      8        0.8        V9 0.3901147
#> 72      8        0.8       V10 0.3726657
#> 73      8        0.8        V1 0.3616566
#> 74      8        0.8        V5 0.4134299
#> 75      8        0.8        V2 0.4031736
#> 76      8        0.8        V3 0.3893412
#> 77      8        0.8        V4 0.3629881
#> 78      8        0.8        V8 0.4024138
#> 79      8        0.8        V6 0.3933983
#> 80      8        0.8        V7 0.4188610
#> 81      9        0.9        V5 0.4128108
#> 82      9        0.9        V7 0.4181765
#> 83      9        0.9        V8 0.4017834
#> 84      9        0.9        V9 0.3895412
#> 85      9        0.9        V6 0.3928632
#> 86      9        0.9       V10 0.3720947
#> 87      9        0.9        V1 0.3611041
#> 88      9        0.9        V2 0.4026095
#> 89      9        0.9        V3 0.3887095
#> 90      9        0.9        V4 0.3624507
#> 91     10        1.0        V1 0.3605524
#> 92     10        1.0        V4 0.3619141
#> 93     10        1.0        V5 0.4121925
#> 94     10        1.0        V9 0.3889685
#> 95     10        1.0        V3 0.3880789
#> 96     10        1.0        V7 0.4174931
#> 97     10        1.0        V8 0.4011540
#> 98     10        1.0        V2 0.4020461
#> 99     10        1.0        V6 0.3923288
#> 100    10        1.0       V10 0.3715245


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
#> 1       0                0          0 0.8757906 0.04787800 0.7846128 0.9456288
#> 2       0               20         20 0.8757906 0.04787800 0.7846128 0.9456288
#> 3       0               40         40 0.8757906 0.04787800 0.7846128 0.9456288
#> 4       0               60         60 0.8757906 0.04787800 0.7846128 0.9456288
#> 5      20                0         20 0.8617131 0.05051333 0.7656542 0.9375907
#> 6      20               20         40 0.8617131 0.05051333 0.7656542 0.9375907
#> 7      20               40         60 0.8617131 0.05051333 0.7656542 0.9375907
#> 8      20               60         80 0.8617131 0.05051333 0.7656542 0.9375907
#> 9      40                0         40 0.8478591 0.05295957 0.7473581 0.9294623
#> 10     40               20         60 0.8478591 0.05295957 0.7473581 0.9294623
#> 11     40               40         80 0.8478591 0.05295957 0.7473581 0.9294623
#> 12     40               60        100 0.8478591 0.05295957 0.7473581 0.9294623
#> 13     60                0         60 0.8342249 0.05524017 0.7296626 0.9212686
#> 14     60               20         80 0.8342249 0.05524017 0.7296626 0.9212686
#> 15     60               40        100 0.8342249 0.05524017 0.7296626 0.9212686
#> 16     60               60        120 0.8342249 0.05524017 0.7296626 0.9212686
#> 17     80                0         80 0.8208071 0.05737349 0.7125188 0.9130295
#> 18     80               20        100 0.8208071 0.05737349 0.7125188 0.9130295
#> 19     80               40        120 0.8208071 0.05737349 0.7125188 0.9130295
#> 20     80               60        140 0.8208071 0.05737349 0.7125188 0.9130295
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11544684 0.180311419 0.5800813
#> 2  0.30574618 0.11169561 0.146859239 0.5194376
#> 3  0.26001915 0.10845178 0.118955659 0.4654498
#> 4  0.22113100 0.10484173 0.095727512 0.4175390
#> 5  0.25589195 0.10148564 0.103736008 0.4510916
#> 6  0.21762106 0.09733695 0.083089021 0.4048136
#> 7  0.18507391 0.09329343 0.065997490 0.3638279
#> 8  0.15739448 0.08902525 0.051917830 0.3275472
#> 9  0.18213629 0.08771956 0.056752787 0.3529491
#> 10 0.15489621 0.08343725 0.044339461 0.3179168
#> 11 0.13173012 0.07917246 0.034221800 0.2868939
#> 12 0.11202872 0.07480996 0.026044610 0.2593957
#> 13 0.12963921 0.07487923 0.028832890 0.2786540
#> 14 0.11025053 0.07068687 0.021725129 0.2520855
#> 15 0.09376159 0.06651997 0.016081762 0.2284897
#> 16 0.07973872 0.06234619 0.011662825 0.2074909
#> 17 0.09227334 0.06328239 0.013152054 0.2222071
#> 18 0.07847305 0.05932458 0.009399648 0.2018911
#> 19 0.06673672 0.05541987 0.006545332 0.1837517
#> 20 0.05675565 0.05156713 0.004422506 0.1675060
```
