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
#> 1         0.1 0.3838513 0.02367252 0.3490701 0.4205064
#> 2         0.2 0.3832760 0.02365863 0.3485804 0.4199671
#> 3         0.3 0.3827015 0.02364485 0.3480914 0.4194284
#> 4         0.4 0.3821279 0.02363116 0.3476030 0.4188905
#> 5         0.5 0.3815551 0.02361758 0.3471154 0.4183532
#> 6         0.6 0.3809832 0.02360411 0.3466284 0.4178167
#> 7         0.7 0.3804122 0.02359073 0.3461421 0.4172808
#> 8         0.8 0.3798420 0.02357746 0.3456565 0.4167456
#> 9         0.9 0.3792726 0.02356428 0.3451716 0.4162111
#> 10        1.0 0.3787041 0.02355121 0.3446873 0.4156773

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3507745
#> 2       1        0.1        V5 0.3726688
#> 3       1        0.1        V9 0.3485753
#> 4       1        0.1        V3 0.4239218
#> 5       1        0.1        V4 0.4087424
#> 6       1        0.1        V8 0.3847888
#> 7       1        0.1        V2 0.3681708
#> 8       1        0.1        V6 0.3724478
#> 9       1        0.1        V7 0.3653454
#> 10      1        0.1       V10 0.3825676
#> 11      2        0.2        V7 0.3648839
#> 12      2        0.2        V8 0.3842335
#> 13      2        0.2        V9 0.3480908
#> 14      2        0.2       V10 0.3819834
#> 15      2        0.2        V1 0.3502668
#> 16      2        0.2        V5 0.3720416
#> 17      2        0.2        V2 0.3676238
#> 18      2        0.2        V3 0.4233851
#> 19      2        0.2        V4 0.4081938
#> 20      2        0.2        V6 0.3718304
#> 21      3        0.3        V4 0.4076460
#> 22      3        0.3        V5 0.3714155
#> 23      3        0.3        V3 0.4228491
#> 24      3        0.3        V7 0.3644231
#> 25      3        0.3        V8 0.3836790
#> 26      3        0.3        V9 0.3476070
#> 27      3        0.3        V6 0.3712140
#> 28      3        0.3       V10 0.3814000
#> 29      3        0.3        V1 0.3497597
#> 30      3        0.3        V2 0.3670777
#> 31      4        0.4        V1 0.3492534
#> 32      4        0.4        V9 0.3471239
#> 33      4        0.4        V3 0.4223138
#> 34      4        0.4        V4 0.4070989
#> 35      4        0.4        V5 0.3707904
#> 36      4        0.4        V2 0.3665323
#> 37      4        0.4        V6 0.3705986
#> 38      4        0.4        V7 0.3639628
#> 39      4        0.4        V8 0.3831253
#> 40      4        0.4       V10 0.3808175
#> 41      5        0.5        V8 0.3825725
#> 42      5        0.5        V9 0.3466414
#> 43      5        0.5       V10 0.3802359
#> 44      5        0.5        V1 0.3487479
#> 45      5        0.5        V5 0.3701664
#> 46      5        0.5        V2 0.3659878
#> 47      5        0.5        V3 0.4217792
#> 48      5        0.5        V4 0.4065526
#> 49      5        0.5        V6 0.3699843
#> 50      5        0.5        V7 0.3635031
#> 51      6        0.6        V4 0.4060070
#> 52      6        0.6        V5 0.3695434
#> 53      6        0.6        V7 0.3630440
#> 54      6        0.6        V8 0.3820204
#> 55      6        0.6        V9 0.3461596
#> 56      6        0.6        V6 0.3693709
#> 57      6        0.6       V10 0.3796552
#> 58      6        0.6        V1 0.3482430
#> 59      6        0.6        V2 0.3654440
#> 60      6        0.6        V3 0.4212453
#> 61      7        0.7        V1 0.3477389
#> 62      7        0.7        V3 0.4207120
#> 63      7        0.7        V4 0.4054621
#> 64      7        0.7        V5 0.3689215
#> 65      7        0.7        V9 0.3456785
#> 66      7        0.7        V6 0.3687586
#> 67      7        0.7        V7 0.3625855
#> 68      7        0.7        V8 0.3814691
#> 69      7        0.7        V2 0.3649011
#> 70      7        0.7       V10 0.3790754
#> 71      8        0.8        V9 0.3451980
#> 72      8        0.8       V10 0.3784965
#> 73      8        0.8        V1 0.3472356
#> 74      8        0.8        V5 0.3683006
#> 75      8        0.8        V2 0.3643590
#> 76      8        0.8        V3 0.4201794
#> 77      8        0.8        V4 0.4049179
#> 78      8        0.8        V8 0.3809186
#> 79      8        0.8        V6 0.3681473
#> 80      8        0.8        V7 0.3621275
#> 81      9        0.9        V5 0.3676808
#> 82      9        0.9        V7 0.3616702
#> 83      9        0.9        V8 0.3803689
#> 84      9        0.9        V9 0.3447183
#> 85      9        0.9        V6 0.3675370
#> 86      9        0.9       V10 0.3779184
#> 87      9        0.9        V1 0.3467329
#> 88      9        0.9        V2 0.3638177
#> 89      9        0.9        V3 0.4196475
#> 90      9        0.9        V4 0.4043745
#> 91     10        1.0        V1 0.3462310
#> 92     10        1.0        V4 0.4038318
#> 93     10        1.0        V5 0.3670620
#> 94     10        1.0        V9 0.3442391
#> 95     10        1.0        V3 0.4191163
#> 96     10        1.0        V7 0.3612134
#> 97     10        1.0        V8 0.3798200
#> 98     10        1.0        V2 0.3632772
#> 99     10        1.0        V6 0.3669277
#> 100    10        1.0       V10 0.3773412

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3507745
#> 2       1        0.1        V5 0.3726688
#> 3       1        0.1        V9 0.3485753
#> 4       1        0.1        V3 0.4239218
#> 5       1        0.1        V4 0.4087424
#> 6       1        0.1        V8 0.3847888
#> 7       1        0.1        V2 0.3681708
#> 8       1        0.1        V6 0.3724478
#> 9       1        0.1        V7 0.3653454
#> 10      1        0.1       V10 0.3825676
#> 11      2        0.2        V7 0.3648839
#> 12      2        0.2        V8 0.3842335
#> 13      2        0.2        V9 0.3480908
#> 14      2        0.2       V10 0.3819834
#> 15      2        0.2        V1 0.3502668
#> 16      2        0.2        V5 0.3720416
#> 17      2        0.2        V2 0.3676238
#> 18      2        0.2        V3 0.4233851
#> 19      2        0.2        V4 0.4081938
#> 20      2        0.2        V6 0.3718304
#> 21      3        0.3        V4 0.4076460
#> 22      3        0.3        V5 0.3714155
#> 23      3        0.3        V3 0.4228491
#> 24      3        0.3        V7 0.3644231
#> 25      3        0.3        V8 0.3836790
#> 26      3        0.3        V9 0.3476070
#> 27      3        0.3        V6 0.3712140
#> 28      3        0.3       V10 0.3814000
#> 29      3        0.3        V1 0.3497597
#> 30      3        0.3        V2 0.3670777
#> 31      4        0.4        V1 0.3492534
#> 32      4        0.4        V9 0.3471239
#> 33      4        0.4        V3 0.4223138
#> 34      4        0.4        V4 0.4070989
#> 35      4        0.4        V5 0.3707904
#> 36      4        0.4        V2 0.3665323
#> 37      4        0.4        V6 0.3705986
#> 38      4        0.4        V7 0.3639628
#> 39      4        0.4        V8 0.3831253
#> 40      4        0.4       V10 0.3808175
#> 41      5        0.5        V8 0.3825725
#> 42      5        0.5        V9 0.3466414
#> 43      5        0.5       V10 0.3802359
#> 44      5        0.5        V1 0.3487479
#> 45      5        0.5        V5 0.3701664
#> 46      5        0.5        V2 0.3659878
#> 47      5        0.5        V3 0.4217792
#> 48      5        0.5        V4 0.4065526
#> 49      5        0.5        V6 0.3699843
#> 50      5        0.5        V7 0.3635031
#> 51      6        0.6        V4 0.4060070
#> 52      6        0.6        V5 0.3695434
#> 53      6        0.6        V7 0.3630440
#> 54      6        0.6        V8 0.3820204
#> 55      6        0.6        V9 0.3461596
#> 56      6        0.6        V6 0.3693709
#> 57      6        0.6       V10 0.3796552
#> 58      6        0.6        V1 0.3482430
#> 59      6        0.6        V2 0.3654440
#> 60      6        0.6        V3 0.4212453
#> 61      7        0.7        V1 0.3477389
#> 62      7        0.7        V3 0.4207120
#> 63      7        0.7        V4 0.4054621
#> 64      7        0.7        V5 0.3689215
#> 65      7        0.7        V9 0.3456785
#> 66      7        0.7        V6 0.3687586
#> 67      7        0.7        V7 0.3625855
#> 68      7        0.7        V8 0.3814691
#> 69      7        0.7        V2 0.3649011
#> 70      7        0.7       V10 0.3790754
#> 71      8        0.8        V9 0.3451980
#> 72      8        0.8       V10 0.3784965
#> 73      8        0.8        V1 0.3472356
#> 74      8        0.8        V5 0.3683006
#> 75      8        0.8        V2 0.3643590
#> 76      8        0.8        V3 0.4201794
#> 77      8        0.8        V4 0.4049179
#> 78      8        0.8        V8 0.3809186
#> 79      8        0.8        V6 0.3681473
#> 80      8        0.8        V7 0.3621275
#> 81      9        0.9        V5 0.3676808
#> 82      9        0.9        V7 0.3616702
#> 83      9        0.9        V8 0.3803689
#> 84      9        0.9        V9 0.3447183
#> 85      9        0.9        V6 0.3675370
#> 86      9        0.9       V10 0.3779184
#> 87      9        0.9        V1 0.3467329
#> 88      9        0.9        V2 0.3638177
#> 89      9        0.9        V3 0.4196475
#> 90      9        0.9        V4 0.4043745
#> 91     10        1.0        V1 0.3462310
#> 92     10        1.0        V4 0.4038318
#> 93     10        1.0        V5 0.3670620
#> 94     10        1.0        V9 0.3442391
#> 95     10        1.0        V3 0.4191163
#> 96     10        1.0        V7 0.3612134
#> 97     10        1.0        V8 0.3798200
#> 98     10        1.0        V2 0.3632772
#> 99     10        1.0        V6 0.3669277
#> 100    10        1.0       V10 0.3773412


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
#> 1       0                0          0 0.8757906 0.04869839 0.7671907 0.9328313
#> 2       0               20         20 0.8757906 0.04869839 0.7671907 0.9328313
#> 3       0               40         40 0.8757906 0.04869839 0.7671907 0.9328313
#> 4       0               60         60 0.8757906 0.04869839 0.7671907 0.9328313
#> 5      20                0         20 0.8617131 0.05009907 0.7526973 0.9256968
#> 6      20               20         40 0.8617131 0.05009907 0.7526973 0.9256968
#> 7      20               40         60 0.8617131 0.05009907 0.7526973 0.9256968
#> 8      20               60         80 0.8617131 0.05009907 0.7526973 0.9256968
#> 9      40                0         40 0.8478591 0.05150970 0.7386195 0.9185208
#> 10     40               20         60 0.8478591 0.05150970 0.7386195 0.9185208
#> 11     40               40         80 0.8478591 0.05150970 0.7386195 0.9185208
#> 12     40               60        100 0.8478591 0.05150970 0.7386195 0.9185208
#> 13     60                0         60 0.8342249 0.05292320 0.7249242 0.9113155
#> 14     60               20         80 0.8342249 0.05292320 0.7249242 0.9113155
#> 15     60               40        100 0.8342249 0.05292320 0.7249242 0.9113155
#> 16     60               60        120 0.8342249 0.05292320 0.7249242 0.9113155
#> 17     80                0         80 0.8208071 0.05433348 0.7115841 0.9040909
#> 18     80               20        100 0.8208071 0.05433348 0.7115841 0.9040909
#> 19     80               40        120 0.8208071 0.05433348 0.7115841 0.9040909
#> 20     80               60        140 0.8208071 0.05433348 0.7115841 0.9040909
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12514330 0.1543440025 0.5817991
#> 2  0.30574618 0.12352916 0.1054411310 0.5173564
#> 3  0.26001915 0.11963682 0.0704076374 0.4604436
#> 4  0.22113100 0.11418284 0.0456345333 0.4103558
#> 5  0.25589195 0.11041621 0.0775186141 0.4458093
#> 6  0.21762106 0.10574699 0.0506280326 0.3974936
#> 7  0.18507391 0.10017071 0.0318865925 0.3550561
#> 8  0.15739448 0.09402879 0.0191612370 0.3177903
#> 9  0.18213629 0.09424469 0.0356362104 0.3441640
#> 10 0.15489621 0.08860831 0.0216735299 0.3082227
#> 11 0.13173012 0.08273077 0.0124464809 0.2766311
#> 12 0.11202872 0.07677950 0.0066326931 0.2488194
#> 13 0.12963921 0.07887101 0.0142430423 0.2685115
#> 14 0.11025053 0.07325228 0.0077371814 0.2416616
#> 15 0.09376159 0.06771796 0.0038341929 0.2179543
#> 16 0.07973872 0.06234725 0.0016848961 0.1969590
#> 17 0.09227334 0.06516276 0.0045569352 0.2118388
#> 18 0.07847305 0.06001094 0.0020652651 0.1915309
#> 19 0.06673672 0.05509245 0.0008076495 0.1734613
#> 20 0.05675565 0.05043627 0.0002607724 0.1573137
```
