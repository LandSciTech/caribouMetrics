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
#> 1         0.1 0.3838513 0.01814506 0.3600739 0.4130462
#> 2         0.2 0.3832760 0.01807532 0.3595950 0.4123285
#> 3         0.3 0.3827015 0.01800599 0.3591168 0.4116121
#> 4         0.4 0.3821279 0.01793707 0.3586391 0.4108969
#> 5         0.5 0.3815551 0.01786856 0.3581622 0.4101829
#> 6         0.6 0.3809832 0.01780047 0.3576858 0.4094703
#> 7         0.7 0.3804122 0.01773278 0.3572101 0.4087588
#> 8         0.8 0.3798420 0.01766550 0.3567350 0.4080486
#> 9         0.9 0.3792726 0.01759864 0.3562606 0.4073396
#> 10        1.0 0.3787041 0.01753218 0.3557868 0.4066318

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3774298
#> 2       1        0.1        V5 0.3695638
#> 3       1        0.1        V9 0.3838601
#> 4       1        0.1        V3 0.4016244
#> 5       1        0.1        V4 0.3573188
#> 6       1        0.1        V8 0.4143516
#> 7       1        0.1        V2 0.3896553
#> 8       1        0.1        V6 0.4085496
#> 9       1        0.1        V7 0.4010009
#> 10      1        0.1       V10 0.3994451
#> 11      2        0.2        V7 0.4004611
#> 12      2        0.2        V8 0.4136397
#> 13      2        0.2        V9 0.3833164
#> 14      2        0.2       V10 0.3987728
#> 15      2        0.2        V1 0.3767880
#> 16      2        0.2        V5 0.3690420
#> 17      2        0.2        V2 0.3891696
#> 18      2        0.2        V3 0.4010007
#> 19      2        0.2        V4 0.3568523
#> 20      2        0.2        V6 0.4078123
#> 21      3        0.3        V4 0.3563865
#> 22      3        0.3        V5 0.3685209
#> 23      3        0.3        V3 0.4003779
#> 24      3        0.3        V7 0.3999219
#> 25      3        0.3        V8 0.4129289
#> 26      3        0.3        V9 0.3827735
#> 27      3        0.3        V6 0.4070764
#> 28      3        0.3       V10 0.3981016
#> 29      3        0.3        V1 0.3761472
#> 30      3        0.3        V2 0.3886845
#> 31      4        0.4        V1 0.3755075
#> 32      4        0.4        V9 0.3822314
#> 33      4        0.4        V3 0.3997561
#> 34      4        0.4        V4 0.3559213
#> 35      4        0.4        V5 0.3680005
#> 36      4        0.4        V2 0.3881999
#> 37      4        0.4        V6 0.4063417
#> 38      4        0.4        V7 0.3993836
#> 39      4        0.4        V8 0.4122194
#> 40      4        0.4       V10 0.3974315
#> 41      5        0.5        V8 0.4115110
#> 42      5        0.5        V9 0.3816900
#> 43      5        0.5       V10 0.3967626
#> 44      5        0.5        V1 0.3748690
#> 45      5        0.5        V5 0.3674809
#> 46      5        0.5        V2 0.3877160
#> 47      5        0.5        V3 0.3991353
#> 48      5        0.5        V4 0.3554567
#> 49      5        0.5        V6 0.4056084
#> 50      5        0.5        V7 0.3988459
#> 51      6        0.6        V4 0.3549927
#> 52      6        0.6        V5 0.3669620
#> 53      6        0.6        V7 0.3983089
#> 54      6        0.6        V8 0.4108039
#> 55      6        0.6        V9 0.3811494
#> 56      6        0.6        V6 0.4048764
#> 57      6        0.6       V10 0.3960948
#> 58      6        0.6        V1 0.3742315
#> 59      6        0.6        V2 0.3872327
#> 60      6        0.6        V3 0.3985154
#> 61      7        0.7        V1 0.3735951
#> 62      7        0.7        V3 0.3978965
#> 63      7        0.7        V4 0.3545293
#> 64      7        0.7        V5 0.3664438
#> 65      7        0.7        V9 0.3806096
#> 66      7        0.7        V6 0.4041458
#> 67      7        0.7        V7 0.3977727
#> 68      7        0.7        V8 0.4100980
#> 69      7        0.7        V2 0.3867500
#> 70      7        0.7       V10 0.3954281
#> 71      8        0.8        V9 0.3800705
#> 72      8        0.8       V10 0.3947625
#> 73      8        0.8        V1 0.3729597
#> 74      8        0.8        V5 0.3659264
#> 75      8        0.8        V2 0.3862679
#> 76      8        0.8        V3 0.3972786
#> 77      8        0.8        V4 0.3540666
#> 78      8        0.8        V8 0.4093934
#> 79      8        0.8        V6 0.4034165
#> 80      8        0.8        V7 0.3972372
#> 81      9        0.9        V5 0.3654097
#> 82      9        0.9        V7 0.3967024
#> 83      9        0.9        V8 0.4086899
#> 84      9        0.9        V9 0.3795322
#> 85      9        0.9        V6 0.4026884
#> 86      9        0.9       V10 0.3940980
#> 87      9        0.9        V1 0.3723255
#> 88      9        0.9        V2 0.3857864
#> 89      9        0.9        V3 0.3966616
#> 90      9        0.9        V4 0.3536044
#> 91     10        1.0        V1 0.3716923
#> 92     10        1.0        V4 0.3531428
#> 93     10        1.0        V5 0.3648937
#> 94     10        1.0        V9 0.3789946
#> 95     10        1.0        V3 0.3960456
#> 96     10        1.0        V7 0.3961684
#> 97     10        1.0        V8 0.4079877
#> 98     10        1.0        V2 0.3853055
#> 99     10        1.0        V6 0.4019617
#> 100    10        1.0       V10 0.3934347

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3774298
#> 2       1        0.1        V5 0.3695638
#> 3       1        0.1        V9 0.3838601
#> 4       1        0.1        V3 0.4016244
#> 5       1        0.1        V4 0.3573188
#> 6       1        0.1        V8 0.4143516
#> 7       1        0.1        V2 0.3896553
#> 8       1        0.1        V6 0.4085496
#> 9       1        0.1        V7 0.4010009
#> 10      1        0.1       V10 0.3994451
#> 11      2        0.2        V7 0.4004611
#> 12      2        0.2        V8 0.4136397
#> 13      2        0.2        V9 0.3833164
#> 14      2        0.2       V10 0.3987728
#> 15      2        0.2        V1 0.3767880
#> 16      2        0.2        V5 0.3690420
#> 17      2        0.2        V2 0.3891696
#> 18      2        0.2        V3 0.4010007
#> 19      2        0.2        V4 0.3568523
#> 20      2        0.2        V6 0.4078123
#> 21      3        0.3        V4 0.3563865
#> 22      3        0.3        V5 0.3685209
#> 23      3        0.3        V3 0.4003779
#> 24      3        0.3        V7 0.3999219
#> 25      3        0.3        V8 0.4129289
#> 26      3        0.3        V9 0.3827735
#> 27      3        0.3        V6 0.4070764
#> 28      3        0.3       V10 0.3981016
#> 29      3        0.3        V1 0.3761472
#> 30      3        0.3        V2 0.3886845
#> 31      4        0.4        V1 0.3755075
#> 32      4        0.4        V9 0.3822314
#> 33      4        0.4        V3 0.3997561
#> 34      4        0.4        V4 0.3559213
#> 35      4        0.4        V5 0.3680005
#> 36      4        0.4        V2 0.3881999
#> 37      4        0.4        V6 0.4063417
#> 38      4        0.4        V7 0.3993836
#> 39      4        0.4        V8 0.4122194
#> 40      4        0.4       V10 0.3974315
#> 41      5        0.5        V8 0.4115110
#> 42      5        0.5        V9 0.3816900
#> 43      5        0.5       V10 0.3967626
#> 44      5        0.5        V1 0.3748690
#> 45      5        0.5        V5 0.3674809
#> 46      5        0.5        V2 0.3877160
#> 47      5        0.5        V3 0.3991353
#> 48      5        0.5        V4 0.3554567
#> 49      5        0.5        V6 0.4056084
#> 50      5        0.5        V7 0.3988459
#> 51      6        0.6        V4 0.3549927
#> 52      6        0.6        V5 0.3669620
#> 53      6        0.6        V7 0.3983089
#> 54      6        0.6        V8 0.4108039
#> 55      6        0.6        V9 0.3811494
#> 56      6        0.6        V6 0.4048764
#> 57      6        0.6       V10 0.3960948
#> 58      6        0.6        V1 0.3742315
#> 59      6        0.6        V2 0.3872327
#> 60      6        0.6        V3 0.3985154
#> 61      7        0.7        V1 0.3735951
#> 62      7        0.7        V3 0.3978965
#> 63      7        0.7        V4 0.3545293
#> 64      7        0.7        V5 0.3664438
#> 65      7        0.7        V9 0.3806096
#> 66      7        0.7        V6 0.4041458
#> 67      7        0.7        V7 0.3977727
#> 68      7        0.7        V8 0.4100980
#> 69      7        0.7        V2 0.3867500
#> 70      7        0.7       V10 0.3954281
#> 71      8        0.8        V9 0.3800705
#> 72      8        0.8       V10 0.3947625
#> 73      8        0.8        V1 0.3729597
#> 74      8        0.8        V5 0.3659264
#> 75      8        0.8        V2 0.3862679
#> 76      8        0.8        V3 0.3972786
#> 77      8        0.8        V4 0.3540666
#> 78      8        0.8        V8 0.4093934
#> 79      8        0.8        V6 0.4034165
#> 80      8        0.8        V7 0.3972372
#> 81      9        0.9        V5 0.3654097
#> 82      9        0.9        V7 0.3967024
#> 83      9        0.9        V8 0.4086899
#> 84      9        0.9        V9 0.3795322
#> 85      9        0.9        V6 0.4026884
#> 86      9        0.9       V10 0.3940980
#> 87      9        0.9        V1 0.3723255
#> 88      9        0.9        V2 0.3857864
#> 89      9        0.9        V3 0.3966616
#> 90      9        0.9        V4 0.3536044
#> 91     10        1.0        V1 0.3716923
#> 92     10        1.0        V4 0.3531428
#> 93     10        1.0        V5 0.3648937
#> 94     10        1.0        V9 0.3789946
#> 95     10        1.0        V3 0.3960456
#> 96     10        1.0        V7 0.3961684
#> 97     10        1.0        V8 0.4079877
#> 98     10        1.0        V2 0.3853055
#> 99     10        1.0        V6 0.4019617
#> 100    10        1.0       V10 0.3934347


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
#> 1       0                0          0 0.8757906 0.04831678 0.7758686 0.9458652
#> 2       0               20         20 0.8757906 0.04831678 0.7758686 0.9458652
#> 3       0               40         40 0.8757906 0.04831678 0.7758686 0.9458652
#> 4       0               60         60 0.8757906 0.04831678 0.7758686 0.9458652
#> 5      20                0         20 0.8617131 0.04964182 0.7635777 0.9364635
#> 6      20               20         40 0.8617131 0.04964182 0.7635777 0.9364635
#> 7      20               40         60 0.8617131 0.04964182 0.7635777 0.9364635
#> 8      20               60         80 0.8617131 0.04964182 0.7635777 0.9364635
#> 9      40                0         40 0.8478591 0.05097531 0.7515681 0.9269477
#> 10     40               20         60 0.8478591 0.05097531 0.7515681 0.9269477
#> 11     40               40         80 0.8478591 0.05097531 0.7515681 0.9269477
#> 12     40               60        100 0.8478591 0.05097531 0.7515681 0.9269477
#> 13     60                0         60 0.8342249 0.05230785 0.7398221 0.9173552
#> 14     60               20         80 0.8342249 0.05230785 0.7398221 0.9173552
#> 15     60               40        100 0.8342249 0.05230785 0.7398221 0.9173552
#> 16     60               60        120 0.8342249 0.05230785 0.7398221 0.9173552
#> 17     80                0         80 0.8208071 0.05363200 0.7283247 0.9077150
#> 18     80               20        100 0.8208071 0.05363200 0.7283247 0.9077150
#> 19     80               40        120 0.8208071 0.05363200 0.7283247 0.9077150
#> 20     80               60        140 0.8208071 0.05363200 0.7283247 0.9077150
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12692828 0.161069270 0.5556504
#> 2  0.30574618 0.11895222 0.130359573 0.5023423
#> 3  0.26001915 0.11145740 0.104845887 0.4544863
#> 4  0.22113100 0.10433013 0.083706780 0.4116226
#> 5  0.25589195 0.11759703 0.093411336 0.4603313
#> 6  0.21762106 0.10916708 0.074259438 0.4168542
#> 7  0.18507391 0.10136807 0.058490283 0.3779561
#> 8  0.15739448 0.09410038 0.045579134 0.3431734
#> 9  0.18213629 0.10588518 0.051484824 0.3827023
#> 10 0.15489621 0.09772548 0.039873090 0.3474173
#> 11 0.13173012 0.09023671 0.030471466 0.3158641
#> 12 0.11202872 0.08332802 0.022930781 0.2876328
#> 13 0.12963921 0.09358879 0.026358616 0.3197145
#> 14 0.11025053 0.08605908 0.019660180 0.2910792
#> 15 0.09376159 0.07917783 0.014384859 0.2654367
#> 16 0.07973872 0.07286401 0.010293176 0.2424429
#> 17 0.09227334 0.08169069 0.012134421 0.2685685
#> 18 0.07847305 0.07493961 0.008571673 0.2452533
#> 19 0.06673672 0.06878457 0.005888244 0.2243144
#> 20 0.05675565 0.06315652 0.003915833 0.2054700
```
