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
#> 1         0.1 0.3838513 0.02286850 0.3506503 0.4161736
#> 2         0.2 0.3832760 0.02284726 0.3501381 0.4155801
#> 3         0.3 0.3827015 0.02282613 0.3496267 0.4149875
#> 4         0.4 0.3821279 0.02280513 0.3491160 0.4143957
#> 5         0.5 0.3815551 0.02278425 0.3486060 0.4138047
#> 6         0.6 0.3809832 0.02276349 0.3480968 0.4132146
#> 7         0.7 0.3804122 0.02274284 0.3475883 0.4126254
#> 8         0.8 0.3798420 0.02272231 0.3470806 0.4120369
#> 9         0.9 0.3792726 0.02270190 0.3465736 0.4114494
#> 10        1.0 0.3787041 0.02268160 0.3460673 0.4108626

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3835849
#> 2       1        0.1        V5 0.3486127
#> 3       1        0.1        V9 0.3605448
#> 4       1        0.1        V3 0.3981327
#> 5       1        0.1        V4 0.3890166
#> 6       1        0.1        V8 0.3753765
#> 7       1        0.1        V2 0.4128088
#> 8       1        0.1        V6 0.3839240
#> 9       1        0.1        V7 0.3576688
#> 10      1        0.1       V10 0.4171505
#> 11      2        0.2        V7 0.3571670
#> 12      2        0.2        V8 0.3747843
#> 13      2        0.2        V9 0.3599743
#> 14      2        0.2       V10 0.4165564
#> 15      2        0.2        V1 0.3830306
#> 16      2        0.2        V5 0.3480975
#> 17      2        0.2        V2 0.4122173
#> 18      2        0.2        V3 0.3975886
#> 19      2        0.2        V4 0.3884831
#> 20      2        0.2        V6 0.3832321
#> 21      3        0.3        V4 0.3879504
#> 22      3        0.3        V5 0.3475830
#> 23      3        0.3        V3 0.3970453
#> 24      3        0.3        V7 0.3566658
#> 25      3        0.3        V8 0.3741931
#> 26      3        0.3        V9 0.3594048
#> 27      3        0.3        V6 0.3825415
#> 28      3        0.3       V10 0.4159632
#> 29      3        0.3        V1 0.3824771
#> 30      3        0.3        V2 0.4116266
#> 31      4        0.4        V1 0.3819244
#> 32      4        0.4        V9 0.3588362
#> 33      4        0.4        V3 0.3965027
#> 34      4        0.4        V4 0.3874184
#> 35      4        0.4        V5 0.3470694
#> 36      4        0.4        V2 0.4110368
#> 37      4        0.4        V6 0.3818521
#> 38      4        0.4        V7 0.3561654
#> 39      4        0.4        V8 0.3736027
#> 40      4        0.4       V10 0.4153709
#> 41      5        0.5        V8 0.3730134
#> 42      5        0.5        V9 0.3582685
#> 43      5        0.5       V10 0.4147793
#> 44      5        0.5        V1 0.3813725
#> 45      5        0.5        V5 0.3465564
#> 46      5        0.5        V2 0.4104478
#> 47      5        0.5        V3 0.3959608
#> 48      5        0.5        V4 0.3868871
#> 49      5        0.5        V6 0.3811639
#> 50      5        0.5        V7 0.3556657
#> 51      6        0.6        V4 0.3863566
#> 52      6        0.6        V5 0.3460442
#> 53      6        0.6        V7 0.3551666
#> 54      6        0.6        V8 0.3724249
#> 55      6        0.6        V9 0.3577017
#> 56      6        0.6        V6 0.3804770
#> 57      6        0.6       V10 0.4141887
#> 58      6        0.6        V1 0.3808214
#> 59      6        0.6        V2 0.4098596
#> 60      6        0.6        V3 0.3954197
#> 61      7        0.7        V1 0.3802711
#> 62      7        0.7        V3 0.3948793
#> 63      7        0.7        V4 0.3858267
#> 64      7        0.7        V5 0.3455328
#> 65      7        0.7        V9 0.3571357
#> 66      7        0.7        V6 0.3797914
#> 67      7        0.7        V7 0.3546683
#> 68      7        0.7        V8 0.3718374
#> 69      7        0.7        V2 0.4092723
#> 70      7        0.7       V10 0.4135988
#> 71      8        0.8        V9 0.3565707
#> 72      8        0.8       V10 0.4130099
#> 73      8        0.8        V1 0.3797216
#> 74      8        0.8        V5 0.3450222
#> 75      8        0.8        V2 0.4086858
#> 76      8        0.8        V3 0.3943397
#> 77      8        0.8        V4 0.3852977
#> 78      8        0.8        V8 0.3712508
#> 79      8        0.8        V6 0.3791069
#> 80      8        0.8        V7 0.3541707
#> 81      9        0.9        V5 0.3445122
#> 82      9        0.9        V7 0.3536738
#> 83      9        0.9        V8 0.3706651
#> 84      9        0.9        V9 0.3560066
#> 85      9        0.9        V6 0.3784237
#> 86      9        0.9       V10 0.4124217
#> 87      9        0.9        V1 0.3791729
#> 88      9        0.9        V2 0.4081002
#> 89      9        0.9        V3 0.3938008
#> 90      9        0.9        V4 0.3847693
#> 91     10        1.0        V1 0.3786250
#> 92     10        1.0        V4 0.3842417
#> 93     10        1.0        V5 0.3440031
#> 94     10        1.0        V9 0.3554433
#> 95     10        1.0        V3 0.3932626
#> 96     10        1.0        V7 0.3531775
#> 97     10        1.0        V8 0.3700804
#> 98     10        1.0        V2 0.4075154
#> 99     10        1.0        V6 0.3777418
#> 100    10        1.0       V10 0.4118344

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3835849
#> 2       1        0.1        V5 0.3486127
#> 3       1        0.1        V9 0.3605448
#> 4       1        0.1        V3 0.3981327
#> 5       1        0.1        V4 0.3890166
#> 6       1        0.1        V8 0.3753765
#> 7       1        0.1        V2 0.4128088
#> 8       1        0.1        V6 0.3839240
#> 9       1        0.1        V7 0.3576688
#> 10      1        0.1       V10 0.4171505
#> 11      2        0.2        V7 0.3571670
#> 12      2        0.2        V8 0.3747843
#> 13      2        0.2        V9 0.3599743
#> 14      2        0.2       V10 0.4165564
#> 15      2        0.2        V1 0.3830306
#> 16      2        0.2        V5 0.3480975
#> 17      2        0.2        V2 0.4122173
#> 18      2        0.2        V3 0.3975886
#> 19      2        0.2        V4 0.3884831
#> 20      2        0.2        V6 0.3832321
#> 21      3        0.3        V4 0.3879504
#> 22      3        0.3        V5 0.3475830
#> 23      3        0.3        V3 0.3970453
#> 24      3        0.3        V7 0.3566658
#> 25      3        0.3        V8 0.3741931
#> 26      3        0.3        V9 0.3594048
#> 27      3        0.3        V6 0.3825415
#> 28      3        0.3       V10 0.4159632
#> 29      3        0.3        V1 0.3824771
#> 30      3        0.3        V2 0.4116266
#> 31      4        0.4        V1 0.3819244
#> 32      4        0.4        V9 0.3588362
#> 33      4        0.4        V3 0.3965027
#> 34      4        0.4        V4 0.3874184
#> 35      4        0.4        V5 0.3470694
#> 36      4        0.4        V2 0.4110368
#> 37      4        0.4        V6 0.3818521
#> 38      4        0.4        V7 0.3561654
#> 39      4        0.4        V8 0.3736027
#> 40      4        0.4       V10 0.4153709
#> 41      5        0.5        V8 0.3730134
#> 42      5        0.5        V9 0.3582685
#> 43      5        0.5       V10 0.4147793
#> 44      5        0.5        V1 0.3813725
#> 45      5        0.5        V5 0.3465564
#> 46      5        0.5        V2 0.4104478
#> 47      5        0.5        V3 0.3959608
#> 48      5        0.5        V4 0.3868871
#> 49      5        0.5        V6 0.3811639
#> 50      5        0.5        V7 0.3556657
#> 51      6        0.6        V4 0.3863566
#> 52      6        0.6        V5 0.3460442
#> 53      6        0.6        V7 0.3551666
#> 54      6        0.6        V8 0.3724249
#> 55      6        0.6        V9 0.3577017
#> 56      6        0.6        V6 0.3804770
#> 57      6        0.6       V10 0.4141887
#> 58      6        0.6        V1 0.3808214
#> 59      6        0.6        V2 0.4098596
#> 60      6        0.6        V3 0.3954197
#> 61      7        0.7        V1 0.3802711
#> 62      7        0.7        V3 0.3948793
#> 63      7        0.7        V4 0.3858267
#> 64      7        0.7        V5 0.3455328
#> 65      7        0.7        V9 0.3571357
#> 66      7        0.7        V6 0.3797914
#> 67      7        0.7        V7 0.3546683
#> 68      7        0.7        V8 0.3718374
#> 69      7        0.7        V2 0.4092723
#> 70      7        0.7       V10 0.4135988
#> 71      8        0.8        V9 0.3565707
#> 72      8        0.8       V10 0.4130099
#> 73      8        0.8        V1 0.3797216
#> 74      8        0.8        V5 0.3450222
#> 75      8        0.8        V2 0.4086858
#> 76      8        0.8        V3 0.3943397
#> 77      8        0.8        V4 0.3852977
#> 78      8        0.8        V8 0.3712508
#> 79      8        0.8        V6 0.3791069
#> 80      8        0.8        V7 0.3541707
#> 81      9        0.9        V5 0.3445122
#> 82      9        0.9        V7 0.3536738
#> 83      9        0.9        V8 0.3706651
#> 84      9        0.9        V9 0.3560066
#> 85      9        0.9        V6 0.3784237
#> 86      9        0.9       V10 0.4124217
#> 87      9        0.9        V1 0.3791729
#> 88      9        0.9        V2 0.4081002
#> 89      9        0.9        V3 0.3938008
#> 90      9        0.9        V4 0.3847693
#> 91     10        1.0        V1 0.3786250
#> 92     10        1.0        V4 0.3842417
#> 93     10        1.0        V5 0.3440031
#> 94     10        1.0        V9 0.3554433
#> 95     10        1.0        V3 0.3932626
#> 96     10        1.0        V7 0.3531775
#> 97     10        1.0        V8 0.3700804
#> 98     10        1.0        V2 0.4075154
#> 99     10        1.0        V6 0.3777418
#> 100    10        1.0       V10 0.4118344


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
#> 1       0                0          0 0.8757906 0.04654429 0.7859138 0.9393896
#> 2       0               20         20 0.8757906 0.04654429 0.7859138 0.9393896
#> 3       0               40         40 0.8757906 0.04654429 0.7859138 0.9393896
#> 4       0               60         60 0.8757906 0.04654429 0.7859138 0.9393896
#> 5      20                0         20 0.8617131 0.04852524 0.7695931 0.9293682
#> 6      20               20         40 0.8617131 0.04852524 0.7695931 0.9293682
#> 7      20               40         60 0.8617131 0.04852524 0.7695931 0.9293682
#> 8      20               60         80 0.8617131 0.04852524 0.7695931 0.9293682
#> 9      40                0         40 0.8478591 0.05044822 0.7537689 0.9192239
#> 10     40               20         60 0.8478591 0.05044822 0.7537689 0.9192239
#> 11     40               40         80 0.8478591 0.05044822 0.7537689 0.9192239
#> 12     40               60        100 0.8478591 0.05044822 0.7537689 0.9192239
#> 13     60                0         60 0.8342249 0.05230956 0.7383997 0.9089966
#> 14     60               20         80 0.8342249 0.05230956 0.7383997 0.9089966
#> 15     60               40        100 0.8342249 0.05230956 0.7383997 0.9089966
#> 16     60               60        120 0.8342249 0.05230956 0.7383997 0.9089966
#> 17     80                0         80 0.8208071 0.05410729 0.7234517 0.8987182
#> 18     80               20        100 0.8208071 0.05410729 0.7234517 0.8987182
#> 19     80               40        120 0.8208071 0.05410729 0.7234517 0.8987182
#> 20     80               60        140 0.8208071 0.05410729 0.7234517 0.8987182
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12433919 0.169548916 0.6098757
#> 2  0.30574618 0.11728401 0.129179280 0.5506996
#> 3  0.26001915 0.11091910 0.097388919 0.4974257
#> 4  0.22113100 0.10492944 0.072476125 0.4496399
#> 5  0.25589195 0.10690950 0.099593297 0.4718378
#> 6  0.21762106 0.09976278 0.074198255 0.4267296
#> 7  0.18507391 0.09340366 0.054429005 0.3863864
#> 8  0.15739448 0.08757069 0.039188357 0.3503334
#> 9  0.18213629 0.09073938 0.055790197 0.3670696
#> 10 0.15489621 0.08424243 0.040231850 0.3330741
#> 11 0.13173012 0.07844971 0.028376634 0.3026947
#> 12 0.11202872 0.07317043 0.019487968 0.2755308
#> 13 0.12963921 0.07650263 0.029182941 0.2881446
#> 14 0.11025053 0.07084755 0.020086971 0.2625114
#> 15 0.09376159 0.06577979 0.013394631 0.2395527
#> 16 0.07973872 0.06116531 0.008596343 0.2189573
#> 17 0.09227334 0.06421442 0.013840795 0.2285305
#> 18 0.07847305 0.05939052 0.008911540 0.2090553
#> 19 0.06673672 0.05504553 0.005481752 0.1915308
#> 20 0.05675565 0.05108588 0.003190756 0.1757225
```
