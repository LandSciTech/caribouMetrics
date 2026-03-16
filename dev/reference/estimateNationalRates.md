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
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
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
#> 1         0.1 0.3838513 0.02004384 0.3563944 0.4187537
#> 2         0.2 0.3832760 0.02001515 0.3558887 0.4181397
#> 3         0.3 0.3827015 0.01998664 0.3553836 0.4175265
#> 4         0.4 0.3821279 0.01995829 0.3548793 0.4169143
#> 5         0.5 0.3815551 0.01993011 0.3543757 0.4163030
#> 6         0.6 0.3809832 0.01990211 0.3538729 0.4156926
#> 7         0.7 0.3804122 0.01987426 0.3533707 0.4150830
#> 8         0.8 0.3798420 0.01984659 0.3528693 0.4144744
#> 9         0.9 0.3792726 0.01981908 0.3523685 0.4138666
#> 10        1.0 0.3787041 0.01979174 0.3518685 0.4132598

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4198458
#> 2       1        0.1        V5 0.3726619
#> 3       1        0.1        V9 0.3938032
#> 4       1        0.1        V3 0.4149920
#> 5       1        0.1        V4 0.3855211
#> 6       1        0.1        V8 0.3958997
#> 7       1        0.1        V2 0.3934441
#> 8       1        0.1        V6 0.3741467
#> 9       1        0.1        V7 0.3866675
#> 10      1        0.1       V10 0.3516716
#> 11      2        0.2        V7 0.3860923
#> 12      2        0.2        V8 0.3952709
#> 13      2        0.2        V9 0.3932479
#> 14      2        0.2       V10 0.3511886
#> 15      2        0.2        V1 0.4192315
#> 16      2        0.2        V5 0.3720778
#> 17      2        0.2        V2 0.3927417
#> 18      2        0.2        V3 0.4143790
#> 19      2        0.2        V4 0.3849361
#> 20      2        0.2        V6 0.3734950
#> 21      3        0.3        V4 0.3843520
#> 22      3        0.3        V5 0.3714946
#> 23      3        0.3        V3 0.4137669
#> 24      3        0.3        V7 0.3855179
#> 25      3        0.3        V8 0.3946431
#> 26      3        0.3        V9 0.3926934
#> 27      3        0.3        V6 0.3728443
#> 28      3        0.3       V10 0.3507063
#> 29      3        0.3        V1 0.4186181
#> 30      3        0.3        V2 0.3920404
#> 31      4        0.4        V1 0.4180055
#> 32      4        0.4        V9 0.3921397
#> 33      4        0.4        V3 0.4131557
#> 34      4        0.4        V4 0.3837688
#> 35      4        0.4        V5 0.3709124
#> 36      4        0.4        V2 0.3913405
#> 37      4        0.4        V6 0.3721949
#> 38      4        0.4        V7 0.3849444
#> 39      4        0.4        V8 0.3940163
#> 40      4        0.4       V10 0.3502246
#> 41      5        0.5        V8 0.3933905
#> 42      5        0.5        V9 0.3915867
#> 43      5        0.5       V10 0.3497436
#> 44      5        0.5        V1 0.4173939
#> 45      5        0.5        V5 0.3703310
#> 46      5        0.5        V2 0.3906418
#> 47      5        0.5        V3 0.4125454
#> 48      5        0.5        V4 0.3831864
#> 49      5        0.5        V6 0.3715465
#> 50      5        0.5        V7 0.3843718
#> 51      6        0.6        V4 0.3826050
#> 52      6        0.6        V5 0.3697506
#> 53      6        0.6        V7 0.3838000
#> 54      6        0.6        V8 0.3927656
#> 55      6        0.6        V9 0.3910346
#> 56      6        0.6        V6 0.3708993
#> 57      6        0.6       V10 0.3492632
#> 58      6        0.6        V1 0.4167832
#> 59      6        0.6        V2 0.3899443
#> 60      6        0.6        V3 0.4119360
#> 61      7        0.7        V1 0.4161734
#> 62      7        0.7        V3 0.4113275
#> 63      7        0.7        V4 0.3820244
#> 64      7        0.7        V5 0.3691711
#> 65      7        0.7        V9 0.3904832
#> 66      7        0.7        V6 0.3702532
#> 67      7        0.7        V7 0.3832290
#> 68      7        0.7        V8 0.3921418
#> 69      7        0.7        V2 0.3892481
#> 70      7        0.7       V10 0.3487835
#> 71      8        0.8        V9 0.3899326
#> 72      8        0.8       V10 0.3483045
#> 73      8        0.8        V1 0.4155644
#> 74      8        0.8        V5 0.3685924
#> 75      8        0.8        V2 0.3885531
#> 76      8        0.8        V3 0.4107199
#> 77      8        0.8        V4 0.3814447
#> 78      8        0.8        V8 0.3915190
#> 79      8        0.8        V6 0.3696082
#> 80      8        0.8        V7 0.3826589
#> 81      9        0.9        V5 0.3680147
#> 82      9        0.9        V7 0.3820897
#> 83      9        0.9        V8 0.3908972
#> 84      9        0.9        V9 0.3893828
#> 85      9        0.9        V6 0.3689643
#> 86      9        0.9       V10 0.3478261
#> 87      9        0.9        V1 0.4149564
#> 88      9        0.9        V2 0.3878594
#> 89      9        0.9        V3 0.4101132
#> 90      9        0.9        V4 0.3808659
#> 91     10        1.0        V1 0.4143492
#> 92     10        1.0        V4 0.3802879
#> 93     10        1.0        V5 0.3674379
#> 94     10        1.0        V9 0.3888337
#> 95     10        1.0        V3 0.4095074
#> 96     10        1.0        V7 0.3815213
#> 97     10        1.0        V8 0.3902763
#> 98     10        1.0        V2 0.3871669
#> 99     10        1.0        V6 0.3683216
#> 100    10        1.0       V10 0.3473484

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4198458
#> 2       1        0.1        V5 0.3726619
#> 3       1        0.1        V9 0.3938032
#> 4       1        0.1        V3 0.4149920
#> 5       1        0.1        V4 0.3855211
#> 6       1        0.1        V8 0.3958997
#> 7       1        0.1        V2 0.3934441
#> 8       1        0.1        V6 0.3741467
#> 9       1        0.1        V7 0.3866675
#> 10      1        0.1       V10 0.3516716
#> 11      2        0.2        V7 0.3860923
#> 12      2        0.2        V8 0.3952709
#> 13      2        0.2        V9 0.3932479
#> 14      2        0.2       V10 0.3511886
#> 15      2        0.2        V1 0.4192315
#> 16      2        0.2        V5 0.3720778
#> 17      2        0.2        V2 0.3927417
#> 18      2        0.2        V3 0.4143790
#> 19      2        0.2        V4 0.3849361
#> 20      2        0.2        V6 0.3734950
#> 21      3        0.3        V4 0.3843520
#> 22      3        0.3        V5 0.3714946
#> 23      3        0.3        V3 0.4137669
#> 24      3        0.3        V7 0.3855179
#> 25      3        0.3        V8 0.3946431
#> 26      3        0.3        V9 0.3926934
#> 27      3        0.3        V6 0.3728443
#> 28      3        0.3       V10 0.3507063
#> 29      3        0.3        V1 0.4186181
#> 30      3        0.3        V2 0.3920404
#> 31      4        0.4        V1 0.4180055
#> 32      4        0.4        V9 0.3921397
#> 33      4        0.4        V3 0.4131557
#> 34      4        0.4        V4 0.3837688
#> 35      4        0.4        V5 0.3709124
#> 36      4        0.4        V2 0.3913405
#> 37      4        0.4        V6 0.3721949
#> 38      4        0.4        V7 0.3849444
#> 39      4        0.4        V8 0.3940163
#> 40      4        0.4       V10 0.3502246
#> 41      5        0.5        V8 0.3933905
#> 42      5        0.5        V9 0.3915867
#> 43      5        0.5       V10 0.3497436
#> 44      5        0.5        V1 0.4173939
#> 45      5        0.5        V5 0.3703310
#> 46      5        0.5        V2 0.3906418
#> 47      5        0.5        V3 0.4125454
#> 48      5        0.5        V4 0.3831864
#> 49      5        0.5        V6 0.3715465
#> 50      5        0.5        V7 0.3843718
#> 51      6        0.6        V4 0.3826050
#> 52      6        0.6        V5 0.3697506
#> 53      6        0.6        V7 0.3838000
#> 54      6        0.6        V8 0.3927656
#> 55      6        0.6        V9 0.3910346
#> 56      6        0.6        V6 0.3708993
#> 57      6        0.6       V10 0.3492632
#> 58      6        0.6        V1 0.4167832
#> 59      6        0.6        V2 0.3899443
#> 60      6        0.6        V3 0.4119360
#> 61      7        0.7        V1 0.4161734
#> 62      7        0.7        V3 0.4113275
#> 63      7        0.7        V4 0.3820244
#> 64      7        0.7        V5 0.3691711
#> 65      7        0.7        V9 0.3904832
#> 66      7        0.7        V6 0.3702532
#> 67      7        0.7        V7 0.3832290
#> 68      7        0.7        V8 0.3921418
#> 69      7        0.7        V2 0.3892481
#> 70      7        0.7       V10 0.3487835
#> 71      8        0.8        V9 0.3899326
#> 72      8        0.8       V10 0.3483045
#> 73      8        0.8        V1 0.4155644
#> 74      8        0.8        V5 0.3685924
#> 75      8        0.8        V2 0.3885531
#> 76      8        0.8        V3 0.4107199
#> 77      8        0.8        V4 0.3814447
#> 78      8        0.8        V8 0.3915190
#> 79      8        0.8        V6 0.3696082
#> 80      8        0.8        V7 0.3826589
#> 81      9        0.9        V5 0.3680147
#> 82      9        0.9        V7 0.3820897
#> 83      9        0.9        V8 0.3908972
#> 84      9        0.9        V9 0.3893828
#> 85      9        0.9        V6 0.3689643
#> 86      9        0.9       V10 0.3478261
#> 87      9        0.9        V1 0.4149564
#> 88      9        0.9        V2 0.3878594
#> 89      9        0.9        V3 0.4101132
#> 90      9        0.9        V4 0.3808659
#> 91     10        1.0        V1 0.4143492
#> 92     10        1.0        V4 0.3802879
#> 93     10        1.0        V5 0.3674379
#> 94     10        1.0        V9 0.3888337
#> 95     10        1.0        V3 0.4095074
#> 96     10        1.0        V7 0.3815213
#> 97     10        1.0        V8 0.3902763
#> 98     10        1.0        V2 0.3871669
#> 99     10        1.0        V6 0.3683216
#> 100    10        1.0       V10 0.3473484


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
#> 1       0                0          0 0.8757906 0.04944015 0.7712459 0.9373200
#> 2       0               20         20 0.8757906 0.04944015 0.7712459 0.9373200
#> 3       0               40         40 0.8757906 0.04944015 0.7712459 0.9373200
#> 4       0               60         60 0.8757906 0.04944015 0.7712459 0.9373200
#> 5      20                0         20 0.8617131 0.05248177 0.7515459 0.9270607
#> 6      20               20         40 0.8617131 0.05248177 0.7515459 0.9270607
#> 7      20               40         60 0.8617131 0.05248177 0.7515459 0.9270607
#> 8      20               60         80 0.8617131 0.05248177 0.7515459 0.9270607
#> 9      40                0         40 0.8478591 0.05533293 0.7325631 0.9166983
#> 10     40               20         60 0.8478591 0.05533293 0.7325631 0.9166983
#> 11     40               40         80 0.8478591 0.05533293 0.7325631 0.9166983
#> 12     40               60        100 0.8478591 0.05533293 0.7325631 0.9166983
#> 13     60                0         60 0.8342249 0.05801029 0.7142302 0.9062707
#> 14     60               20         80 0.8342249 0.05801029 0.7142302 0.9062707
#> 15     60               40        100 0.8342249 0.05801029 0.7142302 0.9062707
#> 16     60               60        120 0.8342249 0.05801029 0.7142302 0.9062707
#> 17     80                0         80 0.8208071 0.06052812 0.6964940 0.8958081
#> 18     80               20        100 0.8208071 0.06052812 0.6964940 0.8958081
#> 19     80               40        120 0.8208071 0.06052812 0.6964940 0.8958081
#> 20     80               60        140 0.8208071 0.06052812 0.6964940 0.8958081
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.13664821 0.1579013268 0.6270161
#> 2  0.30574618 0.14090093 0.1058558731 0.5801266
#> 3  0.26001915 0.14191672 0.0692004058 0.5367902
#> 4  0.22113100 0.14026034 0.0437555330 0.4968590
#> 5  0.25589195 0.12767006 0.0763251887 0.5247591
#> 6  0.21762106 0.12679324 0.0486621030 0.4857893
#> 7  0.18507391 0.12434994 0.0297756211 0.4499704
#> 8  0.15739448 0.12059580 0.0172630648 0.4170847
#> 9  0.18213629 0.11526145 0.0333857635 0.4400521
#> 10 0.15489621 0.11204027 0.0196175043 0.4079829
#> 11 0.13173012 0.10810378 0.0107833440 0.3785609
#> 12 0.11202872 0.10357687 0.0054270527 0.3515701
#> 13 0.12963921 0.10210173 0.0124178986 0.3704193
#> 14 0.11025053 0.09799113 0.0063887628 0.3441004
#> 15 0.09376159 0.09360748 0.0029452523 0.3199477
#> 16 0.07973872 0.08901234 0.0011739298 0.2977692
#> 17 0.09227334 0.08951931 0.0035439641 0.3132606
#> 18 0.07847305 0.08526801 0.0014647757 0.2916255
#> 19 0.06673672 0.08096537 0.0005054121 0.2717358
#> 20 0.05675565 0.07663478 0.0001376191 0.2534290
```
