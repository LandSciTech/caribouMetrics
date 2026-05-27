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
#> 1         0.1 0.3838513 0.02503591 0.3503906 0.4288030
#> 2         0.2 0.3832760 0.02500702 0.3498236 0.4281644
#> 3         0.3 0.3827015 0.02497824 0.3492574 0.4275268
#> 4         0.4 0.3821279 0.02494955 0.3486922 0.4268902
#> 5         0.5 0.3815551 0.02492097 0.3481279 0.4262545
#> 6         0.6 0.3809832 0.02489248 0.3475646 0.4256197
#> 7         0.7 0.3804122 0.02486410 0.3470021 0.4249859
#> 8         0.8 0.3798420 0.02483581 0.3464405 0.4243530
#> 9         0.9 0.3792726 0.02480762 0.3458799 0.4237211
#> 10        1.0 0.3787041 0.02477953 0.3453202 0.4230901

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3845287
#> 2       1        0.1        V5 0.3483139
#> 3       1        0.1        V9 0.3575436
#> 4       1        0.1        V3 0.4364745
#> 5       1        0.1        V4 0.3809664
#> 6       1        0.1        V8 0.4019067
#> 7       1        0.1        V2 0.3787148
#> 8       1        0.1        V6 0.4023789
#> 9       1        0.1        V7 0.3802772
#> 10      1        0.1       V10 0.3696992
#> 11      2        0.2        V7 0.3796861
#> 12      2        0.2        V8 0.4012900
#> 13      2        0.2        V9 0.3569939
#> 14      2        0.2       V10 0.3691956
#> 15      2        0.2        V1 0.3839354
#> 16      2        0.2        V5 0.3477418
#> 17      2        0.2        V2 0.3780685
#> 18      2        0.2        V3 0.4358119
#> 19      2        0.2        V4 0.3804455
#> 20      2        0.2        V6 0.4018232
#> 21      3        0.3        V4 0.3799253
#> 22      3        0.3        V5 0.3471707
#> 23      3        0.3        V3 0.4351503
#> 24      3        0.3        V7 0.3790959
#> 25      3        0.3        V8 0.4006743
#> 26      3        0.3        V9 0.3564451
#> 27      3        0.3        V6 0.4012682
#> 28      3        0.3       V10 0.3686927
#> 29      3        0.3        V1 0.3833431
#> 30      3        0.3        V2 0.3774232
#> 31      4        0.4        V1 0.3827517
#> 32      4        0.4        V9 0.3558971
#> 33      4        0.4        V3 0.4344897
#> 34      4        0.4        V4 0.3794058
#> 35      4        0.4        V5 0.3466005
#> 36      4        0.4        V2 0.3767791
#> 37      4        0.4        V6 0.4007141
#> 38      4        0.4        V7 0.3785066
#> 39      4        0.4        V8 0.4000594
#> 40      4        0.4       V10 0.3681905
#> 41      5        0.5        V8 0.3994456
#> 42      5        0.5        V9 0.3553500
#> 43      5        0.5       V10 0.3676889
#> 44      5        0.5        V1 0.3821611
#> 45      5        0.5        V5 0.3460312
#> 46      5        0.5        V2 0.3761360
#> 47      5        0.5        V3 0.4338301
#> 48      5        0.5        V4 0.3788870
#> 49      5        0.5        V6 0.4001607
#> 50      5        0.5        V7 0.3779182
#> 51      6        0.6        V4 0.3783689
#> 52      6        0.6        V5 0.3454629
#> 53      6        0.6        V7 0.3773308
#> 54      6        0.6        V8 0.3988326
#> 55      6        0.6        V9 0.3548037
#> 56      6        0.6        V6 0.3996081
#> 57      6        0.6       V10 0.3671881
#> 58      6        0.6        V1 0.3815715
#> 59      6        0.6        V2 0.3754941
#> 60      6        0.6        V3 0.4331715
#> 61      7        0.7        V1 0.3809828
#> 62      7        0.7        V3 0.4325139
#> 63      7        0.7        V4 0.3778516
#> 64      7        0.7        V5 0.3448955
#> 65      7        0.7        V9 0.3542582
#> 66      7        0.7        V6 0.3990562
#> 67      7        0.7        V7 0.3767442
#> 68      7        0.7        V8 0.3982206
#> 69      7        0.7        V2 0.3748532
#> 70      7        0.7       V10 0.3666879
#> 71      8        0.8        V9 0.3537136
#> 72      8        0.8       V10 0.3661884
#> 73      8        0.8        V1 0.3803950
#> 74      8        0.8        V5 0.3443290
#> 75      8        0.8        V2 0.3742134
#> 76      8        0.8        V3 0.4318573
#> 77      8        0.8        V4 0.3773349
#> 78      8        0.8        V8 0.3976096
#> 79      8        0.8        V6 0.3985051
#> 80      8        0.8        V7 0.3761586
#> 81      9        0.9        V5 0.3437635
#> 82      9        0.9        V7 0.3755739
#> 83      9        0.9        V8 0.3969995
#> 84      9        0.9        V9 0.3531698
#> 85      9        0.9        V6 0.3979548
#> 86      9        0.9       V10 0.3656896
#> 87      9        0.9        V1 0.3798081
#> 88      9        0.9        V2 0.3735747
#> 89      9        0.9        V3 0.4312016
#> 90      9        0.9        V4 0.3768190
#> 91     10        1.0        V1 0.3792222
#> 92     10        1.0        V4 0.3763037
#> 93     10        1.0        V5 0.3431988
#> 94     10        1.0        V9 0.3526269
#> 95     10        1.0        V3 0.4305470
#> 96     10        1.0        V7 0.3749901
#> 97     10        1.0        V8 0.3963903
#> 98     10        1.0        V2 0.3729371
#> 99     10        1.0        V6 0.3974052
#> 100    10        1.0       V10 0.3651915

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3845287
#> 2       1        0.1        V5 0.3483139
#> 3       1        0.1        V9 0.3575436
#> 4       1        0.1        V3 0.4364745
#> 5       1        0.1        V4 0.3809664
#> 6       1        0.1        V8 0.4019067
#> 7       1        0.1        V2 0.3787148
#> 8       1        0.1        V6 0.4023789
#> 9       1        0.1        V7 0.3802772
#> 10      1        0.1       V10 0.3696992
#> 11      2        0.2        V7 0.3796861
#> 12      2        0.2        V8 0.4012900
#> 13      2        0.2        V9 0.3569939
#> 14      2        0.2       V10 0.3691956
#> 15      2        0.2        V1 0.3839354
#> 16      2        0.2        V5 0.3477418
#> 17      2        0.2        V2 0.3780685
#> 18      2        0.2        V3 0.4358119
#> 19      2        0.2        V4 0.3804455
#> 20      2        0.2        V6 0.4018232
#> 21      3        0.3        V4 0.3799253
#> 22      3        0.3        V5 0.3471707
#> 23      3        0.3        V3 0.4351503
#> 24      3        0.3        V7 0.3790959
#> 25      3        0.3        V8 0.4006743
#> 26      3        0.3        V9 0.3564451
#> 27      3        0.3        V6 0.4012682
#> 28      3        0.3       V10 0.3686927
#> 29      3        0.3        V1 0.3833431
#> 30      3        0.3        V2 0.3774232
#> 31      4        0.4        V1 0.3827517
#> 32      4        0.4        V9 0.3558971
#> 33      4        0.4        V3 0.4344897
#> 34      4        0.4        V4 0.3794058
#> 35      4        0.4        V5 0.3466005
#> 36      4        0.4        V2 0.3767791
#> 37      4        0.4        V6 0.4007141
#> 38      4        0.4        V7 0.3785066
#> 39      4        0.4        V8 0.4000594
#> 40      4        0.4       V10 0.3681905
#> 41      5        0.5        V8 0.3994456
#> 42      5        0.5        V9 0.3553500
#> 43      5        0.5       V10 0.3676889
#> 44      5        0.5        V1 0.3821611
#> 45      5        0.5        V5 0.3460312
#> 46      5        0.5        V2 0.3761360
#> 47      5        0.5        V3 0.4338301
#> 48      5        0.5        V4 0.3788870
#> 49      5        0.5        V6 0.4001607
#> 50      5        0.5        V7 0.3779182
#> 51      6        0.6        V4 0.3783689
#> 52      6        0.6        V5 0.3454629
#> 53      6        0.6        V7 0.3773308
#> 54      6        0.6        V8 0.3988326
#> 55      6        0.6        V9 0.3548037
#> 56      6        0.6        V6 0.3996081
#> 57      6        0.6       V10 0.3671881
#> 58      6        0.6        V1 0.3815715
#> 59      6        0.6        V2 0.3754941
#> 60      6        0.6        V3 0.4331715
#> 61      7        0.7        V1 0.3809828
#> 62      7        0.7        V3 0.4325139
#> 63      7        0.7        V4 0.3778516
#> 64      7        0.7        V5 0.3448955
#> 65      7        0.7        V9 0.3542582
#> 66      7        0.7        V6 0.3990562
#> 67      7        0.7        V7 0.3767442
#> 68      7        0.7        V8 0.3982206
#> 69      7        0.7        V2 0.3748532
#> 70      7        0.7       V10 0.3666879
#> 71      8        0.8        V9 0.3537136
#> 72      8        0.8       V10 0.3661884
#> 73      8        0.8        V1 0.3803950
#> 74      8        0.8        V5 0.3443290
#> 75      8        0.8        V2 0.3742134
#> 76      8        0.8        V3 0.4318573
#> 77      8        0.8        V4 0.3773349
#> 78      8        0.8        V8 0.3976096
#> 79      8        0.8        V6 0.3985051
#> 80      8        0.8        V7 0.3761586
#> 81      9        0.9        V5 0.3437635
#> 82      9        0.9        V7 0.3755739
#> 83      9        0.9        V8 0.3969995
#> 84      9        0.9        V9 0.3531698
#> 85      9        0.9        V6 0.3979548
#> 86      9        0.9       V10 0.3656896
#> 87      9        0.9        V1 0.3798081
#> 88      9        0.9        V2 0.3735747
#> 89      9        0.9        V3 0.4312016
#> 90      9        0.9        V4 0.3768190
#> 91     10        1.0        V1 0.3792222
#> 92     10        1.0        V4 0.3763037
#> 93     10        1.0        V5 0.3431988
#> 94     10        1.0        V9 0.3526269
#> 95     10        1.0        V3 0.4305470
#> 96     10        1.0        V7 0.3749901
#> 97     10        1.0        V8 0.3963903
#> 98     10        1.0        V2 0.3729371
#> 99     10        1.0        V6 0.3974052
#> 100    10        1.0       V10 0.3651915


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
#> 1       0                0          0 0.8757906 0.04742695 0.7806819 0.9369250
#> 2       0               20         20 0.8757906 0.04742695 0.7806819 0.9369250
#> 3       0               40         40 0.8757906 0.04742695 0.7806819 0.9369250
#> 4       0               60         60 0.8757906 0.04742695 0.7806819 0.9369250
#> 5      20                0         20 0.8617131 0.04946684 0.7636778 0.9250538
#> 6      20               20         40 0.8617131 0.04946684 0.7636778 0.9250538
#> 7      20               40         60 0.8617131 0.04946684 0.7636778 0.9250538
#> 8      20               60         80 0.8617131 0.04946684 0.7636778 0.9250538
#> 9      40                0         40 0.8478591 0.05145243 0.7472380 0.9130590
#> 10     40               20         60 0.8478591 0.05145243 0.7472380 0.9130590
#> 11     40               40         80 0.8478591 0.05145243 0.7472380 0.9130590
#> 12     40               60        100 0.8478591 0.05145243 0.7472380 0.9130590
#> 13     60                0         60 0.8342249 0.05338191 0.7313103 0.9009952
#> 14     60               20         80 0.8342249 0.05338191 0.7313103 0.9009952
#> 15     60               40        100 0.8342249 0.05338191 0.7313103 0.9009952
#> 16     60               60        120 0.8342249 0.05338191 0.7313103 0.9009952
#> 17     80                0         80 0.8208071 0.05525410 0.7158528 0.8889036
#> 18     80               20        100 0.8208071 0.05525410 0.7158528 0.8889036
#> 19     80               40        120 0.8208071 0.05525410 0.7158528 0.8889036
#> 20     80               60        140 0.8208071 0.05525410 0.7158528 0.8889036
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12025020 0.1584552030 0.5580056
#> 2  0.30574618 0.12103868 0.1107949248 0.5070414
#> 3  0.26001915 0.11961764 0.0759451020 0.4610006
#> 4  0.22113100 0.11650325 0.0507396518 0.4194973
#> 5  0.25589195 0.10946013 0.0825702623 0.4273706
#> 6  0.21762106 0.10701665 0.0555028689 0.3892180
#> 7  0.18507391 0.10346701 0.0361690442 0.3548884
#> 8  0.15739448 0.09910378 0.0226579043 0.3240056
#> 9  0.18213629 0.09732522 0.0397989936 0.3346101
#> 10 0.15489621 0.09325798 0.0251663197 0.3014877
#> 11 0.13173012 0.08874171 0.0151721324 0.2759466
#> 12 0.11202872 0.08393477 0.0086060361 0.2529362
#> 13 0.12963921 0.08513226 0.0170059409 0.2766377
#> 14 0.11025053 0.08041821 0.0097869626 0.2447033
#> 15 0.09376159 0.07563508 0.0052307104 0.2169984
#> 16 0.07973872 0.07086074 0.0025425832 0.1997013
#> 17 0.09227334 0.07360461 0.0060332504 0.2289866
#> 18 0.07847305 0.06881712 0.0029995918 0.2027080
#> 19 0.06673672 0.06416660 0.0013291077 0.1795206
#> 20 0.05675565 0.05967924 0.0005084438 0.1593007
```
