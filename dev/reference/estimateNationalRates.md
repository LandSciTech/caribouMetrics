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
#> 1         0.1 0.3838513 0.02282505 0.3593194 0.4296010
#> 2         0.2 0.3832760 0.02278208 0.3587550 0.4288850
#> 3         0.3 0.3827015 0.02273934 0.3581837 0.4281701
#> 4         0.4 0.3821279 0.02269684 0.3575894 0.4274564
#> 5         0.5 0.3815551 0.02265457 0.3569960 0.4267440
#> 6         0.6 0.3809832 0.02261254 0.3564036 0.4260327
#> 7         0.7 0.3804122 0.02257074 0.3558122 0.4253226
#> 8         0.8 0.3798420 0.02252918 0.3552218 0.4246138
#> 9         0.9 0.3792726 0.02248785 0.3546324 0.4239061
#> 10        1.0 0.3787041 0.02244675 0.3540439 0.4231996

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3924455
#> 2       1        0.1        V5 0.4225057
#> 3       1        0.1        V9 0.3841077
#> 4       1        0.1        V3 0.3880820
#> 5       1        0.1        V4 0.4316610
#> 6       1        0.1        V8 0.3782452
#> 7       1        0.1        V2 0.3779963
#> 8       1        0.1        V6 0.4072112
#> 9       1        0.1        V7 0.3538971
#> 10      1        0.1       V10 0.3857341
#> 11      2        0.2        V7 0.3533097
#> 12      2        0.2        V8 0.3776179
#> 13      2        0.2        V9 0.3835747
#> 14      2        0.2       V10 0.3851254
#> 15      2        0.2        V1 0.3918014
#> 16      2        0.2        V5 0.4219147
#> 17      2        0.2        V2 0.3775112
#> 18      2        0.2        V3 0.3874353
#> 19      2        0.2        V4 0.4309086
#> 20      2        0.2        V6 0.4065403
#> 21      3        0.3        V4 0.4301575
#> 22      3        0.3        V5 0.4213246
#> 23      3        0.3        V3 0.3867896
#> 24      3        0.3        V7 0.3527234
#> 25      3        0.3        V8 0.3769916
#> 26      3        0.3        V9 0.3830425
#> 27      3        0.3        V6 0.4058706
#> 28      3        0.3       V10 0.3845176
#> 29      3        0.3        V1 0.3911584
#> 30      3        0.3        V2 0.3770267
#> 31      4        0.4        V1 0.3905165
#> 32      4        0.4        V9 0.3825110
#> 33      4        0.4        V3 0.3861450
#> 34      4        0.4        V4 0.4294078
#> 35      4        0.4        V5 0.4207353
#> 36      4        0.4        V2 0.3765428
#> 37      4        0.4        V6 0.4052020
#> 38      4        0.4        V7 0.3521380
#> 39      4        0.4        V8 0.3763664
#> 40      4        0.4       V10 0.3839108
#> 41      5        0.5        V8 0.3757422
#> 42      5        0.5        V9 0.3819802
#> 43      5        0.5       V10 0.3833050
#> 44      5        0.5        V1 0.3898756
#> 45      5        0.5        V5 0.4201468
#> 46      5        0.5        V2 0.3760595
#> 47      5        0.5        V3 0.3855014
#> 48      5        0.5        V4 0.4286593
#> 49      5        0.5        V6 0.4045344
#> 50      5        0.5        V7 0.3515536
#> 51      6        0.6        V4 0.4279121
#> 52      6        0.6        V5 0.4195592
#> 53      6        0.6        V7 0.3509701
#> 54      6        0.6        V8 0.3751190
#> 55      6        0.6        V9 0.3814502
#> 56      6        0.6        V6 0.4038680
#> 57      6        0.6       V10 0.3827001
#> 58      6        0.6        V1 0.3892358
#> 59      6        0.6        V2 0.3755769
#> 60      6        0.6        V3 0.3848589
#> 61      7        0.7        V1 0.3885970
#> 62      7        0.7        V3 0.3842175
#> 63      7        0.7        V4 0.4271663
#> 64      7        0.7        V5 0.4189723
#> 65      7        0.7        V9 0.3809209
#> 66      7        0.7        V6 0.4032027
#> 67      7        0.7        V7 0.3503877
#> 68      7        0.7        V8 0.3744969
#> 69      7        0.7        V2 0.3750949
#> 70      7        0.7       V10 0.3820962
#> 71      8        0.8        V9 0.3803924
#> 72      8        0.8       V10 0.3814932
#> 73      8        0.8        V1 0.3879593
#> 74      8        0.8        V5 0.4183863
#> 75      8        0.8        V2 0.3746135
#> 76      8        0.8        V3 0.3835772
#> 77      8        0.8        V4 0.4264217
#> 78      8        0.8        V8 0.3738758
#> 79      8        0.8        V6 0.4025384
#> 80      8        0.8        V7 0.3498062
#> 81      9        0.9        V5 0.4178011
#> 82      9        0.9        V7 0.3492256
#> 83      9        0.9        V8 0.3732557
#> 84      9        0.9        V9 0.3798646
#> 85      9        0.9        V6 0.4018753
#> 86      9        0.9       V10 0.3808912
#> 87      9        0.9        V1 0.3873226
#> 88      9        0.9        V2 0.3741327
#> 89      9        0.9        V3 0.3829379
#> 90      9        0.9        V4 0.4256785
#> 91     10        1.0        V1 0.3866870
#> 92     10        1.0        V4 0.4249365
#> 93     10        1.0        V5 0.4172167
#> 94     10        1.0        V9 0.3793375
#> 95     10        1.0        V3 0.3822997
#> 96     10        1.0        V7 0.3486460
#> 97     10        1.0        V8 0.3726367
#> 98     10        1.0        V2 0.3736525
#> 99     10        1.0        V6 0.4012133
#> 100    10        1.0       V10 0.3802901

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3924455
#> 2       1        0.1        V5 0.4225057
#> 3       1        0.1        V9 0.3841077
#> 4       1        0.1        V3 0.3880820
#> 5       1        0.1        V4 0.4316610
#> 6       1        0.1        V8 0.3782452
#> 7       1        0.1        V2 0.3779963
#> 8       1        0.1        V6 0.4072112
#> 9       1        0.1        V7 0.3538971
#> 10      1        0.1       V10 0.3857341
#> 11      2        0.2        V7 0.3533097
#> 12      2        0.2        V8 0.3776179
#> 13      2        0.2        V9 0.3835747
#> 14      2        0.2       V10 0.3851254
#> 15      2        0.2        V1 0.3918014
#> 16      2        0.2        V5 0.4219147
#> 17      2        0.2        V2 0.3775112
#> 18      2        0.2        V3 0.3874353
#> 19      2        0.2        V4 0.4309086
#> 20      2        0.2        V6 0.4065403
#> 21      3        0.3        V4 0.4301575
#> 22      3        0.3        V5 0.4213246
#> 23      3        0.3        V3 0.3867896
#> 24      3        0.3        V7 0.3527234
#> 25      3        0.3        V8 0.3769916
#> 26      3        0.3        V9 0.3830425
#> 27      3        0.3        V6 0.4058706
#> 28      3        0.3       V10 0.3845176
#> 29      3        0.3        V1 0.3911584
#> 30      3        0.3        V2 0.3770267
#> 31      4        0.4        V1 0.3905165
#> 32      4        0.4        V9 0.3825110
#> 33      4        0.4        V3 0.3861450
#> 34      4        0.4        V4 0.4294078
#> 35      4        0.4        V5 0.4207353
#> 36      4        0.4        V2 0.3765428
#> 37      4        0.4        V6 0.4052020
#> 38      4        0.4        V7 0.3521380
#> 39      4        0.4        V8 0.3763664
#> 40      4        0.4       V10 0.3839108
#> 41      5        0.5        V8 0.3757422
#> 42      5        0.5        V9 0.3819802
#> 43      5        0.5       V10 0.3833050
#> 44      5        0.5        V1 0.3898756
#> 45      5        0.5        V5 0.4201468
#> 46      5        0.5        V2 0.3760595
#> 47      5        0.5        V3 0.3855014
#> 48      5        0.5        V4 0.4286593
#> 49      5        0.5        V6 0.4045344
#> 50      5        0.5        V7 0.3515536
#> 51      6        0.6        V4 0.4279121
#> 52      6        0.6        V5 0.4195592
#> 53      6        0.6        V7 0.3509701
#> 54      6        0.6        V8 0.3751190
#> 55      6        0.6        V9 0.3814502
#> 56      6        0.6        V6 0.4038680
#> 57      6        0.6       V10 0.3827001
#> 58      6        0.6        V1 0.3892358
#> 59      6        0.6        V2 0.3755769
#> 60      6        0.6        V3 0.3848589
#> 61      7        0.7        V1 0.3885970
#> 62      7        0.7        V3 0.3842175
#> 63      7        0.7        V4 0.4271663
#> 64      7        0.7        V5 0.4189723
#> 65      7        0.7        V9 0.3809209
#> 66      7        0.7        V6 0.4032027
#> 67      7        0.7        V7 0.3503877
#> 68      7        0.7        V8 0.3744969
#> 69      7        0.7        V2 0.3750949
#> 70      7        0.7       V10 0.3820962
#> 71      8        0.8        V9 0.3803924
#> 72      8        0.8       V10 0.3814932
#> 73      8        0.8        V1 0.3879593
#> 74      8        0.8        V5 0.4183863
#> 75      8        0.8        V2 0.3746135
#> 76      8        0.8        V3 0.3835772
#> 77      8        0.8        V4 0.4264217
#> 78      8        0.8        V8 0.3738758
#> 79      8        0.8        V6 0.4025384
#> 80      8        0.8        V7 0.3498062
#> 81      9        0.9        V5 0.4178011
#> 82      9        0.9        V7 0.3492256
#> 83      9        0.9        V8 0.3732557
#> 84      9        0.9        V9 0.3798646
#> 85      9        0.9        V6 0.4018753
#> 86      9        0.9       V10 0.3808912
#> 87      9        0.9        V1 0.3873226
#> 88      9        0.9        V2 0.3741327
#> 89      9        0.9        V3 0.3829379
#> 90      9        0.9        V4 0.4256785
#> 91     10        1.0        V1 0.3866870
#> 92     10        1.0        V4 0.4249365
#> 93     10        1.0        V5 0.4172167
#> 94     10        1.0        V9 0.3793375
#> 95     10        1.0        V3 0.3822997
#> 96     10        1.0        V7 0.3486460
#> 97     10        1.0        V8 0.3726367
#> 98     10        1.0        V2 0.3736525
#> 99     10        1.0        V6 0.4012133
#> 100    10        1.0       V10 0.3802901


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
#> 1       0                0          0 0.8757906 0.04582965 0.7887440 0.9468497
#> 2       0               20         20 0.8757906 0.04582965 0.7887440 0.9468497
#> 3       0               40         40 0.8757906 0.04582965 0.7887440 0.9468497
#> 4       0               60         60 0.8757906 0.04582965 0.7887440 0.9468497
#> 5      20                0         20 0.8617131 0.04853584 0.7693605 0.9378352
#> 6      20               20         40 0.8617131 0.04853584 0.7693605 0.9378352
#> 7      20               40         60 0.8617131 0.04853584 0.7693605 0.9378352
#> 8      20               60         80 0.8617131 0.04853584 0.7693605 0.9378352
#> 9      40                0         40 0.8478591 0.05105246 0.7506736 0.9286876
#> 10     40               20         60 0.8478591 0.05105246 0.7506736 0.9286876
#> 11     40               40         80 0.8478591 0.05105246 0.7506736 0.9286876
#> 12     40               60        100 0.8478591 0.05105246 0.7506736 0.9286876
#> 13     60                0         60 0.8342249 0.05340025 0.7326156 0.9194445
#> 14     60               20         80 0.8342249 0.05340025 0.7326156 0.9194445
#> 15     60               40        100 0.8342249 0.05340025 0.7326156 0.9194445
#> 16     60               60        120 0.8342249 0.05340025 0.7326156 0.9194445
#> 17     80                0         80 0.8208071 0.05559608 0.7151331 0.9101356
#> 18     80               20        100 0.8208071 0.05559608 0.7151331 0.9101356
#> 19     80               40        120 0.8208071 0.05559608 0.7151331 0.9101356
#> 20     80               60        140 0.8208071 0.05559608 0.7151331 0.9101356
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11667160 0.176613658 0.5670748
#> 2  0.30574618 0.10771044 0.144038526 0.4879675
#> 3  0.26001915 0.10105300 0.116842359 0.4205584
#> 4  0.22113100 0.09588005 0.094184160 0.3633411
#> 5  0.25589195 0.10479255 0.105842586 0.4566110
#> 6  0.21762106 0.09562753 0.085039323 0.3939252
#> 7  0.18507391 0.08829818 0.067789557 0.3407639
#> 8  0.15739448 0.08225999 0.053550129 0.2957083
#> 9  0.18213629 0.09263147 0.060858246 0.3691772
#> 10 0.15489621 0.08405299 0.047851250 0.3197901
#> 11 0.13173012 0.07691539 0.037205740 0.2779250
#> 12 0.11202872 0.07084558 0.028557880 0.2423888
#> 13 0.12963921 0.08092731 0.032977167 0.3003040
#> 14 0.11025053 0.07321502 0.025145155 0.2613930
#> 15 0.09376159 0.06663445 0.018870448 0.2283330
#> 16 0.07973872 0.06092295 0.013903060 0.2001640
#> 17 0.09227334 0.07007889 0.016423995 0.2460178
#> 18 0.07847305 0.06328067 0.011986386 0.2152441
#> 19 0.06673672 0.05737546 0.008549948 0.1889815
#> 20 0.05675565 0.05217479 0.005938798 0.1664729
```
