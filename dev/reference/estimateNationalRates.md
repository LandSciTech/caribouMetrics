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
#> 1         0.1 0.3838513 0.02742888 0.3583120 0.4298088
#> 2         0.2 0.3832760 0.02736433 0.3577707 0.4291078
#> 3         0.3 0.3827015 0.02730007 0.3572301 0.4284079
#> 4         0.4 0.3821279 0.02723608 0.3566904 0.4277092
#> 5         0.5 0.3815551 0.02717238 0.3561515 0.4270117
#> 6         0.6 0.3809832 0.02710894 0.3556134 0.4263152
#> 7         0.7 0.3804122 0.02704579 0.3550761 0.4256200
#> 8         0.8 0.3798420 0.02698291 0.3545396 0.4249258
#> 9         0.9 0.3792726 0.02692031 0.3540039 0.4242328
#> 10        1.0 0.3787041 0.02685798 0.3534691 0.4235410

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4282320
#> 2       1        0.1        V5 0.3773062
#> 3       1        0.1        V9 0.3885666
#> 4       1        0.1        V3 0.3564701
#> 5       1        0.1        V4 0.3646565
#> 6       1        0.1        V8 0.3933436
#> 7       1        0.1        V2 0.4299645
#> 8       1        0.1        V6 0.4292725
#> 9       1        0.1        V7 0.3762760
#> 10      1        0.1       V10 0.4075217
#> 11      2        0.2        V7 0.3757379
#> 12      2        0.2        V8 0.3927673
#> 13      2        0.2        V9 0.3880756
#> 14      2        0.2       V10 0.4069606
#> 15      2        0.2        V1 0.4275559
#> 16      2        0.2        V5 0.3766597
#> 17      2        0.2        V2 0.4292886
#> 18      2        0.2        V3 0.3559330
#> 19      2        0.2        V4 0.3641005
#> 20      2        0.2        V6 0.4284848
#> 21      3        0.3        V4 0.3635453
#> 22      3        0.3        V5 0.3760143
#> 23      3        0.3        V3 0.3553967
#> 24      3        0.3        V7 0.3752006
#> 25      3        0.3        V8 0.3921918
#> 26      3        0.3        V9 0.3875852
#> 27      3        0.3        V6 0.4276985
#> 28      3        0.3       V10 0.4064003
#> 29      3        0.3        V1 0.4268810
#> 30      3        0.3        V2 0.4286139
#> 31      4        0.4        V1 0.4262071
#> 32      4        0.4        V9 0.3870954
#> 33      4        0.4        V3 0.3548612
#> 34      4        0.4        V4 0.3629909
#> 35      4        0.4        V5 0.3753701
#> 36      4        0.4        V2 0.4279402
#> 37      4        0.4        V6 0.4269137
#> 38      4        0.4        V7 0.3746640
#> 39      4        0.4        V8 0.3916172
#> 40      4        0.4       V10 0.4058408
#> 41      5        0.5        V8 0.3910434
#> 42      5        0.5        V9 0.3866062
#> 43      5        0.5       V10 0.4052820
#> 44      5        0.5        V1 0.4255342
#> 45      5        0.5        V5 0.3747269
#> 46      5        0.5        V2 0.4272675
#> 47      5        0.5        V3 0.3543265
#> 48      5        0.5        V4 0.3624374
#> 49      5        0.5        V6 0.4261303
#> 50      5        0.5        V7 0.3741282
#> 51      6        0.6        V4 0.3618847
#> 52      6        0.6        V5 0.3740849
#> 53      6        0.6        V7 0.3735932
#> 54      6        0.6        V8 0.3904705
#> 55      6        0.6        V9 0.3861177
#> 56      6        0.6        V6 0.4253483
#> 57      6        0.6       V10 0.4047240
#> 58      6        0.6        V1 0.4248625
#> 59      6        0.6        V2 0.4265960
#> 60      6        0.6        V3 0.3537926
#> 61      7        0.7        V1 0.4241918
#> 62      7        0.7        V3 0.3532596
#> 63      7        0.7        V4 0.3613329
#> 64      7        0.7        V5 0.3734439
#> 65      7        0.7        V9 0.3856298
#> 66      7        0.7        V6 0.4245678
#> 67      7        0.7        V7 0.3730589
#> 68      7        0.7        V8 0.3898984
#> 69      7        0.7        V2 0.4259254
#> 70      7        0.7       V10 0.4041668
#> 71      8        0.8        V9 0.3851424
#> 72      8        0.8       V10 0.4036103
#> 73      8        0.8        V1 0.4235221
#> 74      8        0.8        V5 0.3728040
#> 75      8        0.8        V2 0.4252560
#> 76      8        0.8        V3 0.3527273
#> 77      8        0.8        V4 0.3607819
#> 78      8        0.8        V8 0.3893271
#> 79      8        0.8        V6 0.4237887
#> 80      8        0.8        V7 0.3725254
#> 81      9        0.9        V5 0.3721653
#> 82      9        0.9        V7 0.3719926
#> 83      9        0.9        V8 0.3887567
#> 84      9        0.9        V9 0.3846558
#> 85      9        0.9        V6 0.4230110
#> 86      9        0.9       V10 0.4030546
#> 87      9        0.9        V1 0.4228535
#> 88      9        0.9        V2 0.4245875
#> 89      9        0.9        V3 0.3521958
#> 90      9        0.9        V4 0.3602317
#> 91     10        1.0        V1 0.4221860
#> 92     10        1.0        V4 0.3596824
#> 93     10        1.0        V5 0.3715276
#> 94     10        1.0        V9 0.3841697
#> 95     10        1.0        V3 0.3516652
#> 96     10        1.0        V7 0.3714606
#> 97     10        1.0        V8 0.3881871
#> 98     10        1.0        V2 0.4239202
#> 99     10        1.0        V6 0.4222348
#> 100    10        1.0       V10 0.4024997

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4282320
#> 2       1        0.1        V5 0.3773062
#> 3       1        0.1        V9 0.3885666
#> 4       1        0.1        V3 0.3564701
#> 5       1        0.1        V4 0.3646565
#> 6       1        0.1        V8 0.3933436
#> 7       1        0.1        V2 0.4299645
#> 8       1        0.1        V6 0.4292725
#> 9       1        0.1        V7 0.3762760
#> 10      1        0.1       V10 0.4075217
#> 11      2        0.2        V7 0.3757379
#> 12      2        0.2        V8 0.3927673
#> 13      2        0.2        V9 0.3880756
#> 14      2        0.2       V10 0.4069606
#> 15      2        0.2        V1 0.4275559
#> 16      2        0.2        V5 0.3766597
#> 17      2        0.2        V2 0.4292886
#> 18      2        0.2        V3 0.3559330
#> 19      2        0.2        V4 0.3641005
#> 20      2        0.2        V6 0.4284848
#> 21      3        0.3        V4 0.3635453
#> 22      3        0.3        V5 0.3760143
#> 23      3        0.3        V3 0.3553967
#> 24      3        0.3        V7 0.3752006
#> 25      3        0.3        V8 0.3921918
#> 26      3        0.3        V9 0.3875852
#> 27      3        0.3        V6 0.4276985
#> 28      3        0.3       V10 0.4064003
#> 29      3        0.3        V1 0.4268810
#> 30      3        0.3        V2 0.4286139
#> 31      4        0.4        V1 0.4262071
#> 32      4        0.4        V9 0.3870954
#> 33      4        0.4        V3 0.3548612
#> 34      4        0.4        V4 0.3629909
#> 35      4        0.4        V5 0.3753701
#> 36      4        0.4        V2 0.4279402
#> 37      4        0.4        V6 0.4269137
#> 38      4        0.4        V7 0.3746640
#> 39      4        0.4        V8 0.3916172
#> 40      4        0.4       V10 0.4058408
#> 41      5        0.5        V8 0.3910434
#> 42      5        0.5        V9 0.3866062
#> 43      5        0.5       V10 0.4052820
#> 44      5        0.5        V1 0.4255342
#> 45      5        0.5        V5 0.3747269
#> 46      5        0.5        V2 0.4272675
#> 47      5        0.5        V3 0.3543265
#> 48      5        0.5        V4 0.3624374
#> 49      5        0.5        V6 0.4261303
#> 50      5        0.5        V7 0.3741282
#> 51      6        0.6        V4 0.3618847
#> 52      6        0.6        V5 0.3740849
#> 53      6        0.6        V7 0.3735932
#> 54      6        0.6        V8 0.3904705
#> 55      6        0.6        V9 0.3861177
#> 56      6        0.6        V6 0.4253483
#> 57      6        0.6       V10 0.4047240
#> 58      6        0.6        V1 0.4248625
#> 59      6        0.6        V2 0.4265960
#> 60      6        0.6        V3 0.3537926
#> 61      7        0.7        V1 0.4241918
#> 62      7        0.7        V3 0.3532596
#> 63      7        0.7        V4 0.3613329
#> 64      7        0.7        V5 0.3734439
#> 65      7        0.7        V9 0.3856298
#> 66      7        0.7        V6 0.4245678
#> 67      7        0.7        V7 0.3730589
#> 68      7        0.7        V8 0.3898984
#> 69      7        0.7        V2 0.4259254
#> 70      7        0.7       V10 0.4041668
#> 71      8        0.8        V9 0.3851424
#> 72      8        0.8       V10 0.4036103
#> 73      8        0.8        V1 0.4235221
#> 74      8        0.8        V5 0.3728040
#> 75      8        0.8        V2 0.4252560
#> 76      8        0.8        V3 0.3527273
#> 77      8        0.8        V4 0.3607819
#> 78      8        0.8        V8 0.3893271
#> 79      8        0.8        V6 0.4237887
#> 80      8        0.8        V7 0.3725254
#> 81      9        0.9        V5 0.3721653
#> 82      9        0.9        V7 0.3719926
#> 83      9        0.9        V8 0.3887567
#> 84      9        0.9        V9 0.3846558
#> 85      9        0.9        V6 0.4230110
#> 86      9        0.9       V10 0.4030546
#> 87      9        0.9        V1 0.4228535
#> 88      9        0.9        V2 0.4245875
#> 89      9        0.9        V3 0.3521958
#> 90      9        0.9        V4 0.3602317
#> 91     10        1.0        V1 0.4221860
#> 92     10        1.0        V4 0.3596824
#> 93     10        1.0        V5 0.3715276
#> 94     10        1.0        V9 0.3841697
#> 95     10        1.0        V3 0.3516652
#> 96     10        1.0        V7 0.3714606
#> 97     10        1.0        V8 0.3881871
#> 98     10        1.0        V2 0.4239202
#> 99     10        1.0        V6 0.4222348
#> 100    10        1.0       V10 0.4024997


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
#> 1       0                0          0 0.8757906 0.04372455 0.7894126 0.9447215
#> 2       0               20         20 0.8757906 0.04372455 0.7894126 0.9447215
#> 3       0               40         40 0.8757906 0.04372455 0.7894126 0.9447215
#> 4       0               60         60 0.8757906 0.04372455 0.7894126 0.9447215
#> 5      20                0         20 0.8617131 0.04566677 0.7738461 0.9351069
#> 6      20               20         40 0.8617131 0.04566677 0.7738461 0.9351069
#> 7      20               40         60 0.8617131 0.04566677 0.7738461 0.9351069
#> 8      20               60         80 0.8617131 0.04566677 0.7738461 0.9351069
#> 9      40                0         40 0.8478591 0.04748737 0.7587370 0.9253551
#> 10     40               20         60 0.8478591 0.04748737 0.7587370 0.9253551
#> 11     40               40         80 0.8478591 0.04748737 0.7587370 0.9253551
#> 12     40               60        100 0.8478591 0.04748737 0.7587370 0.9253551
#> 13     60                0         60 0.8342249 0.04919783 0.7440476 0.9155078
#> 14     60               20         80 0.8342249 0.04919783 0.7440476 0.9155078
#> 15     60               40        100 0.8342249 0.04919783 0.7440476 0.9155078
#> 16     60               60        120 0.8342249 0.04919783 0.7440476 0.9155078
#> 17     80                0         80 0.8208071 0.05080767 0.7297469 0.9055978
#> 18     80               20        100 0.8208071 0.05080767 0.7297469 0.9055978
#> 19     80               40        120 0.8208071 0.05080767 0.7297469 0.9055978
#> 20     80               60        140 0.8208071 0.05080767 0.7297469 0.9055978
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12879573 0.167639037 0.5932997
#> 2  0.30574618 0.12766767 0.135327287 0.5480425
#> 3  0.26001915 0.12549142 0.108549019 0.5063577
#> 4  0.22113100 0.12221299 0.086414905 0.4680512
#> 5  0.25589195 0.11324391 0.093318060 0.4691373
#> 6  0.21762106 0.11036900 0.073864863 0.4338991
#> 7  0.18507391 0.10700598 0.057899440 0.4015972
#> 8  0.15739448 0.10310528 0.044874891 0.3720034
#> 9  0.18213629 0.09843221 0.048915635 0.3728416
#> 10 0.15489621 0.09481508 0.037590950 0.3456641
#> 11 0.13173012 0.09101667 0.028475216 0.3207684
#> 12 0.11202872 0.08699554 0.021213948 0.2979556
#> 13 0.12963921 0.08471652 0.023445755 0.2986020
#> 14 0.11025053 0.08093605 0.017249485 0.2776327
#> 15 0.09376159 0.07713439 0.012424042 0.2583929
#> 16 0.07973872 0.07327840 0.008730922 0.2407232
#> 17 0.09227334 0.07232737 0.009848375 0.2412245
#> 18 0.07847305 0.06869214 0.006793615 0.2249380
#> 19 0.06673672 0.06511190 0.004541076 0.2099439
#> 20 0.05675565 0.06156109 0.002926588 0.1961190
```
