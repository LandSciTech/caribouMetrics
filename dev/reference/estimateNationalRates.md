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
#> 1         0.1 0.3838513 0.02984607 0.3535110 0.4391124
#> 2         0.2 0.3832760 0.02980737 0.3529749 0.4385152
#> 3         0.3 0.3827015 0.02976886 0.3524396 0.4379188
#> 4         0.4 0.3821279 0.02973054 0.3519022 0.4373232
#> 5         0.5 0.3815551 0.02969241 0.3513615 0.4367285
#> 6         0.6 0.3809832 0.02965446 0.3508215 0.4361346
#> 7         0.7 0.3804122 0.02961671 0.3502824 0.4355415
#> 8         0.8 0.3798420 0.02957913 0.3497442 0.4349492
#> 9         0.9 0.3792726 0.02954175 0.3492067 0.4343577
#> 10        1.0 0.3787041 0.02950454 0.3486701 0.4337670

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4348217
#> 2       1        0.1        V5 0.3665058
#> 3       1        0.1        V9 0.3497622
#> 4       1        0.1        V3 0.4024042
#> 5       1        0.1        V4 0.3742826
#> 6       1        0.1        V8 0.3664235
#> 7       1        0.1        V2 0.3780071
#> 8       1        0.1        V6 0.3754378
#> 9       1        0.1        V7 0.4403581
#> 10      1        0.1       V10 0.3807055
#> 11      2        0.2        V7 0.4397924
#> 12      2        0.2        V8 0.3658910
#> 13      2        0.2        V9 0.3492250
#> 14      2        0.2       V10 0.3800338
#> 15      2        0.2        V1 0.4341157
#> 16      2        0.2        V5 0.3659415
#> 17      2        0.2        V2 0.3775106
#> 18      2        0.2        V3 0.4017097
#> 19      2        0.2        V4 0.3737457
#> 20      2        0.2        V6 0.3748371
#> 21      3        0.3        V4 0.3732095
#> 22      3        0.3        V5 0.3653781
#> 23      3        0.3        V3 0.4010163
#> 24      3        0.3        V7 0.4392275
#> 25      3        0.3        V8 0.3653593
#> 26      3        0.3        V9 0.3486887
#> 27      3        0.3        V6 0.3742373
#> 28      3        0.3       V10 0.3793633
#> 29      3        0.3        V1 0.4334109
#> 30      3        0.3        V2 0.3770147
#> 31      4        0.4        V1 0.4327071
#> 32      4        0.4        V9 0.3481532
#> 33      4        0.4        V3 0.4003241
#> 34      4        0.4        V4 0.3726740
#> 35      4        0.4        V5 0.3648155
#> 36      4        0.4        V2 0.3765195
#> 37      4        0.4        V6 0.3736385
#> 38      4        0.4        V7 0.4386634
#> 39      4        0.4        V8 0.3648284
#> 40      4        0.4       V10 0.3786940
#> 41      5        0.5        V8 0.3642982
#> 42      5        0.5        V9 0.3476185
#> 43      5        0.5       V10 0.3780259
#> 44      5        0.5        V1 0.4320046
#> 45      5        0.5        V5 0.3642539
#> 46      5        0.5        V2 0.3760249
#> 47      5        0.5        V3 0.3996331
#> 48      5        0.5        V4 0.3721394
#> 49      5        0.5        V6 0.3730407
#> 50      5        0.5        V7 0.4381000
#> 51      6        0.6        V4 0.3716055
#> 52      6        0.6        V5 0.3636931
#> 53      6        0.6        V7 0.4375372
#> 54      6        0.6        V8 0.3637689
#> 55      6        0.6        V9 0.3470846
#> 56      6        0.6        V6 0.3724438
#> 57      6        0.6       V10 0.3773590
#> 58      6        0.6        V1 0.4313031
#> 59      6        0.6        V2 0.3755310
#> 60      6        0.6        V3 0.3989433
#> 61      7        0.7        V1 0.4306028
#> 62      7        0.7        V3 0.3982547
#> 63      7        0.7        V4 0.3710724
#> 64      7        0.7        V5 0.3631331
#> 65      7        0.7        V9 0.3465516
#> 66      7        0.7        V6 0.3718478
#> 67      7        0.7        V7 0.4369752
#> 68      7        0.7        V8 0.3632402
#> 69      7        0.7        V2 0.3750377
#> 70      7        0.7       V10 0.3766933
#> 71      8        0.8        V9 0.3460193
#> 72      8        0.8       V10 0.3760287
#> 73      8        0.8        V1 0.4299037
#> 74      8        0.8        V5 0.3625740
#> 75      8        0.8        V2 0.3745451
#> 76      8        0.8        V3 0.3975673
#> 77      8        0.8        V4 0.3705400
#> 78      8        0.8        V8 0.3627124
#> 79      8        0.8        V6 0.3712528
#> 80      8        0.8        V7 0.4364140
#> 81      9        0.9        V5 0.3620158
#> 82      9        0.9        V7 0.4358534
#> 83      9        0.9        V8 0.3621853
#> 84      9        0.9        V9 0.3454879
#> 85      9        0.9        V6 0.3706588
#> 86      9        0.9       V10 0.3753653
#> 87      9        0.9        V1 0.4292057
#> 88      9        0.9        V2 0.3740531
#> 89      9        0.9        V3 0.3968811
#> 90      9        0.9        V4 0.3700084
#> 91     10        1.0        V1 0.4285088
#> 92     10        1.0        V4 0.3694776
#> 93     10        1.0        V5 0.3614585
#> 94     10        1.0        V9 0.3449573
#> 95     10        1.0        V3 0.3961961
#> 96     10        1.0        V7 0.4352936
#> 97     10        1.0        V8 0.3616590
#> 98     10        1.0        V2 0.3735618
#> 99     10        1.0        V6 0.3700657
#> 100    10        1.0       V10 0.3747030

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4348217
#> 2       1        0.1        V5 0.3665058
#> 3       1        0.1        V9 0.3497622
#> 4       1        0.1        V3 0.4024042
#> 5       1        0.1        V4 0.3742826
#> 6       1        0.1        V8 0.3664235
#> 7       1        0.1        V2 0.3780071
#> 8       1        0.1        V6 0.3754378
#> 9       1        0.1        V7 0.4403581
#> 10      1        0.1       V10 0.3807055
#> 11      2        0.2        V7 0.4397924
#> 12      2        0.2        V8 0.3658910
#> 13      2        0.2        V9 0.3492250
#> 14      2        0.2       V10 0.3800338
#> 15      2        0.2        V1 0.4341157
#> 16      2        0.2        V5 0.3659415
#> 17      2        0.2        V2 0.3775106
#> 18      2        0.2        V3 0.4017097
#> 19      2        0.2        V4 0.3737457
#> 20      2        0.2        V6 0.3748371
#> 21      3        0.3        V4 0.3732095
#> 22      3        0.3        V5 0.3653781
#> 23      3        0.3        V3 0.4010163
#> 24      3        0.3        V7 0.4392275
#> 25      3        0.3        V8 0.3653593
#> 26      3        0.3        V9 0.3486887
#> 27      3        0.3        V6 0.3742373
#> 28      3        0.3       V10 0.3793633
#> 29      3        0.3        V1 0.4334109
#> 30      3        0.3        V2 0.3770147
#> 31      4        0.4        V1 0.4327071
#> 32      4        0.4        V9 0.3481532
#> 33      4        0.4        V3 0.4003241
#> 34      4        0.4        V4 0.3726740
#> 35      4        0.4        V5 0.3648155
#> 36      4        0.4        V2 0.3765195
#> 37      4        0.4        V6 0.3736385
#> 38      4        0.4        V7 0.4386634
#> 39      4        0.4        V8 0.3648284
#> 40      4        0.4       V10 0.3786940
#> 41      5        0.5        V8 0.3642982
#> 42      5        0.5        V9 0.3476185
#> 43      5        0.5       V10 0.3780259
#> 44      5        0.5        V1 0.4320046
#> 45      5        0.5        V5 0.3642539
#> 46      5        0.5        V2 0.3760249
#> 47      5        0.5        V3 0.3996331
#> 48      5        0.5        V4 0.3721394
#> 49      5        0.5        V6 0.3730407
#> 50      5        0.5        V7 0.4381000
#> 51      6        0.6        V4 0.3716055
#> 52      6        0.6        V5 0.3636931
#> 53      6        0.6        V7 0.4375372
#> 54      6        0.6        V8 0.3637689
#> 55      6        0.6        V9 0.3470846
#> 56      6        0.6        V6 0.3724438
#> 57      6        0.6       V10 0.3773590
#> 58      6        0.6        V1 0.4313031
#> 59      6        0.6        V2 0.3755310
#> 60      6        0.6        V3 0.3989433
#> 61      7        0.7        V1 0.4306028
#> 62      7        0.7        V3 0.3982547
#> 63      7        0.7        V4 0.3710724
#> 64      7        0.7        V5 0.3631331
#> 65      7        0.7        V9 0.3465516
#> 66      7        0.7        V6 0.3718478
#> 67      7        0.7        V7 0.4369752
#> 68      7        0.7        V8 0.3632402
#> 69      7        0.7        V2 0.3750377
#> 70      7        0.7       V10 0.3766933
#> 71      8        0.8        V9 0.3460193
#> 72      8        0.8       V10 0.3760287
#> 73      8        0.8        V1 0.4299037
#> 74      8        0.8        V5 0.3625740
#> 75      8        0.8        V2 0.3745451
#> 76      8        0.8        V3 0.3975673
#> 77      8        0.8        V4 0.3705400
#> 78      8        0.8        V8 0.3627124
#> 79      8        0.8        V6 0.3712528
#> 80      8        0.8        V7 0.4364140
#> 81      9        0.9        V5 0.3620158
#> 82      9        0.9        V7 0.4358534
#> 83      9        0.9        V8 0.3621853
#> 84      9        0.9        V9 0.3454879
#> 85      9        0.9        V6 0.3706588
#> 86      9        0.9       V10 0.3753653
#> 87      9        0.9        V1 0.4292057
#> 88      9        0.9        V2 0.3740531
#> 89      9        0.9        V3 0.3968811
#> 90      9        0.9        V4 0.3700084
#> 91     10        1.0        V1 0.4285088
#> 92     10        1.0        V4 0.3694776
#> 93     10        1.0        V5 0.3614585
#> 94     10        1.0        V9 0.3449573
#> 95     10        1.0        V3 0.3961961
#> 96     10        1.0        V7 0.4352936
#> 97     10        1.0        V8 0.3616590
#> 98     10        1.0        V2 0.3735618
#> 99     10        1.0        V6 0.3700657
#> 100    10        1.0       V10 0.3747030


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
#> 1       0                0          0 0.8757906 0.04804776 0.7825237 0.9462538
#> 2       0               20         20 0.8757906 0.04804776 0.7825237 0.9462538
#> 3       0               40         40 0.8757906 0.04804776 0.7825237 0.9462538
#> 4       0               60         60 0.8757906 0.04804776 0.7825237 0.9462538
#> 5      20                0         20 0.8617131 0.04971660 0.7682089 0.9370686
#> 6      20               20         40 0.8617131 0.04971660 0.7682089 0.9370686
#> 7      20               40         60 0.8617131 0.04971660 0.7682089 0.9370686
#> 8      20               60         80 0.8617131 0.04971660 0.7682089 0.9370686
#> 9      40                0         40 0.8478591 0.05131262 0.7542804 0.9277415
#> 10     40               20         60 0.8478591 0.05131262 0.7542804 0.9277415
#> 11     40               40         80 0.8478591 0.05131262 0.7542804 0.9277415
#> 12     40               60        100 0.8478591 0.05131262 0.7542804 0.9277415
#> 13     60                0         60 0.8342249 0.05283905 0.7407091 0.9183124
#> 14     60               20         80 0.8342249 0.05283905 0.7407091 0.9183124
#> 15     60               40        100 0.8342249 0.05283905 0.7407091 0.9183124
#> 16     60               60        120 0.8342249 0.05283905 0.7407091 0.9183124
#> 17     80                0         80 0.8208071 0.05429890 0.7274709 0.9088128
#> 18     80               20        100 0.8208071 0.05429890 0.7274709 0.9088128
#> 19     80               40        120 0.8208071 0.05429890 0.7274709 0.9088128
#> 20     80               60        140 0.8208071 0.05429890 0.7274709 0.9088128
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11798273 0.163120348 0.5859450
#> 2  0.30574618 0.11549516 0.131768343 0.5389983
#> 3  0.26001915 0.11301100 0.105768984 0.4959574
#> 4  0.22113100 0.11006317 0.084267281 0.4565875
#> 5  0.25589195 0.10699068 0.083745378 0.4669402
#> 6  0.21762106 0.10326867 0.066124977 0.4300811
#> 7  0.18507391 0.09950378 0.051684460 0.3964408
#> 8  0.15739448 0.09552975 0.039925450 0.3657567
#> 9  0.18213629 0.09421883 0.039642764 0.3738193
#> 10 0.15489621 0.09003131 0.030197983 0.3451278
#> 11 0.13173012 0.08586391 0.022642806 0.3189629
#> 12 0.11202872 0.08165467 0.016669861 0.2950954
#> 13 0.12963921 0.08137479 0.016528831 0.3013681
#> 14 0.11025053 0.07720991 0.011905637 0.2790387
#> 15 0.09376159 0.07312275 0.008365753 0.2586489
#> 16 0.07973872 0.06909111 0.005711487 0.2400135
#> 17 0.09227334 0.06939353 0.005650810 0.2449152
#> 18 0.07847305 0.06551595 0.003727218 0.2274497
#> 19 0.06673672 0.06175192 0.002364140 0.2114545
#> 20 0.05675565 0.05809522 0.001432949 0.1967848
```
