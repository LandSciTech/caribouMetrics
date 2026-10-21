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
#> 1         0.1 0.3838513 0.01114128 0.3600452 0.3916348
#> 2         0.2 0.3832760 0.01115013 0.3594378 0.3910824
#> 3         0.3 0.3827015 0.01115907 0.3588315 0.3905307
#> 4         0.4 0.3821279 0.01116810 0.3582263 0.3899798
#> 5         0.5 0.3815551 0.01117722 0.3576220 0.3894297
#> 6         0.6 0.3809832 0.01118643 0.3570188 0.3888804
#> 7         0.7 0.3804122 0.01119573 0.3564165 0.3883318
#> 8         0.8 0.3798420 0.01120511 0.3558153 0.3877840
#> 9         0.9 0.3792726 0.01121458 0.3552151 0.3872370
#> 10        1.0 0.3787041 0.01122413 0.3546160 0.3866908

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3819167
#> 2       1        0.1        V5 0.3848238
#> 3       1        0.1        V9 0.3692212
#> 4       1        0.1        V3 0.3576730
#> 5       1        0.1        V4 0.3842554
#> 6       1        0.1        V8 0.3759654
#> 7       1        0.1        V2 0.3894168
#> 8       1        0.1        V6 0.3883133
#> 9       1        0.1        V7 0.3682159
#> 10      1        0.1       V10 0.3922788
#> 11      2        0.2        V7 0.3675828
#> 12      2        0.2        V8 0.3754168
#> 13      2        0.2        V9 0.3686956
#> 14      2        0.2       V10 0.3917341
#> 15      2        0.2        V1 0.3813453
#> 16      2        0.2        V5 0.3841839
#> 17      2        0.2        V2 0.3888374
#> 18      2        0.2        V3 0.3570732
#> 19      2        0.2        V4 0.3836495
#> 20      2        0.2        V6 0.3877870
#> 21      3        0.3        V4 0.3830446
#> 22      3        0.3        V5 0.3835451
#> 23      3        0.3        V3 0.3564743
#> 24      3        0.3        V7 0.3669508
#> 25      3        0.3        V8 0.3748690
#> 26      3        0.3        V9 0.3681708
#> 27      3        0.3        V6 0.3872615
#> 28      3        0.3       V10 0.3911903
#> 29      3        0.3        V1 0.3807747
#> 30      3        0.3        V2 0.3882588
#> 31      4        0.4        V1 0.3802049
#> 32      4        0.4        V9 0.3676467
#> 33      4        0.4        V3 0.3558765
#> 34      4        0.4        V4 0.3824406
#> 35      4        0.4        V5 0.3829074
#> 36      4        0.4        V2 0.3876810
#> 37      4        0.4        V6 0.3867367
#> 38      4        0.4        V7 0.3663199
#> 39      4        0.4        V8 0.3743220
#> 40      4        0.4       V10 0.3906472
#> 41      5        0.5        V8 0.3737758
#> 42      5        0.5        V9 0.3671234
#> 43      5        0.5       V10 0.3901049
#> 44      5        0.5        V1 0.3796360
#> 45      5        0.5        V5 0.3822707
#> 46      5        0.5        V2 0.3871042
#> 47      5        0.5        V3 0.3552797
#> 48      5        0.5        V4 0.3818376
#> 49      5        0.5        V6 0.3862126
#> 50      5        0.5        V7 0.3656901
#> 51      6        0.6        V4 0.3812355
#> 52      6        0.6        V5 0.3816351
#> 53      6        0.6        V7 0.3650613
#> 54      6        0.6        V8 0.3732304
#> 55      6        0.6        V9 0.3666008
#> 56      6        0.6        V6 0.3856892
#> 57      6        0.6       V10 0.3895633
#> 58      6        0.6        V1 0.3790680
#> 59      6        0.6        V2 0.3865282
#> 60      6        0.6        V3 0.3546838
#> 61      7        0.7        V1 0.3785008
#> 62      7        0.7        V3 0.3540890
#> 63      7        0.7        V4 0.3806344
#> 64      7        0.7        V5 0.3810005
#> 65      7        0.7        V9 0.3660789
#> 66      7        0.7        V6 0.3851665
#> 67      7        0.7        V7 0.3644337
#> 68      7        0.7        V8 0.3726857
#> 69      7        0.7        V2 0.3859530
#> 70      7        0.7       V10 0.3890224
#> 71      8        0.8        V9 0.3655578
#> 72      8        0.8       V10 0.3884823
#> 73      8        0.8        V1 0.3779345
#> 74      8        0.8        V5 0.3803670
#> 75      8        0.8        V2 0.3853787
#> 76      8        0.8        V3 0.3534951
#> 77      8        0.8        V4 0.3800342
#> 78      8        0.8        V8 0.3721419
#> 79      8        0.8        V6 0.3846446
#> 80      8        0.8        V7 0.3638071
#> 81      9        0.9        V5 0.3797345
#> 82      9        0.9        V7 0.3631816
#> 83      9        0.9        V8 0.3715989
#> 84      9        0.9        V9 0.3650375
#> 85      9        0.9        V6 0.3841233
#> 86      9        0.9       V10 0.3879430
#> 87      9        0.9        V1 0.3773690
#> 88      9        0.9        V2 0.3848052
#> 89      9        0.9        V3 0.3529023
#> 90      9        0.9        V4 0.3794350
#> 91     10        1.0        V1 0.3768044
#> 92     10        1.0        V4 0.3788367
#> 93     10        1.0        V5 0.3791031
#> 94     10        1.0        V9 0.3645178
#> 95     10        1.0        V3 0.3523104
#> 96     10        1.0        V7 0.3625572
#> 97     10        1.0        V8 0.3710566
#> 98     10        1.0        V2 0.3842327
#> 99     10        1.0        V6 0.3836028
#> 100    10        1.0       V10 0.3874044

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3819167
#> 2       1        0.1        V5 0.3848238
#> 3       1        0.1        V9 0.3692212
#> 4       1        0.1        V3 0.3576730
#> 5       1        0.1        V4 0.3842554
#> 6       1        0.1        V8 0.3759654
#> 7       1        0.1        V2 0.3894168
#> 8       1        0.1        V6 0.3883133
#> 9       1        0.1        V7 0.3682159
#> 10      1        0.1       V10 0.3922788
#> 11      2        0.2        V7 0.3675828
#> 12      2        0.2        V8 0.3754168
#> 13      2        0.2        V9 0.3686956
#> 14      2        0.2       V10 0.3917341
#> 15      2        0.2        V1 0.3813453
#> 16      2        0.2        V5 0.3841839
#> 17      2        0.2        V2 0.3888374
#> 18      2        0.2        V3 0.3570732
#> 19      2        0.2        V4 0.3836495
#> 20      2        0.2        V6 0.3877870
#> 21      3        0.3        V4 0.3830446
#> 22      3        0.3        V5 0.3835451
#> 23      3        0.3        V3 0.3564743
#> 24      3        0.3        V7 0.3669508
#> 25      3        0.3        V8 0.3748690
#> 26      3        0.3        V9 0.3681708
#> 27      3        0.3        V6 0.3872615
#> 28      3        0.3       V10 0.3911903
#> 29      3        0.3        V1 0.3807747
#> 30      3        0.3        V2 0.3882588
#> 31      4        0.4        V1 0.3802049
#> 32      4        0.4        V9 0.3676467
#> 33      4        0.4        V3 0.3558765
#> 34      4        0.4        V4 0.3824406
#> 35      4        0.4        V5 0.3829074
#> 36      4        0.4        V2 0.3876810
#> 37      4        0.4        V6 0.3867367
#> 38      4        0.4        V7 0.3663199
#> 39      4        0.4        V8 0.3743220
#> 40      4        0.4       V10 0.3906472
#> 41      5        0.5        V8 0.3737758
#> 42      5        0.5        V9 0.3671234
#> 43      5        0.5       V10 0.3901049
#> 44      5        0.5        V1 0.3796360
#> 45      5        0.5        V5 0.3822707
#> 46      5        0.5        V2 0.3871042
#> 47      5        0.5        V3 0.3552797
#> 48      5        0.5        V4 0.3818376
#> 49      5        0.5        V6 0.3862126
#> 50      5        0.5        V7 0.3656901
#> 51      6        0.6        V4 0.3812355
#> 52      6        0.6        V5 0.3816351
#> 53      6        0.6        V7 0.3650613
#> 54      6        0.6        V8 0.3732304
#> 55      6        0.6        V9 0.3666008
#> 56      6        0.6        V6 0.3856892
#> 57      6        0.6       V10 0.3895633
#> 58      6        0.6        V1 0.3790680
#> 59      6        0.6        V2 0.3865282
#> 60      6        0.6        V3 0.3546838
#> 61      7        0.7        V1 0.3785008
#> 62      7        0.7        V3 0.3540890
#> 63      7        0.7        V4 0.3806344
#> 64      7        0.7        V5 0.3810005
#> 65      7        0.7        V9 0.3660789
#> 66      7        0.7        V6 0.3851665
#> 67      7        0.7        V7 0.3644337
#> 68      7        0.7        V8 0.3726857
#> 69      7        0.7        V2 0.3859530
#> 70      7        0.7       V10 0.3890224
#> 71      8        0.8        V9 0.3655578
#> 72      8        0.8       V10 0.3884823
#> 73      8        0.8        V1 0.3779345
#> 74      8        0.8        V5 0.3803670
#> 75      8        0.8        V2 0.3853787
#> 76      8        0.8        V3 0.3534951
#> 77      8        0.8        V4 0.3800342
#> 78      8        0.8        V8 0.3721419
#> 79      8        0.8        V6 0.3846446
#> 80      8        0.8        V7 0.3638071
#> 81      9        0.9        V5 0.3797345
#> 82      9        0.9        V7 0.3631816
#> 83      9        0.9        V8 0.3715989
#> 84      9        0.9        V9 0.3650375
#> 85      9        0.9        V6 0.3841233
#> 86      9        0.9       V10 0.3879430
#> 87      9        0.9        V1 0.3773690
#> 88      9        0.9        V2 0.3848052
#> 89      9        0.9        V3 0.3529023
#> 90      9        0.9        V4 0.3794350
#> 91     10        1.0        V1 0.3768044
#> 92     10        1.0        V4 0.3788367
#> 93     10        1.0        V5 0.3791031
#> 94     10        1.0        V9 0.3645178
#> 95     10        1.0        V3 0.3523104
#> 96     10        1.0        V7 0.3625572
#> 97     10        1.0        V8 0.3710566
#> 98     10        1.0        V2 0.3842327
#> 99     10        1.0        V6 0.3836028
#> 100    10        1.0       V10 0.3874044


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
#> 1       0                0          0 0.8757906 0.04499905 0.7862029 0.9410379
#> 2       0               20         20 0.8757906 0.04499905 0.7862029 0.9410379
#> 3       0               40         40 0.8757906 0.04499905 0.7862029 0.9410379
#> 4       0               60         60 0.8757906 0.04499905 0.7862029 0.9410379
#> 5      20                0         20 0.8617131 0.04464353 0.7747324 0.9262589
#> 6      20               20         40 0.8617131 0.04464353 0.7747324 0.9262589
#> 7      20               40         60 0.8617131 0.04464353 0.7747324 0.9262589
#> 8      20               60         80 0.8617131 0.04464353 0.7747324 0.9262589
#> 9      40                0         40 0.8478591 0.04424530 0.7635059 0.9112796
#> 10     40               20         60 0.8478591 0.04424530 0.7635059 0.9112796
#> 11     40               40         80 0.8478591 0.04424530 0.7635059 0.9112796
#> 12     40               60        100 0.8478591 0.04424530 0.7635059 0.9112796
#> 13     60                0         60 0.8342249 0.04382722 0.7525088 0.8962055
#> 14     60               20         80 0.8342249 0.04382722 0.7525088 0.8962055
#> 15     60               40        100 0.8342249 0.04382722 0.7525088 0.8962055
#> 16     60               60        120 0.8342249 0.04382722 0.7525088 0.8962055
#> 17     80                0         80 0.8208071 0.04340596 0.7417285 0.8811114
#> 18     80               20        100 0.8208071 0.04340596 0.7417285 0.8811114
#> 19     80               40        120 0.8208071 0.04340596 0.7417285 0.8811114
#> 20     80               60        140 0.8208071 0.04340596 0.7417285 0.8811114
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11425058 0.163261197 0.5509282
#> 2  0.30574618 0.10785533 0.133819013 0.4810979
#> 3  0.26001915 0.10285627 0.109108678 0.4208093
#> 4  0.22113100 0.09859259 0.088415952 0.3689118
#> 5  0.25589195 0.10307263 0.097872383 0.4414941
#> 6  0.21762106 0.09538639 0.079027859 0.3867081
#> 7  0.18507391 0.08903074 0.063326128 0.3395817
#> 8  0.15739448 0.08355096 0.050300870 0.2990363
#> 9  0.18213629 0.09091166 0.056236125 0.3557413
#> 10 0.15489621 0.08322120 0.044443304 0.3129437
#> 11 0.13173012 0.07666761 0.034745628 0.2760925
#> 12 0.11202872 0.07093608 0.026828718 0.2442906
#> 13 0.12963921 0.07896022 0.030418907 0.2887380
#> 14 0.11025053 0.07182077 0.023319656 0.2552136
#> 15 0.09376159 0.06561235 0.017604567 0.2262253
#> 16 0.07973872 0.06011326 0.013056666 0.2010608
#> 17 0.09227334 0.06782446 0.015103229 0.2361896
#> 18 0.07847305 0.06143305 0.011086869 0.2097235
#> 19 0.06673672 0.05579051 0.007960876 0.1866790
#> 20 0.05675565 0.05073911 0.005572064 0.1665034
```
