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
#> 1         0.1 0.3838513 0.03199635 0.3389177 0.4272130
#> 2         0.2 0.3832760 0.03192286 0.3383983 0.4264816
#> 3         0.3 0.3827015 0.03184961 0.3378796 0.4257515
#> 4         0.4 0.3821279 0.03177659 0.3373618 0.4250226
#> 5         0.5 0.3815551 0.03170380 0.3368447 0.4242949
#> 6         0.6 0.3809832 0.03163124 0.3363285 0.4235685
#> 7         0.7 0.3804122 0.03155892 0.3358130 0.4228433
#> 8         0.8 0.3798420 0.03148682 0.3352983 0.4221194
#> 9         0.9 0.3792726 0.03141495 0.3347845 0.4213967
#> 10        1.0 0.3787041 0.03134331 0.3342714 0.4206753

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3855492
#> 2       1        0.1        V5 0.3621332
#> 3       1        0.1        V9 0.4310413
#> 4       1        0.1        V3 0.3392183
#> 5       1        0.1        V4 0.3432859
#> 6       1        0.1        V8 0.3388304
#> 7       1        0.1        V2 0.3939293
#> 8       1        0.1        V6 0.4140267
#> 9       1        0.1        V7 0.3704120
#> 10      1        0.1       V10 0.3936665
#> 11      2        0.2        V7 0.3699405
#> 12      2        0.2        V8 0.3383000
#> 13      2        0.2        V9 0.4302970
#> 14      2        0.2       V10 0.3930198
#> 15      2        0.2        V1 0.3849721
#> 16      2        0.2        V5 0.3615512
#> 17      2        0.2        V2 0.3932583
#> 18      2        0.2        V3 0.3387368
#> 19      2        0.2        V4 0.3426896
#> 20      2        0.2        V6 0.4133398
#> 21      3        0.3        V4 0.3420944
#> 22      3        0.3        V5 0.3609701
#> 23      3        0.3        V3 0.3382560
#> 24      3        0.3        V7 0.3694696
#> 25      3        0.3        V8 0.3377704
#> 26      3        0.3        V9 0.4295539
#> 27      3        0.3        V6 0.4126541
#> 28      3        0.3       V10 0.3923742
#> 29      3        0.3        V1 0.3843958
#> 30      3        0.3        V2 0.3925885
#> 31      4        0.4        V1 0.3838203
#> 32      4        0.4        V9 0.4288121
#> 33      4        0.4        V3 0.3377758
#> 34      4        0.4        V4 0.3415002
#> 35      4        0.4        V5 0.3603899
#> 36      4        0.4        V2 0.3919198
#> 37      4        0.4        V6 0.4119696
#> 38      4        0.4        V7 0.3689993
#> 39      4        0.4        V8 0.3372416
#> 40      4        0.4       V10 0.3917296
#> 41      5        0.5        V8 0.3367136
#> 42      5        0.5        V9 0.4280716
#> 43      5        0.5       V10 0.3910861
#> 44      5        0.5        V1 0.3832458
#> 45      5        0.5        V5 0.3598107
#> 46      5        0.5        V2 0.3912522
#> 47      5        0.5        V3 0.3372963
#> 48      5        0.5        V4 0.3409070
#> 49      5        0.5        V6 0.4112862
#> 50      5        0.5        V7 0.3685295
#> 51      6        0.6        V4 0.3403148
#> 52      6        0.6        V5 0.3592324
#> 53      6        0.6        V7 0.3680604
#> 54      6        0.6        V8 0.3361865
#> 55      6        0.6        V9 0.4273324
#> 56      6        0.6        V6 0.4106039
#> 57      6        0.6       V10 0.3904437
#> 58      6        0.6        V1 0.3826720
#> 59      6        0.6        V2 0.3905858
#> 60      6        0.6        V3 0.3368175
#> 61      7        0.7        V1 0.3820992
#> 62      7        0.7        V3 0.3363394
#> 63      7        0.7        V4 0.3397237
#> 64      7        0.7        V5 0.3586551
#> 65      7        0.7        V9 0.4265945
#> 66      7        0.7        V6 0.4099227
#> 67      7        0.7        V7 0.3675919
#> 68      7        0.7        V8 0.3356602
#> 69      7        0.7        V2 0.3899205
#> 70      7        0.7       V10 0.3898023
#> 71      8        0.8        V9 0.4258578
#> 72      8        0.8       V10 0.3891619
#> 73      8        0.8        V1 0.3815272
#> 74      8        0.8        V5 0.3580786
#> 75      8        0.8        V2 0.3892564
#> 76      8        0.8        V3 0.3358619
#> 77      8        0.8        V4 0.3391336
#> 78      8        0.8        V8 0.3351347
#> 79      8        0.8        V6 0.4092427
#> 80      8        0.8        V7 0.3671240
#> 81      9        0.9        V5 0.3575031
#> 82      9        0.9        V7 0.3666567
#> 83      9        0.9        V8 0.3346101
#> 84      9        0.9        V9 0.4251224
#> 85      9        0.9        V6 0.4085638
#> 86      9        0.9       V10 0.3885226
#> 87      9        0.9        V1 0.3809561
#> 88      9        0.9        V2 0.3885933
#> 89      9        0.9        V3 0.3353851
#> 90      9        0.9        V4 0.3385446
#> 91     10        1.0        V1 0.3803858
#> 92     10        1.0        V4 0.3379565
#> 93     10        1.0        V5 0.3569285
#> 94     10        1.0        V9 0.4243883
#> 95     10        1.0        V3 0.3349090
#> 96     10        1.0        V7 0.3661899
#> 97     10        1.0        V8 0.3340862
#> 98     10        1.0        V2 0.3879314
#> 99     10        1.0        V6 0.4078860
#> 100    10        1.0       V10 0.3878844

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3855492
#> 2       1        0.1        V5 0.3621332
#> 3       1        0.1        V9 0.4310413
#> 4       1        0.1        V3 0.3392183
#> 5       1        0.1        V4 0.3432859
#> 6       1        0.1        V8 0.3388304
#> 7       1        0.1        V2 0.3939293
#> 8       1        0.1        V6 0.4140267
#> 9       1        0.1        V7 0.3704120
#> 10      1        0.1       V10 0.3936665
#> 11      2        0.2        V7 0.3699405
#> 12      2        0.2        V8 0.3383000
#> 13      2        0.2        V9 0.4302970
#> 14      2        0.2       V10 0.3930198
#> 15      2        0.2        V1 0.3849721
#> 16      2        0.2        V5 0.3615512
#> 17      2        0.2        V2 0.3932583
#> 18      2        0.2        V3 0.3387368
#> 19      2        0.2        V4 0.3426896
#> 20      2        0.2        V6 0.4133398
#> 21      3        0.3        V4 0.3420944
#> 22      3        0.3        V5 0.3609701
#> 23      3        0.3        V3 0.3382560
#> 24      3        0.3        V7 0.3694696
#> 25      3        0.3        V8 0.3377704
#> 26      3        0.3        V9 0.4295539
#> 27      3        0.3        V6 0.4126541
#> 28      3        0.3       V10 0.3923742
#> 29      3        0.3        V1 0.3843958
#> 30      3        0.3        V2 0.3925885
#> 31      4        0.4        V1 0.3838203
#> 32      4        0.4        V9 0.4288121
#> 33      4        0.4        V3 0.3377758
#> 34      4        0.4        V4 0.3415002
#> 35      4        0.4        V5 0.3603899
#> 36      4        0.4        V2 0.3919198
#> 37      4        0.4        V6 0.4119696
#> 38      4        0.4        V7 0.3689993
#> 39      4        0.4        V8 0.3372416
#> 40      4        0.4       V10 0.3917296
#> 41      5        0.5        V8 0.3367136
#> 42      5        0.5        V9 0.4280716
#> 43      5        0.5       V10 0.3910861
#> 44      5        0.5        V1 0.3832458
#> 45      5        0.5        V5 0.3598107
#> 46      5        0.5        V2 0.3912522
#> 47      5        0.5        V3 0.3372963
#> 48      5        0.5        V4 0.3409070
#> 49      5        0.5        V6 0.4112862
#> 50      5        0.5        V7 0.3685295
#> 51      6        0.6        V4 0.3403148
#> 52      6        0.6        V5 0.3592324
#> 53      6        0.6        V7 0.3680604
#> 54      6        0.6        V8 0.3361865
#> 55      6        0.6        V9 0.4273324
#> 56      6        0.6        V6 0.4106039
#> 57      6        0.6       V10 0.3904437
#> 58      6        0.6        V1 0.3826720
#> 59      6        0.6        V2 0.3905858
#> 60      6        0.6        V3 0.3368175
#> 61      7        0.7        V1 0.3820992
#> 62      7        0.7        V3 0.3363394
#> 63      7        0.7        V4 0.3397237
#> 64      7        0.7        V5 0.3586551
#> 65      7        0.7        V9 0.4265945
#> 66      7        0.7        V6 0.4099227
#> 67      7        0.7        V7 0.3675919
#> 68      7        0.7        V8 0.3356602
#> 69      7        0.7        V2 0.3899205
#> 70      7        0.7       V10 0.3898023
#> 71      8        0.8        V9 0.4258578
#> 72      8        0.8       V10 0.3891619
#> 73      8        0.8        V1 0.3815272
#> 74      8        0.8        V5 0.3580786
#> 75      8        0.8        V2 0.3892564
#> 76      8        0.8        V3 0.3358619
#> 77      8        0.8        V4 0.3391336
#> 78      8        0.8        V8 0.3351347
#> 79      8        0.8        V6 0.4092427
#> 80      8        0.8        V7 0.3671240
#> 81      9        0.9        V5 0.3575031
#> 82      9        0.9        V7 0.3666567
#> 83      9        0.9        V8 0.3346101
#> 84      9        0.9        V9 0.4251224
#> 85      9        0.9        V6 0.4085638
#> 86      9        0.9       V10 0.3885226
#> 87      9        0.9        V1 0.3809561
#> 88      9        0.9        V2 0.3885933
#> 89      9        0.9        V3 0.3353851
#> 90      9        0.9        V4 0.3385446
#> 91     10        1.0        V1 0.3803858
#> 92     10        1.0        V4 0.3379565
#> 93     10        1.0        V5 0.3569285
#> 94     10        1.0        V9 0.4243883
#> 95     10        1.0        V3 0.3349090
#> 96     10        1.0        V7 0.3661899
#> 97     10        1.0        V8 0.3340862
#> 98     10        1.0        V2 0.3879314
#> 99     10        1.0        V6 0.4078860
#> 100    10        1.0       V10 0.3878844


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
#> 1       0                0          0 0.8757906 0.04649502 0.7920422 0.9505831
#> 2       0               20         20 0.8757906 0.04649502 0.7920422 0.9505831
#> 3       0               40         40 0.8757906 0.04649502 0.7920422 0.9505831
#> 4       0               60         60 0.8757906 0.04649502 0.7920422 0.9505831
#> 5      20                0         20 0.8617131 0.04864335 0.7718346 0.9385557
#> 6      20               20         40 0.8617131 0.04864335 0.7718346 0.9385557
#> 7      20               40         60 0.8617131 0.04864335 0.7718346 0.9385557
#> 8      20               60         80 0.8617131 0.04864335 0.7718346 0.9385557
#> 9      40                0         40 0.8478591 0.05060321 0.7523979 0.9262906
#> 10     40               20         60 0.8478591 0.05060321 0.7523979 0.9262906
#> 11     40               40         80 0.8478591 0.05060321 0.7523979 0.9262906
#> 12     40               60        100 0.8478591 0.05060321 0.7523979 0.9262906
#> 13     60                0         60 0.8342249 0.05240255 0.7336512 0.9138756
#> 14     60               20         80 0.8342249 0.05240255 0.7336512 0.9138756
#> 15     60               40        100 0.8342249 0.05240255 0.7336512 0.9138756
#> 16     60               60        120 0.8342249 0.05240255 0.7336512 0.9138756
#> 17     80                0         80 0.8208071 0.05406248 0.7155319 0.9013746
#> 18     80               20        100 0.8208071 0.05406248 0.7155319 0.9013746
#> 19     80               40        120 0.8208071 0.05406248 0.7155319 0.9013746
#> 20     80               60        140 0.8208071 0.05406248 0.7155319 0.9013746
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12909261 0.173283316 0.6257995
#> 2  0.30574618 0.12614351 0.115417288 0.5449540
#> 3  0.26001915 0.12197293 0.074952426 0.4749718
#> 4  0.22113100 0.11686656 0.047050085 0.4147665
#> 5  0.25589195 0.11938116 0.101968945 0.5175897
#> 6  0.21762106 0.11302805 0.065625754 0.4514008
#> 7  0.18507391 0.10657460 0.040707596 0.3945325
#> 8  0.15739448 0.10007239 0.024053056 0.3457505
#> 9  0.18213629 0.10755629 0.057278013 0.4291311
#> 10 0.15489621 0.09989481 0.035073442 0.3754265
#> 11 0.13173012 0.09263696 0.020374019 0.3293599
#> 12 0.11202872 0.08575204 0.011039042 0.2898013
#> 13 0.12963921 0.09534673 0.030082656 0.3573855
#> 14 0.11025053 0.08746842 0.017154880 0.3138763
#> 15 0.09376159 0.08020566 0.009070065 0.2764836
#> 16 0.07973872 0.07349909 0.004337295 0.2442470
#> 17 0.09227334 0.08363392 0.014350918 0.2992451
#> 18 0.07847305 0.07609481 0.007388332 0.2638842
#> 19 0.06673672 0.06923816 0.003410283 0.2333535
#> 20 0.05675565 0.06299508 0.001361490 0.2068647
```
