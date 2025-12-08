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
#> 1         0.1 0.3838513 0.01535000 0.3568856 0.4015619
#> 2         0.2 0.3832760 0.01530314 0.3563637 0.4009014
#> 3         0.3 0.3827015 0.01525671 0.3558426 0.4002420
#> 4         0.4 0.3821279 0.01521073 0.3553222 0.3995838
#> 5         0.5 0.3815551 0.01516519 0.3548026 0.3989266
#> 6         0.6 0.3809832 0.01512009 0.3542838 0.3982705
#> 7         0.7 0.3804122 0.01507543 0.3537657 0.3976155
#> 8         0.8 0.3798420 0.01503120 0.3532484 0.3969616
#> 9         0.9 0.3792726 0.01498742 0.3527319 0.3963087
#> 10        1.0 0.3787041 0.01494408 0.3522161 0.3956570

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3852007
#> 2       1        0.1        V5 0.3835769
#> 3       1        0.1        V9 0.3993506
#> 4       1        0.1        V3 0.3667277
#> 5       1        0.1        V4 0.3802343
#> 6       1        0.1        V8 0.3859724
#> 7       1        0.1        V2 0.3540283
#> 8       1        0.1        V6 0.4022039
#> 9       1        0.1        V7 0.3958172
#> 10      1        0.1       V10 0.3688268
#> 11      2        0.2        V7 0.3951577
#> 12      2        0.2        V8 0.3852718
#> 13      2        0.2        V9 0.3988259
#> 14      2        0.2       V10 0.3682799
#> 15      2        0.2        V1 0.3846671
#> 16      2        0.2        V5 0.3829263
#> 17      2        0.2        V2 0.3534871
#> 18      2        0.2        V3 0.3662720
#> 19      2        0.2        V4 0.3796998
#> 20      2        0.2        V6 0.4015040
#> 21      3        0.3        V4 0.3791659
#> 22      3        0.3        V5 0.3822768
#> 23      3        0.3        V3 0.3658169
#> 24      3        0.3        V7 0.3944993
#> 25      3        0.3        V8 0.3845724
#> 26      3        0.3        V9 0.3983018
#> 27      3        0.3        V6 0.4008053
#> 28      3        0.3       V10 0.3677338
#> 29      3        0.3        V1 0.3841342
#> 30      3        0.3        V2 0.3529468
#> 31      4        0.4        V1 0.3836021
#> 32      4        0.4        V9 0.3977784
#> 33      4        0.4        V3 0.3653624
#> 34      4        0.4        V4 0.3786329
#> 35      4        0.4        V5 0.3816284
#> 36      4        0.4        V2 0.3524074
#> 37      4        0.4        V6 0.4001079
#> 38      4        0.4        V7 0.3938421
#> 39      4        0.4        V8 0.3838744
#> 40      4        0.4       V10 0.3671885
#> 41      5        0.5        V8 0.3831775
#> 42      5        0.5        V9 0.3972556
#> 43      5        0.5       V10 0.3666441
#> 44      5        0.5        V1 0.3830707
#> 45      5        0.5        V5 0.3809811
#> 46      5        0.5        V2 0.3518687
#> 47      5        0.5        V3 0.3649084
#> 48      5        0.5        V4 0.3781005
#> 49      5        0.5        V6 0.3994117
#> 50      5        0.5        V7 0.3931859
#> 51      6        0.6        V4 0.3775689
#> 52      6        0.6        V5 0.3803349
#> 53      6        0.6        V7 0.3925308
#> 54      6        0.6        V8 0.3824820
#> 55      6        0.6        V9 0.3967336
#> 56      6        0.6        V6 0.3987167
#> 57      6        0.6       V10 0.3661004
#> 58      6        0.6        V1 0.3825400
#> 59      6        0.6        V2 0.3513309
#> 60      6        0.6        V3 0.3644550
#> 61      7        0.7        V1 0.3820101
#> 62      7        0.7        V3 0.3640022
#> 63      7        0.7        V4 0.3770381
#> 64      7        0.7        V5 0.3796898
#> 65      7        0.7        V9 0.3962123
#> 66      7        0.7        V6 0.3980229
#> 67      7        0.7        V7 0.3918768
#> 68      7        0.7        V8 0.3817877
#> 69      7        0.7        V2 0.3507939
#> 70      7        0.7       V10 0.3655576
#> 71      8        0.8        V9 0.3956916
#> 72      8        0.8       V10 0.3650155
#> 73      8        0.8        V1 0.3814809
#> 74      8        0.8        V5 0.3790458
#> 75      8        0.8        V2 0.3502577
#> 76      8        0.8        V3 0.3635499
#> 77      8        0.8        V4 0.3765080
#> 78      8        0.8        V8 0.3810947
#> 79      8        0.8        V6 0.3973303
#> 80      8        0.8        V7 0.3912239
#> 81      9        0.9        V5 0.3784028
#> 82      9        0.9        V7 0.3905721
#> 83      9        0.9        V8 0.3804029
#> 84      9        0.9        V9 0.3951717
#> 85      9        0.9        V6 0.3966389
#> 86      9        0.9       V10 0.3644743
#> 87      9        0.9        V1 0.3809525
#> 88      9        0.9        V2 0.3497223
#> 89      9        0.9        V3 0.3630982
#> 90      9        0.9        V4 0.3759787
#> 91     10        1.0        V1 0.3804248
#> 92     10        1.0        V4 0.3754501
#> 93     10        1.0        V5 0.3777610
#> 94     10        1.0        V9 0.3946524
#> 95     10        1.0        V3 0.3626470
#> 96     10        1.0        V7 0.3899214
#> 97     10        1.0        V8 0.3797124
#> 98     10        1.0        V2 0.3491878
#> 99     10        1.0        V6 0.3959487
#> 100    10        1.0       V10 0.3639338

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3852007
#> 2       1        0.1        V5 0.3835769
#> 3       1        0.1        V9 0.3993506
#> 4       1        0.1        V3 0.3667277
#> 5       1        0.1        V4 0.3802343
#> 6       1        0.1        V8 0.3859724
#> 7       1        0.1        V2 0.3540283
#> 8       1        0.1        V6 0.4022039
#> 9       1        0.1        V7 0.3958172
#> 10      1        0.1       V10 0.3688268
#> 11      2        0.2        V7 0.3951577
#> 12      2        0.2        V8 0.3852718
#> 13      2        0.2        V9 0.3988259
#> 14      2        0.2       V10 0.3682799
#> 15      2        0.2        V1 0.3846671
#> 16      2        0.2        V5 0.3829263
#> 17      2        0.2        V2 0.3534871
#> 18      2        0.2        V3 0.3662720
#> 19      2        0.2        V4 0.3796998
#> 20      2        0.2        V6 0.4015040
#> 21      3        0.3        V4 0.3791659
#> 22      3        0.3        V5 0.3822768
#> 23      3        0.3        V3 0.3658169
#> 24      3        0.3        V7 0.3944993
#> 25      3        0.3        V8 0.3845724
#> 26      3        0.3        V9 0.3983018
#> 27      3        0.3        V6 0.4008053
#> 28      3        0.3       V10 0.3677338
#> 29      3        0.3        V1 0.3841342
#> 30      3        0.3        V2 0.3529468
#> 31      4        0.4        V1 0.3836021
#> 32      4        0.4        V9 0.3977784
#> 33      4        0.4        V3 0.3653624
#> 34      4        0.4        V4 0.3786329
#> 35      4        0.4        V5 0.3816284
#> 36      4        0.4        V2 0.3524074
#> 37      4        0.4        V6 0.4001079
#> 38      4        0.4        V7 0.3938421
#> 39      4        0.4        V8 0.3838744
#> 40      4        0.4       V10 0.3671885
#> 41      5        0.5        V8 0.3831775
#> 42      5        0.5        V9 0.3972556
#> 43      5        0.5       V10 0.3666441
#> 44      5        0.5        V1 0.3830707
#> 45      5        0.5        V5 0.3809811
#> 46      5        0.5        V2 0.3518687
#> 47      5        0.5        V3 0.3649084
#> 48      5        0.5        V4 0.3781005
#> 49      5        0.5        V6 0.3994117
#> 50      5        0.5        V7 0.3931859
#> 51      6        0.6        V4 0.3775689
#> 52      6        0.6        V5 0.3803349
#> 53      6        0.6        V7 0.3925308
#> 54      6        0.6        V8 0.3824820
#> 55      6        0.6        V9 0.3967336
#> 56      6        0.6        V6 0.3987167
#> 57      6        0.6       V10 0.3661004
#> 58      6        0.6        V1 0.3825400
#> 59      6        0.6        V2 0.3513309
#> 60      6        0.6        V3 0.3644550
#> 61      7        0.7        V1 0.3820101
#> 62      7        0.7        V3 0.3640022
#> 63      7        0.7        V4 0.3770381
#> 64      7        0.7        V5 0.3796898
#> 65      7        0.7        V9 0.3962123
#> 66      7        0.7        V6 0.3980229
#> 67      7        0.7        V7 0.3918768
#> 68      7        0.7        V8 0.3817877
#> 69      7        0.7        V2 0.3507939
#> 70      7        0.7       V10 0.3655576
#> 71      8        0.8        V9 0.3956916
#> 72      8        0.8       V10 0.3650155
#> 73      8        0.8        V1 0.3814809
#> 74      8        0.8        V5 0.3790458
#> 75      8        0.8        V2 0.3502577
#> 76      8        0.8        V3 0.3635499
#> 77      8        0.8        V4 0.3765080
#> 78      8        0.8        V8 0.3810947
#> 79      8        0.8        V6 0.3973303
#> 80      8        0.8        V7 0.3912239
#> 81      9        0.9        V5 0.3784028
#> 82      9        0.9        V7 0.3905721
#> 83      9        0.9        V8 0.3804029
#> 84      9        0.9        V9 0.3951717
#> 85      9        0.9        V6 0.3966389
#> 86      9        0.9       V10 0.3644743
#> 87      9        0.9        V1 0.3809525
#> 88      9        0.9        V2 0.3497223
#> 89      9        0.9        V3 0.3630982
#> 90      9        0.9        V4 0.3759787
#> 91     10        1.0        V1 0.3804248
#> 92     10        1.0        V4 0.3754501
#> 93     10        1.0        V5 0.3777610
#> 94     10        1.0        V9 0.3946524
#> 95     10        1.0        V3 0.3626470
#> 96     10        1.0        V7 0.3899214
#> 97     10        1.0        V8 0.3797124
#> 98     10        1.0        V2 0.3491878
#> 99     10        1.0        V6 0.3959487
#> 100    10        1.0       V10 0.3639338


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
#> 1       0                0          0 0.8757906 0.05298835 0.7666786 0.9435786
#> 2       0               20         20 0.8757906 0.05298835 0.7666786 0.9435786
#> 3       0               40         40 0.8757906 0.05298835 0.7666786 0.9435786
#> 4       0               60         60 0.8757906 0.05298835 0.7666786 0.9435786
#> 5      20                0         20 0.8617131 0.05605742 0.7478409 0.9314673
#> 6      20               20         40 0.8617131 0.05605742 0.7478409 0.9314673
#> 7      20               40         60 0.8617131 0.05605742 0.7478409 0.9314673
#> 8      20               60         80 0.8617131 0.05605742 0.7478409 0.9314673
#> 9      40                0         40 0.8478591 0.05905932 0.7297085 0.9191746
#> 10     40               20         60 0.8478591 0.05905932 0.7297085 0.9191746
#> 11     40               40         80 0.8478591 0.05905932 0.7297085 0.9191746
#> 12     40               60        100 0.8478591 0.05905932 0.7297085 0.9191746
#> 13     60                0         60 0.8342249 0.06198516 0.7122100 0.9067721
#> 14     60               20         80 0.8342249 0.06198516 0.7122100 0.9067721
#> 15     60               40        100 0.8342249 0.06198516 0.7122100 0.9067721
#> 16     60               60        120 0.8342249 0.06198516 0.7122100 0.9067721
#> 17     80                0         80 0.8208071 0.06482915 0.6952891 0.8943131
#> 18     80               20        100 0.8208071 0.06482915 0.6952891 0.8943131
#> 19     80               40        120 0.8208071 0.06482915 0.6952891 0.8943131
#> 20     80               60        140 0.8208071 0.06482915 0.6952891 0.8943131
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12411860 0.1220264590 0.5819015
#> 2  0.30574618 0.11977263 0.0958139725 0.5238618
#> 3  0.26001915 0.11563238 0.0745067289 0.4719343
#> 4  0.22113100 0.11135214 0.0572803733 0.4256176
#> 5  0.25589195 0.11245742 0.0641771762 0.4727616
#> 6  0.21762106 0.10712881 0.0489750557 0.4263549
#> 7  0.18507391 0.10212807 0.0368321632 0.3850292
#> 8  0.15739448 0.09723616 0.0272321213 0.3482509
#> 9  0.18213629 0.09946708 0.0310464277 0.3856867
#> 10 0.15489621 0.09399412 0.0227029559 0.3488360
#> 11 0.13173012 0.08888140 0.0162454482 0.3160345
#> 12 0.11202872 0.08398490 0.0113341868 0.2868130
#> 13 0.12963921 0.08664640 0.0132600696 0.3165566
#> 14 0.11025053 0.08145235 0.0091011044 0.2872784
#> 15 0.09376159 0.07660602 0.0060482656 0.2611594
#> 16 0.07973872 0.07200862 0.0038710009 0.2378136
#> 17 0.09227334 0.07475457 0.0047062074 0.2615757
#> 18 0.07847305 0.07002994 0.0029396703 0.2381861
#> 19 0.06673672 0.06561671 0.0017504433 0.2172298
#> 20 0.05675565 0.06144339 0.0009859507 0.1983985
```
