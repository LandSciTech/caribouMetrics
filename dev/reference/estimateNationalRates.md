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
#> 1         0.1 0.3838513 0.02055157 0.3544563 0.4151907
#> 2         0.2 0.3832760 0.02052603 0.3539046 0.4145432
#> 3         0.3 0.3827015 0.02050065 0.3533537 0.4138967
#> 4         0.4 0.3821279 0.02047543 0.3528038 0.4132512
#> 5         0.5 0.3815551 0.02045038 0.3522547 0.4126068
#> 6         0.6 0.3809832 0.02042549 0.3517064 0.4119633
#> 7         0.7 0.3804122 0.02040075 0.3511590 0.4113209
#> 8         0.8 0.3798420 0.02037618 0.3506124 0.4106794
#> 9         0.9 0.3792726 0.02035176 0.3500667 0.4100390
#> 10        1.0 0.3787041 0.02032751 0.3495219 0.4093995

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4048663
#> 2       1        0.1        V5 0.3867305
#> 3       1        0.1        V9 0.3722622
#> 4       1        0.1        V3 0.3801713
#> 5       1        0.1        V4 0.3858921
#> 6       1        0.1        V8 0.3697280
#> 7       1        0.1        V2 0.3995403
#> 8       1        0.1        V6 0.3533400
#> 9       1        0.1        V7 0.3583010
#> 10      1        0.1       V10 0.4181881
#> 11      2        0.2        V7 0.3577115
#> 12      2        0.2        V8 0.3692235
#> 13      2        0.2        V9 0.3716832
#> 14      2        0.2       V10 0.4175362
#> 15      2        0.2        V1 0.4042340
#> 16      2        0.2        V5 0.3862731
#> 17      2        0.2        V2 0.3989929
#> 18      2        0.2        V3 0.3796451
#> 19      2        0.2        V4 0.3853183
#> 20      2        0.2        V6 0.3527993
#> 21      3        0.3        V4 0.3847454
#> 22      3        0.3        V5 0.3858162
#> 23      3        0.3        V3 0.3791197
#> 24      3        0.3        V7 0.3571230
#> 25      3        0.3        V8 0.3687197
#> 26      3        0.3        V9 0.3711051
#> 27      3        0.3        V6 0.3522595
#> 28      3        0.3       V10 0.4168853
#> 29      3        0.3        V1 0.4036026
#> 30      3        0.3        V2 0.3984462
#> 31      4        0.4        V1 0.4029722
#> 32      4        0.4        V9 0.3705279
#> 33      4        0.4        V3 0.3785949
#> 34      4        0.4        V4 0.3841733
#> 35      4        0.4        V5 0.3853598
#> 36      4        0.4        V2 0.3979003
#> 37      4        0.4        V6 0.3517204
#> 38      4        0.4        V7 0.3565354
#> 39      4        0.4        V8 0.3682166
#> 40      4        0.4       V10 0.4162355
#> 41      5        0.5        V8 0.3677141
#> 42      5        0.5        V9 0.3699516
#> 43      5        0.5       V10 0.4155867
#> 44      5        0.5        V1 0.4023427
#> 45      5        0.5        V5 0.3849040
#> 46      5        0.5        V2 0.3973551
#> 47      5        0.5        V3 0.3780710
#> 48      5        0.5        V4 0.3836021
#> 49      5        0.5        V6 0.3511822
#> 50      5        0.5        V7 0.3559488
#> 51      6        0.6        V4 0.3830317
#> 52      6        0.6        V5 0.3844487
#> 53      6        0.6        V7 0.3553631
#> 54      6        0.6        V8 0.3672124
#> 55      6        0.6        V9 0.3693762
#> 56      6        0.6        V6 0.3506448
#> 57      6        0.6       V10 0.4149388
#> 58      6        0.6        V1 0.4017143
#> 59      6        0.6        V2 0.3968107
#> 60      6        0.6        V3 0.3775477
#> 61      7        0.7        V1 0.4010868
#> 62      7        0.7        V3 0.3770252
#> 63      7        0.7        V4 0.3824621
#> 64      7        0.7        V5 0.3839940
#> 65      7        0.7        V9 0.3688016
#> 66      7        0.7        V6 0.3501082
#> 67      7        0.7        V7 0.3547784
#> 68      7        0.7        V8 0.3667113
#> 69      7        0.7        V2 0.3962670
#> 70      7        0.7       V10 0.4142920
#> 71      8        0.8        V9 0.3682280
#> 72      8        0.8       V10 0.4136462
#> 73      8        0.8        V1 0.4004603
#> 74      8        0.8        V5 0.3835398
#> 75      8        0.8        V2 0.3957240
#> 76      8        0.8        V3 0.3765033
#> 77      8        0.8        V4 0.3818934
#> 78      8        0.8        V8 0.3662109
#> 79      8        0.8        V6 0.3495724
#> 80      8        0.8        V7 0.3541947
#> 81      9        0.9        V5 0.3830861
#> 82      9        0.9        V7 0.3536119
#> 83      9        0.9        V8 0.3657112
#> 84      9        0.9        V9 0.3676553
#> 85      9        0.9        V6 0.3490375
#> 86      9        0.9       V10 0.4130014
#> 87      9        0.9        V1 0.3998348
#> 88      9        0.9        V2 0.3951818
#> 89      9        0.9        V3 0.3759823
#> 90      9        0.9        V4 0.3813256
#> 91     10        1.0        V1 0.3992103
#> 92     10        1.0        V4 0.3807586
#> 93     10        1.0        V5 0.3826330
#> 94     10        1.0        V9 0.3670834
#> 95     10        1.0        V3 0.3754619
#> 96     10        1.0        V7 0.3530301
#> 97     10        1.0        V8 0.3652122
#> 98     10        1.0        V2 0.3946404
#> 99     10        1.0        V6 0.3485034
#> 100    10        1.0       V10 0.4123577

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4048663
#> 2       1        0.1        V5 0.3867305
#> 3       1        0.1        V9 0.3722622
#> 4       1        0.1        V3 0.3801713
#> 5       1        0.1        V4 0.3858921
#> 6       1        0.1        V8 0.3697280
#> 7       1        0.1        V2 0.3995403
#> 8       1        0.1        V6 0.3533400
#> 9       1        0.1        V7 0.3583010
#> 10      1        0.1       V10 0.4181881
#> 11      2        0.2        V7 0.3577115
#> 12      2        0.2        V8 0.3692235
#> 13      2        0.2        V9 0.3716832
#> 14      2        0.2       V10 0.4175362
#> 15      2        0.2        V1 0.4042340
#> 16      2        0.2        V5 0.3862731
#> 17      2        0.2        V2 0.3989929
#> 18      2        0.2        V3 0.3796451
#> 19      2        0.2        V4 0.3853183
#> 20      2        0.2        V6 0.3527993
#> 21      3        0.3        V4 0.3847454
#> 22      3        0.3        V5 0.3858162
#> 23      3        0.3        V3 0.3791197
#> 24      3        0.3        V7 0.3571230
#> 25      3        0.3        V8 0.3687197
#> 26      3        0.3        V9 0.3711051
#> 27      3        0.3        V6 0.3522595
#> 28      3        0.3       V10 0.4168853
#> 29      3        0.3        V1 0.4036026
#> 30      3        0.3        V2 0.3984462
#> 31      4        0.4        V1 0.4029722
#> 32      4        0.4        V9 0.3705279
#> 33      4        0.4        V3 0.3785949
#> 34      4        0.4        V4 0.3841733
#> 35      4        0.4        V5 0.3853598
#> 36      4        0.4        V2 0.3979003
#> 37      4        0.4        V6 0.3517204
#> 38      4        0.4        V7 0.3565354
#> 39      4        0.4        V8 0.3682166
#> 40      4        0.4       V10 0.4162355
#> 41      5        0.5        V8 0.3677141
#> 42      5        0.5        V9 0.3699516
#> 43      5        0.5       V10 0.4155867
#> 44      5        0.5        V1 0.4023427
#> 45      5        0.5        V5 0.3849040
#> 46      5        0.5        V2 0.3973551
#> 47      5        0.5        V3 0.3780710
#> 48      5        0.5        V4 0.3836021
#> 49      5        0.5        V6 0.3511822
#> 50      5        0.5        V7 0.3559488
#> 51      6        0.6        V4 0.3830317
#> 52      6        0.6        V5 0.3844487
#> 53      6        0.6        V7 0.3553631
#> 54      6        0.6        V8 0.3672124
#> 55      6        0.6        V9 0.3693762
#> 56      6        0.6        V6 0.3506448
#> 57      6        0.6       V10 0.4149388
#> 58      6        0.6        V1 0.4017143
#> 59      6        0.6        V2 0.3968107
#> 60      6        0.6        V3 0.3775477
#> 61      7        0.7        V1 0.4010868
#> 62      7        0.7        V3 0.3770252
#> 63      7        0.7        V4 0.3824621
#> 64      7        0.7        V5 0.3839940
#> 65      7        0.7        V9 0.3688016
#> 66      7        0.7        V6 0.3501082
#> 67      7        0.7        V7 0.3547784
#> 68      7        0.7        V8 0.3667113
#> 69      7        0.7        V2 0.3962670
#> 70      7        0.7       V10 0.4142920
#> 71      8        0.8        V9 0.3682280
#> 72      8        0.8       V10 0.4136462
#> 73      8        0.8        V1 0.4004603
#> 74      8        0.8        V5 0.3835398
#> 75      8        0.8        V2 0.3957240
#> 76      8        0.8        V3 0.3765033
#> 77      8        0.8        V4 0.3818934
#> 78      8        0.8        V8 0.3662109
#> 79      8        0.8        V6 0.3495724
#> 80      8        0.8        V7 0.3541947
#> 81      9        0.9        V5 0.3830861
#> 82      9        0.9        V7 0.3536119
#> 83      9        0.9        V8 0.3657112
#> 84      9        0.9        V9 0.3676553
#> 85      9        0.9        V6 0.3490375
#> 86      9        0.9       V10 0.4130014
#> 87      9        0.9        V1 0.3998348
#> 88      9        0.9        V2 0.3951818
#> 89      9        0.9        V3 0.3759823
#> 90      9        0.9        V4 0.3813256
#> 91     10        1.0        V1 0.3992103
#> 92     10        1.0        V4 0.3807586
#> 93     10        1.0        V5 0.3826330
#> 94     10        1.0        V9 0.3670834
#> 95     10        1.0        V3 0.3754619
#> 96     10        1.0        V7 0.3530301
#> 97     10        1.0        V8 0.3652122
#> 98     10        1.0        V2 0.3946404
#> 99     10        1.0        V6 0.3485034
#> 100    10        1.0       V10 0.4123577


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
#> 1       0                0          0 0.8757906 0.05064102 0.7689228 0.9399115
#> 2       0               20         20 0.8757906 0.05064102 0.7689228 0.9399115
#> 3       0               40         40 0.8757906 0.05064102 0.7689228 0.9399115
#> 4       0               60         60 0.8757906 0.05064102 0.7689228 0.9399115
#> 5      20                0         20 0.8617131 0.05315177 0.7492273 0.9306773
#> 6      20               20         40 0.8617131 0.05315177 0.7492273 0.9306773
#> 7      20               40         60 0.8617131 0.05315177 0.7492273 0.9306773
#> 8      20               60         80 0.8617131 0.05315177 0.7492273 0.9306773
#> 9      40                0         40 0.8478591 0.05550188 0.7302406 0.9213509
#> 10     40               20         60 0.8478591 0.05550188 0.7302406 0.9213509
#> 11     40               40         80 0.8478591 0.05550188 0.7302406 0.9213509
#> 12     40               60        100 0.8478591 0.05550188 0.7302406 0.9213509
#> 13     60                0         60 0.8342249 0.05770723 0.7118981 0.9119629
#> 14     60               20         80 0.8342249 0.05770723 0.7118981 0.9119629
#> 15     60               40        100 0.8342249 0.05770723 0.7118981 0.9119629
#> 16     60               60        120 0.8342249 0.05770723 0.7118981 0.9119629
#> 17     80                0         80 0.8208071 0.05978067 0.6941484 0.9025373
#> 18     80               20        100 0.8208071 0.05978067 0.6941484 0.9025373
#> 19     80               40        120 0.8208071 0.05978067 0.6941484 0.9025373
#> 20     80               60        140 0.8208071 0.05978067 0.6941484 0.9025373
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12205002 0.188927925 0.5849587
#> 2  0.30574618 0.11492462 0.135295118 0.5091004
#> 3  0.26001915 0.10765815 0.095421990 0.4436843
#> 4  0.22113100 0.10034255 0.065988819 0.3875144
#> 5  0.25589195 0.11294494 0.102343034 0.4761846
#> 6  0.21762106 0.10393773 0.071076523 0.4154021
#> 7  0.18507391 0.09545700 0.048199992 0.3632647
#> 8  0.15739448 0.08746543 0.031717900 0.3185700
#> 9  0.18213629 0.10115241 0.052135531 0.3891472
#> 10 0.15489621 0.09195563 0.034530583 0.3407593
#> 11 0.13173012 0.08351892 0.022055513 0.2992659
#> 12 0.11202872 0.07576315 0.013454673 0.2636253
#> 13 0.12963921 0.08883574 0.024166078 0.3198693
#> 14 0.11025053 0.08020940 0.014889132 0.2813328
#> 15 0.09376159 0.07238778 0.008676148 0.2481904
#> 16 0.07973872 0.06528163 0.004708694 0.2195882
#> 17 0.09227334 0.07706936 0.009696700 0.2646631
#> 18 0.07847305 0.06931880 0.005344161 0.2338186
#> 19 0.06673672 0.06233270 0.002700914 0.2071444
#> 20 0.05675565 0.05602309 0.001221946 0.1839594
```
