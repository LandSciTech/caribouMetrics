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
#> 1         0.1 0.3838513 0.01664075 0.3628787 0.4071843
#> 2         0.2 0.3832760 0.01660621 0.3623699 0.4065094
#> 3         0.3 0.3827015 0.01657191 0.3618618 0.4058356
#> 4         0.4 0.3821279 0.01653785 0.3613544 0.4051629
#> 5         0.5 0.3815551 0.01650403 0.3608477 0.4044914
#> 6         0.6 0.3809832 0.01647044 0.3603418 0.4038209
#> 7         0.7 0.3804122 0.01643710 0.3598365 0.4031516
#> 8         0.8 0.3798420 0.01640399 0.3593319 0.4024834
#> 9         0.9 0.3792726 0.01637112 0.3588281 0.4018163
#> 10        1.0 0.3787041 0.01633848 0.3583250 0.4011504

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4080658
#> 2       1        0.1        V5 0.3692207
#> 3       1        0.1        V9 0.4028840
#> 4       1        0.1        V3 0.3755903
#> 5       1        0.1        V4 0.4018422
#> 6       1        0.1        V8 0.3767610
#> 7       1        0.1        V2 0.3855508
#> 8       1        0.1        V6 0.3610375
#> 9       1        0.1        V7 0.4041482
#> 10      1        0.1       V10 0.3794140
#> 11      2        0.2        V7 0.4035798
#> 12      2        0.2        V8 0.3761761
#> 13      2        0.2        V9 0.4022196
#> 14      2        0.2       V10 0.3788269
#> 15      2        0.2        V1 0.4073599
#> 16      2        0.2        V5 0.3686808
#> 17      2        0.2        V2 0.3849819
#> 18      2        0.2        V3 0.3749438
#> 19      2        0.2        V4 0.4013035
#> 20      2        0.2        V6 0.3605377
#> 21      3        0.3        V4 0.4007655
#> 22      3        0.3        V5 0.3681417
#> 23      3        0.3        V3 0.3742984
#> 24      3        0.3        V7 0.4030122
#> 25      3        0.3        V8 0.3755920
#> 26      3        0.3        V9 0.4015563
#> 27      3        0.3        V6 0.3600386
#> 28      3        0.3       V10 0.3782408
#> 29      3        0.3        V1 0.4066553
#> 30      3        0.3        V2 0.3844138
#> 31      4        0.4        V1 0.4059519
#> 32      4        0.4        V9 0.4008941
#> 33      4        0.4        V3 0.3736542
#> 34      4        0.4        V4 0.4002282
#> 35      4        0.4        V5 0.3676034
#> 36      4        0.4        V2 0.3838466
#> 37      4        0.4        V6 0.3595402
#> 38      4        0.4        V7 0.4024454
#> 39      4        0.4        V8 0.3750089
#> 40      4        0.4       V10 0.3776556
#> 41      5        0.5        V8 0.3744266
#> 42      5        0.5        V9 0.4002330
#> 43      5        0.5       V10 0.3770713
#> 44      5        0.5        V1 0.4052497
#> 45      5        0.5        V5 0.3670659
#> 46      5        0.5        V2 0.3832802
#> 47      5        0.5        V3 0.3730110
#> 48      5        0.5        V4 0.3996917
#> 49      5        0.5        V6 0.3590424
#> 50      5        0.5        V7 0.4018794
#> 51      6        0.6        V4 0.3991559
#> 52      6        0.6        V5 0.3665292
#> 53      6        0.6        V7 0.4013141
#> 54      6        0.6        V8 0.3738453
#> 55      6        0.6        V9 0.3995730
#> 56      6        0.6        V6 0.3585454
#> 57      6        0.6       V10 0.3764879
#> 58      6        0.6        V1 0.4045487
#> 59      6        0.6        V2 0.3827146
#> 60      6        0.6        V3 0.3723689
#> 61      7        0.7        V1 0.4038489
#> 62      7        0.7        V3 0.3717280
#> 63      7        0.7        V4 0.3986208
#> 64      7        0.7        V5 0.3659933
#> 65      7        0.7        V9 0.3989141
#> 66      7        0.7        V6 0.3580490
#> 67      7        0.7        V7 0.4007497
#> 68      7        0.7        V8 0.3732649
#> 69      7        0.7        V2 0.3821498
#> 70      7        0.7       V10 0.3759054
#> 71      8        0.8        V9 0.3982563
#> 72      8        0.8       V10 0.3753238
#> 73      8        0.8        V1 0.4031504
#> 74      8        0.8        V5 0.3654581
#> 75      8        0.8        V2 0.3815859
#> 76      8        0.8        V3 0.3710881
#> 77      8        0.8        V4 0.3980864
#> 78      8        0.8        V8 0.3726853
#> 79      8        0.8        V6 0.3575534
#> 80      8        0.8        V7 0.4001861
#> 81      9        0.9        V5 0.3649238
#> 82      9        0.9        V7 0.3996232
#> 83      9        0.9        V8 0.3721067
#> 84      9        0.9        V9 0.3975995
#> 85      9        0.9        V6 0.3570584
#> 86      9        0.9       V10 0.3747431
#> 87      9        0.9        V1 0.4024530
#> 88      9        0.9        V2 0.3810229
#> 89      9        0.9        V3 0.3704494
#> 90      9        0.9        V4 0.3975527
#> 91     10        1.0        V1 0.4017569
#> 92     10        1.0        V4 0.3970197
#> 93     10        1.0        V5 0.3643902
#> 94     10        1.0        V9 0.3969439
#> 95     10        1.0        V3 0.3698117
#> 96     10        1.0        V7 0.3990612
#> 97     10        1.0        V8 0.3715290
#> 98     10        1.0        V2 0.3804606
#> 99     10        1.0        V6 0.3565641
#> 100    10        1.0       V10 0.3741633

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4080658
#> 2       1        0.1        V5 0.3692207
#> 3       1        0.1        V9 0.4028840
#> 4       1        0.1        V3 0.3755903
#> 5       1        0.1        V4 0.4018422
#> 6       1        0.1        V8 0.3767610
#> 7       1        0.1        V2 0.3855508
#> 8       1        0.1        V6 0.3610375
#> 9       1        0.1        V7 0.4041482
#> 10      1        0.1       V10 0.3794140
#> 11      2        0.2        V7 0.4035798
#> 12      2        0.2        V8 0.3761761
#> 13      2        0.2        V9 0.4022196
#> 14      2        0.2       V10 0.3788269
#> 15      2        0.2        V1 0.4073599
#> 16      2        0.2        V5 0.3686808
#> 17      2        0.2        V2 0.3849819
#> 18      2        0.2        V3 0.3749438
#> 19      2        0.2        V4 0.4013035
#> 20      2        0.2        V6 0.3605377
#> 21      3        0.3        V4 0.4007655
#> 22      3        0.3        V5 0.3681417
#> 23      3        0.3        V3 0.3742984
#> 24      3        0.3        V7 0.4030122
#> 25      3        0.3        V8 0.3755920
#> 26      3        0.3        V9 0.4015563
#> 27      3        0.3        V6 0.3600386
#> 28      3        0.3       V10 0.3782408
#> 29      3        0.3        V1 0.4066553
#> 30      3        0.3        V2 0.3844138
#> 31      4        0.4        V1 0.4059519
#> 32      4        0.4        V9 0.4008941
#> 33      4        0.4        V3 0.3736542
#> 34      4        0.4        V4 0.4002282
#> 35      4        0.4        V5 0.3676034
#> 36      4        0.4        V2 0.3838466
#> 37      4        0.4        V6 0.3595402
#> 38      4        0.4        V7 0.4024454
#> 39      4        0.4        V8 0.3750089
#> 40      4        0.4       V10 0.3776556
#> 41      5        0.5        V8 0.3744266
#> 42      5        0.5        V9 0.4002330
#> 43      5        0.5       V10 0.3770713
#> 44      5        0.5        V1 0.4052497
#> 45      5        0.5        V5 0.3670659
#> 46      5        0.5        V2 0.3832802
#> 47      5        0.5        V3 0.3730110
#> 48      5        0.5        V4 0.3996917
#> 49      5        0.5        V6 0.3590424
#> 50      5        0.5        V7 0.4018794
#> 51      6        0.6        V4 0.3991559
#> 52      6        0.6        V5 0.3665292
#> 53      6        0.6        V7 0.4013141
#> 54      6        0.6        V8 0.3738453
#> 55      6        0.6        V9 0.3995730
#> 56      6        0.6        V6 0.3585454
#> 57      6        0.6       V10 0.3764879
#> 58      6        0.6        V1 0.4045487
#> 59      6        0.6        V2 0.3827146
#> 60      6        0.6        V3 0.3723689
#> 61      7        0.7        V1 0.4038489
#> 62      7        0.7        V3 0.3717280
#> 63      7        0.7        V4 0.3986208
#> 64      7        0.7        V5 0.3659933
#> 65      7        0.7        V9 0.3989141
#> 66      7        0.7        V6 0.3580490
#> 67      7        0.7        V7 0.4007497
#> 68      7        0.7        V8 0.3732649
#> 69      7        0.7        V2 0.3821498
#> 70      7        0.7       V10 0.3759054
#> 71      8        0.8        V9 0.3982563
#> 72      8        0.8       V10 0.3753238
#> 73      8        0.8        V1 0.4031504
#> 74      8        0.8        V5 0.3654581
#> 75      8        0.8        V2 0.3815859
#> 76      8        0.8        V3 0.3710881
#> 77      8        0.8        V4 0.3980864
#> 78      8        0.8        V8 0.3726853
#> 79      8        0.8        V6 0.3575534
#> 80      8        0.8        V7 0.4001861
#> 81      9        0.9        V5 0.3649238
#> 82      9        0.9        V7 0.3996232
#> 83      9        0.9        V8 0.3721067
#> 84      9        0.9        V9 0.3975995
#> 85      9        0.9        V6 0.3570584
#> 86      9        0.9       V10 0.3747431
#> 87      9        0.9        V1 0.4024530
#> 88      9        0.9        V2 0.3810229
#> 89      9        0.9        V3 0.3704494
#> 90      9        0.9        V4 0.3975527
#> 91     10        1.0        V1 0.4017569
#> 92     10        1.0        V4 0.3970197
#> 93     10        1.0        V5 0.3643902
#> 94     10        1.0        V9 0.3969439
#> 95     10        1.0        V3 0.3698117
#> 96     10        1.0        V7 0.3990612
#> 97     10        1.0        V8 0.3715290
#> 98     10        1.0        V2 0.3804606
#> 99     10        1.0        V6 0.3565641
#> 100    10        1.0       V10 0.3741633


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
#> 1       0                0          0 0.8757906 0.04951606 0.7904905 0.9499376
#> 2       0               20         20 0.8757906 0.04951606 0.7904905 0.9499376
#> 3       0               40         40 0.8757906 0.04951606 0.7904905 0.9499376
#> 4       0               60         60 0.8757906 0.04951606 0.7904905 0.9499376
#> 5      20                0         20 0.8617131 0.04887569 0.7797642 0.9352238
#> 6      20               20         40 0.8617131 0.04887569 0.7797642 0.9352238
#> 7      20               40         60 0.8617131 0.04887569 0.7797642 0.9352238
#> 8      20               60         80 0.8617131 0.04887569 0.7797642 0.9352238
#> 9      40                0         40 0.8478591 0.04814877 0.7692546 0.9202144
#> 10     40               20         60 0.8478591 0.04814877 0.7692546 0.9202144
#> 11     40               40         80 0.8478591 0.04814877 0.7692546 0.9202144
#> 12     40               60        100 0.8478591 0.04814877 0.7692546 0.9202144
#> 13     60                0         60 0.8342249 0.04736728 0.7589488 0.9050466
#> 14     60               20         80 0.8342249 0.04736728 0.7589488 0.9050466
#> 15     60               40        100 0.8342249 0.04736728 0.7589488 0.9050466
#> 16     60               60        120 0.8342249 0.04736728 0.7589488 0.9050466
#> 17     80                0         80 0.8208071 0.04655452 0.7488359 0.8898149
#> 18     80               20        100 0.8208071 0.04655452 0.7488359 0.8898149
#> 19     80               40        120 0.8208071 0.04655452 0.7488359 0.8898149
#> 20     80               60        140 0.8208071 0.04655452 0.7488359 0.8898149
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11350512 0.211870895 0.5311607
#> 2  0.30574618 0.12090798 0.153507282 0.4964568
#> 3  0.26001915 0.12513883 0.109800379 0.4641663
#> 4  0.22113100 0.12682359 0.077238972 0.4341520
#> 5  0.25589195 0.10568517 0.121828097 0.4354060
#> 6  0.21762106 0.10823937 0.086175124 0.4074378
#> 7  0.18507391 0.10898887 0.059768146 0.3814716
#> 8  0.15739448 0.10829198 0.040438503 0.3573707
#> 9  0.18213629 0.09531770 0.066995476 0.3583771
#> 10 0.15489621 0.09513002 0.045700737 0.3359377
#> 11 0.13173012 0.09396260 0.030282004 0.3151120
#> 12 0.11202872 0.09201641 0.019344321 0.2957809
#> 13 0.12963921 0.08414719 0.034459224 0.2965882
#> 14 0.11025053 0.08252567 0.022280649 0.2785819
#> 15 0.09376159 0.08040920 0.013798931 0.2618572
#> 16 0.07973872 0.07790975 0.008091128 0.2463153
#> 17 0.09227334 0.07320215 0.016057302 0.2469648
#> 18 0.07847305 0.07092634 0.009587620 0.2324684
#> 19 0.06673672 0.06843716 0.005365718 0.2189807
#> 20 0.05675565 0.06579715 0.002764325 0.2064219
```
