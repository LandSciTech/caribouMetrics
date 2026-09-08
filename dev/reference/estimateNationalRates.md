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
#> 1         0.1 0.3838513 0.02199526 0.3661382 0.4335835
#> 2         0.2 0.3832760 0.02195501 0.3656091 0.4329358
#> 3         0.3 0.3827015 0.02191496 0.3650808 0.4322891
#> 4         0.4 0.3821279 0.02187513 0.3645532 0.4316434
#> 5         0.5 0.3815551 0.02183551 0.3640264 0.4309986
#> 6         0.6 0.3809832 0.02179609 0.3635003 0.4303548
#> 7         0.7 0.3804122 0.02175689 0.3629750 0.4297119
#> 8         0.8 0.3798420 0.02171790 0.3624505 0.4290700
#> 9         0.9 0.3792726 0.02167912 0.3619267 0.4284291
#> 10        1.0 0.3787041 0.02164055 0.3614037 0.4277891

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3994205
#> 2       1        0.1        V5 0.3675099
#> 3       1        0.1        V9 0.3985736
#> 4       1        0.1        V3 0.3657400
#> 5       1        0.1        V4 0.4383279
#> 6       1        0.1        V8 0.3793531
#> 7       1        0.1        V2 0.4048633
#> 8       1        0.1        V6 0.3969467
#> 9       1        0.1        V7 0.4172419
#> 10      1        0.1       V10 0.3929843
#> 11      2        0.2        V7 0.4165390
#> 12      2        0.2        V8 0.3788072
#> 13      2        0.2        V9 0.3980215
#> 14      2        0.2       V10 0.3923501
#> 15      2        0.2        V1 0.3989226
#> 16      2        0.2        V5 0.3670065
#> 17      2        0.2        V2 0.4043305
#> 18      2        0.2        V3 0.3652034
#> 19      2        0.2        V4 0.4376962
#> 20      2        0.2        V6 0.3963126
#> 21      3        0.3        V4 0.4370655
#> 22      3        0.3        V5 0.3665039
#> 23      3        0.3        V3 0.3646676
#> 24      3        0.3        V7 0.4158372
#> 25      3        0.3        V8 0.3782621
#> 26      3        0.3        V9 0.3974702
#> 27      3        0.3        V6 0.3956795
#> 28      3        0.3       V10 0.3917170
#> 29      3        0.3        V1 0.3984253
#> 30      3        0.3        V2 0.4037985
#> 31      4        0.4        V1 0.3979287
#> 32      4        0.4        V9 0.3969197
#> 33      4        0.4        V3 0.3641326
#> 34      4        0.4        V4 0.4364356
#> 35      4        0.4        V5 0.3660019
#> 36      4        0.4        V2 0.4032671
#> 37      4        0.4        V6 0.3950475
#> 38      4        0.4        V7 0.4151366
#> 39      4        0.4        V8 0.3777178
#> 40      4        0.4       V10 0.3910848
#> 41      5        0.5        V8 0.3771743
#> 42      5        0.5        V9 0.3963699
#> 43      5        0.5       V10 0.3904537
#> 44      5        0.5        V1 0.3974326
#> 45      5        0.5        V5 0.3655007
#> 46      5        0.5        V2 0.4027364
#> 47      5        0.5        V3 0.3635983
#> 48      5        0.5        V4 0.4358067
#> 49      5        0.5        V6 0.3944165
#> 50      5        0.5        V7 0.4144372
#> 51      6        0.6        V4 0.4351787
#> 52      6        0.6        V5 0.3650001
#> 53      6        0.6        V7 0.4137390
#> 54      6        0.6        V8 0.3766315
#> 55      6        0.6        V9 0.3958208
#> 56      6        0.6        V6 0.3937865
#> 57      6        0.6       V10 0.3898236
#> 58      6        0.6        V1 0.3969372
#> 59      6        0.6        V2 0.4022065
#> 60      6        0.6        V3 0.3630649
#> 61      7        0.7        V1 0.3964423
#> 62      7        0.7        V3 0.3625322
#> 63      7        0.7        V4 0.4345516
#> 64      7        0.7        V5 0.3645002
#> 65      7        0.7        V9 0.3952726
#> 66      7        0.7        V6 0.3931574
#> 67      7        0.7        V7 0.4130419
#> 68      7        0.7        V8 0.3760895
#> 69      7        0.7        V2 0.4016772
#> 70      7        0.7       V10 0.3891946
#> 71      8        0.8        V9 0.3947251
#> 72      8        0.8       V10 0.3885665
#> 73      8        0.8        V1 0.3959481
#> 74      8        0.8        V5 0.3640010
#> 75      8        0.8        V2 0.4011486
#> 76      8        0.8        V3 0.3620003
#> 77      8        0.8        V4 0.4339254
#> 78      8        0.8        V8 0.3755483
#> 79      8        0.8        V6 0.3925294
#> 80      8        0.8        V7 0.4123461
#> 81      9        0.9        V5 0.3635025
#> 82      9        0.9        V7 0.4116514
#> 83      9        0.9        V8 0.3750079
#> 84      9        0.9        V9 0.3941783
#> 85      9        0.9        V6 0.3919024
#> 86      9        0.9       V10 0.3879395
#> 87      9        0.9        V1 0.3954546
#> 88      9        0.9        V2 0.4006207
#> 89      9        0.9        V3 0.3614692
#> 90      9        0.9        V4 0.4333000
#> 91     10        1.0        V1 0.3949616
#> 92     10        1.0        V4 0.4326756
#> 93     10        1.0        V5 0.3630046
#> 94     10        1.0        V9 0.3936323
#> 95     10        1.0        V3 0.3609389
#> 96     10        1.0        V7 0.4109578
#> 97     10        1.0        V8 0.3744683
#> 98     10        1.0        V2 0.4000935
#> 99     10        1.0        V6 0.3912764
#> 100    10        1.0       V10 0.3873134

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3994205
#> 2       1        0.1        V5 0.3675099
#> 3       1        0.1        V9 0.3985736
#> 4       1        0.1        V3 0.3657400
#> 5       1        0.1        V4 0.4383279
#> 6       1        0.1        V8 0.3793531
#> 7       1        0.1        V2 0.4048633
#> 8       1        0.1        V6 0.3969467
#> 9       1        0.1        V7 0.4172419
#> 10      1        0.1       V10 0.3929843
#> 11      2        0.2        V7 0.4165390
#> 12      2        0.2        V8 0.3788072
#> 13      2        0.2        V9 0.3980215
#> 14      2        0.2       V10 0.3923501
#> 15      2        0.2        V1 0.3989226
#> 16      2        0.2        V5 0.3670065
#> 17      2        0.2        V2 0.4043305
#> 18      2        0.2        V3 0.3652034
#> 19      2        0.2        V4 0.4376962
#> 20      2        0.2        V6 0.3963126
#> 21      3        0.3        V4 0.4370655
#> 22      3        0.3        V5 0.3665039
#> 23      3        0.3        V3 0.3646676
#> 24      3        0.3        V7 0.4158372
#> 25      3        0.3        V8 0.3782621
#> 26      3        0.3        V9 0.3974702
#> 27      3        0.3        V6 0.3956795
#> 28      3        0.3       V10 0.3917170
#> 29      3        0.3        V1 0.3984253
#> 30      3        0.3        V2 0.4037985
#> 31      4        0.4        V1 0.3979287
#> 32      4        0.4        V9 0.3969197
#> 33      4        0.4        V3 0.3641326
#> 34      4        0.4        V4 0.4364356
#> 35      4        0.4        V5 0.3660019
#> 36      4        0.4        V2 0.4032671
#> 37      4        0.4        V6 0.3950475
#> 38      4        0.4        V7 0.4151366
#> 39      4        0.4        V8 0.3777178
#> 40      4        0.4       V10 0.3910848
#> 41      5        0.5        V8 0.3771743
#> 42      5        0.5        V9 0.3963699
#> 43      5        0.5       V10 0.3904537
#> 44      5        0.5        V1 0.3974326
#> 45      5        0.5        V5 0.3655007
#> 46      5        0.5        V2 0.4027364
#> 47      5        0.5        V3 0.3635983
#> 48      5        0.5        V4 0.4358067
#> 49      5        0.5        V6 0.3944165
#> 50      5        0.5        V7 0.4144372
#> 51      6        0.6        V4 0.4351787
#> 52      6        0.6        V5 0.3650001
#> 53      6        0.6        V7 0.4137390
#> 54      6        0.6        V8 0.3766315
#> 55      6        0.6        V9 0.3958208
#> 56      6        0.6        V6 0.3937865
#> 57      6        0.6       V10 0.3898236
#> 58      6        0.6        V1 0.3969372
#> 59      6        0.6        V2 0.4022065
#> 60      6        0.6        V3 0.3630649
#> 61      7        0.7        V1 0.3964423
#> 62      7        0.7        V3 0.3625322
#> 63      7        0.7        V4 0.4345516
#> 64      7        0.7        V5 0.3645002
#> 65      7        0.7        V9 0.3952726
#> 66      7        0.7        V6 0.3931574
#> 67      7        0.7        V7 0.4130419
#> 68      7        0.7        V8 0.3760895
#> 69      7        0.7        V2 0.4016772
#> 70      7        0.7       V10 0.3891946
#> 71      8        0.8        V9 0.3947251
#> 72      8        0.8       V10 0.3885665
#> 73      8        0.8        V1 0.3959481
#> 74      8        0.8        V5 0.3640010
#> 75      8        0.8        V2 0.4011486
#> 76      8        0.8        V3 0.3620003
#> 77      8        0.8        V4 0.4339254
#> 78      8        0.8        V8 0.3755483
#> 79      8        0.8        V6 0.3925294
#> 80      8        0.8        V7 0.4123461
#> 81      9        0.9        V5 0.3635025
#> 82      9        0.9        V7 0.4116514
#> 83      9        0.9        V8 0.3750079
#> 84      9        0.9        V9 0.3941783
#> 85      9        0.9        V6 0.3919024
#> 86      9        0.9       V10 0.3879395
#> 87      9        0.9        V1 0.3954546
#> 88      9        0.9        V2 0.4006207
#> 89      9        0.9        V3 0.3614692
#> 90      9        0.9        V4 0.4333000
#> 91     10        1.0        V1 0.3949616
#> 92     10        1.0        V4 0.4326756
#> 93     10        1.0        V5 0.3630046
#> 94     10        1.0        V9 0.3936323
#> 95     10        1.0        V3 0.3609389
#> 96     10        1.0        V7 0.4109578
#> 97     10        1.0        V8 0.3744683
#> 98     10        1.0        V2 0.4000935
#> 99     10        1.0        V6 0.3912764
#> 100    10        1.0       V10 0.3873134


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
#> 1       0                0          0 0.8757906 0.04980189 0.7709334 0.9373458
#> 2       0               20         20 0.8757906 0.04980189 0.7709334 0.9373458
#> 3       0               40         40 0.8757906 0.04980189 0.7709334 0.9373458
#> 4       0               60         60 0.8757906 0.04980189 0.7709334 0.9373458
#> 5      20                0         20 0.8617131 0.05233807 0.7520741 0.9272252
#> 6      20               20         40 0.8617131 0.05233807 0.7520741 0.9272252
#> 7      20               40         60 0.8617131 0.05233807 0.7520741 0.9272252
#> 8      20               60         80 0.8617131 0.05233807 0.7520741 0.9272252
#> 9      40                0         40 0.8478591 0.05469413 0.7338591 0.9169983
#> 10     40               20         60 0.8478591 0.05469413 0.7338591 0.9169983
#> 11     40               40         80 0.8478591 0.05469413 0.7338591 0.9169983
#> 12     40               60        100 0.8478591 0.05469413 0.7338591 0.9169983
#> 13     60                0         60 0.8342249 0.05688901 0.7162327 0.9067027
#> 14     60               20         80 0.8342249 0.05688901 0.7162327 0.9067027
#> 15     60               40        100 0.8342249 0.05688901 0.7162327 0.9067027
#> 16     60               60        120 0.8342249 0.05688901 0.7162327 0.9067027
#> 17     80                0         80 0.8208071 0.05893821 0.6991497 0.8963680
#> 18     80               20        100 0.8208071 0.05893821 0.6991497 0.8963680
#> 19     80               40        120 0.8208071 0.05893821 0.6991497 0.8963680
#> 20     80               60        140 0.8208071 0.05893821 0.6991497 0.8963680
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11386373 0.182807374 0.5665299
#> 2  0.30574618 0.11016334 0.141095025 0.5105048
#> 3  0.26001915 0.10588789 0.107923118 0.4603492
#> 4  0.22113100 0.10117567 0.081641637 0.4155658
#> 5  0.25589195 0.09921070 0.101406944 0.4452134
#> 6  0.21762106 0.09480779 0.076496299 0.4020661
#> 7  0.18507391 0.09019833 0.056906647 0.3636056
#> 8  0.15739448 0.08545943 0.041632891 0.3293341
#> 9  0.18213629 0.08527657 0.053094350 0.3520181
#> 10 0.15489621 0.08082039 0.038680152 0.3190076
#> 11 0.13173012 0.07634822 0.027599180 0.2895778
#> 12 0.11202872 0.07190356 0.019208309 0.2633148
#> 13 0.12963921 0.07269147 0.025479951 0.2807039
#> 14 0.11025053 0.06848591 0.017621871 0.2553887
#> 15 0.09376159 0.06436184 0.011813529 0.2327527
#> 16 0.07973872 0.06034220 0.007628529 0.2124722
#> 17 0.09227334 0.06161927 0.010735364 0.2259102
#> 18 0.07847305 0.05779554 0.006866686 0.2063321
#> 19 0.06673672 0.05409792 0.004190258 0.1887338
#> 20 0.05675565 0.05053581 0.002416085 0.1728688
```
