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
[`simTrajectory()`](https://landscitech.github.io/caribouMetrics/dev/reference/simTrajectory.md),
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
#> 1         0.1 0.3838513 0.02177558 0.3574263 0.4262336
#> 2         0.2 0.3832760 0.02175528 0.3568879 0.4256019
#> 3         0.3 0.3827015 0.02173511 0.3563503 0.4249712
#> 4         0.4 0.3821279 0.02171509 0.3558135 0.4243414
#> 5         0.5 0.3815551 0.02169521 0.3552775 0.4237125
#> 6         0.6 0.3809832 0.02167546 0.3547423 0.4230846
#> 7         0.7 0.3804122 0.02165585 0.3542079 0.4224576
#> 8         0.8 0.3798420 0.02163638 0.3536743 0.4218315
#> 9         0.9 0.3792726 0.02161704 0.3531415 0.4212064
#> 10        1.0 0.3787041 0.02159784 0.3526095 0.4205822

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3895838
#> 2       1        0.1        V5 0.4037353
#> 3       1        0.1        V9 0.3821947
#> 4       1        0.1        V3 0.3982661
#> 5       1        0.1        V4 0.4043742
#> 6       1        0.1        V8 0.4072203
#> 7       1        0.1        V2 0.4015372
#> 8       1        0.1        V6 0.4237234
#> 9       1        0.1        V7 0.3502355
#> 10      1        0.1       V10 0.4269624
#> 11      2        0.2        V7 0.3497139
#> 12      2        0.2        V8 0.4066961
#> 13      2        0.2        V9 0.3815984
#> 14      2        0.2       V10 0.4263172
#> 15      2        0.2        V1 0.3889035
#> 16      2        0.2        V5 0.4031204
#> 17      2        0.2        V2 0.4009021
#> 18      2        0.2        V3 0.3977059
#> 19      2        0.2        V4 0.4037007
#> 20      2        0.2        V6 0.4231383
#> 21      3        0.3        V4 0.4030284
#> 22      3        0.3        V5 0.4025065
#> 23      3        0.3        V3 0.3971464
#> 24      3        0.3        V7 0.3491930
#> 25      3        0.3        V8 0.4061726
#> 26      3        0.3        V9 0.3810031
#> 27      3        0.3        V6 0.4225539
#> 28      3        0.3       V10 0.4256730
#> 29      3        0.3        V1 0.3882244
#> 30      3        0.3        V2 0.4002681
#> 31      4        0.4        V1 0.3875465
#> 32      4        0.4        V9 0.3804087
#> 33      4        0.4        V3 0.3965877
#> 34      4        0.4        V4 0.4023571
#> 35      4        0.4        V5 0.4018935
#> 36      4        0.4        V2 0.3996350
#> 37      4        0.4        V6 0.4219704
#> 38      4        0.4        V7 0.3486729
#> 39      4        0.4        V8 0.4056497
#> 40      4        0.4       V10 0.4250297
#> 41      5        0.5        V8 0.4051275
#> 42      5        0.5        V9 0.3798152
#> 43      5        0.5       V10 0.4243875
#> 44      5        0.5        V1 0.3868697
#> 45      5        0.5        V5 0.4012814
#> 46      5        0.5        V2 0.3990030
#> 47      5        0.5        V3 0.3960299
#> 48      5        0.5        V4 0.4016870
#> 49      5        0.5        V6 0.4213877
#> 50      5        0.5        V7 0.3481536
#> 51      6        0.6        V4 0.4010180
#> 52      6        0.6        V5 0.4006703
#> 53      6        0.6        V7 0.3476351
#> 54      6        0.6        V8 0.4046060
#> 55      6        0.6        V9 0.3792226
#> 56      6        0.6        V6 0.4208058
#> 57      6        0.6       V10 0.4237462
#> 58      6        0.6        V1 0.3861942
#> 59      6        0.6        V2 0.3983719
#> 60      6        0.6        V3 0.3954728
#> 61      7        0.7        V1 0.3855198
#> 62      7        0.7        V3 0.3949165
#> 63      7        0.7        V4 0.4003501
#> 64      7        0.7        V5 0.4000601
#> 65      7        0.7        V9 0.3786310
#> 66      7        0.7        V6 0.4202247
#> 67      7        0.7        V7 0.3471173
#> 68      7        0.7        V8 0.4040852
#> 69      7        0.7        V2 0.3977419
#> 70      7        0.7       V10 0.4231059
#> 71      8        0.8        V9 0.3780403
#> 72      8        0.8       V10 0.4224665
#> 73      8        0.8        V1 0.3848466
#> 74      8        0.8        V5 0.3994508
#> 75      8        0.8        V2 0.3971128
#> 76      8        0.8        V3 0.3943609
#> 77      8        0.8        V4 0.3996834
#> 78      8        0.8        V8 0.4035650
#> 79      8        0.8        V6 0.4196444
#> 80      8        0.8        V7 0.3466003
#> 81      9        0.9        V5 0.3988425
#> 82      9        0.9        V7 0.3460841
#> 83      9        0.9        V8 0.4030455
#> 84      9        0.9        V9 0.3774505
#> 85      9        0.9        V6 0.4190649
#> 86      9        0.9       V10 0.4218281
#> 87      9        0.9        V1 0.3841746
#> 88      9        0.9        V2 0.3964848
#> 89      9        0.9        V3 0.3938062
#> 90      9        0.9        V4 0.3990177
#> 91     10        1.0        V1 0.3835038
#> 92     10        1.0        V4 0.3983532
#> 93     10        1.0        V5 0.3982351
#> 94     10        1.0        V9 0.3768616
#> 95     10        1.0        V3 0.3932522
#> 96     10        1.0        V7 0.3455686
#> 97     10        1.0        V8 0.4025267
#> 98     10        1.0        V2 0.3958577
#> 99     10        1.0        V6 0.4184862
#> 100    10        1.0       V10 0.4211907

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3895838
#> 2       1        0.1        V5 0.4037353
#> 3       1        0.1        V9 0.3821947
#> 4       1        0.1        V3 0.3982661
#> 5       1        0.1        V4 0.4043742
#> 6       1        0.1        V8 0.4072203
#> 7       1        0.1        V2 0.4015372
#> 8       1        0.1        V6 0.4237234
#> 9       1        0.1        V7 0.3502355
#> 10      1        0.1       V10 0.4269624
#> 11      2        0.2        V7 0.3497139
#> 12      2        0.2        V8 0.4066961
#> 13      2        0.2        V9 0.3815984
#> 14      2        0.2       V10 0.4263172
#> 15      2        0.2        V1 0.3889035
#> 16      2        0.2        V5 0.4031204
#> 17      2        0.2        V2 0.4009021
#> 18      2        0.2        V3 0.3977059
#> 19      2        0.2        V4 0.4037007
#> 20      2        0.2        V6 0.4231383
#> 21      3        0.3        V4 0.4030284
#> 22      3        0.3        V5 0.4025065
#> 23      3        0.3        V3 0.3971464
#> 24      3        0.3        V7 0.3491930
#> 25      3        0.3        V8 0.4061726
#> 26      3        0.3        V9 0.3810031
#> 27      3        0.3        V6 0.4225539
#> 28      3        0.3       V10 0.4256730
#> 29      3        0.3        V1 0.3882244
#> 30      3        0.3        V2 0.4002681
#> 31      4        0.4        V1 0.3875465
#> 32      4        0.4        V9 0.3804087
#> 33      4        0.4        V3 0.3965877
#> 34      4        0.4        V4 0.4023571
#> 35      4        0.4        V5 0.4018935
#> 36      4        0.4        V2 0.3996350
#> 37      4        0.4        V6 0.4219704
#> 38      4        0.4        V7 0.3486729
#> 39      4        0.4        V8 0.4056497
#> 40      4        0.4       V10 0.4250297
#> 41      5        0.5        V8 0.4051275
#> 42      5        0.5        V9 0.3798152
#> 43      5        0.5       V10 0.4243875
#> 44      5        0.5        V1 0.3868697
#> 45      5        0.5        V5 0.4012814
#> 46      5        0.5        V2 0.3990030
#> 47      5        0.5        V3 0.3960299
#> 48      5        0.5        V4 0.4016870
#> 49      5        0.5        V6 0.4213877
#> 50      5        0.5        V7 0.3481536
#> 51      6        0.6        V4 0.4010180
#> 52      6        0.6        V5 0.4006703
#> 53      6        0.6        V7 0.3476351
#> 54      6        0.6        V8 0.4046060
#> 55      6        0.6        V9 0.3792226
#> 56      6        0.6        V6 0.4208058
#> 57      6        0.6       V10 0.4237462
#> 58      6        0.6        V1 0.3861942
#> 59      6        0.6        V2 0.3983719
#> 60      6        0.6        V3 0.3954728
#> 61      7        0.7        V1 0.3855198
#> 62      7        0.7        V3 0.3949165
#> 63      7        0.7        V4 0.4003501
#> 64      7        0.7        V5 0.4000601
#> 65      7        0.7        V9 0.3786310
#> 66      7        0.7        V6 0.4202247
#> 67      7        0.7        V7 0.3471173
#> 68      7        0.7        V8 0.4040852
#> 69      7        0.7        V2 0.3977419
#> 70      7        0.7       V10 0.4231059
#> 71      8        0.8        V9 0.3780403
#> 72      8        0.8       V10 0.4224665
#> 73      8        0.8        V1 0.3848466
#> 74      8        0.8        V5 0.3994508
#> 75      8        0.8        V2 0.3971128
#> 76      8        0.8        V3 0.3943609
#> 77      8        0.8        V4 0.3996834
#> 78      8        0.8        V8 0.4035650
#> 79      8        0.8        V6 0.4196444
#> 80      8        0.8        V7 0.3466003
#> 81      9        0.9        V5 0.3988425
#> 82      9        0.9        V7 0.3460841
#> 83      9        0.9        V8 0.4030455
#> 84      9        0.9        V9 0.3774505
#> 85      9        0.9        V6 0.4190649
#> 86      9        0.9       V10 0.4218281
#> 87      9        0.9        V1 0.3841746
#> 88      9        0.9        V2 0.3964848
#> 89      9        0.9        V3 0.3938062
#> 90      9        0.9        V4 0.3990177
#> 91     10        1.0        V1 0.3835038
#> 92     10        1.0        V4 0.3983532
#> 93     10        1.0        V5 0.3982351
#> 94     10        1.0        V9 0.3768616
#> 95     10        1.0        V3 0.3932522
#> 96     10        1.0        V7 0.3455686
#> 97     10        1.0        V8 0.4025267
#> 98     10        1.0        V2 0.3958577
#> 99     10        1.0        V6 0.4184862
#> 100    10        1.0       V10 0.4211907


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
#> 1       0                0          0 0.8757906 0.04917779 0.7857728 0.9474711
#> 2       0               20         20 0.8757906 0.04917779 0.7857728 0.9474711
#> 3       0               40         40 0.8757906 0.04917779 0.7857728 0.9474711
#> 4       0               60         60 0.8757906 0.04917779 0.7857728 0.9474711
#> 5      20                0         20 0.8617131 0.05131296 0.7680159 0.9367887
#> 6      20               20         40 0.8617131 0.05131296 0.7680159 0.9367887
#> 7      20               40         60 0.8617131 0.05131296 0.7680159 0.9367887
#> 8      20               60         80 0.8617131 0.05131296 0.7680159 0.9367887
#> 9      40                0         40 0.8478591 0.05328910 0.7508393 0.9259278
#> 10     40               20         60 0.8478591 0.05328910 0.7508393 0.9259278
#> 11     40               40         80 0.8478591 0.05328910 0.7508393 0.9259278
#> 12     40               60        100 0.8478591 0.05328910 0.7508393 0.9259278
#> 13     60                0         60 0.8342249 0.05512444 0.7341920 0.9149479
#> 14     60               20         80 0.8342249 0.05512444 0.7341920 0.9149479
#> 15     60               40        100 0.8342249 0.05512444 0.7341920 0.9149479
#> 16     60               60        120 0.8342249 0.05512444 0.7341920 0.9149479
#> 17     80                0         80 0.8208071 0.05683363 0.7180329 0.9038938
#> 18     80               20        100 0.8208071 0.05683363 0.7180329 0.9038938
#> 19     80               40        120 0.8208071 0.05683363 0.7180329 0.9038938
#> 20     80               60        140 0.8208071 0.05683363 0.7180329 0.9038938
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12334607 0.181250168 0.6036328
#> 2  0.30574618 0.11128994 0.128890524 0.5202398
#> 3  0.26001915 0.10060462 0.090156322 0.4489291
#> 4  0.22113100 0.09104631 0.061732458 0.3882799
#> 5  0.25589195 0.11382728 0.106621941 0.4883231
#> 6  0.21762106 0.10149727 0.073778529 0.4217579
#> 7  0.18507391 0.09075297 0.049828039 0.3652142
#> 8  0.15739448 0.08130461 0.032633859 0.3172382
#> 9  0.18213629 0.10183704 0.059955741 0.3964213
#> 10 0.15489621 0.09035481 0.039865190 0.3437159
#> 11 0.13173012 0.08038179 0.025598145 0.2989919
#> 12 0.11202872 0.07165411 0.015726110 0.2609894
#> 13 0.12963921 0.08935827 0.031575971 0.3236772
#> 14 0.11025053 0.07912334 0.019825611 0.2819744
#> 15 0.09376159 0.07022773 0.011837332 0.2465036
#> 16 0.07973872 0.06244866 0.006628827 0.2162359
#> 17 0.09227334 0.07740408 0.015134348 0.2660961
#> 18 0.07847305 0.06849280 0.008747418 0.2329687
#> 19 0.06673672 0.06073156 0.004699295 0.2046511
#> 20 0.05675565 0.05393658 0.002299116 0.1803256
```
