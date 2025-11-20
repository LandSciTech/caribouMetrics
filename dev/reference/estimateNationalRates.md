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
#> 1         0.1 0.3838513 0.02865579 0.3273277 0.4149186
#> 2         0.2 0.3832760 0.02862649 0.3268327 0.4143993
#> 3         0.3 0.3827015 0.02859731 0.3263384 0.4138807
#> 4         0.4 0.3821279 0.02856825 0.3258448 0.4133627
#> 5         0.5 0.3815551 0.02853932 0.3253520 0.4128454
#> 6         0.6 0.3809832 0.02851051 0.3248600 0.4123288
#> 7         0.7 0.3804122 0.02848183 0.3243686 0.4118128
#> 8         0.8 0.3798420 0.02845327 0.3238781 0.4112975
#> 9         0.9 0.3792726 0.02842483 0.3233882 0.4107828
#> 10        1.0 0.3787041 0.02839651 0.3228991 0.4102688

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4124190
#> 2       1        0.1        V5 0.3792603
#> 3       1        0.1        V9 0.3148093
#> 4       1        0.1        V3 0.3765767
#> 5       1        0.1        V4 0.3944593
#> 6       1        0.1        V8 0.3999652
#> 7       1        0.1        V2 0.4156443
#> 8       1        0.1        V6 0.3787813
#> 9       1        0.1        V7 0.3704467
#> 10      1        0.1       V10 0.3958979
#> 11      2        0.2        V7 0.3698709
#> 12      2        0.2        V8 0.3994005
#> 13      2        0.2        V9 0.3143377
#> 14      2        0.2       V10 0.3952522
#> 15      2        0.2        V1 0.4117671
#> 16      2        0.2        V5 0.3787000
#> 17      2        0.2        V2 0.4151635
#> 18      2        0.2        V3 0.3760525
#> 19      2        0.2        V4 0.3938955
#> 20      2        0.2        V6 0.3782540
#> 21      3        0.3        V4 0.3933324
#> 22      3        0.3        V5 0.3781405
#> 23      3        0.3        V3 0.3755290
#> 24      3        0.3        V7 0.3692960
#> 25      3        0.3        V8 0.3988367
#> 26      3        0.3        V9 0.3138668
#> 27      3        0.3        V6 0.3777274
#> 28      3        0.3       V10 0.3946076
#> 29      3        0.3        V1 0.4111162
#> 30      3        0.3        V2 0.4146833
#> 31      4        0.4        V1 0.4104663
#> 32      4        0.4        V9 0.3133966
#> 33      4        0.4        V3 0.3750063
#> 34      4        0.4        V4 0.3927702
#> 35      4        0.4        V5 0.3775818
#> 36      4        0.4        V2 0.4142036
#> 37      4        0.4        V6 0.3772016
#> 38      4        0.4        V7 0.3687220
#> 39      4        0.4        V8 0.3982737
#> 40      4        0.4       V10 0.3939640
#> 41      5        0.5        V8 0.3977115
#> 42      5        0.5        V9 0.3129271
#> 43      5        0.5       V10 0.3933215
#> 44      5        0.5        V1 0.4098174
#> 45      5        0.5        V5 0.3770240
#> 46      5        0.5        V2 0.4137245
#> 47      5        0.5        V3 0.3744843
#> 48      5        0.5        V4 0.3922088
#> 49      5        0.5        V6 0.3766765
#> 50      5        0.5        V7 0.3681489
#> 51      6        0.6        V4 0.3916481
#> 52      6        0.6        V5 0.3764670
#> 53      6        0.6        V7 0.3675767
#> 54      6        0.6        V8 0.3971500
#> 55      6        0.6        V9 0.3124583
#> 56      6        0.6        V6 0.3761521
#> 57      6        0.6       V10 0.3926800
#> 58      6        0.6        V1 0.4091695
#> 59      6        0.6        V2 0.4132460
#> 60      6        0.6        V3 0.3739630
#> 61      7        0.7        V1 0.4085227
#> 62      7        0.7        V3 0.3734424
#> 63      7        0.7        V4 0.3910883
#> 64      7        0.7        V5 0.3759108
#> 65      7        0.7        V9 0.3119902
#> 66      7        0.7        V6 0.3756285
#> 67      7        0.7        V7 0.3670053
#> 68      7        0.7        V8 0.3965894
#> 69      7        0.7        V2 0.4127680
#> 70      7        0.7       V10 0.3920396
#> 71      8        0.8        V9 0.3115229
#> 72      8        0.8       V10 0.3914002
#> 73      8        0.8        V1 0.4078769
#> 74      8        0.8        V5 0.3753555
#> 75      8        0.8        V2 0.4122905
#> 76      8        0.8        V3 0.3729225
#> 77      8        0.8        V4 0.3905293
#> 78      8        0.8        V8 0.3960295
#> 79      8        0.8        V6 0.3751056
#> 80      8        0.8        V7 0.3664349
#> 81      9        0.9        V5 0.3748009
#> 82      9        0.9        V7 0.3658653
#> 83      9        0.9        V8 0.3954704
#> 84      9        0.9        V9 0.3110562
#> 85      9        0.9        V6 0.3745834
#> 86      9        0.9       V10 0.3907618
#> 87      9        0.9        V1 0.4072321
#> 88      9        0.9        V2 0.4118136
#> 89      9        0.9        V3 0.3724034
#> 90      9        0.9        V4 0.3899711
#> 91     10        1.0        V1 0.4065884
#> 92     10        1.0        V4 0.3894136
#> 93     10        1.0        V5 0.3742472
#> 94     10        1.0        V9 0.3105902
#> 95     10        1.0        V3 0.3718850
#> 96     10        1.0        V7 0.3652966
#> 97     10        1.0        V8 0.3949122
#> 98     10        1.0        V2 0.4113373
#> 99     10        1.0        V6 0.3740619
#> 100    10        1.0       V10 0.3901245

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4124190
#> 2       1        0.1        V5 0.3792603
#> 3       1        0.1        V9 0.3148093
#> 4       1        0.1        V3 0.3765767
#> 5       1        0.1        V4 0.3944593
#> 6       1        0.1        V8 0.3999652
#> 7       1        0.1        V2 0.4156443
#> 8       1        0.1        V6 0.3787813
#> 9       1        0.1        V7 0.3704467
#> 10      1        0.1       V10 0.3958979
#> 11      2        0.2        V7 0.3698709
#> 12      2        0.2        V8 0.3994005
#> 13      2        0.2        V9 0.3143377
#> 14      2        0.2       V10 0.3952522
#> 15      2        0.2        V1 0.4117671
#> 16      2        0.2        V5 0.3787000
#> 17      2        0.2        V2 0.4151635
#> 18      2        0.2        V3 0.3760525
#> 19      2        0.2        V4 0.3938955
#> 20      2        0.2        V6 0.3782540
#> 21      3        0.3        V4 0.3933324
#> 22      3        0.3        V5 0.3781405
#> 23      3        0.3        V3 0.3755290
#> 24      3        0.3        V7 0.3692960
#> 25      3        0.3        V8 0.3988367
#> 26      3        0.3        V9 0.3138668
#> 27      3        0.3        V6 0.3777274
#> 28      3        0.3       V10 0.3946076
#> 29      3        0.3        V1 0.4111162
#> 30      3        0.3        V2 0.4146833
#> 31      4        0.4        V1 0.4104663
#> 32      4        0.4        V9 0.3133966
#> 33      4        0.4        V3 0.3750063
#> 34      4        0.4        V4 0.3927702
#> 35      4        0.4        V5 0.3775818
#> 36      4        0.4        V2 0.4142036
#> 37      4        0.4        V6 0.3772016
#> 38      4        0.4        V7 0.3687220
#> 39      4        0.4        V8 0.3982737
#> 40      4        0.4       V10 0.3939640
#> 41      5        0.5        V8 0.3977115
#> 42      5        0.5        V9 0.3129271
#> 43      5        0.5       V10 0.3933215
#> 44      5        0.5        V1 0.4098174
#> 45      5        0.5        V5 0.3770240
#> 46      5        0.5        V2 0.4137245
#> 47      5        0.5        V3 0.3744843
#> 48      5        0.5        V4 0.3922088
#> 49      5        0.5        V6 0.3766765
#> 50      5        0.5        V7 0.3681489
#> 51      6        0.6        V4 0.3916481
#> 52      6        0.6        V5 0.3764670
#> 53      6        0.6        V7 0.3675767
#> 54      6        0.6        V8 0.3971500
#> 55      6        0.6        V9 0.3124583
#> 56      6        0.6        V6 0.3761521
#> 57      6        0.6       V10 0.3926800
#> 58      6        0.6        V1 0.4091695
#> 59      6        0.6        V2 0.4132460
#> 60      6        0.6        V3 0.3739630
#> 61      7        0.7        V1 0.4085227
#> 62      7        0.7        V3 0.3734424
#> 63      7        0.7        V4 0.3910883
#> 64      7        0.7        V5 0.3759108
#> 65      7        0.7        V9 0.3119902
#> 66      7        0.7        V6 0.3756285
#> 67      7        0.7        V7 0.3670053
#> 68      7        0.7        V8 0.3965894
#> 69      7        0.7        V2 0.4127680
#> 70      7        0.7       V10 0.3920396
#> 71      8        0.8        V9 0.3115229
#> 72      8        0.8       V10 0.3914002
#> 73      8        0.8        V1 0.4078769
#> 74      8        0.8        V5 0.3753555
#> 75      8        0.8        V2 0.4122905
#> 76      8        0.8        V3 0.3729225
#> 77      8        0.8        V4 0.3905293
#> 78      8        0.8        V8 0.3960295
#> 79      8        0.8        V6 0.3751056
#> 80      8        0.8        V7 0.3664349
#> 81      9        0.9        V5 0.3748009
#> 82      9        0.9        V7 0.3658653
#> 83      9        0.9        V8 0.3954704
#> 84      9        0.9        V9 0.3110562
#> 85      9        0.9        V6 0.3745834
#> 86      9        0.9       V10 0.3907618
#> 87      9        0.9        V1 0.4072321
#> 88      9        0.9        V2 0.4118136
#> 89      9        0.9        V3 0.3724034
#> 90      9        0.9        V4 0.3899711
#> 91     10        1.0        V1 0.4065884
#> 92     10        1.0        V4 0.3894136
#> 93     10        1.0        V5 0.3742472
#> 94     10        1.0        V9 0.3105902
#> 95     10        1.0        V3 0.3718850
#> 96     10        1.0        V7 0.3652966
#> 97     10        1.0        V8 0.3949122
#> 98     10        1.0        V2 0.4113373
#> 99     10        1.0        V6 0.3740619
#> 100    10        1.0       V10 0.3901245


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
#> 1       0                0          0 0.8757906 0.05076932 0.7761190 0.9452819
#> 2       0               20         20 0.8757906 0.05076932 0.7761190 0.9452819
#> 3       0               40         40 0.8757906 0.05076932 0.7761190 0.9452819
#> 4       0               60         60 0.8757906 0.05076932 0.7761190 0.9452819
#> 5      20                0         20 0.8617131 0.05257199 0.7588838 0.9342184
#> 6      20               20         40 0.8617131 0.05257199 0.7588838 0.9342184
#> 7      20               40         60 0.8617131 0.05257199 0.7588838 0.9342184
#> 8      20               60         80 0.8617131 0.05257199 0.7588838 0.9342184
#> 9      40                0         40 0.8478591 0.05423998 0.7421960 0.9229941
#> 10     40               20         60 0.8478591 0.05423998 0.7421960 0.9229941
#> 11     40               40         80 0.8478591 0.05423998 0.7421960 0.9229941
#> 12     40               60        100 0.8478591 0.05423998 0.7421960 0.9229941
#> 13     60                0         60 0.8342249 0.05579035 0.7260097 0.9116676
#> 14     60               20         80 0.8342249 0.05579035 0.7260097 0.9116676
#> 15     60               40        100 0.8342249 0.05579035 0.7260097 0.9116676
#> 16     60               60        120 0.8342249 0.05579035 0.7260097 0.9116676
#> 17     80                0         80 0.8208071 0.05723646 0.7102878 0.9002831
#> 18     80               20        100 0.8208071 0.05723646 0.7102878 0.9002831
#> 19     80               40        120 0.8208071 0.05723646 0.7102878 0.9002831
#> 20     80               60        140 0.8208071 0.05723646 0.7102878 0.9002831
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12535531 0.1608206050 0.5576781
#> 2  0.30574618 0.12060506 0.1135370801 0.5246214
#> 3  0.26001915 0.11716035 0.0786825233 0.4936236
#> 4  0.22113100 0.11413939 0.0532453141 0.4645891
#> 5  0.25589195 0.11145002 0.0757929432 0.4385101
#> 6  0.21762106 0.10615249 0.0511528961 0.4130234
#> 7  0.18507391 0.10177769 0.0334729423 0.3891938
#> 8  0.15739448 0.09784766 0.0210638961 0.3669190
#> 9  0.18213629 0.09624076 0.0320367124 0.3469375
#> 10 0.15489621 0.09103839 0.0200718552 0.3274234
#> 11 0.13173012 0.08657066 0.0119531294 0.3091818
#> 12 0.11202872 0.08255572 0.0066730499 0.2921257
#> 13 0.12963921 0.08154859 0.0113197836 0.2768153
#> 14 0.11025053 0.07681725 0.0062736880 0.2618475
#> 15 0.09376159 0.07270148 0.0031957186 0.2478358
#> 16 0.07973872 0.06902001 0.0014618766 0.2347116
#> 17 0.09227334 0.06828261 0.0029736188 0.2229066
#> 18 0.07847305 0.06419861 0.0013439244 0.2113385
#> 19 0.06673672 0.06063068 0.0005281190 0.2004800
#> 20 0.05675565 0.05744699 0.0001739595 0.1902784
```
