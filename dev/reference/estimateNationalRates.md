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
#> 1         0.1 0.3838513 0.01892671 0.3351746 0.3939381
#> 2         0.2 0.3832760 0.01890947 0.3346689 0.3933925
#> 3         0.3 0.3827015 0.01889238 0.3341639 0.3928477
#> 4         0.4 0.3821279 0.01887542 0.3336596 0.3923036
#> 5         0.5 0.3815551 0.01885859 0.3331562 0.3917603
#> 6         0.6 0.3809832 0.01884190 0.3326535 0.3912178
#> 7         0.7 0.3804122 0.01882534 0.3321516 0.3906760
#> 8         0.8 0.3798420 0.01880892 0.3316504 0.3901349
#> 9         0.9 0.3792726 0.01879263 0.3311501 0.3895947
#> 10        1.0 0.3787041 0.01877647 0.3306505 0.3890551

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3731246
#> 2       1        0.1        V5 0.3276257
#> 3       1        0.1        V9 0.3899541
#> 4       1        0.1        V3 0.3660632
#> 5       1        0.1        V4 0.3655798
#> 6       1        0.1        V8 0.3611766
#> 7       1        0.1        V2 0.3950947
#> 8       1        0.1        V6 0.3800189
#> 9       1        0.1        V7 0.3830990
#> 10      1        0.1       V10 0.3664532
#> 11      2        0.2        V7 0.3825187
#> 12      2        0.2        V8 0.3607504
#> 13      2        0.2        V9 0.3893410
#> 14      2        0.2       V10 0.3658577
#> 15      2        0.2        V1 0.3725419
#> 16      2        0.2        V5 0.3270968
#> 17      2        0.2        V2 0.3945688
#> 18      2        0.2        V3 0.3654953
#> 19      2        0.2        V4 0.3650359
#> 20      2        0.2        V6 0.3794666
#> 21      3        0.3        V4 0.3644929
#> 22      3        0.3        V5 0.3265688
#> 23      3        0.3        V3 0.3649282
#> 24      3        0.3        V7 0.3819393
#> 25      3        0.3        V8 0.3603248
#> 26      3        0.3        V9 0.3887288
#> 27      3        0.3        V6 0.3789152
#> 28      3        0.3       V10 0.3652631
#> 29      3        0.3        V1 0.3719601
#> 30      3        0.3        V2 0.3940435
#> 31      4        0.4        V1 0.3713792
#> 32      4        0.4        V9 0.3881176
#> 33      4        0.4        V3 0.3643621
#> 34      4        0.4        V4 0.3639507
#> 35      4        0.4        V5 0.3260416
#> 36      4        0.4        V2 0.3935189
#> 37      4        0.4        V6 0.3783645
#> 38      4        0.4        V7 0.3813608
#> 39      4        0.4        V8 0.3598996
#> 40      4        0.4       V10 0.3646696
#> 41      5        0.5        V8 0.3594749
#> 42      5        0.5        V9 0.3875074
#> 43      5        0.5       V10 0.3640770
#> 44      5        0.5        V1 0.3707992
#> 45      5        0.5        V5 0.3255153
#> 46      5        0.5        V2 0.3929950
#> 47      5        0.5        V3 0.3637968
#> 48      5        0.5        V4 0.3634092
#> 49      5        0.5        V6 0.3778147
#> 50      5        0.5        V7 0.3807832
#> 51      6        0.6        V4 0.3628686
#> 52      6        0.6        V5 0.3249898
#> 53      6        0.6        V7 0.3802064
#> 54      6        0.6        V8 0.3590507
#> 55      6        0.6        V9 0.3868981
#> 56      6        0.6        V6 0.3772657
#> 57      6        0.6       V10 0.3634853
#> 58      6        0.6        V1 0.3702202
#> 59      6        0.6        V2 0.3924719
#> 60      6        0.6        V3 0.3632324
#> 61      7        0.7        V1 0.3696420
#> 62      7        0.7        V3 0.3626689
#> 63      7        0.7        V4 0.3623288
#> 64      7        0.7        V5 0.3244652
#> 65      7        0.7        V9 0.3862898
#> 66      7        0.7        V6 0.3767174
#> 67      7        0.7        V7 0.3796305
#> 68      7        0.7        V8 0.3586270
#> 69      7        0.7        V2 0.3919494
#> 70      7        0.7       V10 0.3628946
#> 71      8        0.8        V9 0.3856825
#> 72      8        0.8       V10 0.3623049
#> 73      8        0.8        V1 0.3690648
#> 74      8        0.8        V5 0.3239414
#> 75      8        0.8        V2 0.3914276
#> 76      8        0.8        V3 0.3621062
#> 77      8        0.8        V4 0.3617898
#> 78      8        0.8        V8 0.3582038
#> 79      8        0.8        V6 0.3761700
#> 80      8        0.8        V7 0.3790555
#> 81      9        0.9        V5 0.3234184
#> 82      9        0.9        V7 0.3784814
#> 83      9        0.9        V8 0.3577812
#> 84      9        0.9        V9 0.3850761
#> 85      9        0.9        V6 0.3756233
#> 86      9        0.9       V10 0.3617161
#> 87      9        0.9        V1 0.3684884
#> 88      9        0.9        V2 0.3909065
#> 89      9        0.9        V3 0.3615445
#> 90      9        0.9        V4 0.3612516
#> 91     10        1.0        V1 0.3679129
#> 92     10        1.0        V4 0.3607142
#> 93     10        1.0        V5 0.3228964
#> 94     10        1.0        V9 0.3844707
#> 95     10        1.0        V3 0.3609835
#> 96     10        1.0        V7 0.3779081
#> 97     10        1.0        V8 0.3573590
#> 98     10        1.0        V2 0.3903861
#> 99     10        1.0        V6 0.3750775
#> 100    10        1.0       V10 0.3611283

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3731246
#> 2       1        0.1        V5 0.3276257
#> 3       1        0.1        V9 0.3899541
#> 4       1        0.1        V3 0.3660632
#> 5       1        0.1        V4 0.3655798
#> 6       1        0.1        V8 0.3611766
#> 7       1        0.1        V2 0.3950947
#> 8       1        0.1        V6 0.3800189
#> 9       1        0.1        V7 0.3830990
#> 10      1        0.1       V10 0.3664532
#> 11      2        0.2        V7 0.3825187
#> 12      2        0.2        V8 0.3607504
#> 13      2        0.2        V9 0.3893410
#> 14      2        0.2       V10 0.3658577
#> 15      2        0.2        V1 0.3725419
#> 16      2        0.2        V5 0.3270968
#> 17      2        0.2        V2 0.3945688
#> 18      2        0.2        V3 0.3654953
#> 19      2        0.2        V4 0.3650359
#> 20      2        0.2        V6 0.3794666
#> 21      3        0.3        V4 0.3644929
#> 22      3        0.3        V5 0.3265688
#> 23      3        0.3        V3 0.3649282
#> 24      3        0.3        V7 0.3819393
#> 25      3        0.3        V8 0.3603248
#> 26      3        0.3        V9 0.3887288
#> 27      3        0.3        V6 0.3789152
#> 28      3        0.3       V10 0.3652631
#> 29      3        0.3        V1 0.3719601
#> 30      3        0.3        V2 0.3940435
#> 31      4        0.4        V1 0.3713792
#> 32      4        0.4        V9 0.3881176
#> 33      4        0.4        V3 0.3643621
#> 34      4        0.4        V4 0.3639507
#> 35      4        0.4        V5 0.3260416
#> 36      4        0.4        V2 0.3935189
#> 37      4        0.4        V6 0.3783645
#> 38      4        0.4        V7 0.3813608
#> 39      4        0.4        V8 0.3598996
#> 40      4        0.4       V10 0.3646696
#> 41      5        0.5        V8 0.3594749
#> 42      5        0.5        V9 0.3875074
#> 43      5        0.5       V10 0.3640770
#> 44      5        0.5        V1 0.3707992
#> 45      5        0.5        V5 0.3255153
#> 46      5        0.5        V2 0.3929950
#> 47      5        0.5        V3 0.3637968
#> 48      5        0.5        V4 0.3634092
#> 49      5        0.5        V6 0.3778147
#> 50      5        0.5        V7 0.3807832
#> 51      6        0.6        V4 0.3628686
#> 52      6        0.6        V5 0.3249898
#> 53      6        0.6        V7 0.3802064
#> 54      6        0.6        V8 0.3590507
#> 55      6        0.6        V9 0.3868981
#> 56      6        0.6        V6 0.3772657
#> 57      6        0.6       V10 0.3634853
#> 58      6        0.6        V1 0.3702202
#> 59      6        0.6        V2 0.3924719
#> 60      6        0.6        V3 0.3632324
#> 61      7        0.7        V1 0.3696420
#> 62      7        0.7        V3 0.3626689
#> 63      7        0.7        V4 0.3623288
#> 64      7        0.7        V5 0.3244652
#> 65      7        0.7        V9 0.3862898
#> 66      7        0.7        V6 0.3767174
#> 67      7        0.7        V7 0.3796305
#> 68      7        0.7        V8 0.3586270
#> 69      7        0.7        V2 0.3919494
#> 70      7        0.7       V10 0.3628946
#> 71      8        0.8        V9 0.3856825
#> 72      8        0.8       V10 0.3623049
#> 73      8        0.8        V1 0.3690648
#> 74      8        0.8        V5 0.3239414
#> 75      8        0.8        V2 0.3914276
#> 76      8        0.8        V3 0.3621062
#> 77      8        0.8        V4 0.3617898
#> 78      8        0.8        V8 0.3582038
#> 79      8        0.8        V6 0.3761700
#> 80      8        0.8        V7 0.3790555
#> 81      9        0.9        V5 0.3234184
#> 82      9        0.9        V7 0.3784814
#> 83      9        0.9        V8 0.3577812
#> 84      9        0.9        V9 0.3850761
#> 85      9        0.9        V6 0.3756233
#> 86      9        0.9       V10 0.3617161
#> 87      9        0.9        V1 0.3684884
#> 88      9        0.9        V2 0.3909065
#> 89      9        0.9        V3 0.3615445
#> 90      9        0.9        V4 0.3612516
#> 91     10        1.0        V1 0.3679129
#> 92     10        1.0        V4 0.3607142
#> 93     10        1.0        V5 0.3228964
#> 94     10        1.0        V9 0.3844707
#> 95     10        1.0        V3 0.3609835
#> 96     10        1.0        V7 0.3779081
#> 97     10        1.0        V8 0.3573590
#> 98     10        1.0        V2 0.3903861
#> 99     10        1.0        V6 0.3750775
#> 100    10        1.0       V10 0.3611283


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
#> 1       0                0          0 0.8757906 0.04799267 0.7821376 0.9533484
#> 2       0               20         20 0.8757906 0.04799267 0.7821376 0.9533484
#> 3       0               40         40 0.8757906 0.04799267 0.7821376 0.9533484
#> 4       0               60         60 0.8757906 0.04799267 0.7821376 0.9533484
#> 5      20                0         20 0.8617131 0.05002039 0.7657783 0.9443228
#> 6      20               20         40 0.8617131 0.05002039 0.7657783 0.9443228
#> 7      20               40         60 0.8617131 0.05002039 0.7657783 0.9443228
#> 8      20               60         80 0.8617131 0.05002039 0.7657783 0.9443228
#> 9      40                0         40 0.8478591 0.05196399 0.7499238 0.9351178
#> 10     40               20         60 0.8478591 0.05196399 0.7499238 0.9351178
#> 11     40               40         80 0.8478591 0.05196399 0.7499238 0.9351178
#> 12     40               60        100 0.8478591 0.05196399 0.7499238 0.9351178
#> 13     60                0         60 0.8342249 0.05382773 0.7345313 0.9257818
#> 14     60               20         80 0.8342249 0.05382773 0.7345313 0.9257818
#> 15     60               40        100 0.8342249 0.05382773 0.7345313 0.9257818
#> 16     60               60        120 0.8342249 0.05382773 0.7345313 0.9257818
#> 17     80                0         80 0.8208071 0.05561557 0.7195660 0.9163523
#> 18     80               20        100 0.8208071 0.05561557 0.7195660 0.9163523
#> 19     80               40        120 0.8208071 0.05561557 0.7195660 0.9163523
#> 20     80               60        140 0.8208071 0.05561557 0.7195660 0.9163523
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.13034652 0.1300953431 0.6095551
#> 2  0.30574618 0.12516558 0.0889072218 0.5534380
#> 3  0.26001915 0.11943656 0.0592488571 0.5026672
#> 4  0.22113100 0.11315410 0.0382256249 0.4568950
#> 5  0.25589195 0.11415421 0.0677707334 0.4856537
#> 6  0.21762106 0.10796197 0.0442231263 0.4415810
#> 7  0.18507391 0.10182325 0.0277748084 0.4019493
#> 8  0.15739448 0.09565553 0.0166137639 0.3663388
#> 9  0.18213629 0.09839646 0.0324365720 0.3887033
#> 10 0.15489621 0.09225965 0.0197366582 0.3544385
#> 11 0.13173012 0.08642295 0.0113376343 0.3236489
#> 12 0.11202872 0.08079479 0.0060511557 0.2959631
#> 13 0.12963921 0.08406114 0.0136610802 0.3133549
#> 14 0.11025053 0.07844962 0.0074816067 0.2867000
#> 15 0.09376159 0.07322044 0.0037566566 0.2626915
#> 16 0.07973872 0.06828878 0.0016884221 0.2410287
#> 17 0.09227334 0.07146654 0.0047451400 0.2546469
#> 18 0.07847305 0.06655501 0.0022170828 0.2337596
#> 19 0.06673672 0.06202157 0.0009092145 0.2148557
#> 20 0.05675565 0.05778715 0.0003161082 0.1977008
```
