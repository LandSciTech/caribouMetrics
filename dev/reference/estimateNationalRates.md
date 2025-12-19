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
#> 1         0.1 0.3838513 0.02336296 0.3419085 0.4067104
#> 2         0.2 0.3832760 0.02333747 0.3413747 0.4061283
#> 3         0.3 0.3827015 0.02331212 0.3408417 0.4055470
#> 4         0.4 0.3821279 0.02328689 0.3403095 0.4049666
#> 5         0.5 0.3815551 0.02326179 0.3397782 0.4043869
#> 6         0.6 0.3809832 0.02323682 0.3392477 0.4038082
#> 7         0.7 0.3804122 0.02321198 0.3387181 0.4032302
#> 8         0.8 0.3798420 0.02318726 0.3381892 0.4026531
#> 9         0.9 0.3792726 0.02316267 0.3376612 0.4020768
#> 10        1.0 0.3787041 0.02313821 0.3371340 0.4015013

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4059331
#> 2       1        0.1        V5 0.3672968
#> 3       1        0.1        V9 0.3606096
#> 4       1        0.1        V3 0.3396519
#> 5       1        0.1        V4 0.3630986
#> 6       1        0.1        V8 0.3803084
#> 7       1        0.1        V2 0.3945035
#> 8       1        0.1        V6 0.3903469
#> 9       1        0.1        V7 0.4069360
#> 10      1        0.1       V10 0.3496814
#> 11      2        0.2        V7 0.4063484
#> 12      2        0.2        V8 0.3798249
#> 13      2        0.2        V9 0.3601420
#> 14      2        0.2       V10 0.3491165
#> 15      2        0.2        V1 0.4053699
#> 16      2        0.2        V5 0.3668249
#> 17      2        0.2        V2 0.3939748
#> 18      2        0.2        V3 0.3391271
#> 19      2        0.2        V4 0.3626116
#> 20      2        0.2        V6 0.3897105
#> 21      3        0.3        V4 0.3621252
#> 22      3        0.3        V5 0.3663537
#> 23      3        0.3        V3 0.3386031
#> 24      3        0.3        V7 0.4057617
#> 25      3        0.3        V8 0.3793421
#> 26      3        0.3        V9 0.3596749
#> 27      3        0.3        V6 0.3890752
#> 28      3        0.3       V10 0.3485525
#> 29      3        0.3        V1 0.4048075
#> 30      3        0.3        V2 0.3934468
#> 31      4        0.4        V1 0.4042459
#> 32      4        0.4        V9 0.3592085
#> 33      4        0.4        V3 0.3380799
#> 34      4        0.4        V4 0.3616395
#> 35      4        0.4        V5 0.3658830
#> 36      4        0.4        V2 0.3929194
#> 37      4        0.4        V6 0.3884409
#> 38      4        0.4        V7 0.4051758
#> 39      4        0.4        V8 0.3788599
#> 40      4        0.4       V10 0.3479894
#> 41      5        0.5        V8 0.3783783
#> 42      5        0.5        V9 0.3587427
#> 43      5        0.5       V10 0.3474272
#> 44      5        0.5        V1 0.4036851
#> 45      5        0.5        V5 0.3654129
#> 46      5        0.5        V2 0.3923928
#> 47      5        0.5        V3 0.3375576
#> 48      5        0.5        V4 0.3611545
#> 49      5        0.5        V6 0.3878076
#> 50      5        0.5        V7 0.4045907
#> 51      6        0.6        V4 0.3606700
#> 52      6        0.6        V5 0.3649434
#> 53      6        0.6        V7 0.4040065
#> 54      6        0.6        V8 0.3778973
#> 55      6        0.6        V9 0.3582774
#> 56      6        0.6        V6 0.3871754
#> 57      6        0.6       V10 0.3468659
#> 58      6        0.6        V1 0.4031251
#> 59      6        0.6        V2 0.3918669
#> 60      6        0.6        V3 0.3370360
#> 61      7        0.7        V1 0.4025658
#> 62      7        0.7        V3 0.3365152
#> 63      7        0.7        V4 0.3601863
#> 64      7        0.7        V5 0.3644745
#> 65      7        0.7        V9 0.3578128
#> 66      7        0.7        V6 0.3865441
#> 67      7        0.7        V7 0.4034231
#> 68      7        0.7        V8 0.3774170
#> 69      7        0.7        V2 0.3913417
#> 70      7        0.7       V10 0.3463056
#> 71      8        0.8        V9 0.3573488
#> 72      8        0.8       V10 0.3457461
#> 73      8        0.8        V1 0.4020073
#> 74      8        0.8        V5 0.3640063
#> 75      8        0.8        V2 0.3908172
#> 76      8        0.8        V3 0.3359953
#> 77      8        0.8        V4 0.3597032
#> 78      8        0.8        V8 0.3769372
#> 79      8        0.8        V6 0.3859140
#> 80      8        0.8        V7 0.4028406
#> 81      9        0.9        V5 0.3635386
#> 82      9        0.9        V7 0.4022589
#> 83      9        0.9        V8 0.3764580
#> 84      9        0.9        V9 0.3568854
#> 85      9        0.9        V6 0.3852848
#> 86      9        0.9       V10 0.3451876
#> 87      9        0.9        V1 0.4014496
#> 88      9        0.9        V2 0.3902934
#> 89      9        0.9        V3 0.3354762
#> 90      9        0.9        V4 0.3592207
#> 91     10        1.0        V1 0.4008926
#> 92     10        1.0        V4 0.3587389
#> 93     10        1.0        V5 0.3630715
#> 94     10        1.0        V9 0.3564226
#> 95     10        1.0        V3 0.3349578
#> 96     10        1.0        V7 0.4016780
#> 97     10        1.0        V8 0.3759795
#> 98     10        1.0        V2 0.3897703
#> 99     10        1.0        V6 0.3846567
#> 100    10        1.0       V10 0.3446299

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4059331
#> 2       1        0.1        V5 0.3672968
#> 3       1        0.1        V9 0.3606096
#> 4       1        0.1        V3 0.3396519
#> 5       1        0.1        V4 0.3630986
#> 6       1        0.1        V8 0.3803084
#> 7       1        0.1        V2 0.3945035
#> 8       1        0.1        V6 0.3903469
#> 9       1        0.1        V7 0.4069360
#> 10      1        0.1       V10 0.3496814
#> 11      2        0.2        V7 0.4063484
#> 12      2        0.2        V8 0.3798249
#> 13      2        0.2        V9 0.3601420
#> 14      2        0.2       V10 0.3491165
#> 15      2        0.2        V1 0.4053699
#> 16      2        0.2        V5 0.3668249
#> 17      2        0.2        V2 0.3939748
#> 18      2        0.2        V3 0.3391271
#> 19      2        0.2        V4 0.3626116
#> 20      2        0.2        V6 0.3897105
#> 21      3        0.3        V4 0.3621252
#> 22      3        0.3        V5 0.3663537
#> 23      3        0.3        V3 0.3386031
#> 24      3        0.3        V7 0.4057617
#> 25      3        0.3        V8 0.3793421
#> 26      3        0.3        V9 0.3596749
#> 27      3        0.3        V6 0.3890752
#> 28      3        0.3       V10 0.3485525
#> 29      3        0.3        V1 0.4048075
#> 30      3        0.3        V2 0.3934468
#> 31      4        0.4        V1 0.4042459
#> 32      4        0.4        V9 0.3592085
#> 33      4        0.4        V3 0.3380799
#> 34      4        0.4        V4 0.3616395
#> 35      4        0.4        V5 0.3658830
#> 36      4        0.4        V2 0.3929194
#> 37      4        0.4        V6 0.3884409
#> 38      4        0.4        V7 0.4051758
#> 39      4        0.4        V8 0.3788599
#> 40      4        0.4       V10 0.3479894
#> 41      5        0.5        V8 0.3783783
#> 42      5        0.5        V9 0.3587427
#> 43      5        0.5       V10 0.3474272
#> 44      5        0.5        V1 0.4036851
#> 45      5        0.5        V5 0.3654129
#> 46      5        0.5        V2 0.3923928
#> 47      5        0.5        V3 0.3375576
#> 48      5        0.5        V4 0.3611545
#> 49      5        0.5        V6 0.3878076
#> 50      5        0.5        V7 0.4045907
#> 51      6        0.6        V4 0.3606700
#> 52      6        0.6        V5 0.3649434
#> 53      6        0.6        V7 0.4040065
#> 54      6        0.6        V8 0.3778973
#> 55      6        0.6        V9 0.3582774
#> 56      6        0.6        V6 0.3871754
#> 57      6        0.6       V10 0.3468659
#> 58      6        0.6        V1 0.4031251
#> 59      6        0.6        V2 0.3918669
#> 60      6        0.6        V3 0.3370360
#> 61      7        0.7        V1 0.4025658
#> 62      7        0.7        V3 0.3365152
#> 63      7        0.7        V4 0.3601863
#> 64      7        0.7        V5 0.3644745
#> 65      7        0.7        V9 0.3578128
#> 66      7        0.7        V6 0.3865441
#> 67      7        0.7        V7 0.4034231
#> 68      7        0.7        V8 0.3774170
#> 69      7        0.7        V2 0.3913417
#> 70      7        0.7       V10 0.3463056
#> 71      8        0.8        V9 0.3573488
#> 72      8        0.8       V10 0.3457461
#> 73      8        0.8        V1 0.4020073
#> 74      8        0.8        V5 0.3640063
#> 75      8        0.8        V2 0.3908172
#> 76      8        0.8        V3 0.3359953
#> 77      8        0.8        V4 0.3597032
#> 78      8        0.8        V8 0.3769372
#> 79      8        0.8        V6 0.3859140
#> 80      8        0.8        V7 0.4028406
#> 81      9        0.9        V5 0.3635386
#> 82      9        0.9        V7 0.4022589
#> 83      9        0.9        V8 0.3764580
#> 84      9        0.9        V9 0.3568854
#> 85      9        0.9        V6 0.3852848
#> 86      9        0.9       V10 0.3451876
#> 87      9        0.9        V1 0.4014496
#> 88      9        0.9        V2 0.3902934
#> 89      9        0.9        V3 0.3354762
#> 90      9        0.9        V4 0.3592207
#> 91     10        1.0        V1 0.4008926
#> 92     10        1.0        V4 0.3587389
#> 93     10        1.0        V5 0.3630715
#> 94     10        1.0        V9 0.3564226
#> 95     10        1.0        V3 0.3349578
#> 96     10        1.0        V7 0.4016780
#> 97     10        1.0        V8 0.3759795
#> 98     10        1.0        V2 0.3897703
#> 99     10        1.0        V6 0.3846567
#> 100    10        1.0       V10 0.3446299


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
#> 1       0                0          0 0.8757906 0.04093389 0.8054063 0.9395164
#> 2       0               20         20 0.8757906 0.04093389 0.8054063 0.9395164
#> 3       0               40         40 0.8757906 0.04093389 0.8054063 0.9395164
#> 4       0               60         60 0.8757906 0.04093389 0.8054063 0.9395164
#> 5      20                0         20 0.8617131 0.04348685 0.7866715 0.9291522
#> 6      20               20         40 0.8617131 0.04348685 0.7866715 0.9291522
#> 7      20               40         60 0.8617131 0.04348685 0.7866715 0.9291522
#> 8      20               60         80 0.8617131 0.04348685 0.7866715 0.9291522
#> 9      40                0         40 0.8478591 0.04589467 0.7685980 0.9186793
#> 10     40               20         60 0.8478591 0.04589467 0.7685980 0.9186793
#> 11     40               40         80 0.8478591 0.04589467 0.7685980 0.9186793
#> 12     40               60        100 0.8478591 0.04589467 0.7685980 0.9186793
#> 13     60                0         60 0.8342249 0.04816804 0.7511185 0.9081385
#> 14     60               20         80 0.8342249 0.04816804 0.7511185 0.9081385
#> 15     60               40        100 0.8342249 0.04816804 0.7511185 0.9081385
#> 16     60               60        120 0.8342249 0.04816804 0.7511185 0.9081385
#> 17     80                0         80 0.8208071 0.05031644 0.7341809 0.8975614
#> 18     80               20        100 0.8208071 0.05031644 0.7341809 0.8975614
#> 19     80               40        120 0.8208071 0.05031644 0.7341809 0.8975614
#> 20     80               60        140 0.8208071 0.05031644 0.7341809 0.8975614
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11514662 0.191053001 0.5820029
#> 2  0.30574618 0.10661061 0.157898337 0.5146835
#> 3  0.26001915 0.09927522 0.129911757 0.4555732
#> 4  0.22113100 0.09276359 0.106321373 0.4038613
#> 5  0.25589195 0.11079136 0.095119693 0.4847685
#> 6  0.21762106 0.10143811 0.077083524 0.4293858
#> 7  0.18507391 0.09324395 0.061993647 0.3809862
#> 8  0.15739448 0.08594865 0.049425486 0.3387345
#> 9  0.18213629 0.10151605 0.043525832 0.4048718
#> 10 0.15489621 0.09237559 0.034153091 0.3595837
#> 11 0.13173012 0.08430582 0.026473014 0.3200523
#> 12 0.11202872 0.07710976 0.020234997 0.2855285
#> 13 0.12963921 0.09041595 0.017373230 0.3395599
#> 14 0.11025053 0.08199179 0.012942847 0.3025690
#> 15 0.09376159 0.07452490 0.009453310 0.2702462
#> 16 0.07973872 0.06786057 0.006749422 0.2419547
#> 17 0.09227334 0.07923167 0.005561142 0.2862033
#> 18 0.07847305 0.07173262 0.003807568 0.2559289
#> 19 0.06673672 0.06507208 0.002523890 0.2294015
#> 20 0.05675565 0.05912685 0.001611929 0.2060945
```
