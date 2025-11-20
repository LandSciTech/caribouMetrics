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
#> 1         0.1 0.3838513 0.01160040 0.3587514 0.3919403
#> 2         0.2 0.3832760 0.01156135 0.3582188 0.3912853
#> 3         0.3 0.3827015 0.01152260 0.3576870 0.3906313
#> 4         0.4 0.3821279 0.01148416 0.3571560 0.3899785
#> 5         0.5 0.3815551 0.01144602 0.3566258 0.3893267
#> 6         0.6 0.3809832 0.01140819 0.3560964 0.3886760
#> 7         0.7 0.3804122 0.01137066 0.3555677 0.3880265
#> 8         0.8 0.3798420 0.01133344 0.3550399 0.3873780
#> 9         0.9 0.3792726 0.01129652 0.3545128 0.3867306
#> 10        1.0 0.3787041 0.01125991 0.3539865 0.3860843

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3568557
#> 2       1        0.1        V5 0.3903802
#> 3       1        0.1        V9 0.3839780
#> 4       1        0.1        V3 0.3679182
#> 5       1        0.1        V4 0.3923933
#> 6       1        0.1        V8 0.3761516
#> 7       1        0.1        V2 0.3777828
#> 8       1        0.1        V6 0.3826105
#> 9       1        0.1        V7 0.3869300
#> 10      1        0.1       V10 0.3652811
#> 11      2        0.2        V7 0.3862962
#> 12      2        0.2        V8 0.3757018
#> 13      2        0.2        V9 0.3834235
#> 14      2        0.2       V10 0.3647565
#> 15      2        0.2        V1 0.3563208
#> 16      2        0.2        V5 0.3898070
#> 17      2        0.2        V2 0.3772008
#> 18      2        0.2        V3 0.3673838
#> 19      2        0.2        V4 0.3917144
#> 20      2        0.2        V6 0.3820175
#> 21      3        0.3        V4 0.3910368
#> 22      3        0.3        V5 0.3892347
#> 23      3        0.3        V3 0.3668503
#> 24      3        0.3        V7 0.3856635
#> 25      3        0.3        V8 0.3752526
#> 26      3        0.3        V9 0.3828698
#> 27      3        0.3        V6 0.3814255
#> 28      3        0.3       V10 0.3642326
#> 29      3        0.3        V1 0.3557867
#> 30      3        0.3        V2 0.3766196
#> 31      4        0.4        V1 0.3552534
#> 32      4        0.4        V9 0.3823168
#> 33      4        0.4        V3 0.3663175
#> 34      4        0.4        V4 0.3903603
#> 35      4        0.4        V5 0.3886632
#> 36      4        0.4        V2 0.3760394
#> 37      4        0.4        V6 0.3808344
#> 38      4        0.4        V7 0.3850318
#> 39      4        0.4        V8 0.3748038
#> 40      4        0.4       V10 0.3637095
#> 41      5        0.5        V8 0.3743556
#> 42      5        0.5        V9 0.3817647
#> 43      5        0.5       V10 0.3631872
#> 44      5        0.5        V1 0.3547209
#> 45      5        0.5        V5 0.3880926
#> 46      5        0.5        V2 0.3754600
#> 47      5        0.5        V3 0.3657855
#> 48      5        0.5        V4 0.3896850
#> 49      5        0.5        V6 0.3802442
#> 50      5        0.5        V7 0.3844011
#> 51      6        0.6        V4 0.3890108
#> 52      6        0.6        V5 0.3875228
#> 53      6        0.6        V7 0.3837715
#> 54      6        0.6        V8 0.3739080
#> 55      6        0.6        V9 0.3812134
#> 56      6        0.6        V6 0.3796549
#> 57      6        0.6       V10 0.3626656
#> 58      6        0.6        V1 0.3541892
#> 59      6        0.6        V2 0.3748816
#> 60      6        0.6        V3 0.3652543
#> 61      7        0.7        V1 0.3536583
#> 62      7        0.7        V3 0.3647238
#> 63      7        0.7        V4 0.3883379
#> 64      7        0.7        V5 0.3869538
#> 65      7        0.7        V9 0.3806628
#> 66      7        0.7        V6 0.3790665
#> 67      7        0.7        V7 0.3831429
#> 68      7        0.7        V8 0.3734609
#> 69      7        0.7        V2 0.3743040
#> 70      7        0.7       V10 0.3621447
#> 71      8        0.8        V9 0.3801131
#> 72      8        0.8       V10 0.3616246
#> 73      8        0.8        V1 0.3531282
#> 74      8        0.8        V5 0.3863857
#> 75      8        0.8        V2 0.3737274
#> 76      8        0.8        V3 0.3641941
#> 77      8        0.8        V4 0.3876661
#> 78      8        0.8        V8 0.3730143
#> 79      8        0.8        V6 0.3784790
#> 80      8        0.8        V7 0.3825153
#> 81      9        0.9        V5 0.3858184
#> 82      9        0.9        V7 0.3818887
#> 83      9        0.9        V8 0.3725682
#> 84      9        0.9        V9 0.3795641
#> 85      9        0.9        V6 0.3778925
#> 86      9        0.9       V10 0.3611053
#> 87      9        0.9        V1 0.3525988
#> 88      9        0.9        V2 0.3731516
#> 89      9        0.9        V3 0.3636652
#> 90      9        0.9        V4 0.3869954
#> 91     10        1.0        V1 0.3520703
#> 92     10        1.0        V4 0.3863259
#> 93     10        1.0        V5 0.3852520
#> 94     10        1.0        V9 0.3790160
#> 95     10        1.0        V3 0.3631371
#> 96     10        1.0        V7 0.3812632
#> 97     10        1.0        V8 0.3721227
#> 98     10        1.0        V2 0.3725767
#> 99     10        1.0        V6 0.3773068
#> 100    10        1.0       V10 0.3605867

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3568557
#> 2       1        0.1        V5 0.3903802
#> 3       1        0.1        V9 0.3839780
#> 4       1        0.1        V3 0.3679182
#> 5       1        0.1        V4 0.3923933
#> 6       1        0.1        V8 0.3761516
#> 7       1        0.1        V2 0.3777828
#> 8       1        0.1        V6 0.3826105
#> 9       1        0.1        V7 0.3869300
#> 10      1        0.1       V10 0.3652811
#> 11      2        0.2        V7 0.3862962
#> 12      2        0.2        V8 0.3757018
#> 13      2        0.2        V9 0.3834235
#> 14      2        0.2       V10 0.3647565
#> 15      2        0.2        V1 0.3563208
#> 16      2        0.2        V5 0.3898070
#> 17      2        0.2        V2 0.3772008
#> 18      2        0.2        V3 0.3673838
#> 19      2        0.2        V4 0.3917144
#> 20      2        0.2        V6 0.3820175
#> 21      3        0.3        V4 0.3910368
#> 22      3        0.3        V5 0.3892347
#> 23      3        0.3        V3 0.3668503
#> 24      3        0.3        V7 0.3856635
#> 25      3        0.3        V8 0.3752526
#> 26      3        0.3        V9 0.3828698
#> 27      3        0.3        V6 0.3814255
#> 28      3        0.3       V10 0.3642326
#> 29      3        0.3        V1 0.3557867
#> 30      3        0.3        V2 0.3766196
#> 31      4        0.4        V1 0.3552534
#> 32      4        0.4        V9 0.3823168
#> 33      4        0.4        V3 0.3663175
#> 34      4        0.4        V4 0.3903603
#> 35      4        0.4        V5 0.3886632
#> 36      4        0.4        V2 0.3760394
#> 37      4        0.4        V6 0.3808344
#> 38      4        0.4        V7 0.3850318
#> 39      4        0.4        V8 0.3748038
#> 40      4        0.4       V10 0.3637095
#> 41      5        0.5        V8 0.3743556
#> 42      5        0.5        V9 0.3817647
#> 43      5        0.5       V10 0.3631872
#> 44      5        0.5        V1 0.3547209
#> 45      5        0.5        V5 0.3880926
#> 46      5        0.5        V2 0.3754600
#> 47      5        0.5        V3 0.3657855
#> 48      5        0.5        V4 0.3896850
#> 49      5        0.5        V6 0.3802442
#> 50      5        0.5        V7 0.3844011
#> 51      6        0.6        V4 0.3890108
#> 52      6        0.6        V5 0.3875228
#> 53      6        0.6        V7 0.3837715
#> 54      6        0.6        V8 0.3739080
#> 55      6        0.6        V9 0.3812134
#> 56      6        0.6        V6 0.3796549
#> 57      6        0.6       V10 0.3626656
#> 58      6        0.6        V1 0.3541892
#> 59      6        0.6        V2 0.3748816
#> 60      6        0.6        V3 0.3652543
#> 61      7        0.7        V1 0.3536583
#> 62      7        0.7        V3 0.3647238
#> 63      7        0.7        V4 0.3883379
#> 64      7        0.7        V5 0.3869538
#> 65      7        0.7        V9 0.3806628
#> 66      7        0.7        V6 0.3790665
#> 67      7        0.7        V7 0.3831429
#> 68      7        0.7        V8 0.3734609
#> 69      7        0.7        V2 0.3743040
#> 70      7        0.7       V10 0.3621447
#> 71      8        0.8        V9 0.3801131
#> 72      8        0.8       V10 0.3616246
#> 73      8        0.8        V1 0.3531282
#> 74      8        0.8        V5 0.3863857
#> 75      8        0.8        V2 0.3737274
#> 76      8        0.8        V3 0.3641941
#> 77      8        0.8        V4 0.3876661
#> 78      8        0.8        V8 0.3730143
#> 79      8        0.8        V6 0.3784790
#> 80      8        0.8        V7 0.3825153
#> 81      9        0.9        V5 0.3858184
#> 82      9        0.9        V7 0.3818887
#> 83      9        0.9        V8 0.3725682
#> 84      9        0.9        V9 0.3795641
#> 85      9        0.9        V6 0.3778925
#> 86      9        0.9       V10 0.3611053
#> 87      9        0.9        V1 0.3525988
#> 88      9        0.9        V2 0.3731516
#> 89      9        0.9        V3 0.3636652
#> 90      9        0.9        V4 0.3869954
#> 91     10        1.0        V1 0.3520703
#> 92     10        1.0        V4 0.3863259
#> 93     10        1.0        V5 0.3852520
#> 94     10        1.0        V9 0.3790160
#> 95     10        1.0        V3 0.3631371
#> 96     10        1.0        V7 0.3812632
#> 97     10        1.0        V8 0.3721227
#> 98     10        1.0        V2 0.3725767
#> 99     10        1.0        V6 0.3773068
#> 100    10        1.0       V10 0.3605867


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
#> 1       0                0          0 0.8757906 0.04421810 0.7894016 0.9411546
#> 2       0               20         20 0.8757906 0.04421810 0.7894016 0.9411546
#> 3       0               40         40 0.8757906 0.04421810 0.7894016 0.9411546
#> 4       0               60         60 0.8757906 0.04421810 0.7894016 0.9411546
#> 5      20                0         20 0.8617131 0.04480763 0.7739805 0.9290959
#> 6      20               20         40 0.8617131 0.04480763 0.7739805 0.9290959
#> 7      20               40         60 0.8617131 0.04480763 0.7739805 0.9290959
#> 8      20               60         80 0.8617131 0.04480763 0.7739805 0.9290959
#> 9      40                0         40 0.8478591 0.04527524 0.7590124 0.9168838
#> 10     40               20         60 0.8478591 0.04527524 0.7590124 0.9168838
#> 11     40               40         80 0.8478591 0.04527524 0.7590124 0.9168838
#> 12     40               60        100 0.8478591 0.04527524 0.7590124 0.9168838
#> 13     60                0         60 0.8342249 0.04564436 0.7444595 0.9045826
#> 14     60               20         80 0.8342249 0.04564436 0.7444595 0.9045826
#> 15     60               40        100 0.8342249 0.04564436 0.7444595 0.9045826
#> 16     60               60        120 0.8342249 0.04564436 0.7444595 0.9045826
#> 17     80                0         80 0.8208071 0.04593297 0.7302911 0.8922403
#> 18     80               20        100 0.8208071 0.04593297 0.7302911 0.8922403
#> 19     80               40        120 0.8208071 0.04593297 0.7302911 0.8922403
#> 20     80               60        140 0.8208071 0.04593297 0.7302911 0.8922403
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12168546 0.1603861272 0.5712200
#> 2  0.30574618 0.11257392 0.1118990313 0.4970146
#> 3  0.26001915 0.10390463 0.0765038644 0.4330581
#> 4  0.22113100 0.09567368 0.0509516436 0.3781431
#> 5  0.25589195 0.11084589 0.0810182278 0.4615094
#> 6  0.21762106 0.10090311 0.0541888280 0.4025561
#> 7  0.18507391 0.09186444 0.0350931084 0.3519880
#> 8  0.15739448 0.08361058 0.0218126870 0.3086354
#> 9  0.18213629 0.09803166 0.0374928158 0.3744657
#> 10 0.15489621 0.08848326 0.0234600589 0.3279071
#> 11 0.13173012 0.07995047 0.0139577513 0.2879798
#> 12 0.11202872 0.07228891 0.0077872424 0.2536891
#> 13 0.12963921 0.08509993 0.0151189757 0.3057320
#> 14 0.11025053 0.07648283 0.0085237342 0.2689439
#> 15 0.09376159 0.06884395 0.0044395797 0.2373105
#> 16 0.07973872 0.06204253 0.0020899504 0.2100271
#> 17 0.09227334 0.07303933 0.0049140970 0.2513889
#> 18 0.07847305 0.06553868 0.0023516657 0.2221815
#> 19 0.06673672 0.05891329 0.0009938533 0.1969396
#> 20 0.05675565 0.05303107 0.0003587947 0.1750272
```
