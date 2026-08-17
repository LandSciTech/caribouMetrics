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
#> 1         0.1 0.3838513 0.01727924 0.3490694 0.3969994
#> 2         0.2 0.3832760 0.01723227 0.3486200 0.3964039
#> 3         0.3 0.3827015 0.01718554 0.3481711 0.3958094
#> 4         0.4 0.3821279 0.01713904 0.3477229 0.3952158
#> 5         0.5 0.3815551 0.01709279 0.3472752 0.3946230
#> 6         0.6 0.3809832 0.01704679 0.3468280 0.3940312
#> 7         0.7 0.3804122 0.01700102 0.3463815 0.3934402
#> 8         0.8 0.3798420 0.01695550 0.3459356 0.3928501
#> 9         0.9 0.3792726 0.01691021 0.3454902 0.3922609
#> 10        1.0 0.3787041 0.01686517 0.3450454 0.3916726

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3659886
#> 2       1        0.1        V5 0.3678029
#> 3       1        0.1        V9 0.3906599
#> 4       1        0.1        V3 0.3940345
#> 5       1        0.1        V4 0.3689362
#> 6       1        0.1        V8 0.3793059
#> 7       1        0.1        V2 0.3949887
#> 8       1        0.1        V6 0.3975831
#> 9       1        0.1        V7 0.3441574
#> 10      1        0.1       V10 0.3665927
#> 11      2        0.2        V7 0.3437336
#> 12      2        0.2        V8 0.3786575
#> 13      2        0.2        V9 0.3901128
#> 14      2        0.2       V10 0.3661138
#> 15      2        0.2        V1 0.3654509
#> 16      2        0.2        V5 0.3671972
#> 17      2        0.2        V2 0.3943522
#> 18      2        0.2        V3 0.3934678
#> 19      2        0.2        V4 0.3683976
#> 20      2        0.2        V6 0.3969996
#> 21      3        0.3        V4 0.3678597
#> 22      3        0.3        V5 0.3665924
#> 23      3        0.3        V3 0.3929018
#> 24      3        0.3        V7 0.3433103
#> 25      3        0.3        V8 0.3780102
#> 26      3        0.3        V9 0.3895665
#> 27      3        0.3        V6 0.3964170
#> 28      3        0.3       V10 0.3656355
#> 29      3        0.3        V1 0.3649140
#> 30      3        0.3        V2 0.3937166
#> 31      4        0.4        V1 0.3643779
#> 32      4        0.4        V9 0.3890209
#> 33      4        0.4        V3 0.3923367
#> 34      4        0.4        V4 0.3673226
#> 35      4        0.4        V5 0.3659886
#> 36      4        0.4        V2 0.3930821
#> 37      4        0.4        V6 0.3958352
#> 38      4        0.4        V7 0.3428875
#> 39      4        0.4        V8 0.3773640
#> 40      4        0.4       V10 0.3651579
#> 41      5        0.5        V8 0.3767189
#> 42      5        0.5        V9 0.3884761
#> 43      5        0.5       V10 0.3646809
#> 44      5        0.5        V1 0.3638425
#> 45      5        0.5        V5 0.3653859
#> 46      5        0.5        V2 0.3924486
#> 47      5        0.5        V3 0.3917724
#> 48      5        0.5        V4 0.3667863
#> 49      5        0.5        V6 0.3952543
#> 50      5        0.5        V7 0.3424653
#> 51      6        0.6        V4 0.3662508
#> 52      6        0.6        V5 0.3647841
#> 53      6        0.6        V7 0.3420435
#> 54      6        0.6        V8 0.3760748
#> 55      6        0.6        V9 0.3879321
#> 56      6        0.6        V6 0.3946742
#> 57      6        0.6       V10 0.3642045
#> 58      6        0.6        V1 0.3633080
#> 59      6        0.6        V2 0.3918161
#> 60      6        0.6        V3 0.3912090
#> 61      7        0.7        V1 0.3627742
#> 62      7        0.7        V3 0.3906463
#> 63      7        0.7        V4 0.3657161
#> 64      7        0.7        V5 0.3641833
#> 65      7        0.7        V9 0.3873888
#> 66      7        0.7        V6 0.3940950
#> 67      7        0.7        V7 0.3416223
#> 68      7        0.7        V8 0.3754319
#> 69      7        0.7        V2 0.3911847
#> 70      7        0.7       V10 0.3637287
#> 71      8        0.8        V9 0.3868463
#> 72      8        0.8       V10 0.3632536
#> 73      8        0.8        V1 0.3622413
#> 74      8        0.8        V5 0.3635835
#> 75      8        0.8        V2 0.3905542
#> 76      8        0.8        V3 0.3900844
#> 77      8        0.8        V4 0.3651821
#> 78      8        0.8        V8 0.3747901
#> 79      8        0.8        V6 0.3935167
#> 80      8        0.8        V7 0.3412016
#> 81      9        0.9        V5 0.3629846
#> 82      9        0.9        V7 0.3407815
#> 83      9        0.9        V8 0.3741494
#> 84      9        0.9        V9 0.3863046
#> 85      9        0.9        V6 0.3929392
#> 86      9        0.9       V10 0.3627790
#> 87      9        0.9        V1 0.3617091
#> 88      9        0.9        V2 0.3899248
#> 89      9        0.9        V3 0.3895234
#> 90      9        0.9        V4 0.3646490
#> 91     10        1.0        V1 0.3611777
#> 92     10        1.0        V4 0.3641166
#> 93     10        1.0        V5 0.3623868
#> 94     10        1.0        V9 0.3857636
#> 95     10        1.0        V3 0.3889631
#> 96     10        1.0        V7 0.3403618
#> 97     10        1.0        V8 0.3735098
#> 98     10        1.0        V2 0.3892964
#> 99     10        1.0        V6 0.3923625
#> 100    10        1.0       V10 0.3623051

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3659886
#> 2       1        0.1        V5 0.3678029
#> 3       1        0.1        V9 0.3906599
#> 4       1        0.1        V3 0.3940345
#> 5       1        0.1        V4 0.3689362
#> 6       1        0.1        V8 0.3793059
#> 7       1        0.1        V2 0.3949887
#> 8       1        0.1        V6 0.3975831
#> 9       1        0.1        V7 0.3441574
#> 10      1        0.1       V10 0.3665927
#> 11      2        0.2        V7 0.3437336
#> 12      2        0.2        V8 0.3786575
#> 13      2        0.2        V9 0.3901128
#> 14      2        0.2       V10 0.3661138
#> 15      2        0.2        V1 0.3654509
#> 16      2        0.2        V5 0.3671972
#> 17      2        0.2        V2 0.3943522
#> 18      2        0.2        V3 0.3934678
#> 19      2        0.2        V4 0.3683976
#> 20      2        0.2        V6 0.3969996
#> 21      3        0.3        V4 0.3678597
#> 22      3        0.3        V5 0.3665924
#> 23      3        0.3        V3 0.3929018
#> 24      3        0.3        V7 0.3433103
#> 25      3        0.3        V8 0.3780102
#> 26      3        0.3        V9 0.3895665
#> 27      3        0.3        V6 0.3964170
#> 28      3        0.3       V10 0.3656355
#> 29      3        0.3        V1 0.3649140
#> 30      3        0.3        V2 0.3937166
#> 31      4        0.4        V1 0.3643779
#> 32      4        0.4        V9 0.3890209
#> 33      4        0.4        V3 0.3923367
#> 34      4        0.4        V4 0.3673226
#> 35      4        0.4        V5 0.3659886
#> 36      4        0.4        V2 0.3930821
#> 37      4        0.4        V6 0.3958352
#> 38      4        0.4        V7 0.3428875
#> 39      4        0.4        V8 0.3773640
#> 40      4        0.4       V10 0.3651579
#> 41      5        0.5        V8 0.3767189
#> 42      5        0.5        V9 0.3884761
#> 43      5        0.5       V10 0.3646809
#> 44      5        0.5        V1 0.3638425
#> 45      5        0.5        V5 0.3653859
#> 46      5        0.5        V2 0.3924486
#> 47      5        0.5        V3 0.3917724
#> 48      5        0.5        V4 0.3667863
#> 49      5        0.5        V6 0.3952543
#> 50      5        0.5        V7 0.3424653
#> 51      6        0.6        V4 0.3662508
#> 52      6        0.6        V5 0.3647841
#> 53      6        0.6        V7 0.3420435
#> 54      6        0.6        V8 0.3760748
#> 55      6        0.6        V9 0.3879321
#> 56      6        0.6        V6 0.3946742
#> 57      6        0.6       V10 0.3642045
#> 58      6        0.6        V1 0.3633080
#> 59      6        0.6        V2 0.3918161
#> 60      6        0.6        V3 0.3912090
#> 61      7        0.7        V1 0.3627742
#> 62      7        0.7        V3 0.3906463
#> 63      7        0.7        V4 0.3657161
#> 64      7        0.7        V5 0.3641833
#> 65      7        0.7        V9 0.3873888
#> 66      7        0.7        V6 0.3940950
#> 67      7        0.7        V7 0.3416223
#> 68      7        0.7        V8 0.3754319
#> 69      7        0.7        V2 0.3911847
#> 70      7        0.7       V10 0.3637287
#> 71      8        0.8        V9 0.3868463
#> 72      8        0.8       V10 0.3632536
#> 73      8        0.8        V1 0.3622413
#> 74      8        0.8        V5 0.3635835
#> 75      8        0.8        V2 0.3905542
#> 76      8        0.8        V3 0.3900844
#> 77      8        0.8        V4 0.3651821
#> 78      8        0.8        V8 0.3747901
#> 79      8        0.8        V6 0.3935167
#> 80      8        0.8        V7 0.3412016
#> 81      9        0.9        V5 0.3629846
#> 82      9        0.9        V7 0.3407815
#> 83      9        0.9        V8 0.3741494
#> 84      9        0.9        V9 0.3863046
#> 85      9        0.9        V6 0.3929392
#> 86      9        0.9       V10 0.3627790
#> 87      9        0.9        V1 0.3617091
#> 88      9        0.9        V2 0.3899248
#> 89      9        0.9        V3 0.3895234
#> 90      9        0.9        V4 0.3646490
#> 91     10        1.0        V1 0.3611777
#> 92     10        1.0        V4 0.3641166
#> 93     10        1.0        V5 0.3623868
#> 94     10        1.0        V9 0.3857636
#> 95     10        1.0        V3 0.3889631
#> 96     10        1.0        V7 0.3403618
#> 97     10        1.0        V8 0.3735098
#> 98     10        1.0        V2 0.3892964
#> 99     10        1.0        V6 0.3923625
#> 100    10        1.0       V10 0.3623051


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
#> 1       0                0          0 0.8757906 0.04629076 0.7937337 0.9535915
#> 2       0               20         20 0.8757906 0.04629076 0.7937337 0.9535915
#> 3       0               40         40 0.8757906 0.04629076 0.7937337 0.9535915
#> 4       0               60         60 0.8757906 0.04629076 0.7937337 0.9535915
#> 5      20                0         20 0.8617131 0.04958665 0.7753291 0.9468446
#> 6      20               20         40 0.8617131 0.04958665 0.7753291 0.9468446
#> 7      20               40         60 0.8617131 0.04958665 0.7753291 0.9468446
#> 8      20               60         80 0.8617131 0.04958665 0.7753291 0.9468446
#> 9      40                0         40 0.8478591 0.05271780 0.7575491 0.9399899
#> 10     40               20         60 0.8478591 0.05271780 0.7575491 0.9399899
#> 11     40               40         80 0.8478591 0.05271780 0.7575491 0.9399899
#> 12     40               60        100 0.8478591 0.05271780 0.7575491 0.9399899
#> 13     60                0         60 0.8342249 0.05569729 0.7403355 0.9330497
#> 14     60               20         80 0.8342249 0.05569729 0.7403355 0.9330497
#> 15     60               40        100 0.8342249 0.05569729 0.7403355 0.9330497
#> 16     60               60        120 0.8342249 0.05569729 0.7403355 0.9330497
#> 17     80                0         80 0.8208071 0.05853624 0.7236421 0.9260425
#> 18     80               20        100 0.8208071 0.05853624 0.7236421 0.9260425
#> 19     80               40        120 0.8208071 0.05853624 0.7236421 0.9260425
#> 20     80               60        140 0.8208071 0.05853624 0.7236421 0.9260425
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.13114512 0.178089051 0.6389652
#> 2  0.30574618 0.13591677 0.136183698 0.6085033
#> 3  0.26001915 0.13901688 0.103124770 0.5794714
#> 4  0.22113100 0.14018199 0.077156503 0.5518480
#> 5  0.25589195 0.12025595 0.098250065 0.5228559
#> 6  0.21762106 0.12162020 0.073341045 0.4980766
#> 7  0.18507391 0.12194980 0.053928035 0.4745725
#> 8  0.15739448 0.12114734 0.038939668 0.4522906
#> 9  0.18213629 0.10728913 0.051092800 0.4289725
#> 10 0.15489621 0.10681851 0.036765764 0.4090852
#> 11 0.13173012 0.10571874 0.025865443 0.3902477
#> 12 0.11202872 0.10396493 0.017707794 0.3724053
#> 13 0.12963921 0.09417704 0.024301339 0.3537414
#> 14 0.11025053 0.09278053 0.016551219 0.3378234
#> 15 0.09376159 0.09101512 0.010895488 0.3227400
#> 16 0.07973872 0.08888562 0.006882295 0.3084433
#> 17 0.09227334 0.08184197 0.010108261 0.2934721
#> 18 0.07847305 0.08004731 0.006335072 0.2806857
#> 19 0.06673672 0.07805000 0.003769695 0.2685501
#> 20 0.05675565 0.07586483 0.002106134 0.2570258
```
