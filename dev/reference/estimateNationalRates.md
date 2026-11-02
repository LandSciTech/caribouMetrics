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
[`addMissingYears()`](https://landscitech.github.io/caribouMetrics/dev/reference/addMissingYears.md),
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
[`getCaribouYear()`](https://landscitech.github.io/caribouMetrics/dev/reference/getCaribouYear.md),
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
#> 1         0.1 0.3838513 0.01852023 0.3385657 0.3919999
#> 2         0.2 0.3832760 0.01848545 0.3380510 0.3914335
#> 3         0.3 0.3827015 0.01845096 0.3375370 0.3908680
#> 4         0.4 0.3821279 0.01841676 0.3370239 0.3903032
#> 5         0.5 0.3815551 0.01838284 0.3365115 0.3897393
#> 6         0.6 0.3809832 0.01834920 0.3359999 0.3891763
#> 7         0.7 0.3804122 0.01831584 0.3354891 0.3886140
#> 8         0.8 0.3798420 0.01828276 0.3349790 0.3880526
#> 9         0.9 0.3792726 0.01824996 0.3344698 0.3874921
#> 10        1.0 0.3787041 0.01821744 0.3339613 0.3869323

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3906940
#> 2       1        0.1        V5 0.3647000
#> 3       1        0.1        V9 0.3923791
#> 4       1        0.1        V3 0.3852854
#> 5       1        0.1        V4 0.3585292
#> 6       1        0.1        V8 0.3793103
#> 7       1        0.1        V2 0.3566114
#> 8       1        0.1        V6 0.3333266
#> 9       1        0.1        V7 0.3703142
#> 10      1        0.1       V10 0.3570558
#> 11      2        0.2        V7 0.3698328
#> 12      2        0.2        V8 0.3786963
#> 13      2        0.2        V9 0.3918591
#> 14      2        0.2       V10 0.3564955
#> 15      2        0.2        V1 0.3899677
#> 16      2        0.2        V5 0.3641503
#> 17      2        0.2        V2 0.3560021
#> 18      2        0.2        V3 0.3847000
#> 19      2        0.2        V4 0.3579921
#> 20      2        0.2        V6 0.3328394
#> 21      3        0.3        V4 0.3574557
#> 22      3        0.3        V5 0.3636014
#> 23      3        0.3        V3 0.3841155
#> 24      3        0.3        V7 0.3693520
#> 25      3        0.3        V8 0.3780833
#> 26      3        0.3        V9 0.3913398
#> 27      3        0.3        V6 0.3323528
#> 28      3        0.3       V10 0.3559360
#> 29      3        0.3        V1 0.3892427
#> 30      3        0.3        V2 0.3553938
#> 31      4        0.4        V1 0.3885191
#> 32      4        0.4        V9 0.3908212
#> 33      4        0.4        V3 0.3835319
#> 34      4        0.4        V4 0.3569201
#> 35      4        0.4        V5 0.3630533
#> 36      4        0.4        V2 0.3547865
#> 37      4        0.4        V6 0.3318670
#> 38      4        0.4        V7 0.3688719
#> 39      4        0.4        V8 0.3774712
#> 40      4        0.4       V10 0.3553774
#> 41      5        0.5        V8 0.3768602
#> 42      5        0.5        V9 0.3903033
#> 43      5        0.5       V10 0.3548197
#> 44      5        0.5        V1 0.3877969
#> 45      5        0.5        V5 0.3625060
#> 46      5        0.5        V2 0.3541803
#> 47      5        0.5        V3 0.3829492
#> 48      5        0.5        V4 0.3563853
#> 49      5        0.5        V6 0.3313818
#> 50      5        0.5        V7 0.3683923
#> 51      6        0.6        V4 0.3558514
#> 52      6        0.6        V5 0.3619596
#> 53      6        0.6        V7 0.3679134
#> 54      6        0.6        V8 0.3762501
#> 55      6        0.6        V9 0.3897860
#> 56      6        0.6        V6 0.3308974
#> 57      6        0.6       V10 0.3542628
#> 58      6        0.6        V1 0.3870759
#> 59      6        0.6        V2 0.3535751
#> 60      6        0.6        V3 0.3823673
#> 61      7        0.7        V1 0.3863564
#> 62      7        0.7        V3 0.3817863
#> 63      7        0.7        V4 0.3553182
#> 64      7        0.7        V5 0.3614140
#> 65      7        0.7        V9 0.3892695
#> 66      7        0.7        V6 0.3304137
#> 67      7        0.7        V7 0.3674351
#> 68      7        0.7        V8 0.3756410
#> 69      7        0.7        V2 0.3529710
#> 70      7        0.7       V10 0.3537068
#> 71      8        0.8        V9 0.3887536
#> 72      8        0.8       V10 0.3531517
#> 73      8        0.8        V1 0.3856381
#> 74      8        0.8        V5 0.3608692
#> 75      8        0.8        V2 0.3523679
#> 76      8        0.8        V3 0.3812063
#> 77      8        0.8        V4 0.3547858
#> 78      8        0.8        V8 0.3750330
#> 79      8        0.8        V6 0.3299307
#> 80      8        0.8        V7 0.3669575
#> 81      9        0.9        V5 0.3603252
#> 82      9        0.9        V7 0.3664804
#> 83      9        0.9        V8 0.3744258
#> 84      9        0.9        V9 0.3882384
#> 85      9        0.9        V6 0.3294484
#> 86      9        0.9       V10 0.3525975
#> 87      9        0.9        V1 0.3849212
#> 88      9        0.9        V2 0.3517658
#> 89      9        0.9        V3 0.3806271
#> 90      9        0.9        V4 0.3542543
#> 91     10        1.0        V1 0.3842057
#> 92     10        1.0        V4 0.3537235
#> 93     10        1.0        V5 0.3597821
#> 94     10        1.0        V9 0.3877239
#> 95     10        1.0        V3 0.3800487
#> 96     10        1.0        V7 0.3660040
#> 97     10        1.0        V8 0.3738197
#> 98     10        1.0        V2 0.3511647
#> 99     10        1.0        V6 0.3289668
#> 100    10        1.0       V10 0.3520441

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3906940
#> 2       1        0.1        V5 0.3647000
#> 3       1        0.1        V9 0.3923791
#> 4       1        0.1        V3 0.3852854
#> 5       1        0.1        V4 0.3585292
#> 6       1        0.1        V8 0.3793103
#> 7       1        0.1        V2 0.3566114
#> 8       1        0.1        V6 0.3333266
#> 9       1        0.1        V7 0.3703142
#> 10      1        0.1       V10 0.3570558
#> 11      2        0.2        V7 0.3698328
#> 12      2        0.2        V8 0.3786963
#> 13      2        0.2        V9 0.3918591
#> 14      2        0.2       V10 0.3564955
#> 15      2        0.2        V1 0.3899677
#> 16      2        0.2        V5 0.3641503
#> 17      2        0.2        V2 0.3560021
#> 18      2        0.2        V3 0.3847000
#> 19      2        0.2        V4 0.3579921
#> 20      2        0.2        V6 0.3328394
#> 21      3        0.3        V4 0.3574557
#> 22      3        0.3        V5 0.3636014
#> 23      3        0.3        V3 0.3841155
#> 24      3        0.3        V7 0.3693520
#> 25      3        0.3        V8 0.3780833
#> 26      3        0.3        V9 0.3913398
#> 27      3        0.3        V6 0.3323528
#> 28      3        0.3       V10 0.3559360
#> 29      3        0.3        V1 0.3892427
#> 30      3        0.3        V2 0.3553938
#> 31      4        0.4        V1 0.3885191
#> 32      4        0.4        V9 0.3908212
#> 33      4        0.4        V3 0.3835319
#> 34      4        0.4        V4 0.3569201
#> 35      4        0.4        V5 0.3630533
#> 36      4        0.4        V2 0.3547865
#> 37      4        0.4        V6 0.3318670
#> 38      4        0.4        V7 0.3688719
#> 39      4        0.4        V8 0.3774712
#> 40      4        0.4       V10 0.3553774
#> 41      5        0.5        V8 0.3768602
#> 42      5        0.5        V9 0.3903033
#> 43      5        0.5       V10 0.3548197
#> 44      5        0.5        V1 0.3877969
#> 45      5        0.5        V5 0.3625060
#> 46      5        0.5        V2 0.3541803
#> 47      5        0.5        V3 0.3829492
#> 48      5        0.5        V4 0.3563853
#> 49      5        0.5        V6 0.3313818
#> 50      5        0.5        V7 0.3683923
#> 51      6        0.6        V4 0.3558514
#> 52      6        0.6        V5 0.3619596
#> 53      6        0.6        V7 0.3679134
#> 54      6        0.6        V8 0.3762501
#> 55      6        0.6        V9 0.3897860
#> 56      6        0.6        V6 0.3308974
#> 57      6        0.6       V10 0.3542628
#> 58      6        0.6        V1 0.3870759
#> 59      6        0.6        V2 0.3535751
#> 60      6        0.6        V3 0.3823673
#> 61      7        0.7        V1 0.3863564
#> 62      7        0.7        V3 0.3817863
#> 63      7        0.7        V4 0.3553182
#> 64      7        0.7        V5 0.3614140
#> 65      7        0.7        V9 0.3892695
#> 66      7        0.7        V6 0.3304137
#> 67      7        0.7        V7 0.3674351
#> 68      7        0.7        V8 0.3756410
#> 69      7        0.7        V2 0.3529710
#> 70      7        0.7       V10 0.3537068
#> 71      8        0.8        V9 0.3887536
#> 72      8        0.8       V10 0.3531517
#> 73      8        0.8        V1 0.3856381
#> 74      8        0.8        V5 0.3608692
#> 75      8        0.8        V2 0.3523679
#> 76      8        0.8        V3 0.3812063
#> 77      8        0.8        V4 0.3547858
#> 78      8        0.8        V8 0.3750330
#> 79      8        0.8        V6 0.3299307
#> 80      8        0.8        V7 0.3669575
#> 81      9        0.9        V5 0.3603252
#> 82      9        0.9        V7 0.3664804
#> 83      9        0.9        V8 0.3744258
#> 84      9        0.9        V9 0.3882384
#> 85      9        0.9        V6 0.3294484
#> 86      9        0.9       V10 0.3525975
#> 87      9        0.9        V1 0.3849212
#> 88      9        0.9        V2 0.3517658
#> 89      9        0.9        V3 0.3806271
#> 90      9        0.9        V4 0.3542543
#> 91     10        1.0        V1 0.3842057
#> 92     10        1.0        V4 0.3537235
#> 93     10        1.0        V5 0.3597821
#> 94     10        1.0        V9 0.3877239
#> 95     10        1.0        V3 0.3800487
#> 96     10        1.0        V7 0.3660040
#> 97     10        1.0        V8 0.3738197
#> 98     10        1.0        V2 0.3511647
#> 99     10        1.0        V6 0.3289668
#> 100    10        1.0       V10 0.3520441


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
#> 1       0                0          0 0.8757906 0.04748345 0.7872640 0.9523004
#> 2       0               20         20 0.8757906 0.04748345 0.7872640 0.9523004
#> 3       0               40         40 0.8757906 0.04748345 0.7872640 0.9523004
#> 4       0               60         60 0.8757906 0.04748345 0.7872640 0.9523004
#> 5      20                0         20 0.8617131 0.04946439 0.7690891 0.9397499
#> 6      20               20         40 0.8617131 0.04946439 0.7690891 0.9397499
#> 7      20               40         60 0.8617131 0.04946439 0.7690891 0.9397499
#> 8      20               60         80 0.8617131 0.04946439 0.7690891 0.9397499
#> 9      40                0         40 0.8478591 0.05129415 0.7515252 0.9269129
#> 10     40               20         60 0.8478591 0.05129415 0.7515252 0.9269129
#> 11     40               40         80 0.8478591 0.05129415 0.7515252 0.9269129
#> 12     40               60        100 0.8478591 0.05129415 0.7515252 0.9269129
#> 13     60                0         60 0.8342249 0.05299647 0.7345165 0.9138955
#> 14     60               20         80 0.8342249 0.05299647 0.7345165 0.9138955
#> 15     60               40        100 0.8342249 0.05299647 0.7345165 0.9138955
#> 16     60               60        120 0.8342249 0.05299647 0.7345165 0.9138955
#> 17     80                0         80 0.8208071 0.05458890 0.7180188 0.9007741
#> 18     80               20        100 0.8208071 0.05458890 0.7180188 0.9007741
#> 19     80               40        120 0.8208071 0.05458890 0.7180188 0.9007741
#> 20     80               60        140 0.8208071 0.05458890 0.7180188 0.9007741
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11447693 0.165055221 0.5642412
#> 2  0.30574618 0.10473439 0.133860882 0.5111462
#> 3  0.26001915 0.09769369 0.107907517 0.4633779
#> 4  0.22113100 0.09219758 0.086369623 0.4205051
#> 5  0.25589195 0.10467462 0.089348905 0.4551078
#> 6  0.21762106 0.09498423 0.071019466 0.4130894
#> 7  0.18507391 0.08750727 0.055923283 0.3754306
#> 8  0.15739448 0.08144334 0.043561502 0.3416891
#> 9  0.18213629 0.09304978 0.045259845 0.3689193
#> 10 0.15489621 0.08423490 0.034886575 0.3358546
#> 11 0.13173012 0.07717100 0.026513583 0.3062159
#> 12 0.11202872 0.07129965 0.019823955 0.2796223
#> 13 0.12963921 0.08130057 0.020731784 0.3010879
#> 14 0.11025053 0.07361116 0.015256951 0.2750175
#> 15 0.09376159 0.06728980 0.010990312 0.2515844
#> 16 0.07973872 0.06194286 0.007722950 0.2304802
#> 17 0.09227334 0.07026425 0.008157032 0.2475216
#> 18 0.07847305 0.06370432 0.005595079 0.2268161
#> 19 0.06673672 0.05820875 0.003714485 0.2081146
#> 20 0.05675565 0.05349738 0.002374301 0.1911746
```
