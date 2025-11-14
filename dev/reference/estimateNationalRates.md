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
#> 1         0.1 0.3838513 0.01888476 0.3522928 0.4122606
#> 2         0.2 0.3832760 0.01884111 0.3517551 0.4115563
#> 3         0.3 0.3827015 0.01879770 0.3512182 0.4108532
#> 4         0.4 0.3821279 0.01875452 0.3506821 0.4101513
#> 5         0.5 0.3815551 0.01871157 0.3501468 0.4094507
#> 6         0.6 0.3809832 0.01866885 0.3496124 0.4087512
#> 7         0.7 0.3804122 0.01862637 0.3490787 0.4080529
#> 8         0.8 0.3798420 0.01858411 0.3485459 0.4073558
#> 9         0.9 0.3792726 0.01854208 0.3480139 0.4066599
#> 10        1.0 0.3787041 0.01850028 0.3474827 0.4059653

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3824236
#> 2       1        0.1        V5 0.3778068
#> 3       1        0.1        V9 0.4145054
#> 4       1        0.1        V3 0.3884498
#> 5       1        0.1        V4 0.3965051
#> 6       1        0.1        V8 0.3448855
#> 7       1        0.1        V2 0.4045286
#> 8       1        0.1        V6 0.3832453
#> 9       1        0.1        V7 0.3998937
#> 10      1        0.1       V10 0.3900034
#> 11      2        0.2        V7 0.3992492
#> 12      2        0.2        V8 0.3443590
#> 13      2        0.2        V9 0.4137722
#> 14      2        0.2       V10 0.3894466
#> 15      2        0.2        V1 0.3818263
#> 16      2        0.2        V5 0.3772305
#> 17      2        0.2        V2 0.4039238
#> 18      2        0.2        V3 0.3879612
#> 19      2        0.2        V4 0.3958905
#> 20      2        0.2        V6 0.3826184
#> 21      3        0.3        V4 0.3952770
#> 22      3        0.3        V5 0.3766551
#> 23      3        0.3        V3 0.3874731
#> 24      3        0.3        V7 0.3986058
#> 25      3        0.3        V8 0.3438332
#> 26      3        0.3        V9 0.4130403
#> 27      3        0.3        V6 0.3819926
#> 28      3        0.3       V10 0.3888906
#> 29      3        0.3        V1 0.3812300
#> 30      3        0.3        V2 0.4033199
#> 31      4        0.4        V1 0.3806345
#> 32      4        0.4        V9 0.4123098
#> 33      4        0.4        V3 0.3869857
#> 34      4        0.4        V4 0.3946643
#> 35      4        0.4        V5 0.3760805
#> 36      4        0.4        V2 0.4027168
#> 37      4        0.4        V6 0.3813677
#> 38      4        0.4        V7 0.3979633
#> 39      4        0.4        V8 0.3433083
#> 40      4        0.4       V10 0.3883354
#> 41      5        0.5        V8 0.3427842
#> 42      5        0.5        V9 0.4115805
#> 43      5        0.5       V10 0.3877810
#> 44      5        0.5        V1 0.3800400
#> 45      5        0.5        V5 0.3755069
#> 46      5        0.5        V2 0.4021147
#> 47      5        0.5        V3 0.3864989
#> 48      5        0.5        V4 0.3940527
#> 49      5        0.5        V6 0.3807439
#> 50      5        0.5        V7 0.3973219
#> 51      6        0.6        V4 0.3934419
#> 52      6        0.6        V5 0.3749341
#> 53      6        0.6        V7 0.3966816
#> 54      6        0.6        V8 0.3422609
#> 55      6        0.6        V9 0.4108525
#> 56      6        0.6        V6 0.3801211
#> 57      6        0.6       V10 0.3872274
#> 58      6        0.6        V1 0.3794465
#> 59      6        0.6        V2 0.4015135
#> 60      6        0.6        V3 0.3860127
#> 61      7        0.7        V1 0.3788538
#> 62      7        0.7        V3 0.3855271
#> 63      7        0.7        V4 0.3928321
#> 64      7        0.7        V5 0.3743621
#> 65      7        0.7        V9 0.4101257
#> 66      7        0.7        V6 0.3794993
#> 67      7        0.7        V7 0.3960422
#> 68      7        0.7        V8 0.3417384
#> 69      7        0.7        V2 0.4009132
#> 70      7        0.7       V10 0.3866745
#> 71      8        0.8        V9 0.4094003
#> 72      8        0.8       V10 0.3861225
#> 73      8        0.8        V1 0.3782621
#> 74      8        0.8        V5 0.3737911
#> 75      8        0.8        V2 0.4003137
#> 76      8        0.8        V3 0.3850421
#> 77      8        0.8        V4 0.3922233
#> 78      8        0.8        V8 0.3412167
#> 79      8        0.8        V6 0.3788786
#> 80      8        0.8        V7 0.3954039
#> 81      9        0.9        V5 0.3732209
#> 82      9        0.9        V7 0.3947667
#> 83      9        0.9        V8 0.3406958
#> 84      9        0.9        V9 0.4086762
#> 85      9        0.9        V6 0.3782588
#> 86      9        0.9       V10 0.3855712
#> 87      9        0.9        V1 0.3776714
#> 88      9        0.9        V2 0.3997152
#> 89      9        0.9        V3 0.3845578
#> 90      9        0.9        V4 0.3916154
#> 91     10        1.0        V1 0.3770815
#> 92     10        1.0        V4 0.3910085
#> 93     10        1.0        V5 0.3726516
#> 94     10        1.0        V9 0.4079533
#> 95     10        1.0        V3 0.3840740
#> 96     10        1.0        V7 0.3941304
#> 97     10        1.0        V8 0.3401756
#> 98     10        1.0        V2 0.3991176
#> 99     10        1.0        V6 0.3776401
#> 100    10        1.0       V10 0.3850208

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3824236
#> 2       1        0.1        V5 0.3778068
#> 3       1        0.1        V9 0.4145054
#> 4       1        0.1        V3 0.3884498
#> 5       1        0.1        V4 0.3965051
#> 6       1        0.1        V8 0.3448855
#> 7       1        0.1        V2 0.4045286
#> 8       1        0.1        V6 0.3832453
#> 9       1        0.1        V7 0.3998937
#> 10      1        0.1       V10 0.3900034
#> 11      2        0.2        V7 0.3992492
#> 12      2        0.2        V8 0.3443590
#> 13      2        0.2        V9 0.4137722
#> 14      2        0.2       V10 0.3894466
#> 15      2        0.2        V1 0.3818263
#> 16      2        0.2        V5 0.3772305
#> 17      2        0.2        V2 0.4039238
#> 18      2        0.2        V3 0.3879612
#> 19      2        0.2        V4 0.3958905
#> 20      2        0.2        V6 0.3826184
#> 21      3        0.3        V4 0.3952770
#> 22      3        0.3        V5 0.3766551
#> 23      3        0.3        V3 0.3874731
#> 24      3        0.3        V7 0.3986058
#> 25      3        0.3        V8 0.3438332
#> 26      3        0.3        V9 0.4130403
#> 27      3        0.3        V6 0.3819926
#> 28      3        0.3       V10 0.3888906
#> 29      3        0.3        V1 0.3812300
#> 30      3        0.3        V2 0.4033199
#> 31      4        0.4        V1 0.3806345
#> 32      4        0.4        V9 0.4123098
#> 33      4        0.4        V3 0.3869857
#> 34      4        0.4        V4 0.3946643
#> 35      4        0.4        V5 0.3760805
#> 36      4        0.4        V2 0.4027168
#> 37      4        0.4        V6 0.3813677
#> 38      4        0.4        V7 0.3979633
#> 39      4        0.4        V8 0.3433083
#> 40      4        0.4       V10 0.3883354
#> 41      5        0.5        V8 0.3427842
#> 42      5        0.5        V9 0.4115805
#> 43      5        0.5       V10 0.3877810
#> 44      5        0.5        V1 0.3800400
#> 45      5        0.5        V5 0.3755069
#> 46      5        0.5        V2 0.4021147
#> 47      5        0.5        V3 0.3864989
#> 48      5        0.5        V4 0.3940527
#> 49      5        0.5        V6 0.3807439
#> 50      5        0.5        V7 0.3973219
#> 51      6        0.6        V4 0.3934419
#> 52      6        0.6        V5 0.3749341
#> 53      6        0.6        V7 0.3966816
#> 54      6        0.6        V8 0.3422609
#> 55      6        0.6        V9 0.4108525
#> 56      6        0.6        V6 0.3801211
#> 57      6        0.6       V10 0.3872274
#> 58      6        0.6        V1 0.3794465
#> 59      6        0.6        V2 0.4015135
#> 60      6        0.6        V3 0.3860127
#> 61      7        0.7        V1 0.3788538
#> 62      7        0.7        V3 0.3855271
#> 63      7        0.7        V4 0.3928321
#> 64      7        0.7        V5 0.3743621
#> 65      7        0.7        V9 0.4101257
#> 66      7        0.7        V6 0.3794993
#> 67      7        0.7        V7 0.3960422
#> 68      7        0.7        V8 0.3417384
#> 69      7        0.7        V2 0.4009132
#> 70      7        0.7       V10 0.3866745
#> 71      8        0.8        V9 0.4094003
#> 72      8        0.8       V10 0.3861225
#> 73      8        0.8        V1 0.3782621
#> 74      8        0.8        V5 0.3737911
#> 75      8        0.8        V2 0.4003137
#> 76      8        0.8        V3 0.3850421
#> 77      8        0.8        V4 0.3922233
#> 78      8        0.8        V8 0.3412167
#> 79      8        0.8        V6 0.3788786
#> 80      8        0.8        V7 0.3954039
#> 81      9        0.9        V5 0.3732209
#> 82      9        0.9        V7 0.3947667
#> 83      9        0.9        V8 0.3406958
#> 84      9        0.9        V9 0.4086762
#> 85      9        0.9        V6 0.3782588
#> 86      9        0.9       V10 0.3855712
#> 87      9        0.9        V1 0.3776714
#> 88      9        0.9        V2 0.3997152
#> 89      9        0.9        V3 0.3845578
#> 90      9        0.9        V4 0.3916154
#> 91     10        1.0        V1 0.3770815
#> 92     10        1.0        V4 0.3910085
#> 93     10        1.0        V5 0.3726516
#> 94     10        1.0        V9 0.4079533
#> 95     10        1.0        V3 0.3840740
#> 96     10        1.0        V7 0.3941304
#> 97     10        1.0        V8 0.3401756
#> 98     10        1.0        V2 0.3991176
#> 99     10        1.0        V6 0.3776401
#> 100    10        1.0       V10 0.3850208


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
#> 1       0                0          0 0.8757906 0.04705034 0.7795994 0.9382899
#> 2       0               20         20 0.8757906 0.04705034 0.7795994 0.9382899
#> 3       0               40         40 0.8757906 0.04705034 0.7795994 0.9382899
#> 4       0               60         60 0.8757906 0.04705034 0.7795994 0.9382899
#> 5      20                0         20 0.8617131 0.04880015 0.7608653 0.9265026
#> 6      20               20         40 0.8617131 0.04880015 0.7608653 0.9265026
#> 7      20               40         60 0.8617131 0.04880015 0.7608653 0.9265026
#> 8      20               60         80 0.8617131 0.04880015 0.7608653 0.9265026
#> 9      40                0         40 0.8478591 0.05042216 0.7427584 0.9145841
#> 10     40               20         60 0.8478591 0.05042216 0.7427584 0.9145841
#> 11     40               40         80 0.8478591 0.05042216 0.7427584 0.9145841
#> 12     40               60        100 0.8478591 0.05042216 0.7427584 0.9145841
#> 13     60                0         60 0.8342249 0.05193343 0.7252247 0.9025904
#> 14     60               20         80 0.8342249 0.05193343 0.7252247 0.9025904
#> 15     60               40        100 0.8342249 0.05193343 0.7252247 0.9025904
#> 16     60               60        120 0.8342249 0.05193343 0.7252247 0.9025904
#> 17     80                0         80 0.8208071 0.05334706 0.7082209 0.8905636
#> 18     80               20        100 0.8208071 0.05334706 0.7082209 0.8905636
#> 19     80               40        120 0.8208071 0.05334706 0.7082209 0.8905636
#> 20     80               60        140 0.8208071 0.05334706 0.7082209 0.8905636
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.12745317 0.178015079 0.6117568
#> 2  0.30574618 0.12277475 0.144318920 0.5474456
#> 3  0.26001915 0.11764054 0.116346476 0.4901497
#> 4  0.22113100 0.11212039 0.092733335 0.4393217
#> 5  0.25589195 0.11370398 0.105106364 0.4816695
#> 6  0.21762106 0.10764115 0.083886961 0.4318107
#> 7  0.18507391 0.10156301 0.066398470 0.3876918
#> 8  0.15739448 0.09549748 0.049631317 0.3486823
#> 9  0.18213629 0.09910613 0.059419466 0.3811775
#> 10 0.15489621 0.09281143 0.046352231 0.3429225
#> 11 0.13173012 0.08669334 0.035730205 0.3090878
#> 12 0.11202872 0.08076232 0.024283459 0.2791301
#> 13 0.12963921 0.08514100 0.031542910 0.3040893
#> 14 0.11025053 0.07913726 0.023811209 0.2747002
#> 15 0.09376159 0.07338909 0.016892640 0.2486241
#> 16 0.07973872 0.06789859 0.010347607 0.2254301
#> 17 0.09227334 0.07239661 0.015294716 0.2447622
#> 18 0.07847305 0.06690615 0.011009701 0.2219892
#> 19 0.06673672 0.06169712 0.006610065 0.2016600
#> 20 0.05675565 0.05676878 0.003553233 0.1834451
```
