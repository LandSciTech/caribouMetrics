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
#> 1         0.1 0.3838513 0.02388632 0.3460071 0.4238489
#> 2         0.2 0.3832760 0.02381575 0.3454955 0.4230989
#> 3         0.3 0.3827015 0.02374545 0.3449846 0.4223502
#> 4         0.4 0.3821279 0.02367541 0.3444745 0.4216029
#> 5         0.5 0.3815551 0.02360564 0.3439652 0.4208569
#> 6         0.6 0.3809832 0.02353614 0.3434566 0.4201122
#> 7         0.7 0.3804122 0.02346689 0.3429487 0.4193688
#> 8         0.8 0.3798420 0.02339792 0.3424416 0.4186267
#> 9         0.9 0.3792726 0.02332921 0.3419353 0.4178860
#> 10        1.0 0.3787041 0.02326076 0.3414297 0.4171466

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3694363
#> 2       1        0.1        V5 0.3423097
#> 3       1        0.1        V9 0.3734287
#> 4       1        0.1        V3 0.3799546
#> 5       1        0.1        V4 0.3944225
#> 6       1        0.1        V8 0.4323920
#> 7       1        0.1        V2 0.3587425
#> 8       1        0.1        V6 0.3868033
#> 9       1        0.1        V7 0.3843607
#> 10      1        0.1       V10 0.3702089
#> 11      2        0.2        V7 0.3838301
#> 12      2        0.2        V8 0.4316389
#> 13      2        0.2        V9 0.3728791
#> 14      2        0.2       V10 0.3695743
#> 15      2        0.2        V1 0.3688546
#> 16      2        0.2        V5 0.3418065
#> 17      2        0.2        V2 0.3582020
#> 18      2        0.2        V3 0.3793599
#> 19      2        0.2        V4 0.3936834
#> 20      2        0.2        V6 0.3861593
#> 21      3        0.3        V4 0.3929456
#> 22      3        0.3        V5 0.3413040
#> 23      3        0.3        V3 0.3787661
#> 24      3        0.3        V7 0.3833002
#> 25      3        0.3        V8 0.4308870
#> 26      3        0.3        V9 0.3723303
#> 27      3        0.3        V6 0.3855164
#> 28      3        0.3       V10 0.3689408
#> 29      3        0.3        V1 0.3682738
#> 30      3        0.3        V2 0.3576623
#> 31      4        0.4        V1 0.3676939
#> 32      4        0.4        V9 0.3717823
#> 33      4        0.4        V3 0.3781732
#> 34      4        0.4        V4 0.3922092
#> 35      4        0.4        V5 0.3408023
#> 36      4        0.4        V2 0.3571233
#> 37      4        0.4        V6 0.3848746
#> 38      4        0.4        V7 0.3827710
#> 39      4        0.4        V8 0.4301365
#> 40      4        0.4       V10 0.3683084
#> 41      5        0.5        V8 0.4293873
#> 42      5        0.5        V9 0.3712351
#> 43      5        0.5       V10 0.3676771
#> 44      5        0.5        V1 0.3671149
#> 45      5        0.5        V5 0.3403013
#> 46      5        0.5        V2 0.3565852
#> 47      5        0.5        V3 0.3775812
#> 48      5        0.5        V4 0.3914742
#> 49      5        0.5        V6 0.3842338
#> 50      5        0.5        V7 0.3822425
#> 51      6        0.6        V4 0.3907405
#> 52      6        0.6        V5 0.3398010
#> 53      6        0.6        V7 0.3817148
#> 54      6        0.6        V8 0.4286394
#> 55      6        0.6        V9 0.3706887
#> 56      6        0.6        V6 0.3835942
#> 57      6        0.6       V10 0.3670469
#> 58      6        0.6        V1 0.3665369
#> 59      6        0.6        V2 0.3560479
#> 60      6        0.6        V3 0.3769902
#> 61      7        0.7        V1 0.3659597
#> 62      7        0.7        V3 0.3764001
#> 63      7        0.7        V4 0.3900083
#> 64      7        0.7        V5 0.3393015
#> 65      7        0.7        V9 0.3701431
#> 66      7        0.7        V6 0.3829555
#> 67      7        0.7        V7 0.3811878
#> 68      7        0.7        V8 0.4278928
#> 69      7        0.7        V2 0.3555115
#> 70      7        0.7       V10 0.3664177
#> 71      8        0.8        V9 0.3695983
#> 72      8        0.8       V10 0.3657896
#> 73      8        0.8        V1 0.3653835
#> 74      8        0.8        V5 0.3388027
#> 75      8        0.8        V2 0.3549758
#> 76      8        0.8        V3 0.3758109
#> 77      8        0.8        V4 0.3892774
#> 78      8        0.8        V8 0.4271475
#> 79      8        0.8        V6 0.3823180
#> 80      8        0.8        V7 0.3806615
#> 81      9        0.9        V5 0.3383046
#> 82      9        0.9        V7 0.3801360
#> 83      9        0.9        V8 0.4264035
#> 84      9        0.9        V9 0.3690544
#> 85      9        0.9        V6 0.3816815
#> 86      9        0.9       V10 0.3651626
#> 87      9        0.9        V1 0.3648082
#> 88      9        0.9        V2 0.3544409
#> 89      9        0.9        V3 0.3752227
#> 90      9        0.9        V4 0.3885478
#> 91     10        1.0        V1 0.3642338
#> 92     10        1.0        V4 0.3878197
#> 93     10        1.0        V5 0.3378073
#> 94     10        1.0        V9 0.3685112
#> 95     10        1.0        V3 0.3746353
#> 96     10        1.0        V7 0.3796112
#> 97     10        1.0        V8 0.4256608
#> 98     10        1.0        V2 0.3539069
#> 99     10        1.0        V6 0.3810460
#> 100    10        1.0       V10 0.3645367

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3694363
#> 2       1        0.1        V5 0.3423097
#> 3       1        0.1        V9 0.3734287
#> 4       1        0.1        V3 0.3799546
#> 5       1        0.1        V4 0.3944225
#> 6       1        0.1        V8 0.4323920
#> 7       1        0.1        V2 0.3587425
#> 8       1        0.1        V6 0.3868033
#> 9       1        0.1        V7 0.3843607
#> 10      1        0.1       V10 0.3702089
#> 11      2        0.2        V7 0.3838301
#> 12      2        0.2        V8 0.4316389
#> 13      2        0.2        V9 0.3728791
#> 14      2        0.2       V10 0.3695743
#> 15      2        0.2        V1 0.3688546
#> 16      2        0.2        V5 0.3418065
#> 17      2        0.2        V2 0.3582020
#> 18      2        0.2        V3 0.3793599
#> 19      2        0.2        V4 0.3936834
#> 20      2        0.2        V6 0.3861593
#> 21      3        0.3        V4 0.3929456
#> 22      3        0.3        V5 0.3413040
#> 23      3        0.3        V3 0.3787661
#> 24      3        0.3        V7 0.3833002
#> 25      3        0.3        V8 0.4308870
#> 26      3        0.3        V9 0.3723303
#> 27      3        0.3        V6 0.3855164
#> 28      3        0.3       V10 0.3689408
#> 29      3        0.3        V1 0.3682738
#> 30      3        0.3        V2 0.3576623
#> 31      4        0.4        V1 0.3676939
#> 32      4        0.4        V9 0.3717823
#> 33      4        0.4        V3 0.3781732
#> 34      4        0.4        V4 0.3922092
#> 35      4        0.4        V5 0.3408023
#> 36      4        0.4        V2 0.3571233
#> 37      4        0.4        V6 0.3848746
#> 38      4        0.4        V7 0.3827710
#> 39      4        0.4        V8 0.4301365
#> 40      4        0.4       V10 0.3683084
#> 41      5        0.5        V8 0.4293873
#> 42      5        0.5        V9 0.3712351
#> 43      5        0.5       V10 0.3676771
#> 44      5        0.5        V1 0.3671149
#> 45      5        0.5        V5 0.3403013
#> 46      5        0.5        V2 0.3565852
#> 47      5        0.5        V3 0.3775812
#> 48      5        0.5        V4 0.3914742
#> 49      5        0.5        V6 0.3842338
#> 50      5        0.5        V7 0.3822425
#> 51      6        0.6        V4 0.3907405
#> 52      6        0.6        V5 0.3398010
#> 53      6        0.6        V7 0.3817148
#> 54      6        0.6        V8 0.4286394
#> 55      6        0.6        V9 0.3706887
#> 56      6        0.6        V6 0.3835942
#> 57      6        0.6       V10 0.3670469
#> 58      6        0.6        V1 0.3665369
#> 59      6        0.6        V2 0.3560479
#> 60      6        0.6        V3 0.3769902
#> 61      7        0.7        V1 0.3659597
#> 62      7        0.7        V3 0.3764001
#> 63      7        0.7        V4 0.3900083
#> 64      7        0.7        V5 0.3393015
#> 65      7        0.7        V9 0.3701431
#> 66      7        0.7        V6 0.3829555
#> 67      7        0.7        V7 0.3811878
#> 68      7        0.7        V8 0.4278928
#> 69      7        0.7        V2 0.3555115
#> 70      7        0.7       V10 0.3664177
#> 71      8        0.8        V9 0.3695983
#> 72      8        0.8       V10 0.3657896
#> 73      8        0.8        V1 0.3653835
#> 74      8        0.8        V5 0.3388027
#> 75      8        0.8        V2 0.3549758
#> 76      8        0.8        V3 0.3758109
#> 77      8        0.8        V4 0.3892774
#> 78      8        0.8        V8 0.4271475
#> 79      8        0.8        V6 0.3823180
#> 80      8        0.8        V7 0.3806615
#> 81      9        0.9        V5 0.3383046
#> 82      9        0.9        V7 0.3801360
#> 83      9        0.9        V8 0.4264035
#> 84      9        0.9        V9 0.3690544
#> 85      9        0.9        V6 0.3816815
#> 86      9        0.9       V10 0.3651626
#> 87      9        0.9        V1 0.3648082
#> 88      9        0.9        V2 0.3544409
#> 89      9        0.9        V3 0.3752227
#> 90      9        0.9        V4 0.3885478
#> 91     10        1.0        V1 0.3642338
#> 92     10        1.0        V4 0.3878197
#> 93     10        1.0        V5 0.3378073
#> 94     10        1.0        V9 0.3685112
#> 95     10        1.0        V3 0.3746353
#> 96     10        1.0        V7 0.3796112
#> 97     10        1.0        V8 0.4256608
#> 98     10        1.0        V2 0.3539069
#> 99     10        1.0        V6 0.3810460
#> 100    10        1.0       V10 0.3645367


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
#> 1       0                0          0 0.8757906 0.04948040 0.7728714 0.9337690
#> 2       0               20         20 0.8757906 0.04948040 0.7728714 0.9337690
#> 3       0               40         40 0.8757906 0.04948040 0.7728714 0.9337690
#> 4       0               60         60 0.8757906 0.04948040 0.7728714 0.9337690
#> 5      20                0         20 0.8617131 0.05070310 0.7580835 0.9218309
#> 6      20               20         40 0.8617131 0.05070310 0.7580835 0.9218309
#> 7      20               40         60 0.8617131 0.05070310 0.7580835 0.9218309
#> 8      20               60         80 0.8617131 0.05070310 0.7580835 0.9218309
#> 9      40                0         40 0.8478591 0.05184847 0.7437187 0.9097907
#> 10     40               20         60 0.8478591 0.05184847 0.7437187 0.9097907
#> 11     40               40         80 0.8478591 0.05184847 0.7437187 0.9097907
#> 12     40               60        100 0.8478591 0.05184847 0.7437187 0.9097907
#> 13     60                0         60 0.8342249 0.05292939 0.7297437 0.8976985
#> 14     60               20         80 0.8342249 0.05292939 0.7297437 0.8976985
#> 15     60               40        100 0.8342249 0.05292939 0.7297437 0.8976985
#> 16     60               60        120 0.8342249 0.05292939 0.7297437 0.8976985
#> 17     80                0         80 0.8208071 0.05395511 0.7161308 0.8855920
#> 18     80               20        100 0.8208071 0.05395511 0.7161308 0.8855920
#> 19     80               40        120 0.8208071 0.05395511 0.7161308 0.8855920
#> 20     80               60        140 0.8208071 0.05395511 0.7161308 0.8855920
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11360342 0.174056810 0.5632972
#> 2  0.30574618 0.10753737 0.133960494 0.5091432
#> 3  0.26001915 0.10228628 0.102125788 0.4605270
#> 4  0.22113100 0.09735173 0.076956434 0.4169897
#> 5  0.25589195 0.10277325 0.095840031 0.4446648
#> 6  0.21762106 0.09593493 0.072004956 0.4027993
#> 7  0.18507391 0.08998755 0.053312196 0.3653665
#> 8  0.15739448 0.08457085 0.038786397 0.3319054
#> 9  0.18213629 0.09034070 0.049658516 0.3531716
#> 10 0.15489621 0.08370812 0.035967276 0.3210025
#> 11 0.13173012 0.07790167 0.025487444 0.2922288
#> 12 0.11202872 0.07266223 0.017593935 0.2664652
#> 13 0.12963921 0.07806258 0.023476865 0.2828458
#> 14 0.11025053 0.07203265 0.016098012 0.2580557
#> 15 0.09376159 0.06672379 0.010681825 0.2358118
#> 16 0.07973872 0.06194788 0.006813222 0.2158111
#> 17 0.09227334 0.06669381 0.009675400 0.2285383
#> 18 0.07847305 0.06139606 0.006109211 0.2092605
#> 19 0.06673672 0.05670924 0.003670349 0.1918671
#> 20 0.05675565 0.05249439 0.002076838 0.1761269
```
