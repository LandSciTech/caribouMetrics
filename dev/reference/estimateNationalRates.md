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
#> 1         0.1 0.3838513 0.02596350 0.3384268 0.4225437
#> 2         0.2 0.3832760 0.02592779 0.3379286 0.4219440
#> 3         0.3 0.3827015 0.02589221 0.3374312 0.4213452
#> 4         0.4 0.3821279 0.02585679 0.3369344 0.4207472
#> 5         0.5 0.3815551 0.02582151 0.3364385 0.4201501
#> 6         0.6 0.3809832 0.02578637 0.3359432 0.4195538
#> 7         0.7 0.3804122 0.02575138 0.3354487 0.4189584
#> 8         0.8 0.3798420 0.02571653 0.3349549 0.4183639
#> 9         0.9 0.3792726 0.02568183 0.3344618 0.4177702
#> 10        1.0 0.3787041 0.02564727 0.3339695 0.4171773

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3725111
#> 2       1        0.1        V5 0.3961503
#> 3       1        0.1        V9 0.3842881
#> 4       1        0.1        V3 0.3730025
#> 5       1        0.1        V4 0.4059921
#> 6       1        0.1        V8 0.3954571
#> 7       1        0.1        V2 0.3744026
#> 8       1        0.1        V6 0.4273490
#> 9       1        0.1        V7 0.3285313
#> 10      1        0.1       V10 0.3862033
#> 11      2        0.2        V7 0.3280375
#> 12      2        0.2        V8 0.3948840
#> 13      2        0.2        V9 0.3837195
#> 14      2        0.2       V10 0.3856841
#> 15      2        0.2        V1 0.3719979
#> 16      2        0.2        V5 0.3956054
#> 17      2        0.2        V2 0.3738364
#> 18      2        0.2        V3 0.3725098
#> 19      2        0.2        V4 0.4052890
#> 20      2        0.2        V6 0.4267793
#> 21      3        0.3        V4 0.4045870
#> 22      3        0.3        V5 0.3950612
#> 23      3        0.3        V3 0.3720178
#> 24      3        0.3        V7 0.3275445
#> 25      3        0.3        V8 0.3943119
#> 26      3        0.3        V9 0.3831517
#> 27      3        0.3        V6 0.4262104
#> 28      3        0.3       V10 0.3851657
#> 29      3        0.3        V1 0.3714853
#> 30      3        0.3        V2 0.3732712
#> 31      4        0.4        V1 0.3709735
#> 32      4        0.4        V9 0.3825848
#> 33      4        0.4        V3 0.3715265
#> 34      4        0.4        V4 0.4038864
#> 35      4        0.4        V5 0.3945178
#> 36      4        0.4        V2 0.3727068
#> 37      4        0.4        V6 0.4256423
#> 38      4        0.4        V7 0.3270521
#> 39      4        0.4        V8 0.3937405
#> 40      4        0.4       V10 0.3846480
#> 41      5        0.5        V8 0.3931700
#> 42      5        0.5        V9 0.3820187
#> 43      5        0.5       V10 0.3841310
#> 44      5        0.5        V1 0.3704624
#> 45      5        0.5        V5 0.3939751
#> 46      5        0.5        V2 0.3721432
#> 47      5        0.5        V3 0.3710358
#> 48      5        0.5        V4 0.4031869
#> 49      5        0.5        V6 0.4250749
#> 50      5        0.5        V7 0.3265606
#> 51      6        0.6        V4 0.4024886
#> 52      6        0.6        V5 0.3934331
#> 53      6        0.6        V7 0.3260697
#> 54      6        0.6        V8 0.3926003
#> 55      6        0.6        V9 0.3814535
#> 56      6        0.6        V6 0.4245082
#> 57      6        0.6       V10 0.3836146
#> 58      6        0.6        V1 0.3699519
#> 59      6        0.6        V2 0.3715805
#> 60      6        0.6        V3 0.3705457
#> 61      7        0.7        V1 0.3694422
#> 62      7        0.7        V3 0.3700563
#> 63      7        0.7        V4 0.4017916
#> 64      7        0.7        V5 0.3928920
#> 65      7        0.7        V9 0.3808891
#> 66      7        0.7        V6 0.4239423
#> 67      7        0.7        V7 0.3255796
#> 68      7        0.7        V8 0.3920314
#> 69      7        0.7        V2 0.3710186
#> 70      7        0.7       V10 0.3830990
#> 71      8        0.8        V9 0.3803255
#> 72      8        0.8       V10 0.3825840
#> 73      8        0.8        V1 0.3689332
#> 74      8        0.8        V5 0.3923515
#> 75      8        0.8        V2 0.3704576
#> 76      8        0.8        V3 0.3695675
#> 77      8        0.8        V4 0.4010957
#> 78      8        0.8        V8 0.3914634
#> 79      8        0.8        V6 0.4233772
#> 80      8        0.8        V7 0.3250902
#> 81      9        0.9        V5 0.3918118
#> 82      9        0.9        V7 0.3246016
#> 83      9        0.9        V8 0.3908962
#> 84      9        0.9        V9 0.3797627
#> 85      9        0.9        V6 0.4228128
#> 86      9        0.9       V10 0.3820698
#> 87      9        0.9        V1 0.3684248
#> 88      9        0.9        V2 0.3698975
#> 89      9        0.9        V3 0.3690794
#> 90      9        0.9        V4 0.4004011
#> 91     10        1.0        V1 0.3679172
#> 92     10        1.0        V4 0.3997076
#> 93     10        1.0        V5 0.3912728
#> 94     10        1.0        V9 0.3792008
#> 95     10        1.0        V3 0.3685920
#> 96     10        1.0        V7 0.3241137
#> 97     10        1.0        V8 0.3903298
#> 98     10        1.0        V2 0.3693382
#> 99     10        1.0        V6 0.4222492
#> 100    10        1.0       V10 0.3815562

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3725111
#> 2       1        0.1        V5 0.3961503
#> 3       1        0.1        V9 0.3842881
#> 4       1        0.1        V3 0.3730025
#> 5       1        0.1        V4 0.4059921
#> 6       1        0.1        V8 0.3954571
#> 7       1        0.1        V2 0.3744026
#> 8       1        0.1        V6 0.4273490
#> 9       1        0.1        V7 0.3285313
#> 10      1        0.1       V10 0.3862033
#> 11      2        0.2        V7 0.3280375
#> 12      2        0.2        V8 0.3948840
#> 13      2        0.2        V9 0.3837195
#> 14      2        0.2       V10 0.3856841
#> 15      2        0.2        V1 0.3719979
#> 16      2        0.2        V5 0.3956054
#> 17      2        0.2        V2 0.3738364
#> 18      2        0.2        V3 0.3725098
#> 19      2        0.2        V4 0.4052890
#> 20      2        0.2        V6 0.4267793
#> 21      3        0.3        V4 0.4045870
#> 22      3        0.3        V5 0.3950612
#> 23      3        0.3        V3 0.3720178
#> 24      3        0.3        V7 0.3275445
#> 25      3        0.3        V8 0.3943119
#> 26      3        0.3        V9 0.3831517
#> 27      3        0.3        V6 0.4262104
#> 28      3        0.3       V10 0.3851657
#> 29      3        0.3        V1 0.3714853
#> 30      3        0.3        V2 0.3732712
#> 31      4        0.4        V1 0.3709735
#> 32      4        0.4        V9 0.3825848
#> 33      4        0.4        V3 0.3715265
#> 34      4        0.4        V4 0.4038864
#> 35      4        0.4        V5 0.3945178
#> 36      4        0.4        V2 0.3727068
#> 37      4        0.4        V6 0.4256423
#> 38      4        0.4        V7 0.3270521
#> 39      4        0.4        V8 0.3937405
#> 40      4        0.4       V10 0.3846480
#> 41      5        0.5        V8 0.3931700
#> 42      5        0.5        V9 0.3820187
#> 43      5        0.5       V10 0.3841310
#> 44      5        0.5        V1 0.3704624
#> 45      5        0.5        V5 0.3939751
#> 46      5        0.5        V2 0.3721432
#> 47      5        0.5        V3 0.3710358
#> 48      5        0.5        V4 0.4031869
#> 49      5        0.5        V6 0.4250749
#> 50      5        0.5        V7 0.3265606
#> 51      6        0.6        V4 0.4024886
#> 52      6        0.6        V5 0.3934331
#> 53      6        0.6        V7 0.3260697
#> 54      6        0.6        V8 0.3926003
#> 55      6        0.6        V9 0.3814535
#> 56      6        0.6        V6 0.4245082
#> 57      6        0.6       V10 0.3836146
#> 58      6        0.6        V1 0.3699519
#> 59      6        0.6        V2 0.3715805
#> 60      6        0.6        V3 0.3705457
#> 61      7        0.7        V1 0.3694422
#> 62      7        0.7        V3 0.3700563
#> 63      7        0.7        V4 0.4017916
#> 64      7        0.7        V5 0.3928920
#> 65      7        0.7        V9 0.3808891
#> 66      7        0.7        V6 0.4239423
#> 67      7        0.7        V7 0.3255796
#> 68      7        0.7        V8 0.3920314
#> 69      7        0.7        V2 0.3710186
#> 70      7        0.7       V10 0.3830990
#> 71      8        0.8        V9 0.3803255
#> 72      8        0.8       V10 0.3825840
#> 73      8        0.8        V1 0.3689332
#> 74      8        0.8        V5 0.3923515
#> 75      8        0.8        V2 0.3704576
#> 76      8        0.8        V3 0.3695675
#> 77      8        0.8        V4 0.4010957
#> 78      8        0.8        V8 0.3914634
#> 79      8        0.8        V6 0.4233772
#> 80      8        0.8        V7 0.3250902
#> 81      9        0.9        V5 0.3918118
#> 82      9        0.9        V7 0.3246016
#> 83      9        0.9        V8 0.3908962
#> 84      9        0.9        V9 0.3797627
#> 85      9        0.9        V6 0.4228128
#> 86      9        0.9       V10 0.3820698
#> 87      9        0.9        V1 0.3684248
#> 88      9        0.9        V2 0.3698975
#> 89      9        0.9        V3 0.3690794
#> 90      9        0.9        V4 0.4004011
#> 91     10        1.0        V1 0.3679172
#> 92     10        1.0        V4 0.3997076
#> 93     10        1.0        V5 0.3912728
#> 94     10        1.0        V9 0.3792008
#> 95     10        1.0        V3 0.3685920
#> 96     10        1.0        V7 0.3241137
#> 97     10        1.0        V8 0.3903298
#> 98     10        1.0        V2 0.3693382
#> 99     10        1.0        V6 0.4222492
#> 100    10        1.0       V10 0.3815562


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
#> 1       0                0          0 0.8757906 0.04639415 0.7717261 0.9325841
#> 2       0               20         20 0.8757906 0.04639415 0.7717261 0.9325841
#> 3       0               40         40 0.8757906 0.04639415 0.7717261 0.9325841
#> 4       0               60         60 0.8757906 0.04639415 0.7717261 0.9325841
#> 5      20                0         20 0.8617131 0.04769377 0.7546550 0.9230695
#> 6      20               20         40 0.8617131 0.04769377 0.7546550 0.9230695
#> 7      20               40         60 0.8617131 0.04769377 0.7546550 0.9230695
#> 8      20               60         80 0.8617131 0.04769377 0.7546550 0.9230695
#> 9      40                0         40 0.8478591 0.04893172 0.7381115 0.9134902
#> 10     40               20         60 0.8478591 0.04893172 0.7381115 0.9134902
#> 11     40               40         80 0.8478591 0.04893172 0.7381115 0.9134902
#> 12     40               60        100 0.8478591 0.04893172 0.7381115 0.9134902
#> 13     60                0         60 0.8342249 0.05012325 0.7220537 0.9038724
#> 14     60               20         80 0.8342249 0.05012325 0.7220537 0.9038724
#> 15     60               40        100 0.8342249 0.05012325 0.7220537 0.9038724
#> 16     60               60        120 0.8342249 0.05012325 0.7220537 0.9038724
#> 17     80                0         80 0.8208071 0.05127765 0.7064473 0.8942369
#> 18     80               20        100 0.8208071 0.05127765 0.7064473 0.8942369
#> 19     80               40        120 0.8208071 0.05127765 0.7064473 0.8942369
#> 20     80               60        140 0.8208071 0.05127765 0.7064473 0.8942369
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11509341 0.173382170 0.5521989
#> 2  0.30574618 0.10687885 0.143425258 0.4937503
#> 3  0.26001915 0.09923436 0.118081708 0.4418874
#> 4  0.22113100 0.09212743 0.096679465 0.3959826
#> 5  0.25589195 0.10195183 0.097040335 0.4362015
#> 6  0.21762106 0.09416724 0.078956138 0.3909544
#> 7  0.18507391 0.08703821 0.063775150 0.3509595
#> 8  0.15739448 0.08049272 0.051085760 0.3156174
#> 9  0.18213629 0.08869860 0.051298428 0.3465802
#> 10 0.15489621 0.08170248 0.040710753 0.3117471
#> 11 0.13173012 0.07532868 0.031961366 0.2809524
#> 12 0.11202872 0.06950146 0.024785368 0.2537007
#> 13 0.12963921 0.07620165 0.024904337 0.2775783
#> 14 0.11025053 0.07009104 0.019048624 0.2507125
#> 15 0.09376159 0.06452787 0.014338372 0.2268975
#> 16 0.07973872 0.05944465 0.010597062 0.2057444
#> 17 0.09227334 0.06487772 0.010657963 0.2242831
#> 18 0.07847305 0.05962607 0.007716805 0.2034191
#> 19 0.06673672 0.05483833 0.005455103 0.1848358
#> 20 0.05675565 0.05045841 0.003751940 0.1682347
```
