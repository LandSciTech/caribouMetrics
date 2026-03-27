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
#> 1         0.1 0.3838513 0.01397921 0.3535634 0.3981293
#> 2         0.2 0.3832760 0.01396504 0.3530677 0.3975603
#> 3         0.3 0.3827015 0.01395110 0.3525726 0.3969920
#> 4         0.4 0.3821279 0.01393738 0.3520783 0.3964245
#> 5         0.5 0.3815551 0.01392388 0.3515847 0.3958579
#> 6         0.6 0.3809832 0.01391060 0.3510917 0.3952921
#> 7         0.7 0.3804122 0.01389755 0.3505995 0.3947271
#> 8         0.8 0.3798420 0.01388471 0.3501079 0.3941629
#> 9         0.9 0.3792726 0.01387209 0.3496171 0.3935995
#> 10        1.0 0.3787041 0.01385969 0.3491269 0.3930369

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3631315
#> 2       1        0.1        V5 0.3635580
#> 3       1        0.1        V9 0.3509073
#> 4       1        0.1        V3 0.3726477
#> 5       1        0.1        V4 0.3845795
#> 6       1        0.1        V8 0.3692476
#> 7       1        0.1        V2 0.3680885
#> 8       1        0.1        V6 0.3733421
#> 9       1        0.1        V7 0.4020631
#> 10      1        0.1       V10 0.3627122
#> 11      2        0.2        V7 0.4014975
#> 12      2        0.2        V8 0.3686715
#> 13      2        0.2        V9 0.3504380
#> 14      2        0.2       V10 0.3621253
#> 15      2        0.2        V1 0.3626524
#> 16      2        0.2        V5 0.3629137
#> 17      2        0.2        V2 0.3675178
#> 18      2        0.2        V3 0.3721631
#> 19      2        0.2        V4 0.3839986
#> 20      2        0.2        V6 0.3728216
#> 21      3        0.3        V4 0.3834185
#> 22      3        0.3        V5 0.3622705
#> 23      3        0.3        V3 0.3716790
#> 24      3        0.3        V7 0.4009327
#> 25      3        0.3        V8 0.3680963
#> 26      3        0.3        V9 0.3499694
#> 27      3        0.3        V6 0.3723019
#> 28      3        0.3       V10 0.3615394
#> 29      3        0.3        V1 0.3621740
#> 30      3        0.3        V2 0.3669479
#> 31      4        0.4        V1 0.3616963
#> 32      4        0.4        V9 0.3495014
#> 33      4        0.4        V3 0.3711956
#> 34      4        0.4        V4 0.3828394
#> 35      4        0.4        V5 0.3616285
#> 36      4        0.4        V2 0.3663790
#> 37      4        0.4        V6 0.3717828
#> 38      4        0.4        V7 0.4003686
#> 39      4        0.4        V8 0.3675220
#> 40      4        0.4       V10 0.3609544
#> 41      5        0.5        V8 0.3669486
#> 42      5        0.5        V9 0.3490340
#> 43      5        0.5       V10 0.3603704
#> 44      5        0.5        V1 0.3612192
#> 45      5        0.5        V5 0.3609875
#> 46      5        0.5        V2 0.3658109
#> 47      5        0.5        V3 0.3707129
#> 48      5        0.5        V4 0.3822610
#> 49      5        0.5        V6 0.3712645
#> 50      5        0.5        V7 0.3998054
#> 51      6        0.6        V4 0.3816836
#> 52      6        0.6        V5 0.3603478
#> 53      6        0.6        V7 0.3992429
#> 54      6        0.6        V8 0.3663761
#> 55      6        0.6        V9 0.3485672
#> 56      6        0.6        V6 0.3707469
#> 57      6        0.6       V10 0.3597873
#> 58      6        0.6        V1 0.3607427
#> 59      6        0.6        V2 0.3652438
#> 60      6        0.6        V3 0.3702307
#> 61      7        0.7        V1 0.3602668
#> 62      7        0.7        V3 0.3697492
#> 63      7        0.7        V4 0.3811070
#> 64      7        0.7        V5 0.3597091
#> 65      7        0.7        V9 0.3481011
#> 66      7        0.7        V6 0.3702300
#> 67      7        0.7        V7 0.3986813
#> 68      7        0.7        V8 0.3658045
#> 69      7        0.7        V2 0.3646775
#> 70      7        0.7       V10 0.3592052
#> 71      8        0.8        V9 0.3476355
#> 72      8        0.8       V10 0.3586240
#> 73      8        0.8        V1 0.3597916
#> 74      8        0.8        V5 0.3590716
#> 75      8        0.8        V2 0.3641120
#> 76      8        0.8        V3 0.3692684
#> 77      8        0.8        V4 0.3805314
#> 78      8        0.8        V8 0.3652338
#> 79      8        0.8        V6 0.3697139
#> 80      8        0.8        V7 0.3981204
#> 81      9        0.9        V5 0.3584352
#> 82      9        0.9        V7 0.3975603
#> 83      9        0.9        V8 0.3646640
#> 84      9        0.9        V9 0.3471706
#> 85      9        0.9        V6 0.3691985
#> 86      9        0.9       V10 0.3580437
#> 87      9        0.9        V1 0.3593169
#> 88      9        0.9        V2 0.3635475
#> 89      9        0.9        V3 0.3687881
#> 90      9        0.9        V4 0.3799565
#> 91     10        1.0        V1 0.3588430
#> 92     10        1.0        V4 0.3793826
#> 93     10        1.0        V5 0.3578000
#> 94     10        1.0        V9 0.3467064
#> 95     10        1.0        V3 0.3683085
#> 96     10        1.0        V7 0.3970010
#> 97     10        1.0        V8 0.3640950
#> 98     10        1.0        V2 0.3629838
#> 99     10        1.0        V6 0.3686837
#> 100    10        1.0       V10 0.3574644

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3631315
#> 2       1        0.1        V5 0.3635580
#> 3       1        0.1        V9 0.3509073
#> 4       1        0.1        V3 0.3726477
#> 5       1        0.1        V4 0.3845795
#> 6       1        0.1        V8 0.3692476
#> 7       1        0.1        V2 0.3680885
#> 8       1        0.1        V6 0.3733421
#> 9       1        0.1        V7 0.4020631
#> 10      1        0.1       V10 0.3627122
#> 11      2        0.2        V7 0.4014975
#> 12      2        0.2        V8 0.3686715
#> 13      2        0.2        V9 0.3504380
#> 14      2        0.2       V10 0.3621253
#> 15      2        0.2        V1 0.3626524
#> 16      2        0.2        V5 0.3629137
#> 17      2        0.2        V2 0.3675178
#> 18      2        0.2        V3 0.3721631
#> 19      2        0.2        V4 0.3839986
#> 20      2        0.2        V6 0.3728216
#> 21      3        0.3        V4 0.3834185
#> 22      3        0.3        V5 0.3622705
#> 23      3        0.3        V3 0.3716790
#> 24      3        0.3        V7 0.4009327
#> 25      3        0.3        V8 0.3680963
#> 26      3        0.3        V9 0.3499694
#> 27      3        0.3        V6 0.3723019
#> 28      3        0.3       V10 0.3615394
#> 29      3        0.3        V1 0.3621740
#> 30      3        0.3        V2 0.3669479
#> 31      4        0.4        V1 0.3616963
#> 32      4        0.4        V9 0.3495014
#> 33      4        0.4        V3 0.3711956
#> 34      4        0.4        V4 0.3828394
#> 35      4        0.4        V5 0.3616285
#> 36      4        0.4        V2 0.3663790
#> 37      4        0.4        V6 0.3717828
#> 38      4        0.4        V7 0.4003686
#> 39      4        0.4        V8 0.3675220
#> 40      4        0.4       V10 0.3609544
#> 41      5        0.5        V8 0.3669486
#> 42      5        0.5        V9 0.3490340
#> 43      5        0.5       V10 0.3603704
#> 44      5        0.5        V1 0.3612192
#> 45      5        0.5        V5 0.3609875
#> 46      5        0.5        V2 0.3658109
#> 47      5        0.5        V3 0.3707129
#> 48      5        0.5        V4 0.3822610
#> 49      5        0.5        V6 0.3712645
#> 50      5        0.5        V7 0.3998054
#> 51      6        0.6        V4 0.3816836
#> 52      6        0.6        V5 0.3603478
#> 53      6        0.6        V7 0.3992429
#> 54      6        0.6        V8 0.3663761
#> 55      6        0.6        V9 0.3485672
#> 56      6        0.6        V6 0.3707469
#> 57      6        0.6       V10 0.3597873
#> 58      6        0.6        V1 0.3607427
#> 59      6        0.6        V2 0.3652438
#> 60      6        0.6        V3 0.3702307
#> 61      7        0.7        V1 0.3602668
#> 62      7        0.7        V3 0.3697492
#> 63      7        0.7        V4 0.3811070
#> 64      7        0.7        V5 0.3597091
#> 65      7        0.7        V9 0.3481011
#> 66      7        0.7        V6 0.3702300
#> 67      7        0.7        V7 0.3986813
#> 68      7        0.7        V8 0.3658045
#> 69      7        0.7        V2 0.3646775
#> 70      7        0.7       V10 0.3592052
#> 71      8        0.8        V9 0.3476355
#> 72      8        0.8       V10 0.3586240
#> 73      8        0.8        V1 0.3597916
#> 74      8        0.8        V5 0.3590716
#> 75      8        0.8        V2 0.3641120
#> 76      8        0.8        V3 0.3692684
#> 77      8        0.8        V4 0.3805314
#> 78      8        0.8        V8 0.3652338
#> 79      8        0.8        V6 0.3697139
#> 80      8        0.8        V7 0.3981204
#> 81      9        0.9        V5 0.3584352
#> 82      9        0.9        V7 0.3975603
#> 83      9        0.9        V8 0.3646640
#> 84      9        0.9        V9 0.3471706
#> 85      9        0.9        V6 0.3691985
#> 86      9        0.9       V10 0.3580437
#> 87      9        0.9        V1 0.3593169
#> 88      9        0.9        V2 0.3635475
#> 89      9        0.9        V3 0.3687881
#> 90      9        0.9        V4 0.3799565
#> 91     10        1.0        V1 0.3588430
#> 92     10        1.0        V4 0.3793826
#> 93     10        1.0        V5 0.3578000
#> 94     10        1.0        V9 0.3467064
#> 95     10        1.0        V3 0.3683085
#> 96     10        1.0        V7 0.3970010
#> 97     10        1.0        V8 0.3640950
#> 98     10        1.0        V2 0.3629838
#> 99     10        1.0        V6 0.3686837
#> 100    10        1.0       V10 0.3574644


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
#> 1       0                0          0 0.8757906 0.04575959 0.7953619 0.9415010
#> 2       0               20         20 0.8757906 0.04575959 0.7953619 0.9415010
#> 3       0               40         40 0.8757906 0.04575959 0.7953619 0.9415010
#> 4       0               60         60 0.8757906 0.04575959 0.7953619 0.9415010
#> 5      20                0         20 0.8617131 0.04856489 0.7764063 0.9330811
#> 6      20               20         40 0.8617131 0.04856489 0.7764063 0.9330811
#> 7      20               40         60 0.8617131 0.04856489 0.7764063 0.9330811
#> 8      20               60         80 0.8617131 0.04856489 0.7764063 0.9330811
#> 9      40                0         40 0.8478591 0.05128241 0.7580993 0.9245739
#> 10     40               20         60 0.8478591 0.05128241 0.7580993 0.9245739
#> 11     40               40         80 0.8478591 0.05128241 0.7580993 0.9245739
#> 12     40               60        100 0.8478591 0.05128241 0.7580993 0.9245739
#> 13     60                0         60 0.8342249 0.05391083 0.7403810 0.9160049
#> 14     60               20         80 0.8342249 0.05391083 0.7403810 0.9160049
#> 15     60               40        100 0.8342249 0.05391083 0.7403810 0.9160049
#> 16     60               60        120 0.8342249 0.05391083 0.7403810 0.9160049
#> 17     80                0         80 0.8208071 0.05645027 0.7232037 0.9073945
#> 18     80               20        100 0.8208071 0.05645027 0.7232037 0.9073945
#> 19     80               40        120 0.8208071 0.05645027 0.7232037 0.9073945
#> 20     80               60        140 0.8208071 0.05645027 0.7232037 0.9073945
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.13044948 0.1365878761 0.5863966
#> 2  0.30574618 0.12841969 0.0978019069 0.5562075
#> 3  0.26001915 0.12760230 0.0687619137 0.5276285
#> 4  0.22113100 0.12666798 0.0472478183 0.5006058
#> 5  0.25589195 0.11931510 0.0741686266 0.4830964
#> 6  0.21762106 0.11573965 0.0512311412 0.4585468
#> 7  0.18507391 0.11317037 0.0344331961 0.4353788
#> 8  0.15739448 0.11076541 0.0223681621 0.4135231
#> 9  0.18213629 0.10651985 0.0375250505 0.3993818
#> 10 0.15489621 0.10241864 0.0245671675 0.3795743
#> 11 0.13173012 0.09914499 0.0154445812 0.3608953
#> 12 0.11202872 0.09615083 0.0092295504 0.3432794
#> 13 0.12963921 0.09365393 0.0170905703 0.3318808
#> 14 0.11025053 0.08956696 0.0103324081 0.3159098
#> 15 0.09376159 0.08615767 0.0058795508 0.3008383
#> 16 0.07973872 0.08306037 0.0031010918 0.2866105
#> 17 0.09227334 0.08151488 0.0066563404 0.2773945
#> 18 0.07847305 0.07771307 0.0035723873 0.2644650
#> 19 0.06673672 0.07446330 0.0017509940 0.2522425
#> 20 0.05675565 0.07151506 0.0007662781 0.2406807
```
