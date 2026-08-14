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
#> 1         0.1 0.3838513 0.02828642 0.3315460 0.4094790
#> 2         0.2 0.3832760 0.02823141 0.3311173 0.4089484
#> 3         0.3 0.3827015 0.02817656 0.3306890 0.4084185
#> 4         0.4 0.3821279 0.02812188 0.3302614 0.4078893
#> 5         0.5 0.3815551 0.02806737 0.3298342 0.4073608
#> 6         0.6 0.3809832 0.02801302 0.3294077 0.4068330
#> 7         0.7 0.3804122 0.02795885 0.3289817 0.4063058
#> 8         0.8 0.3798420 0.02790483 0.3285562 0.4057794
#> 9         0.9 0.3792726 0.02785099 0.3281313 0.4052536
#> 10        1.0 0.3787041 0.02779731 0.3277069 0.4047285

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3932495
#> 2       1        0.1        V5 0.4048150
#> 3       1        0.1        V9 0.3283700
#> 4       1        0.1        V3 0.4103054
#> 5       1        0.1        V4 0.3817663
#> 6       1        0.1        V8 0.4017509
#> 7       1        0.1        V2 0.3424857
#> 8       1        0.1        V6 0.3983260
#> 9       1        0.1        V7 0.4066325
#> 10      1        0.1       V10 0.3982031
#> 11      2        0.2        V7 0.4061281
#> 12      2        0.2        V8 0.4011712
#> 13      2        0.2        V9 0.3279471
#> 14      2        0.2       V10 0.3976405
#> 15      2        0.2        V1 0.3926982
#> 16      2        0.2        V5 0.4041810
#> 17      2        0.2        V2 0.3420366
#> 18      2        0.2        V3 0.4097672
#> 19      2        0.2        V4 0.3812781
#> 20      2        0.2        V6 0.3976974
#> 21      3        0.3        V4 0.3807906
#> 22      3        0.3        V5 0.4035479
#> 23      3        0.3        V3 0.4092298
#> 24      3        0.3        V7 0.4056242
#> 25      3        0.3        V8 0.4005924
#> 26      3        0.3        V9 0.3275248
#> 27      3        0.3        V6 0.3970697
#> 28      3        0.3       V10 0.3970786
#> 29      3        0.3        V1 0.3921477
#> 30      3        0.3        V2 0.3415880
#> 31      4        0.4        V1 0.3915979
#> 32      4        0.4        V9 0.3271030
#> 33      4        0.4        V3 0.4086930
#> 34      4        0.4        V4 0.3803037
#> 35      4        0.4        V5 0.4029158
#> 36      4        0.4        V2 0.3411400
#> 37      4        0.4        V6 0.3964431
#> 38      4        0.4        V7 0.4051210
#> 39      4        0.4        V8 0.4000144
#> 40      4        0.4       V10 0.3965176
#> 41      5        0.5        V8 0.3994372
#> 42      5        0.5        V9 0.3266818
#> 43      5        0.5       V10 0.3959573
#> 44      5        0.5        V1 0.3910489
#> 45      5        0.5        V5 0.4022847
#> 46      5        0.5        V2 0.3406926
#> 47      5        0.5        V3 0.4081570
#> 48      5        0.5        V4 0.3798174
#> 49      5        0.5        V6 0.3958174
#> 50      5        0.5        V7 0.4046184
#> 51      6        0.6        V4 0.3793318
#> 52      6        0.6        V5 0.4016546
#> 53      6        0.6        V7 0.4041165
#> 54      6        0.6        V8 0.3988609
#> 55      6        0.6        V9 0.3262611
#> 56      6        0.6        V6 0.3951928
#> 57      6        0.6       V10 0.3953979
#> 58      6        0.6        V1 0.3905007
#> 59      6        0.6        V2 0.3402458
#> 60      6        0.6        V3 0.4076217
#> 61      7        0.7        V1 0.3899532
#> 62      7        0.7        V3 0.4070870
#> 63      7        0.7        V4 0.3788467
#> 64      7        0.7        V5 0.4010255
#> 65      7        0.7        V9 0.3258410
#> 66      7        0.7        V6 0.3945691
#> 67      7        0.7        V7 0.4036151
#> 68      7        0.7        V8 0.3982853
#> 69      7        0.7        V2 0.3397996
#> 70      7        0.7       V10 0.3948392
#> 71      8        0.8        V9 0.3254214
#> 72      8        0.8       V10 0.3942813
#> 73      8        0.8        V1 0.3894065
#> 74      8        0.8        V5 0.4003974
#> 75      8        0.8        V2 0.3393539
#> 76      8        0.8        V3 0.4065531
#> 77      8        0.8        V4 0.3783623
#> 78      8        0.8        V8 0.3977107
#> 79      8        0.8        V6 0.3939464
#> 80      8        0.8        V7 0.4031144
#> 81      9        0.9        V5 0.3997702
#> 82      9        0.9        V7 0.4026143
#> 83      9        0.9        V8 0.3971368
#> 84      9        0.9        V9 0.3250023
#> 85      9        0.9        V6 0.3933247
#> 86      9        0.9       V10 0.3937242
#> 87      9        0.9        V1 0.3888606
#> 88      9        0.9        V2 0.3389089
#> 89      9        0.9        V3 0.4060198
#> 90      9        0.9        V4 0.3778785
#> 91     10        1.0        V1 0.3883154
#> 92     10        1.0        V4 0.3773953
#> 93     10        1.0        V5 0.3991441
#> 94     10        1.0        V9 0.3245838
#> 95     10        1.0        V3 0.4054873
#> 96     10        1.0        V7 0.4021148
#> 97     10        1.0        V8 0.3965638
#> 98     10        1.0        V2 0.3384644
#> 99     10        1.0        V6 0.3927039
#> 100    10        1.0       V10 0.3931679

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3932495
#> 2       1        0.1        V5 0.4048150
#> 3       1        0.1        V9 0.3283700
#> 4       1        0.1        V3 0.4103054
#> 5       1        0.1        V4 0.3817663
#> 6       1        0.1        V8 0.4017509
#> 7       1        0.1        V2 0.3424857
#> 8       1        0.1        V6 0.3983260
#> 9       1        0.1        V7 0.4066325
#> 10      1        0.1       V10 0.3982031
#> 11      2        0.2        V7 0.4061281
#> 12      2        0.2        V8 0.4011712
#> 13      2        0.2        V9 0.3279471
#> 14      2        0.2       V10 0.3976405
#> 15      2        0.2        V1 0.3926982
#> 16      2        0.2        V5 0.4041810
#> 17      2        0.2        V2 0.3420366
#> 18      2        0.2        V3 0.4097672
#> 19      2        0.2        V4 0.3812781
#> 20      2        0.2        V6 0.3976974
#> 21      3        0.3        V4 0.3807906
#> 22      3        0.3        V5 0.4035479
#> 23      3        0.3        V3 0.4092298
#> 24      3        0.3        V7 0.4056242
#> 25      3        0.3        V8 0.4005924
#> 26      3        0.3        V9 0.3275248
#> 27      3        0.3        V6 0.3970697
#> 28      3        0.3       V10 0.3970786
#> 29      3        0.3        V1 0.3921477
#> 30      3        0.3        V2 0.3415880
#> 31      4        0.4        V1 0.3915979
#> 32      4        0.4        V9 0.3271030
#> 33      4        0.4        V3 0.4086930
#> 34      4        0.4        V4 0.3803037
#> 35      4        0.4        V5 0.4029158
#> 36      4        0.4        V2 0.3411400
#> 37      4        0.4        V6 0.3964431
#> 38      4        0.4        V7 0.4051210
#> 39      4        0.4        V8 0.4000144
#> 40      4        0.4       V10 0.3965176
#> 41      5        0.5        V8 0.3994372
#> 42      5        0.5        V9 0.3266818
#> 43      5        0.5       V10 0.3959573
#> 44      5        0.5        V1 0.3910489
#> 45      5        0.5        V5 0.4022847
#> 46      5        0.5        V2 0.3406926
#> 47      5        0.5        V3 0.4081570
#> 48      5        0.5        V4 0.3798174
#> 49      5        0.5        V6 0.3958174
#> 50      5        0.5        V7 0.4046184
#> 51      6        0.6        V4 0.3793318
#> 52      6        0.6        V5 0.4016546
#> 53      6        0.6        V7 0.4041165
#> 54      6        0.6        V8 0.3988609
#> 55      6        0.6        V9 0.3262611
#> 56      6        0.6        V6 0.3951928
#> 57      6        0.6       V10 0.3953979
#> 58      6        0.6        V1 0.3905007
#> 59      6        0.6        V2 0.3402458
#> 60      6        0.6        V3 0.4076217
#> 61      7        0.7        V1 0.3899532
#> 62      7        0.7        V3 0.4070870
#> 63      7        0.7        V4 0.3788467
#> 64      7        0.7        V5 0.4010255
#> 65      7        0.7        V9 0.3258410
#> 66      7        0.7        V6 0.3945691
#> 67      7        0.7        V7 0.4036151
#> 68      7        0.7        V8 0.3982853
#> 69      7        0.7        V2 0.3397996
#> 70      7        0.7       V10 0.3948392
#> 71      8        0.8        V9 0.3254214
#> 72      8        0.8       V10 0.3942813
#> 73      8        0.8        V1 0.3894065
#> 74      8        0.8        V5 0.4003974
#> 75      8        0.8        V2 0.3393539
#> 76      8        0.8        V3 0.4065531
#> 77      8        0.8        V4 0.3783623
#> 78      8        0.8        V8 0.3977107
#> 79      8        0.8        V6 0.3939464
#> 80      8        0.8        V7 0.4031144
#> 81      9        0.9        V5 0.3997702
#> 82      9        0.9        V7 0.4026143
#> 83      9        0.9        V8 0.3971368
#> 84      9        0.9        V9 0.3250023
#> 85      9        0.9        V6 0.3933247
#> 86      9        0.9       V10 0.3937242
#> 87      9        0.9        V1 0.3888606
#> 88      9        0.9        V2 0.3389089
#> 89      9        0.9        V3 0.4060198
#> 90      9        0.9        V4 0.3778785
#> 91     10        1.0        V1 0.3883154
#> 92     10        1.0        V4 0.3773953
#> 93     10        1.0        V5 0.3991441
#> 94     10        1.0        V9 0.3245838
#> 95     10        1.0        V3 0.4054873
#> 96     10        1.0        V7 0.4021148
#> 97     10        1.0        V8 0.3965638
#> 98     10        1.0        V2 0.3384644
#> 99     10        1.0        V6 0.3927039
#> 100    10        1.0       V10 0.3931679


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
#> 1       0                0          0 0.8757906 0.04970312 0.7738566 0.9479748
#> 2       0               20         20 0.8757906 0.04970312 0.7738566 0.9479748
#> 3       0               40         40 0.8757906 0.04970312 0.7738566 0.9479748
#> 4       0               60         60 0.8757906 0.04970312 0.7738566 0.9479748
#> 5      20                0         20 0.8617131 0.05244071 0.7543969 0.9401705
#> 6      20               20         40 0.8617131 0.05244071 0.7543969 0.9401705
#> 7      20               40         60 0.8617131 0.05244071 0.7543969 0.9401705
#> 8      20               60         80 0.8617131 0.05244071 0.7543969 0.9401705
#> 9      40                0         40 0.8478591 0.05504156 0.7356612 0.9322554
#> 10     40               20         60 0.8478591 0.05504156 0.7356612 0.9322554
#> 11     40               40         80 0.8478591 0.05504156 0.7356612 0.9322554
#> 12     40               60        100 0.8478591 0.05504156 0.7356612 0.9322554
#> 13     60                0         60 0.8342249 0.05751690 0.7175776 0.9242562
#> 14     60               20         80 0.8342249 0.05751690 0.7175776 0.9242562
#> 15     60               40        100 0.8342249 0.05751690 0.7175776 0.9242562
#> 16     60               60        120 0.8342249 0.05751690 0.7175776 0.9242562
#> 17     80                0         80 0.8208071 0.05987578 0.7000897 0.9161947
#> 18     80               20        100 0.8208071 0.05987578 0.7000897 0.9161947
#> 19     80               40        120 0.8208071 0.05987578 0.7000897 0.9161947
#> 20     80               60        140 0.8208071 0.05987578 0.7000897 0.9161947
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.12832402 0.1520989328 0.5641961
#> 2  0.30574618 0.12476679 0.1189216493 0.5241949
#> 3  0.26001915 0.12088350 0.0921440741 0.4871841
#> 4  0.22113100 0.11667866 0.0706264454 0.4529948
#> 5  0.25589195 0.11491230 0.0751591060 0.4521610
#> 6  0.21762106 0.11027058 0.0570514098 0.4206770
#> 7  0.18507391 0.10561963 0.0426781622 0.3916448
#> 8  0.15739448 0.10094708 0.0313807993 0.3648829
#> 9  0.18213629 0.10017472 0.0337360328 0.3642308
#> 10 0.15489621 0.09536011 0.0244282410 0.3396148
#> 11 0.13173012 0.09068614 0.0172853303 0.3169225
#> 12 0.11202872 0.08613511 0.0119026202 0.2959962
#> 13 0.12963921 0.08590517 0.0130031172 0.2954860
#> 14 0.11025053 0.08136580 0.0087390597 0.2762174
#> 15 0.09376159 0.07703090 0.0056622324 0.2584263
#> 16 0.07973872 0.07288104 0.0035133003 0.2419851
#> 17 0.09227334 0.07293449 0.0039373229 0.2415837
#> 18 0.07847305 0.06887059 0.0023504815 0.2264044
#> 19 0.06673672 0.06502676 0.0013239961 0.2123455
#> 20 0.05675565 0.06138262 0.0006964962 0.1993069
```
