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
#> 1         0.1 0.3838513 0.02715407 0.3510690 0.4330433
#> 2         0.2 0.3832760 0.02712238 0.3505098 0.4323775
#> 3         0.3 0.3827015 0.02709077 0.3499514 0.4317127
#> 4         0.4 0.3821279 0.02705924 0.3493940 0.4310489
#> 5         0.5 0.3815551 0.02702777 0.3488374 0.4303862
#> 6         0.6 0.3809832 0.02699639 0.3482817 0.4297245
#> 7         0.7 0.3804122 0.02696507 0.3477269 0.4290638
#> 8         0.8 0.3798420 0.02693382 0.3471730 0.4284042
#> 9         0.9 0.3792726 0.02690265 0.3466199 0.4277455
#> 10        1.0 0.3787041 0.02687155 0.3460678 0.4270879

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3798739
#> 2       1        0.1        V5 0.4378343
#> 3       1        0.1        V9 0.3503513
#> 4       1        0.1        V3 0.3936041
#> 5       1        0.1        V4 0.3873739
#> 6       1        0.1        V8 0.3647233
#> 7       1        0.1        V2 0.3827747
#> 8       1        0.1        V6 0.3535412
#> 9       1        0.1        V7 0.4165410
#> 10      1        0.1       V10 0.3728033
#> 11      2        0.2        V7 0.4159544
#> 12      2        0.2        V8 0.3641739
#> 13      2        0.2        V9 0.3497906
#> 14      2        0.2       V10 0.3722111
#> 15      2        0.2        V1 0.3792883
#> 16      2        0.2        V5 0.4371455
#> 17      2        0.2        V2 0.3822375
#> 18      2        0.2        V3 0.3930296
#> 19      2        0.2        V4 0.3868490
#> 20      2        0.2        V6 0.3529870
#> 21      3        0.3        V4 0.3863249
#> 22      3        0.3        V5 0.4364577
#> 23      3        0.3        V3 0.3924559
#> 24      3        0.3        V7 0.4153687
#> 25      3        0.3        V8 0.3636253
#> 26      3        0.3        V9 0.3492307
#> 27      3        0.3        V6 0.3524338
#> 28      3        0.3       V10 0.3716198
#> 29      3        0.3        V1 0.3787036
#> 30      3        0.3        V2 0.3817010
#> 31      4        0.4        V1 0.3781198
#> 32      4        0.4        V9 0.3486718
#> 33      4        0.4        V3 0.3918830
#> 34      4        0.4        V4 0.3858014
#> 35      4        0.4        V5 0.4357711
#> 36      4        0.4        V2 0.3811654
#> 37      4        0.4        V6 0.3518814
#> 38      4        0.4        V7 0.4147838
#> 39      4        0.4        V8 0.3630775
#> 40      4        0.4       V10 0.3710294
#> 41      5        0.5        V8 0.3625305
#> 42      5        0.5        V9 0.3481138
#> 43      5        0.5       V10 0.3704400
#> 44      5        0.5        V1 0.3775369
#> 45      5        0.5        V5 0.4350855
#> 46      5        0.5        V2 0.3806304
#> 47      5        0.5        V3 0.3913110
#> 48      5        0.5        V4 0.3852786
#> 49      5        0.5        V6 0.3513298
#> 50      5        0.5        V7 0.4141997
#> 51      6        0.6        V4 0.3847566
#> 52      6        0.6        V5 0.4344011
#> 53      6        0.6        V7 0.4136164
#> 54      6        0.6        V8 0.3619844
#> 55      6        0.6        V9 0.3475566
#> 56      6        0.6        V6 0.3507792
#> 57      6        0.6       V10 0.3698516
#> 58      6        0.6        V1 0.3769548
#> 59      6        0.6        V2 0.3800962
#> 60      6        0.6        V3 0.3907398
#> 61      7        0.7        V1 0.3763737
#> 62      7        0.7        V3 0.3901694
#> 63      7        0.7        V4 0.3842352
#> 64      7        0.7        V5 0.4337177
#> 65      7        0.7        V9 0.3470004
#> 66      7        0.7        V6 0.3502294
#> 67      7        0.7        V7 0.4130339
#> 68      7        0.7        V8 0.3614390
#> 69      7        0.7        V2 0.3795628
#> 70      7        0.7       V10 0.3692640
#> 71      8        0.8        V9 0.3464450
#> 72      8        0.8       V10 0.3686774
#> 73      8        0.8        V1 0.3757935
#> 74      8        0.8        V5 0.4330353
#> 75      8        0.8        V2 0.3790301
#> 76      8        0.8        V3 0.3895999
#> 77      8        0.8        V4 0.3837146
#> 78      8        0.8        V8 0.3608945
#> 79      8        0.8        V6 0.3496804
#> 80      8        0.8        V7 0.4124523
#> 81      9        0.9        V5 0.4323541
#> 82      9        0.9        V7 0.4118715
#> 83      9        0.9        V8 0.3603508
#> 84      9        0.9        V9 0.3458905
#> 85      9        0.9        V6 0.3491324
#> 86      9        0.9       V10 0.3680918
#> 87      9        0.9        V1 0.3752142
#> 88      9        0.9        V2 0.3784982
#> 89      9        0.9        V3 0.3890312
#> 90      9        0.9        V4 0.3831946
#> 91     10        1.0        V1 0.3746358
#> 92     10        1.0        V4 0.3826754
#> 93     10        1.0        V5 0.4316739
#> 94     10        1.0        V9 0.3453369
#> 95     10        1.0        V3 0.3884633
#> 96     10        1.0        V7 0.4112915
#> 97     10        1.0        V8 0.3598080
#> 98     10        1.0        V2 0.3779670
#> 99     10        1.0        V6 0.3485851
#> 100    10        1.0       V10 0.3675070

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3798739
#> 2       1        0.1        V5 0.4378343
#> 3       1        0.1        V9 0.3503513
#> 4       1        0.1        V3 0.3936041
#> 5       1        0.1        V4 0.3873739
#> 6       1        0.1        V8 0.3647233
#> 7       1        0.1        V2 0.3827747
#> 8       1        0.1        V6 0.3535412
#> 9       1        0.1        V7 0.4165410
#> 10      1        0.1       V10 0.3728033
#> 11      2        0.2        V7 0.4159544
#> 12      2        0.2        V8 0.3641739
#> 13      2        0.2        V9 0.3497906
#> 14      2        0.2       V10 0.3722111
#> 15      2        0.2        V1 0.3792883
#> 16      2        0.2        V5 0.4371455
#> 17      2        0.2        V2 0.3822375
#> 18      2        0.2        V3 0.3930296
#> 19      2        0.2        V4 0.3868490
#> 20      2        0.2        V6 0.3529870
#> 21      3        0.3        V4 0.3863249
#> 22      3        0.3        V5 0.4364577
#> 23      3        0.3        V3 0.3924559
#> 24      3        0.3        V7 0.4153687
#> 25      3        0.3        V8 0.3636253
#> 26      3        0.3        V9 0.3492307
#> 27      3        0.3        V6 0.3524338
#> 28      3        0.3       V10 0.3716198
#> 29      3        0.3        V1 0.3787036
#> 30      3        0.3        V2 0.3817010
#> 31      4        0.4        V1 0.3781198
#> 32      4        0.4        V9 0.3486718
#> 33      4        0.4        V3 0.3918830
#> 34      4        0.4        V4 0.3858014
#> 35      4        0.4        V5 0.4357711
#> 36      4        0.4        V2 0.3811654
#> 37      4        0.4        V6 0.3518814
#> 38      4        0.4        V7 0.4147838
#> 39      4        0.4        V8 0.3630775
#> 40      4        0.4       V10 0.3710294
#> 41      5        0.5        V8 0.3625305
#> 42      5        0.5        V9 0.3481138
#> 43      5        0.5       V10 0.3704400
#> 44      5        0.5        V1 0.3775369
#> 45      5        0.5        V5 0.4350855
#> 46      5        0.5        V2 0.3806304
#> 47      5        0.5        V3 0.3913110
#> 48      5        0.5        V4 0.3852786
#> 49      5        0.5        V6 0.3513298
#> 50      5        0.5        V7 0.4141997
#> 51      6        0.6        V4 0.3847566
#> 52      6        0.6        V5 0.4344011
#> 53      6        0.6        V7 0.4136164
#> 54      6        0.6        V8 0.3619844
#> 55      6        0.6        V9 0.3475566
#> 56      6        0.6        V6 0.3507792
#> 57      6        0.6       V10 0.3698516
#> 58      6        0.6        V1 0.3769548
#> 59      6        0.6        V2 0.3800962
#> 60      6        0.6        V3 0.3907398
#> 61      7        0.7        V1 0.3763737
#> 62      7        0.7        V3 0.3901694
#> 63      7        0.7        V4 0.3842352
#> 64      7        0.7        V5 0.4337177
#> 65      7        0.7        V9 0.3470004
#> 66      7        0.7        V6 0.3502294
#> 67      7        0.7        V7 0.4130339
#> 68      7        0.7        V8 0.3614390
#> 69      7        0.7        V2 0.3795628
#> 70      7        0.7       V10 0.3692640
#> 71      8        0.8        V9 0.3464450
#> 72      8        0.8       V10 0.3686774
#> 73      8        0.8        V1 0.3757935
#> 74      8        0.8        V5 0.4330353
#> 75      8        0.8        V2 0.3790301
#> 76      8        0.8        V3 0.3895999
#> 77      8        0.8        V4 0.3837146
#> 78      8        0.8        V8 0.3608945
#> 79      8        0.8        V6 0.3496804
#> 80      8        0.8        V7 0.4124523
#> 81      9        0.9        V5 0.4323541
#> 82      9        0.9        V7 0.4118715
#> 83      9        0.9        V8 0.3603508
#> 84      9        0.9        V9 0.3458905
#> 85      9        0.9        V6 0.3491324
#> 86      9        0.9       V10 0.3680918
#> 87      9        0.9        V1 0.3752142
#> 88      9        0.9        V2 0.3784982
#> 89      9        0.9        V3 0.3890312
#> 90      9        0.9        V4 0.3831946
#> 91     10        1.0        V1 0.3746358
#> 92     10        1.0        V4 0.3826754
#> 93     10        1.0        V5 0.4316739
#> 94     10        1.0        V9 0.3453369
#> 95     10        1.0        V3 0.3884633
#> 96     10        1.0        V7 0.4112915
#> 97     10        1.0        V8 0.3598080
#> 98     10        1.0        V2 0.3779670
#> 99     10        1.0        V6 0.3485851
#> 100    10        1.0       V10 0.3675070


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
#> 1       0                0          0 0.8757906 0.05322235 0.7652580 0.9538323
#> 2       0               20         20 0.8757906 0.05322235 0.7652580 0.9538323
#> 3       0               40         40 0.8757906 0.05322235 0.7652580 0.9538323
#> 4       0               60         60 0.8757906 0.05322235 0.7652580 0.9538323
#> 5      20                0         20 0.8617131 0.05424978 0.7472841 0.9389913
#> 6      20               20         40 0.8617131 0.05424978 0.7472841 0.9389913
#> 7      20               40         60 0.8617131 0.05424978 0.7472841 0.9389913
#> 8      20               60         80 0.8617131 0.05424978 0.7472841 0.9389913
#> 9      40                0         40 0.8478591 0.05512948 0.7299464 0.9237783
#> 10     40               20         60 0.8478591 0.05512948 0.7299464 0.9237783
#> 11     40               40         80 0.8478591 0.05512948 0.7299464 0.9237783
#> 12     40               60        100 0.8478591 0.05512948 0.7299464 0.9237783
#> 13     60                0         60 0.8342249 0.05590205 0.7131842 0.9083569
#> 14     60               20         80 0.8342249 0.05590205 0.7131842 0.9083569
#> 15     60               40        100 0.8342249 0.05590205 0.7131842 0.9083569
#> 16     60               60        120 0.8342249 0.05590205 0.7131842 0.9083569
#> 17     80                0         80 0.8208071 0.05659548 0.6969493 0.8928383
#> 18     80               20        100 0.8208071 0.05659548 0.6969493 0.8928383
#> 19     80               40        120 0.8208071 0.05659548 0.6969493 0.8928383
#> 20     80               60        140 0.8208071 0.05659548 0.6969493 0.8928383
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11703928 0.198933309 0.5636801
#> 2  0.30574618 0.11774197 0.157295134 0.5413947
#> 3  0.26001915 0.11845347 0.123539481 0.5200276
#> 4  0.22113100 0.11878523 0.096238574 0.4995535
#> 5  0.25589195 0.10768965 0.114903572 0.4566236
#> 6  0.21762106 0.10648134 0.089270270 0.4388539
#> 7  0.18507391 0.10543361 0.068647196 0.4218529
#> 8  0.15739448 0.10429383 0.052151162 0.4055904
#> 9  0.18213629 0.09627524 0.063406897 0.3715582
#> 10 0.15489621 0.09420768 0.047979257 0.3574899
#> 11 0.13173012 0.09234899 0.035762797 0.3440361
#> 12 0.11202872 0.09053415 0.026186946 0.3311694
#> 13 0.12963921 0.08455254 0.032697924 0.3042405
#> 14 0.11025053 0.08218763 0.023803780 0.2931026
#> 15 0.09376159 0.08003763 0.016949667 0.2824452
#> 16 0.07973872 0.07799282 0.011756171 0.2722455
#> 17 0.09227334 0.07338824 0.015266639 0.2508665
#> 18 0.07847305 0.07102543 0.010497864 0.2420077
#> 19 0.06673672 0.06886579 0.006984863 0.2335198
#> 20 0.05675565 0.06683461 0.004468391 0.2253846
```
