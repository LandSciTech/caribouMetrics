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
[`dataFromSheets()`](https://landscitech.github.io/caribouMetrics/dev/reference/dataFromSheets.md),
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
#> 1         0.1 0.3838513 0.02214075 0.3340439 0.3995338
#> 2         0.2 0.3832760 0.02211494 0.3335123 0.3989257
#> 3         0.3 0.3827015 0.02208933 0.3329816 0.3983186
#> 4         0.4 0.3821279 0.02206394 0.3324518 0.3977124
#> 5         0.5 0.3815551 0.02203875 0.3319227 0.3971072
#> 6         0.6 0.3809832 0.02201377 0.3313946 0.3965028
#> 7         0.7 0.3804122 0.02198900 0.3308672 0.3958994
#> 8         0.8 0.3798420 0.02196443 0.3303407 0.3952969
#> 9         0.9 0.3792726 0.02194006 0.3298151 0.3946953
#> 10        1.0 0.3787041 0.02191590 0.3292903 0.3940946

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3806073
#> 2       1        0.1        V5 0.3771326
#> 3       1        0.1        V9 0.4007173
#> 4       1        0.1        V3 0.3627189
#> 5       1        0.1        V4 0.3913068
#> 6       1        0.1        V8 0.3954572
#> 7       1        0.1        V2 0.3282914
#> 8       1        0.1        V6 0.3538580
#> 9       1        0.1        V7 0.3909125
#> 10      1        0.1       V10 0.3744057
#> 11      2        0.2        V7 0.3903838
#> 12      2        0.2        V8 0.3948446
#> 13      2        0.2        V9 0.4001106
#> 14      2        0.2       V10 0.3739054
#> 15      2        0.2        V1 0.3801715
#> 16      2        0.2        V5 0.3766394
#> 17      2        0.2        V2 0.3277490
#> 18      2        0.2        V3 0.3621728
#> 19      2        0.2        V4 0.3906421
#> 20      2        0.2        V6 0.3533638
#> 21      3        0.3        V4 0.3899785
#> 22      3        0.3        V5 0.3761469
#> 23      3        0.3        V3 0.3616276
#> 24      3        0.3        V7 0.3898558
#> 25      3        0.3        V8 0.3942330
#> 26      3        0.3        V9 0.3995048
#> 27      3        0.3        V6 0.3528703
#> 28      3        0.3       V10 0.3734059
#> 29      3        0.3        V1 0.3797361
#> 30      3        0.3        V2 0.3272075
#> 31      4        0.4        V1 0.3793012
#> 32      4        0.4        V9 0.3988999
#> 33      4        0.4        V3 0.3610831
#> 34      4        0.4        V4 0.3893161
#> 35      4        0.4        V5 0.3756550
#> 36      4        0.4        V2 0.3266669
#> 37      4        0.4        V6 0.3523775
#> 38      4        0.4        V7 0.3893285
#> 39      4        0.4        V8 0.3936224
#> 40      4        0.4       V10 0.3729070
#> 41      5        0.5        V8 0.3930127
#> 42      5        0.5        V9 0.3982959
#> 43      5        0.5       V10 0.3724088
#> 44      5        0.5        V1 0.3788669
#> 45      5        0.5        V5 0.3751637
#> 46      5        0.5        V2 0.3261271
#> 47      5        0.5        V3 0.3605395
#> 48      5        0.5        V4 0.3886548
#> 49      5        0.5        V6 0.3518854
#> 50      5        0.5        V7 0.3888020
#> 51      6        0.6        V4 0.3879946
#> 52      6        0.6        V5 0.3746731
#> 53      6        0.6        V7 0.3882761
#> 54      6        0.6        V8 0.3924039
#> 55      6        0.6        V9 0.3976928
#> 56      6        0.6        V6 0.3513940
#> 57      6        0.6       V10 0.3719113
#> 58      6        0.6        V1 0.3784330
#> 59      6        0.6        V2 0.3255883
#> 60      6        0.6        V3 0.3599967
#> 61      7        0.7        V1 0.3779997
#> 62      7        0.7        V3 0.3594547
#> 63      7        0.7        V4 0.3873355
#> 64      7        0.7        V5 0.3741831
#> 65      7        0.7        V9 0.3970906
#> 66      7        0.7        V6 0.3509033
#> 67      7        0.7        V7 0.3877510
#> 68      7        0.7        V8 0.3917961
#> 69      7        0.7        V2 0.3250503
#> 70      7        0.7       V10 0.3714144
#> 71      8        0.8        V9 0.3964894
#> 72      8        0.8       V10 0.3709182
#> 73      8        0.8        V1 0.3775668
#> 74      8        0.8        V5 0.3736938
#> 75      8        0.8        V2 0.3245132
#> 76      8        0.8        V3 0.3589136
#> 77      8        0.8        V4 0.3866776
#> 78      8        0.8        V8 0.3911892
#> 79      8        0.8        V6 0.3504132
#> 80      8        0.8        V7 0.3872265
#> 81      9        0.9        V5 0.3732051
#> 82      9        0.9        V7 0.3867028
#> 83      9        0.9        V8 0.3905833
#> 84      9        0.9        V9 0.3958891
#> 85      9        0.9        V6 0.3499239
#> 86      9        0.9       V10 0.3704226
#> 87      9        0.9        V1 0.3771344
#> 88      9        0.9        V2 0.3239771
#> 89      9        0.9        V3 0.3583732
#> 90      9        0.9        V4 0.3860208
#> 91     10        1.0        V1 0.3767025
#> 92     10        1.0        V4 0.3853651
#> 93     10        1.0        V5 0.3727171
#> 94     10        1.0        V9 0.3952896
#> 95     10        1.0        V3 0.3578337
#> 96     10        1.0        V7 0.3861798
#> 97     10        1.0        V8 0.3899783
#> 98     10        1.0        V2 0.3234418
#> 99     10        1.0        V6 0.3494352
#> 100    10        1.0       V10 0.3699277

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.3806073
#> 2       1        0.1        V5 0.3771326
#> 3       1        0.1        V9 0.4007173
#> 4       1        0.1        V3 0.3627189
#> 5       1        0.1        V4 0.3913068
#> 6       1        0.1        V8 0.3954572
#> 7       1        0.1        V2 0.3282914
#> 8       1        0.1        V6 0.3538580
#> 9       1        0.1        V7 0.3909125
#> 10      1        0.1       V10 0.3744057
#> 11      2        0.2        V7 0.3903838
#> 12      2        0.2        V8 0.3948446
#> 13      2        0.2        V9 0.4001106
#> 14      2        0.2       V10 0.3739054
#> 15      2        0.2        V1 0.3801715
#> 16      2        0.2        V5 0.3766394
#> 17      2        0.2        V2 0.3277490
#> 18      2        0.2        V3 0.3621728
#> 19      2        0.2        V4 0.3906421
#> 20      2        0.2        V6 0.3533638
#> 21      3        0.3        V4 0.3899785
#> 22      3        0.3        V5 0.3761469
#> 23      3        0.3        V3 0.3616276
#> 24      3        0.3        V7 0.3898558
#> 25      3        0.3        V8 0.3942330
#> 26      3        0.3        V9 0.3995048
#> 27      3        0.3        V6 0.3528703
#> 28      3        0.3       V10 0.3734059
#> 29      3        0.3        V1 0.3797361
#> 30      3        0.3        V2 0.3272075
#> 31      4        0.4        V1 0.3793012
#> 32      4        0.4        V9 0.3988999
#> 33      4        0.4        V3 0.3610831
#> 34      4        0.4        V4 0.3893161
#> 35      4        0.4        V5 0.3756550
#> 36      4        0.4        V2 0.3266669
#> 37      4        0.4        V6 0.3523775
#> 38      4        0.4        V7 0.3893285
#> 39      4        0.4        V8 0.3936224
#> 40      4        0.4       V10 0.3729070
#> 41      5        0.5        V8 0.3930127
#> 42      5        0.5        V9 0.3982959
#> 43      5        0.5       V10 0.3724088
#> 44      5        0.5        V1 0.3788669
#> 45      5        0.5        V5 0.3751637
#> 46      5        0.5        V2 0.3261271
#> 47      5        0.5        V3 0.3605395
#> 48      5        0.5        V4 0.3886548
#> 49      5        0.5        V6 0.3518854
#> 50      5        0.5        V7 0.3888020
#> 51      6        0.6        V4 0.3879946
#> 52      6        0.6        V5 0.3746731
#> 53      6        0.6        V7 0.3882761
#> 54      6        0.6        V8 0.3924039
#> 55      6        0.6        V9 0.3976928
#> 56      6        0.6        V6 0.3513940
#> 57      6        0.6       V10 0.3719113
#> 58      6        0.6        V1 0.3784330
#> 59      6        0.6        V2 0.3255883
#> 60      6        0.6        V3 0.3599967
#> 61      7        0.7        V1 0.3779997
#> 62      7        0.7        V3 0.3594547
#> 63      7        0.7        V4 0.3873355
#> 64      7        0.7        V5 0.3741831
#> 65      7        0.7        V9 0.3970906
#> 66      7        0.7        V6 0.3509033
#> 67      7        0.7        V7 0.3877510
#> 68      7        0.7        V8 0.3917961
#> 69      7        0.7        V2 0.3250503
#> 70      7        0.7       V10 0.3714144
#> 71      8        0.8        V9 0.3964894
#> 72      8        0.8       V10 0.3709182
#> 73      8        0.8        V1 0.3775668
#> 74      8        0.8        V5 0.3736938
#> 75      8        0.8        V2 0.3245132
#> 76      8        0.8        V3 0.3589136
#> 77      8        0.8        V4 0.3866776
#> 78      8        0.8        V8 0.3911892
#> 79      8        0.8        V6 0.3504132
#> 80      8        0.8        V7 0.3872265
#> 81      9        0.9        V5 0.3732051
#> 82      9        0.9        V7 0.3867028
#> 83      9        0.9        V8 0.3905833
#> 84      9        0.9        V9 0.3958891
#> 85      9        0.9        V6 0.3499239
#> 86      9        0.9       V10 0.3704226
#> 87      9        0.9        V1 0.3771344
#> 88      9        0.9        V2 0.3239771
#> 89      9        0.9        V3 0.3583732
#> 90      9        0.9        V4 0.3860208
#> 91     10        1.0        V1 0.3767025
#> 92     10        1.0        V4 0.3853651
#> 93     10        1.0        V5 0.3727171
#> 94     10        1.0        V9 0.3952896
#> 95     10        1.0        V3 0.3578337
#> 96     10        1.0        V7 0.3861798
#> 97     10        1.0        V8 0.3899783
#> 98     10        1.0        V2 0.3234418
#> 99     10        1.0        V6 0.3494352
#> 100    10        1.0       V10 0.3699277


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
#> 1       0                0          0 0.8757906 0.04541221 0.7934586 0.9435692
#> 2       0               20         20 0.8757906 0.04541221 0.7934586 0.9435692
#> 3       0               40         40 0.8757906 0.04541221 0.7934586 0.9435692
#> 4       0               60         60 0.8757906 0.04541221 0.7934586 0.9435692
#> 5      20                0         20 0.8617131 0.04837267 0.7744947 0.9346474
#> 6      20               20         40 0.8617131 0.04837267 0.7744947 0.9346474
#> 7      20               40         60 0.8617131 0.04837267 0.7744947 0.9346474
#> 8      20               60         80 0.8617131 0.04837267 0.7744947 0.9346474
#> 9      40                0         40 0.8478591 0.05111291 0.7562171 0.9256112
#> 10     40               20         60 0.8478591 0.05111291 0.7562171 0.9256112
#> 11     40               40         80 0.8478591 0.05111291 0.7562171 0.9256112
#> 12     40               60        100 0.8478591 0.05111291 0.7562171 0.9256112
#> 13     60                0         60 0.8342249 0.05365957 0.7385561 0.9164936
#> 14     60               20         80 0.8342249 0.05365957 0.7385561 0.9164936
#> 15     60               40        100 0.8342249 0.05365957 0.7385561 0.9164936
#> 16     60               60        120 0.8342249 0.05365957 0.7385561 0.9164936
#> 17     80                0         80 0.8208071 0.05603372 0.7214574 0.9073207
#> 18     80               20        100 0.8208071 0.05603372 0.7214574 0.9073207
#> 19     80               40        120 0.8208071 0.05603372 0.7214574 0.9073207
#> 20     80               60        140 0.8208071 0.05603372 0.7214574 0.9073207
#>         R_bar   R_stdErr     R_PIlow  R_PIhigh
#> 1  0.35951478 0.11584962 0.192064003 0.5811132
#> 2  0.30574618 0.11977648 0.146458355 0.5499457
#> 3  0.26001915 0.12347382 0.110587416 0.5204908
#> 4  0.22113100 0.12599350 0.082480801 0.4926852
#> 5  0.25589195 0.10239038 0.117828222 0.4561759
#> 6  0.21762106 0.10409488 0.088141751 0.4320497
#> 7  0.18507391 0.10553883 0.064992313 0.4093239
#> 8  0.15739448 0.10620195 0.047094430 0.3879242
#> 9  0.18213629 0.08924785 0.069643326 0.3598863
#> 10 0.15489621 0.08954613 0.050675196 0.3413851
#> 11 0.13173012 0.08963018 0.036137113 0.3239704
#> 12 0.11202872 0.08921148 0.025151784 0.3075775
#> 13 0.12963921 0.07701159 0.039033496 0.2860994
#> 14 0.11025053 0.07646496 0.027325492 0.2719210
#> 15 0.09376159 0.07576879 0.018600262 0.2585667
#> 16 0.07973872 0.07476141 0.012240358 0.2459852
#> 17 0.09227334 0.06592109 0.020315359 0.2294784
#> 18 0.07847305 0.06490765 0.013477310 0.2185637
#> 19 0.06673672 0.06380498 0.008598097 0.2082669
#> 20 0.05675565 0.06252097 0.005232483 0.1985484
```
