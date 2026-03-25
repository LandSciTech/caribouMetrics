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
#> 1         0.1 0.3838513 0.01690935 0.3660162 0.4092285
#> 2         0.2 0.3832760 0.01690145 0.3654427 0.4086010
#> 3         0.3 0.3827015 0.01689373 0.3648702 0.4079744
#> 4         0.4 0.3821279 0.01688619 0.3642985 0.4073488
#> 5         0.5 0.3815551 0.01687883 0.3637278 0.4067242
#> 6         0.6 0.3809832 0.01687164 0.3631579 0.4061005
#> 7         0.7 0.3804122 0.01686462 0.3625890 0.4054778
#> 8         0.8 0.3798420 0.01685778 0.3620210 0.4048561
#> 9         0.9 0.3792726 0.01685111 0.3614538 0.4042353
#> 10        1.0 0.3787041 0.01684461 0.3608875 0.4036155

# return one row per replicate * scenario
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE)
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4021340
#> 2       1        0.1        V5 0.3977104
#> 3       1        0.1        V9 0.4107341
#> 4       1        0.1        V3 0.3758944
#> 5       1        0.1        V4 0.4040427
#> 6       1        0.1        V8 0.3654231
#> 7       1        0.1        V2 0.4038891
#> 8       1        0.1        V6 0.3680588
#> 9       1        0.1        V7 0.3959103
#> 10      1        0.1       V10 0.3749178
#> 11      2        0.2        V7 0.3952284
#> 12      2        0.2        V8 0.3648700
#> 13      2        0.2        V9 0.4100781
#> 14      2        0.2       V10 0.3743130
#> 15      2        0.2        V1 0.4016009
#> 16      2        0.2        V5 0.3971433
#> 17      2        0.2        V2 0.4032353
#> 18      2        0.2        V3 0.3753482
#> 19      2        0.2        V4 0.4035132
#> 20      2        0.2        V6 0.3674153
#> 21      3        0.3        V4 0.4029845
#> 22      3        0.3        V5 0.3965771
#> 23      3        0.3        V3 0.3748028
#> 24      3        0.3        V7 0.3945476
#> 25      3        0.3        V8 0.3643178
#> 26      3        0.3        V9 0.4094231
#> 27      3        0.3        V6 0.3667729
#> 28      3        0.3       V10 0.3737092
#> 29      3        0.3        V1 0.4010686
#> 30      3        0.3        V2 0.4025826
#> 31      4        0.4        V1 0.4005369
#> 32      4        0.4        V9 0.4087692
#> 33      4        0.4        V3 0.3742582
#> 34      4        0.4        V4 0.4024564
#> 35      4        0.4        V5 0.3960116
#> 36      4        0.4        V2 0.4019309
#> 37      4        0.4        V6 0.3661317
#> 38      4        0.4        V7 0.3938679
#> 39      4        0.4        V8 0.3637663
#> 40      4        0.4       V10 0.3731064
#> 41      5        0.5        V8 0.3632157
#> 42      5        0.5        V9 0.4081163
#> 43      5        0.5       V10 0.3725046
#> 44      5        0.5        V1 0.4000059
#> 45      5        0.5        V5 0.3954470
#> 46      5        0.5        V2 0.4012803
#> 47      5        0.5        V3 0.3737143
#> 48      5        0.5        V4 0.4019291
#> 49      5        0.5        V6 0.3654915
#> 50      5        0.5        V7 0.3931895
#> 51      6        0.6        V4 0.4014024
#> 52      6        0.6        V5 0.3948831
#> 53      6        0.6        V7 0.3925122
#> 54      6        0.6        V8 0.3626660
#> 55      6        0.6        V9 0.4074645
#> 56      6        0.6        V6 0.3648525
#> 57      6        0.6       V10 0.3719037
#> 58      6        0.6        V1 0.3994756
#> 59      6        0.6        V2 0.4006307
#> 60      6        0.6        V3 0.3731713
#> 61      7        0.7        V1 0.3989461
#> 62      7        0.7        V3 0.3726291
#> 63      7        0.7        V4 0.4008764
#> 64      7        0.7        V5 0.3943201
#> 65      7        0.7        V9 0.4068137
#> 66      7        0.7        V6 0.3642146
#> 67      7        0.7        V7 0.3918360
#> 68      7        0.7        V8 0.3621170
#> 69      7        0.7        V2 0.3999822
#> 70      7        0.7       V10 0.3713038
#> 71      8        0.8        V9 0.4061640
#> 72      8        0.8       V10 0.3707048
#> 73      8        0.8        V1 0.3984172
#> 74      8        0.8        V5 0.3937578
#> 75      8        0.8        V2 0.3993347
#> 76      8        0.8        V3 0.3720876
#> 77      8        0.8        V4 0.4003511
#> 78      8        0.8        V8 0.3615689
#> 79      8        0.8        V6 0.3635779
#> 80      8        0.8        V7 0.3911611
#> 81      9        0.9        V5 0.3931964
#> 82      9        0.9        V7 0.3904873
#> 83      9        0.9        V8 0.3610217
#> 84      9        0.9        V9 0.4055153
#> 85      9        0.9        V6 0.3629422
#> 86      9        0.9       V10 0.3701068
#> 87      9        0.9        V1 0.3978890
#> 88      9        0.9        V2 0.3986883
#> 89      9        0.9        V3 0.3715469
#> 90      9        0.9        V4 0.3998265
#> 91     10        1.0        V1 0.3973616
#> 92     10        1.0        V4 0.3993026
#> 93     10        1.0        V5 0.3926357
#> 94     10        1.0        V9 0.4048676
#> 95     10        1.0        V3 0.3710070
#> 96     10        1.0        V7 0.3898146
#> 97     10        1.0        V8 0.3604752
#> 98     10        1.0        V2 0.3980430
#> 99     10        1.0        V6 0.3623077
#> 100    10        1.0       V10 0.3695098

# return one row per replicate * scenario with replicates assigned to a quantile
estimateNationalRate(distScen, cfSamps$coefSamples, cfSamps$coefValues,
            "Johnson", "recruitment", ignorePrecision = TRUE, 
            returnSample = TRUE, 
            quantilesToUse = quantile(x = c(0, 1),
                                      probs = seq(0.025, 0.975, length.out = 10)))
#>     scnID Total_dist replicate     value
#> 1       1        0.1        V1 0.4021340
#> 2       1        0.1        V5 0.3977104
#> 3       1        0.1        V9 0.4107341
#> 4       1        0.1        V3 0.3758944
#> 5       1        0.1        V4 0.4040427
#> 6       1        0.1        V8 0.3654231
#> 7       1        0.1        V2 0.4038891
#> 8       1        0.1        V6 0.3680588
#> 9       1        0.1        V7 0.3959103
#> 10      1        0.1       V10 0.3749178
#> 11      2        0.2        V7 0.3952284
#> 12      2        0.2        V8 0.3648700
#> 13      2        0.2        V9 0.4100781
#> 14      2        0.2       V10 0.3743130
#> 15      2        0.2        V1 0.4016009
#> 16      2        0.2        V5 0.3971433
#> 17      2        0.2        V2 0.4032353
#> 18      2        0.2        V3 0.3753482
#> 19      2        0.2        V4 0.4035132
#> 20      2        0.2        V6 0.3674153
#> 21      3        0.3        V4 0.4029845
#> 22      3        0.3        V5 0.3965771
#> 23      3        0.3        V3 0.3748028
#> 24      3        0.3        V7 0.3945476
#> 25      3        0.3        V8 0.3643178
#> 26      3        0.3        V9 0.4094231
#> 27      3        0.3        V6 0.3667729
#> 28      3        0.3       V10 0.3737092
#> 29      3        0.3        V1 0.4010686
#> 30      3        0.3        V2 0.4025826
#> 31      4        0.4        V1 0.4005369
#> 32      4        0.4        V9 0.4087692
#> 33      4        0.4        V3 0.3742582
#> 34      4        0.4        V4 0.4024564
#> 35      4        0.4        V5 0.3960116
#> 36      4        0.4        V2 0.4019309
#> 37      4        0.4        V6 0.3661317
#> 38      4        0.4        V7 0.3938679
#> 39      4        0.4        V8 0.3637663
#> 40      4        0.4       V10 0.3731064
#> 41      5        0.5        V8 0.3632157
#> 42      5        0.5        V9 0.4081163
#> 43      5        0.5       V10 0.3725046
#> 44      5        0.5        V1 0.4000059
#> 45      5        0.5        V5 0.3954470
#> 46      5        0.5        V2 0.4012803
#> 47      5        0.5        V3 0.3737143
#> 48      5        0.5        V4 0.4019291
#> 49      5        0.5        V6 0.3654915
#> 50      5        0.5        V7 0.3931895
#> 51      6        0.6        V4 0.4014024
#> 52      6        0.6        V5 0.3948831
#> 53      6        0.6        V7 0.3925122
#> 54      6        0.6        V8 0.3626660
#> 55      6        0.6        V9 0.4074645
#> 56      6        0.6        V6 0.3648525
#> 57      6        0.6       V10 0.3719037
#> 58      6        0.6        V1 0.3994756
#> 59      6        0.6        V2 0.4006307
#> 60      6        0.6        V3 0.3731713
#> 61      7        0.7        V1 0.3989461
#> 62      7        0.7        V3 0.3726291
#> 63      7        0.7        V4 0.4008764
#> 64      7        0.7        V5 0.3943201
#> 65      7        0.7        V9 0.4068137
#> 66      7        0.7        V6 0.3642146
#> 67      7        0.7        V7 0.3918360
#> 68      7        0.7        V8 0.3621170
#> 69      7        0.7        V2 0.3999822
#> 70      7        0.7       V10 0.3713038
#> 71      8        0.8        V9 0.4061640
#> 72      8        0.8       V10 0.3707048
#> 73      8        0.8        V1 0.3984172
#> 74      8        0.8        V5 0.3937578
#> 75      8        0.8        V2 0.3993347
#> 76      8        0.8        V3 0.3720876
#> 77      8        0.8        V4 0.4003511
#> 78      8        0.8        V8 0.3615689
#> 79      8        0.8        V6 0.3635779
#> 80      8        0.8        V7 0.3911611
#> 81      9        0.9        V5 0.3931964
#> 82      9        0.9        V7 0.3904873
#> 83      9        0.9        V8 0.3610217
#> 84      9        0.9        V9 0.4055153
#> 85      9        0.9        V6 0.3629422
#> 86      9        0.9       V10 0.3701068
#> 87      9        0.9        V1 0.3978890
#> 88      9        0.9        V2 0.3986883
#> 89      9        0.9        V3 0.3715469
#> 90      9        0.9        V4 0.3998265
#> 91     10        1.0        V1 0.3973616
#> 92     10        1.0        V4 0.3993026
#> 93     10        1.0        V5 0.3926357
#> 94     10        1.0        V9 0.4048676
#> 95     10        1.0        V3 0.3710070
#> 96     10        1.0        V7 0.3898146
#> 97     10        1.0        V8 0.3604752
#> 98     10        1.0        V2 0.3980430
#> 99     10        1.0        V6 0.3623077
#> 100    10        1.0       V10 0.3695098


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
#> 1       0                0          0 0.8757906 0.04369744 0.7994686 0.9450013
#> 2       0               20         20 0.8757906 0.04369744 0.7994686 0.9450013
#> 3       0               40         40 0.8757906 0.04369744 0.7994686 0.9450013
#> 4       0               60         60 0.8757906 0.04369744 0.7994686 0.9450013
#> 5      20                0         20 0.8617131 0.04572794 0.7853844 0.9364669
#> 6      20               20         40 0.8617131 0.04572794 0.7853844 0.9364669
#> 7      20               40         60 0.8617131 0.04572794 0.7853844 0.9364669
#> 8      20               60         80 0.8617131 0.04572794 0.7853844 0.9364669
#> 9      40                0         40 0.8478591 0.04776427 0.7716698 0.9278137
#> 10     40               20         60 0.8478591 0.04776427 0.7716698 0.9278137
#> 11     40               40         80 0.8478591 0.04776427 0.7716698 0.9278137
#> 12     40               60        100 0.8478591 0.04776427 0.7716698 0.9278137
#> 13     60                0         60 0.8342249 0.04978991 0.7582969 0.9190732
#> 14     60               20         80 0.8342249 0.04978991 0.7582969 0.9190732
#> 15     60               40        100 0.8342249 0.04978991 0.7582969 0.9190732
#> 16     60               60        120 0.8342249 0.04978991 0.7582969 0.9190732
#> 17     80                0         80 0.8208071 0.05179228 0.7452423 0.9102706
#> 18     80               20        100 0.8208071 0.05179228 0.7452423 0.9102706
#> 19     80               40        120 0.8208071 0.05179228 0.7452423 0.9102706
#> 20     80               60        140 0.8208071 0.05179228 0.7452423 0.9102706
#>         R_bar   R_stdErr      R_PIlow  R_PIhigh
#> 1  0.35951478 0.11796188 0.1360894649 0.5180345
#> 2  0.30574618 0.12501733 0.0907557804 0.4831640
#> 3  0.26001915 0.12878766 0.0588554917 0.4508210
#> 4  0.22113100 0.12984040 0.0367937600 0.4208508
#> 5  0.25589195 0.10474395 0.0752409832 0.4174192
#> 6  0.21762106 0.10753920 0.0480673020 0.3899197
#> 7  0.18507391 0.10832242 0.0294667568 0.3644638
#> 8  0.15739448 0.10747006 0.0171177012 0.3409029
#> 9  0.18213629 0.09146945 0.0389430842 0.3382062
#> 10 0.15489621 0.09178992 0.0233533462 0.3165987
#> 11 0.13173012 0.09090550 0.0131814485 0.2965941
#> 12 0.11202872 0.08906811 0.0068723918 0.2780669
#> 13 0.12963921 0.07894274 0.0182917215 0.2759452
#> 14 0.11025053 0.07791293 0.0099950165 0.2589329
#> 15 0.09376159 0.07618913 0.0049938118 0.2431584
#> 16 0.07973872 0.07393949 0.0022181403 0.2285205
#> 17 0.09227334 0.06756839 0.0074481637 0.2268421
#> 18 0.07847305 0.06588406 0.0035467162 0.2133660
#> 19 0.06673672 0.06382588 0.0014790524 0.2008364
#> 20 0.05675565 0.06150132 0.0005185175 0.1891738
```
