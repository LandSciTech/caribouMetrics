# Sample demographic regression model coefficients

Select the regression coefficient values and standard errors for the
desired model version (see `popGrowthTableJohnsonECCC` for options) and
then sample from the Gaussian distribution for each replicate
population. `getNationalCoefficients` is a wrapper around
`subsetNationalCoefs()`, which selects coefficients and
`sampleNationalCoefs()`, which samples coefficients, for both the
survival and recruitment models.

## Usage

``` r
getNationalCoefficients(
  replicates,
  modelVersion = "Johnson",
  survivalModelNumber = "M1",
  recruitmentModelNumber = "M4",
  useQuantiles = TRUE,
  populationGrowthTable = popGrowthTableJohnsonECCC
)

sampleNationalCoefs(coefTable, replicates)

subsetNationalCoefs(populationGrowthTable, resVar, modelVersion, modNum)
```

## Arguments

- replicates:

  integer. Number of replicate populations.

- modelVersion:

  character. Which model version to use. Currently the only option is
  "Johnson" for the model used in Johnson et. al. (2020), but additional
  options may be added in the future.

- survivalModelNumber, recruitmentModelNumber:

  character. Which model number to use see
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  for options.

- useQuantiles:

  logical or numeric. If it is a numeric vector it must be length 2 and
  give the low and high limits of the quantiles to use. If
  `useQuantiles != FALSE`, each replicate population is assigned to a
  quantile of the distribution of variation around the expected values,
  and remains in that quantile as covariates change. If
  `useQuantiles = TRUE`, replicate populations will be assigned to
  quantiles in the default range of 0.025 and 0.975.

- populationGrowthTable:

  data.frame.[popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  is included in the package and should be used in most cases. A custom
  table of model coefficients and standard errors or confidence
  intervals can be provided but it must match the column names of
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md).
  If the table does not contain the standard error it is calculated from
  the confidence interval.

- coefTable:

  data.table. Table must have columns "Coefficient" for the name of the
  coefficient, "Value" for the value of the coefficient and "StdErr" for
  the standard error of coefficients. Typically created with
  `subsetNationalCoefs()`

- resVar:

  character. Response variable, typically "femaleSurvival" or
  "recruitment"

- modNum:

  character vector. Which model number(s) to use see
  [popGrowthTableJohnsonECCC](https://landscitech.github.io/caribouMetrics/dev/reference/popGrowthTableJohnsonECCC.md)
  for typical options.

## Value

For `getNationalCoefficients` a list with elements:

- "modelVersion": The name of the model version

- "coefSamples_Survival" and"coefSamples_Recruitment": lists with
  elements:

  - "coefSamples": Bootstrapped coefficients with `replicates` rows

  - "coefValues": Coefficient values taken from `populationGrowthTable`

  - "quantiles": A vector of randomly selected quantiles between 0.025
    and 0.975 with length `replicates`

For `sampleNationalCoefs` a list with elements:

- "coefSamples": Bootstrapped coefficients with `replicates` rows

- "coefValues": Coefficient values taken from `populationGrowthTable`

For `subsetNationalCoefs`: a named list with one element per model
version. The names are `modelVersion_modNum_Type`. Each element contains
a data.frame that is a subset of `populationGrowthTable` for the
selected model

## Details

Each population is optionally assigned to quantiles of the error
distributions for survival and recruitment. Using quantiles means that
the population will stay in these quantiles as disturbance changes over
time, so there is persistent variation in recruitment and survival among
example populations. See
[`estimateNationalRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md)
for more details.

## References

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
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
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
# sample coefficients for default models
getNationalCoefficients(10)
#> $modelVersion
#> [1] "Johnson"
#> 
#> $coefSamples_Survival
#> $coefSamples_Survival$coefSamples
#>        Intercept        Anthro Precision
#>  [1,] -0.1392748 -0.0007621126  68.83647
#>  [2,] -0.1393385 -0.0010255801  60.36325
#>  [3,] -0.1441847 -0.0007709322  57.79665
#>  [4,] -0.1426833 -0.0008070736  56.62803
#>  [5,] -0.1394975 -0.0009195023  57.07883
#>  [6,] -0.1407351 -0.0010538044  65.53183
#>  [7,] -0.1491469 -0.0007292385  73.02784
#>  [8,] -0.1468109 -0.0008512666  55.72058
#>  [9,] -0.1458068 -0.0008377470  75.22492
#> [10,] -0.1379537 -0.0009721341  71.22320
#> 
#> $coefSamples_Survival$coefValues
#>    Intercept Anthro Precision
#>        <num>  <num>     <num>
#> 1:    -0.142 -8e-04  63.43724
#> 
#> $coefSamples_Survival$coefStdErrs
#>      Intercept      Anthro Precision
#>          <num>       <num>     <num>
#> 1: 0.007908163 0.000127551  8.272731
#> 
#> $coefSamples_Survival$quantiles
#>  [1] 0.5527778 0.9750000 0.6583333 0.8694444 0.2361111 0.7638889 0.1305556
#>  [8] 0.3416667 0.4472222 0.0250000
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept      Anthro fire_excl_anthro Precision
#>  [1,] -0.9807072 -0.01537949     -0.006732436  17.47412
#>  [2,] -1.0966455 -0.01545213     -0.009292509  21.38604
#>  [3,] -0.9590578 -0.01439435     -0.006279100  21.94885
#>  [4,] -1.0068339 -0.01518384     -0.009792803  23.46515
#>  [5,] -1.0220588 -0.01925490     -0.010436944  18.33165
#>  [6,] -0.7821580 -0.01428498     -0.009031424  20.15949
#>  [7,] -0.9509732 -0.01609041     -0.009460336  22.38481
#>  [8,] -1.0123316 -0.01814691     -0.008819888  24.96329
#>  [9,] -1.0228535 -0.01740828     -0.009686527  21.05482
#> [10,] -1.0424280 -0.01662056     -0.008219040  16.38331
#> 
#> $coefSamples_Recruitment$coefValues
#>    Intercept Anthro fire_excl_anthro Precision
#>        <num>  <num>            <num>     <num>
#> 1:    -1.023 -0.017          -0.0081  19.86189
#> 
#> $coefSamples_Recruitment$coefStdErrs
#>     Intercept      Anthro fire_excl_anthro Precision
#>         <num>       <num>            <num>     <num>
#> 1: 0.06122449 0.001530612      0.002040816  2.228655
#> 
#> $coefSamples_Recruitment$quantiles
#>  [1] 0.0250000 0.4472222 0.1305556 0.5527778 0.3416667 0.8694444 0.7638889
#>  [8] 0.9750000 0.6583333 0.2361111
#> 
#> 

# try a different model
getNationalCoefficients(10, modelVersion = "Johnson", survivalModelNumber = "M1",
                        recruitmentModelNumber = "M3")
#> $modelVersion
#> [1] "Johnson"
#> 
#> $coefSamples_Survival
#> $coefSamples_Survival$coefSamples
#>        Intercept        Anthro Precision
#>  [1,] -0.1569708 -0.0008582988  69.30472
#>  [2,] -0.1399055 -0.0007301496  57.57007
#>  [3,] -0.1466777 -0.0008493100  55.80301
#>  [4,] -0.1365231 -0.0008122076  63.59664
#>  [5,] -0.1329864 -0.0006809292  55.76730
#>  [6,] -0.1307850 -0.0006051873  70.38670
#>  [7,] -0.1409838 -0.0007936979  64.69578
#>  [8,] -0.1412547 -0.0008977646  53.37333
#>  [9,] -0.1425925 -0.0009630006  55.98207
#> [10,] -0.1400507 -0.0006219704  62.11808
#> 
#> $coefSamples_Survival$coefValues
#>    Intercept Anthro Precision
#>        <num>  <num>     <num>
#> 1:    -0.142 -8e-04  63.43724
#> 
#> $coefSamples_Survival$coefStdErrs
#>      Intercept      Anthro Precision
#>          <num>       <num>     <num>
#> 1: 0.007908163 0.000127551  8.272731
#> 
#> $coefSamples_Survival$quantiles
#>  [1] 0.0250000 0.9750000 0.3416667 0.1305556 0.7638889 0.8694444 0.4472222
#>  [8] 0.5527778 0.2361111 0.6583333
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.8593744 -0.01655265
#>  [2,] -0.9862915 -0.01699925
#>  [3,] -0.9765410 -0.01379829
#>  [4,] -0.9809084 -0.01671771
#>  [5,] -0.9295554 -0.01524384
#>  [6,] -1.0201937 -0.01690986
#>  [7,] -1.0399455 -0.01211002
#>  [8,] -0.9505392 -0.01477119
#>  [9,] -0.9523207 -0.01291817
#> [10,] -0.9120106 -0.01283717
#> 
#> $coefSamples_Recruitment$coefValues
#>    Intercept Total_dist
#>        <num>      <num>
#> 1:    -0.956     -0.015
#> 
#> $coefSamples_Recruitment$coefStdErrs
#>    Intercept  Total_dist
#>        <num>       <num>
#> 1: 0.0619898 0.001530612
#> 
#> $coefSamples_Recruitment$quantiles
#>  [1] 0.3416667 0.2361111 0.6583333 0.0250000 0.8694444 0.9750000 0.5527778
#>  [8] 0.1305556 0.7638889 0.4472222
#> 
#> 

cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

sampleNationalCoefs(cfs[[1]], 10)
#> $coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9681444 -0.01398264
#>  [2,] -0.9573189 -0.01475347
#>  [3,] -0.9713287 -0.01776771
#>  [4,] -0.9459989 -0.01451957
#>  [5,] -0.8822363 -0.01149836
#>  [6,] -0.9055004 -0.01565601
#>  [7,] -0.9501705 -0.01289928
#>  [8,] -0.8478571 -0.01705309
#>  [9,] -1.0058990 -0.01257477
#> [10,] -0.8743282 -0.01479708
#> 
#> $coefValues
#>    Intercept Total_dist
#>        <num>      <num>
#> 1:    -0.956     -0.015
#> 
#> $coefStdErrs
#>    Intercept  Total_dist
#>        <num>       <num>
#> 1: 0.0619898 0.001530612
#> 

subsetNationalCoefs(popGrowthTableJohnsonECCC, "femaleSurvival", "Johnson", "M1")
#> $Johnson_M1_National
#>    modelVersion responseVariable ModelNumber     Type Coefficient    Value
#>          <char>           <char>      <char>   <char>      <char>    <num>
#> 1:      Johnson   femaleSurvival          M1 National   Intercept -0.14200
#> 2:      Johnson   femaleSurvival          M1 National      Anthro -0.00080
#> 3:      Johnson   femaleSurvival          M1 National   Precision 63.43724
#>         StdErr lowerCI upperCI
#>          <num>   <num>   <num>
#> 1: 0.007908163  -0.158 -0.1270
#> 2: 0.000127551  -0.001 -0.0005
#> 3: 8.272730950      NA      NA
#> 
```
