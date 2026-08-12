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
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
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
# sample coefficients for default models
getNationalCoefficients(10)
#> $modelVersion
#> [1] "Johnson"
#> 
#> $coefSamples_Survival
#> $coefSamples_Survival$coefSamples
#>        Intercept        Anthro Precision
#>  [1,] -0.1267497 -0.0007297408  61.36137
#>  [2,] -0.1434377 -0.0009339489  62.82264
#>  [3,] -0.1498346 -0.0006632000  57.69875
#>  [4,] -0.1524195 -0.0006620696  88.06442
#>  [5,] -0.1520567 -0.0007539970  65.09325
#>  [6,] -0.1431232 -0.0006641330  54.59174
#>  [7,] -0.1321505 -0.0006718344  67.11186
#>  [8,] -0.1519888 -0.0006906389  77.15316
#>  [9,] -0.1503150 -0.0007770357  59.70923
#> [10,] -0.1404523 -0.0010567603  62.95751
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
#>  [1] 0.3416667 0.1305556 0.4472222 0.8694444 0.7638889 0.0250000 0.5527778
#>  [8] 0.6583333 0.2361111 0.9750000
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept      Anthro Fire_excl_anthro Precision
#>  [1,] -1.0847249 -0.01603358     -0.010618909  23.90156
#>  [2,] -0.9998196 -0.01517987     -0.006717466  20.67937
#>  [3,] -1.0402687 -0.01723239     -0.006465975  18.53268
#>  [4,] -0.8931000 -0.01774963     -0.011040135  18.59799
#>  [5,] -1.0087908 -0.01486880     -0.008488228  20.86175
#>  [6,] -1.0445856 -0.01759048     -0.007447013  18.12417
#>  [7,] -1.1557038 -0.01607495     -0.007827396  20.52868
#>  [8,] -1.0048343 -0.01834319     -0.006450545  19.97543
#>  [9,] -1.0293981 -0.01853872     -0.006073476  22.16698
#> [10,] -0.9680017 -0.01718680     -0.009460592  24.77689
#> 
#> $coefSamples_Recruitment$coefValues
#>    Intercept Anthro Fire_excl_anthro Precision
#>        <num>  <num>            <num>     <num>
#> 1:    -1.023 -0.017          -0.0081  19.86189
#> 
#> $coefSamples_Recruitment$coefStdErrs
#>     Intercept      Anthro Fire_excl_anthro Precision
#>         <num>       <num>            <num>     <num>
#> 1: 0.06122449 0.001530612      0.002040816  2.228655
#> 
#> $coefSamples_Recruitment$quantiles
#>  [1] 0.1305556 0.8694444 0.5527778 0.7638889 0.6583333 0.4472222 0.0250000
#>  [8] 0.3416667 0.9750000 0.2361111
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
#>  [1,] -0.1394212 -0.0009143944  70.99067
#>  [2,] -0.1479595 -0.0008182044  64.00527
#>  [3,] -0.1291905 -0.0006890285  72.83197
#>  [4,] -0.1402579 -0.0007963141  67.40686
#>  [5,] -0.1382653 -0.0007643340  57.28996
#>  [6,] -0.1336491 -0.0008524628  56.33201
#>  [7,] -0.1246191 -0.0009280135  74.78511
#>  [8,] -0.1409499 -0.0008831935  72.74880
#>  [9,] -0.1410947 -0.0008077740  44.25386
#> [10,] -0.1283713 -0.0007741672  58.35136
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
#>  [1] 0.4472222 0.2361111 0.5527778 0.0250000 0.7638889 0.9750000 0.1305556
#>  [8] 0.6583333 0.3416667 0.8694444
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9735633 -0.01633396
#>  [2,] -0.9960999 -0.01425738
#>  [3,] -0.9277828 -0.01509315
#>  [4,] -1.0408267 -0.01886687
#>  [5,] -0.9629287 -0.01652454
#>  [6,] -0.9458981 -0.01747052
#>  [7,] -0.8996886 -0.01639462
#>  [8,] -0.9116151 -0.01299306
#>  [9,] -0.8469247 -0.01465965
#> [10,] -0.9684036 -0.01479619
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
#>  [1] 0.9750000 0.1305556 0.8694444 0.6583333 0.4472222 0.0250000 0.7638889
#>  [8] 0.3416667 0.5527778 0.2361111
#> 
#> 

cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

sampleNationalCoefs(cfs[[1]], 10)
#> $coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.8824942 -0.01700552
#>  [2,] -0.8138640 -0.01344958
#>  [3,] -0.9683592 -0.01597249
#>  [4,] -0.9393436 -0.01723292
#>  [5,] -1.0253017 -0.01207913
#>  [6,] -0.9623051 -0.01555135
#>  [7,] -0.9968580 -0.01422298
#>  [8,] -0.8952682 -0.01433137
#>  [9,] -0.8529288 -0.01231108
#> [10,] -0.9908045 -0.01369081
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
