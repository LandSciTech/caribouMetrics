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
#>  [1,] -0.1375893 -0.0008214780  66.39165
#>  [2,] -0.1503580 -0.0006914206  54.24947
#>  [3,] -0.1218574 -0.0007674631  57.26584
#>  [4,] -0.1482105 -0.0008157216  68.17130
#>  [5,] -0.1491151 -0.0007580918  52.06607
#>  [6,] -0.1599969 -0.0007487956  67.75248
#>  [7,] -0.1351537 -0.0008658241  58.66916
#>  [8,] -0.1331251 -0.0010378301  67.73499
#>  [9,] -0.1492377 -0.0008264594  55.38230
#> [10,] -0.1529072 -0.0007537961  47.33583
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
#>  [1] 0.0250000 0.5527778 0.7638889 0.8694444 0.6583333 0.3416667 0.4472222
#>  [8] 0.1305556 0.9750000 0.2361111
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept      Anthro fire_excl_anthro Precision
#>  [1,] -1.0404962 -0.01540243     -0.006966033  19.12575
#>  [2,] -1.0573687 -0.01691326     -0.010830153  16.33570
#>  [3,] -1.0530327 -0.01737363     -0.009945861  15.73955
#>  [4,] -1.0875604 -0.01772134     -0.006544805  19.36795
#>  [5,] -1.0322034 -0.01792035     -0.007450377  16.09558
#>  [6,] -0.9395177 -0.01614165     -0.006679011  16.55621
#>  [7,] -0.9997499 -0.01600013     -0.008019057  17.44081
#>  [8,] -1.0720669 -0.01574681     -0.007152510  21.85707
#>  [9,] -0.8860691 -0.01798761     -0.008534226  18.16164
#> [10,] -0.9961209 -0.01862720     -0.011257839  22.06328
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
#>  [1] 0.2361111 0.1305556 0.9750000 0.4472222 0.5527778 0.0250000 0.8694444
#>  [8] 0.7638889 0.3416667 0.6583333
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
#>  [1,] -0.1288397 -0.0008848235  64.50029
#>  [2,] -0.1410068 -0.0008226110  64.21693
#>  [3,] -0.1333665 -0.0009005080  62.81746
#>  [4,] -0.1458198 -0.0007836418  65.47641
#>  [5,] -0.1344640 -0.0010414647  59.65608
#>  [6,] -0.1456477 -0.0007662172  67.96761
#>  [7,] -0.1529465 -0.0008754472  60.23908
#>  [8,] -0.1540877 -0.0007116626  62.64547
#>  [9,] -0.1387897 -0.0006546201  71.15996
#> [10,] -0.1586024 -0.0006191131  76.07244
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
#>  [1] 0.7638889 0.2361111 0.3416667 0.8694444 0.6583333 0.4472222 0.9750000
#>  [8] 0.5527778 0.1305556 0.0250000
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9529372 -0.01476715
#>  [2,] -1.0035136 -0.01686201
#>  [3,] -1.0352183 -0.01637935
#>  [4,] -0.8694776 -0.01524407
#>  [5,] -0.9120334 -0.01261418
#>  [6,] -0.9999643 -0.01574794
#>  [7,] -1.0132053 -0.01550719
#>  [8,] -0.9548055 -0.01561502
#>  [9,] -1.0134729 -0.01434705
#> [10,] -0.9039258 -0.01658503
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
#>  [1] 0.0250000 0.6583333 0.9750000 0.2361111 0.5527778 0.1305556 0.4472222
#>  [8] 0.8694444 0.7638889 0.3416667
#> 
#> 

cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

sampleNationalCoefs(cfs[[1]], 10)
#> $coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.8510095 -0.01392656
#>  [2,] -0.9163447 -0.01547171
#>  [3,] -0.8909582 -0.01237876
#>  [4,] -0.9664610 -0.01699347
#>  [5,] -0.9859341 -0.01645362
#>  [6,] -0.9305850 -0.01157922
#>  [7,] -0.9578151 -0.01313435
#>  [8,] -0.9193104 -0.01607172
#>  [9,] -0.9142721 -0.01634361
#> [10,] -0.9983274 -0.01510461
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
