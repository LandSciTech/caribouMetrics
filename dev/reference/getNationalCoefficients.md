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
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
[`getScenarioDefaults()`](https://landscitech.github.io/caribouMetrics/dev/reference/getScenarioDefaults.md),
[`plotCompareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotCompareTrajectories.md),
[`plotSurvivalSeries()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotSurvivalSeries.md),
[`plotTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/plotTrajectories.md),
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
#>  [1,] -0.1480462 -0.0010438585  58.03616
#>  [2,] -0.1430032 -0.0009454453  59.33736
#>  [3,] -0.1489133 -0.0006387924  51.67832
#>  [4,] -0.1277558 -0.0010612015  68.16593
#>  [5,] -0.1389300 -0.0005665318  73.33024
#>  [6,] -0.1289576 -0.0007483172  72.31326
#>  [7,] -0.1363444 -0.0007572366  51.84803
#>  [8,] -0.1458662 -0.0005932405  38.26306
#>  [9,] -0.1537381 -0.0007009531  68.11915
#> [10,] -0.1505659 -0.0006804294  64.74954
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
#>  [1] 0.8694444 0.5527778 0.4472222 0.6583333 0.0250000 0.2361111 0.7638889
#>  [8] 0.9750000 0.3416667 0.1305556
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept      Anthro fire_excl_anthro Precision
#>  [1,] -1.1863976 -0.01765815     -0.005496848  21.78181
#>  [2,] -1.0776251 -0.01571174     -0.006672900  19.46246
#>  [3,] -1.0826485 -0.01586546     -0.006026198  21.63728
#>  [4,] -0.9957944 -0.01717795     -0.008396980  22.25531
#>  [5,] -1.0013989 -0.01898397     -0.007321524  21.02224
#>  [6,] -1.0788348 -0.01839223     -0.009030325  22.42766
#>  [7,] -1.0554040 -0.01879328     -0.007762457  17.97300
#>  [8,] -1.0626236 -0.01619433     -0.007952730  21.47961
#>  [9,] -1.1299056 -0.01548108     -0.009833262  18.55831
#> [10,] -0.9719267 -0.01555126     -0.009716902  19.85290
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
#>  [1] 0.2361111 0.5527778 0.3416667 0.0250000 0.8694444 0.1305556 0.7638889
#>  [8] 0.9750000 0.6583333 0.4472222
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
#>  [1,] -0.1356785 -0.0009017750  74.34194
#>  [2,] -0.1344399 -0.0006684816  66.72762
#>  [3,] -0.1325078 -0.0007554310  60.84536
#>  [4,] -0.1378479 -0.0005618761  61.78707
#>  [5,] -0.1397351 -0.0008536154  56.36560
#>  [6,] -0.1302149 -0.0008946659  76.43922
#>  [7,] -0.1416202 -0.0009844278  63.44977
#>  [8,] -0.1406308 -0.0007403799  72.98995
#>  [9,] -0.1517740 -0.0007090066  68.65632
#> [10,] -0.1404511 -0.0006406354  54.97784
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
#>  [1] 0.0250000 0.9750000 0.3416667 0.7638889 0.8694444 0.6583333 0.2361111
#>  [8] 0.1305556 0.4472222 0.5527778
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.8646352 -0.01603464
#>  [2,] -0.9369916 -0.01371640
#>  [3,] -0.9418396 -0.01400316
#>  [4,] -0.9622292 -0.01435096
#>  [5,] -0.8909215 -0.01696744
#>  [6,] -0.9320926 -0.01452815
#>  [7,] -0.9017412 -0.01693457
#>  [8,] -1.0190819 -0.01639863
#>  [9,] -1.0354865 -0.01332016
#> [10,] -0.9700190 -0.01636591
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
#>  [1] 0.8694444 0.3416667 0.0250000 0.2361111 0.7638889 0.4472222 0.5527778
#>  [8] 0.9750000 0.6583333 0.1305556
#> 
#> 

cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

sampleNationalCoefs(cfs[[1]], 10)
#> $coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9735333 -0.01528894
#>  [2,] -0.8303363 -0.01489409
#>  [3,] -0.8892294 -0.01371337
#>  [4,] -0.9528565 -0.01583717
#>  [5,] -0.8973713 -0.01420578
#>  [6,] -0.9853699 -0.01458906
#>  [7,] -0.9713435 -0.01289383
#>  [8,] -1.0521418 -0.01462253
#>  [9,] -0.9867926 -0.01675020
#> [10,] -0.8718322 -0.01749507
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
