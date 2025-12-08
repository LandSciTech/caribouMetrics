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
#>  [1,] -0.1439429 -0.0009774425  56.09693
#>  [2,] -0.1457064 -0.0008004473  61.48855
#>  [3,] -0.1422335 -0.0007659977  73.04347
#>  [4,] -0.1409010 -0.0008591432  64.20958
#>  [5,] -0.1367324 -0.0008372087  66.41723
#>  [6,] -0.1316502 -0.0011589137  60.57292
#>  [7,] -0.1257453 -0.0006491190  60.94579
#>  [8,] -0.1528124 -0.0007745746  60.14864
#>  [9,] -0.1393779 -0.0007738779  55.14949
#> [10,] -0.1487257 -0.0005469178  46.64238
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
#>  [1] 0.3416667 0.4472222 0.0250000 0.8694444 0.6583333 0.7638889 0.2361111
#>  [8] 0.9750000 0.1305556 0.5527778
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept      Anthro fire_excl_anthro Precision
#>  [1,] -1.0886241 -0.01711594     -0.007456206  22.75677
#>  [2,] -0.9739582 -0.01540915     -0.005941910  19.11029
#>  [3,] -1.0004546 -0.02085116     -0.013454155  21.11714
#>  [4,] -1.0143668 -0.01769882     -0.007502312  19.03260
#>  [5,] -1.1369351 -0.01789710     -0.011098517  20.17579
#>  [6,] -0.9993384 -0.01626498     -0.007808392  17.59848
#>  [7,] -0.9486961 -0.01673297     -0.010302625  16.99404
#>  [8,] -1.0232608 -0.01733236     -0.010129449  21.93867
#>  [9,] -1.0817212 -0.01823328     -0.009319808  19.87137
#> [10,] -1.0402215 -0.01794420     -0.009549139  21.13804
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
#>  [1] 0.3416667 0.6583333 0.8694444 0.9750000 0.5527778 0.0250000 0.7638889
#>  [8] 0.1305556 0.4472222 0.2361111
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
#>  [1,] -0.1496414 -0.0009832711  62.49886
#>  [2,] -0.1358365 -0.0005234796  59.34979
#>  [3,] -0.1310439 -0.0008973436  74.66406
#>  [4,] -0.1393536 -0.0006621852  61.54484
#>  [5,] -0.1431382 -0.0007348578  61.65064
#>  [6,] -0.1416100 -0.0009417600  58.19539
#>  [7,] -0.1532636 -0.0007634347  52.50444
#>  [8,] -0.1510141 -0.0008534174  66.37139
#>  [9,] -0.1316167 -0.0009233339  65.85342
#> [10,] -0.1425711 -0.0007484273  80.88584
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
#>  [1] 0.3416667 0.1305556 0.7638889 0.5527778 0.9750000 0.4472222 0.6583333
#>  [8] 0.8694444 0.0250000 0.2361111
#> 
#> 
#> $coefSamples_Recruitment
#> $coefSamples_Recruitment$coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9704581 -0.01564570
#>  [2,] -0.9845122 -0.01324768
#>  [3,] -0.9942972 -0.01686176
#>  [4,] -0.9827596 -0.01899211
#>  [5,] -0.9467670 -0.01399213
#>  [6,] -0.9419558 -0.01518320
#>  [7,] -0.9378449 -0.01353034
#>  [8,] -0.8604821 -0.01399528
#>  [9,] -0.9624754 -0.01427855
#> [10,] -0.9985000 -0.01656535
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
#>  [1] 0.2361111 0.3416667 0.8694444 0.5527778 0.7638889 0.1305556 0.6583333
#>  [8] 0.0250000 0.9750000 0.4472222
#> 
#> 

cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")

sampleNationalCoefs(cfs[[1]], 10)
#> $coefSamples
#>        Intercept  Total_dist
#>  [1,] -0.9929092 -0.01434329
#>  [2,] -0.9828987 -0.01555055
#>  [3,] -1.0648389 -0.01242014
#>  [4,] -1.0378557 -0.01532705
#>  [5,] -0.9390160 -0.01576930
#>  [6,] -1.0501705 -0.01682715
#>  [7,] -0.9090130 -0.01558927
#>  [8,] -1.0157227 -0.01632534
#>  [9,] -1.0050548 -0.01434024
#> [10,] -0.8843576 -0.01340417
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
