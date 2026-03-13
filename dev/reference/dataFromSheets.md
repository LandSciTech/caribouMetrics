# Demographic data from Google sheet

Download data from Google sheets stored in bboutools format. See the
[bboutools
website](https://poissonconsulting.github.io/bboutools/articles/bboutools.html#providing-data)
and this [template
sheet](https://docs.google.com/spreadsheets/d/1i53nQrJXgrq3B6jO0ATHhSIbibtLq5TmmFL-PxGQNm8/edit?usp=sharing)
for more details on the format of the data.

## Usage

``` r
dataFromSheets(survey_url, shiny_progress = FALSE, i18n = NULL)
```

## Arguments

- survey_url:

  character. Google Sheet url.

- shiny_progress:

  logical. Is this inside a shiny app and called with
  [`shiny::withProgress`](https://rdrr.io/pkg/shiny/man/withProgress.html)?

- i18n:

  shiny.i18n translator, or NULL

## Value

a list containing:

- survey_surv: survey data for survival

- survey_recruit: survey data for recruitment

- N0: Initial population data

- pops_run: a vector of population names that are present in all the
  sheets

- dat_desc: description of the survey data for use in the app.

## See also

Caribou demography functions:
[`bayesianScenariosWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianScenariosWorkflow.md),
[`bayesianTrajectoryWorkflow()`](https://landscitech.github.io/caribouMetrics/dev/reference/bayesianTrajectoryWorkflow.md),
[`betaNationalPriors()`](https://landscitech.github.io/caribouMetrics/dev/reference/betaNationalPriors.md),
[`caribouPopGrowth()`](https://landscitech.github.io/caribouMetrics/dev/reference/caribouPopGrowth.md),
[`compareTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/compareTrajectories.md),
[`compositionBiasCorrection()`](https://landscitech.github.io/caribouMetrics/dev/reference/compositionBiasCorrection.md),
[`convertTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md),
[`demographicProjectionApp()`](https://landscitech.github.io/caribouMetrics/dev/reference/demographicProjectionApp.md),
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md),
[`estimateNationalRate()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateNationalRates.md),
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
 
dataFromSheets("https://docs.google.com/spreadsheets/d/1i53nQrJXgrq3B6jO0ATHhSIbibtLq5TmmFL-PxGQNm8/edit?usp=sharing")
#> ✔ Reading from bbou_example_data.
#> ✔ Range ''recruit_data''.
#> ✔ Reading from bbou_example_data.
#> ✔ Range ''surv_data''.
#> ✔ Reading from bbou_example_data.
#> ✔ Range ''population_estimates''.
#> ✔ Reading from bbou_example_data.
#> ✔ Range ''data_description''.
#> $survey_surv
#> # A tibble: 570 × 6
#>    PopulationName  Year Month StartTotal MortalitiesCertain MortalitiesUncertain
#>    <chr>          <dbl> <dbl>      <dbl>              <dbl>                <dbl>
#>  1 A               1986     1          0                  0                    0
#>  2 A               1986     2          8                  0                    0
#>  3 A               1986     3          8                  0                    0
#>  4 A               1986     4          8                  0                    0
#>  5 A               1986     5          8                  0                    0
#>  6 A               1986     6          8                  0                    0
#>  7 A               1986     7          8                  0                    0
#>  8 A               1986     8          8                  0                    0
#>  9 A               1986     9          8                  0                    0
#> 10 A               1986    10          8                  0                    0
#> # ℹ 560 more rows
#> 
#> $pops_run
#> [1] "A" "B"
#> 
#> $survey_recruit
#> # A tibble: 1,177 × 9
#>    PopulationName  Year Month   Day  Cows Bulls UnknownAdults Yearlings Calves
#>    <chr>          <dbl> <dbl> <dbl> <dbl> <dbl>         <dbl>     <dbl>  <dbl>
#>  1 A               1990     3     9     1     1             0         0      0
#>  2 A               1990     3     9     5     1             0         0      0
#>  3 A               1990     3     9     4     1             0         0      0
#>  4 A               1990     3     9     2     0             0         0      0
#>  5 A               1990     3     9     6     0             0         0      0
#>  6 A               1990     3     9     4     1             0         0      0
#>  7 A               1990     3     9     5     0             0         0      0
#>  8 A               1990     3     9     2     0             0         0      0
#>  9 A               1990     3     9     3     2             0         0      1
#> 10 A               1990     3     9     4     0             0         0      1
#> # ℹ 1,167 more rows
#> 
#> $N0
#> # A tibble: 2 × 4
#> # Groups:   PopulationName [2]
#>   PopulationName  Year FemalePopulationLower FemalePopulationUpper
#>   <chr>          <dbl>                 <dbl>                 <dbl>
#> 1 A               2020                   450                   550
#> 2 B               2021                   975                  1025
#> 
#> $dat_desc
#> # A tibble: 1 × 4
#>   Description_en                  Description_en2 Description_fr Description_fr2
#>   <chr>                           <chr>           <chr>          <chr>          
#> 1 "This is simulated example dat… "This is simul… "Il s'agit d'… "Il s'agit d'u…
#> 
```
