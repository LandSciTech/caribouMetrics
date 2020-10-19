
<!-- README.md is generated from README.Rmd. Please edit that file -->

# caribouMetrics

<!-- badges: start -->

<!-- badges: end -->

caribouMetrics is built off of two different models, one calculates
predictors described in Table 52 of Environment Canada (2011) Scientific
Assessment to Inform the Identification of Critical Habitat for Woodland
Caribou and the other implements the caribou resource selection
probability functions described in Rempel, R. and M. Hornseth. 2018.
Range-specific seasonal resource selection probability functions for 13
caribou ranges in Northern Ontario. The second part of the package is
the actively developed part and the rest of the README focuses on it.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LandSciTech/caribouMetrics")
```

## Example

This is a basic example which shows you how to use the primary function:

``` r
#library(caribouMetrics)
devtools::load_all()
#> Loading caribouMetrics
#> Skipping missing files: lstools.R
#> Adding files missing in collate: caribouMetrics.R
pthBase <- "tests/testthat/data/"

# load example data
plcD = raster(paste0(pthBase, "plc", ".tif"))
eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE)
friD = raster(paste0(pthBase, "fri", ".tif"))
ageD = raster(paste0(pthBase, "age", ".tif"))
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
harvD = raster(paste0(pthBase, "harv", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE)

# Calaulate habitat use
carHab1 <- caribouHabitat(
  plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = projectPolyD,
  friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE) %>%
    mutate(RFU = toupper(RFU) %>% stringr::str_replace("HRDMW", "HRDMX")),
  caribouRange = "Churchill", 
  winArea = 500
)
#> projectPoly being transformed to have crs matching plc
#> cropping plc to extent of projectPoly
#> extending plc to extent of projectPoly
#> cropping esker to extent of plc
#> extending esker to extent of plc
#> cropping linFeat to extent of plc
#> extending linFeat to extent of plc
#> Warning in .checkLevels(levels(x)[[1]], value): the number of rows in the raster
#> attributes (factors) data.frame is unexpected
#> Warning in sort(newv[, 1]) == sort(old[, 1]): longer object length is not a
#> multiple of shorter object length
#> Warning in .checkLevels(levels(x)[[1]], value): the values in the "ID" column in
#> the raster attributes (factors) data.frame have changed

# plot the results
plot(carHab1)
```

<img src="man/figures/README-example-1.png" width="100%" />

More detailed examples are provided in the
vignettes/overall\_vignette.Rmd and in scripts/notesOnHandRModel.Rmd
