
<!-- README.md is generated from README.Rmd. Please edit that file -->

# caribouMetrics

<!-- badges: start -->
<!-- badges: end -->

caribouMetrics is built off of two different models, one calculates
predictors described in Table 52 of Environment Canada (2011) Scientific
Assessment to Inform the Identification of Critical Habitat for Woodland
Caribou and the other implements the caribou resource selection
probability functions described in Hornseth and Rempel (2016) Seasonal
resource selection of woodland caribou (*Rangifer tarandus caribou*)
across a gradient of anthropogenic disturbance.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LandSciTech/caribouMetrics")
```

## Example

This is a basic example which shows you how to use the main functions:

``` r
library(caribouMetrics)

pthBase <- system.file("extdata", package = "caribouMetrics")

# load example data
landCoverD <- raster::raster(file.path(pthBase, "landCover.tif")) 
  # convert PLC classes to resource types used in the model 
landCoverD <- reclassPLC(landCoverD)
eskerDras <- raster::raster(file.path(pthBase, "eskerTif.tif"))
eskerDshp <- sf::read_sf(file.path(pthBase, "esker.shp"))
natDistD <- raster::raster(file.path(pthBase, "natDist.tif"))
anthroDistD <-raster::raster(file.path(pthBase, "anthroDist.tif"))
linFeatDras <- raster::raster(file.path(pthBase, "linFeatTif.tif"))
projectPolyD <- sf::read_sf(file.path(pthBase, "projectPoly.shp"))

# Calculate habitat use
carHab1 <- caribouHabitat(
  landCover = landCoverD,
  esker = eskerDras, 
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  caribouRange = "Churchill"
)
#> cropping landCover to extent of projectPoly
#> cropping natDist to extent of projectPoly
#> cropping anthroDist to extent of projectPoly
#> cropping esker to extent of projectPoly
#> cropping linFeat to extent of projectPoly
#> resampling linFeat to match landCover resolution
#> resampling esker to match landCover resolution
#> Applying moving window.

# plot the results
plot(carHab1)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
# calculate disturbance 
disturb <- disturbanceMetrics(landCover = landCoverD,
                              linFeat = linFeatDras,  
                              natDist = natDistD,
                              projectPoly = projectPolyD)
#> cropping landCover to extent of projectPoly
#> cropping natDist to extent of projectPoly
#> cropping linFeat to extent of projectPoly
#> buffering anthropogenic disturbance
#> calculating disturbance metrics

results(disturb)
#>   zone    Anthro       Fire Total_dist fire_excl_anthro FID
#> 1    1 0.3997933 0.01734581  0.4056555      0.005862182   0
```

More detailed examples are provided in the vignettes.
