# Load Spatial Input Data

Load spatial input data from file and then align the inputs to the
`projectPoly` and `refRast`. Inputs can be file names or spatial
objects.

## Usage

``` r
loadSpatialInputs(
  projectPoly,
  refRast,
  inputsList,
  convertToRast = NULL,
  convertToRastDens = NULL,
  altTemplate = NULL,
  useTemplate = NULL,
  reclassOptions = NULL,
  bufferWidth = NULL,
  ptDensity = 1,
  rastOut = "terra"
)
```

## Arguments

- projectPoly:

  character or sf. A polygon delineating the study area.

- refRast:

  character or raster. A raster which will be used as the template that
  all other rasters must align to (but see `alttemplate`).

- inputsList:

  list. A named list of inputs that are either spatial objects or file
  names of spatial objects. ".shp" is the only extension accepted for
  vector data and all other extensions will be passed to
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).
  If an element is a list these are assumed to be linear features and
  they are combined.

- convertToRast:

  character. Optional. Names of elements of `inputsList` that should be
  converted to raster after loading.

- convertToRastDens:

  character. Optional. Names of elements of `inputsList` that should be
  converted to raster line density after loading.

- altTemplate:

  raster. Optional template raster for raster inputs that can have a
  different resolution from the `refRast`.

- useTemplate:

  character. Optional. Names of elements of `inputsList` that use
  `altTemplate`.

- reclassOptions:

  list. An optional named list containing a function, a list where the
  first element is a function that takes the corresponding `inputsList`
  element as its first argument and the subsequent elements are named
  arguments for that function, or a matrix that will be passed to
  [`terra::classify()`](https://rspatial.github.io/terra/reference/classify.html).

- bufferWidth:

  numeric. The width of a moving window that will be applied to the
  data. If it is supplied a buffer of 3\*`bufferWidth` around the
  `projectPoly` is used so that rasters will be cropped to a larger
  area. This can be used to avoid edge effects in moving window
  calculations

- ptDensity:

  number. Only used if an element of `inputsList` is a list that
  contains a mixture of rasters and lines and is included in
  convertToRast. See
  [`rasterizeLineDensity()`](https://landscitech.github.io/caribouMetrics/reference/rasterizeLineDensity.md).

- rastOut:

  character. The format that rasters should be output with. "raster" for
  RasterLayers and "terra" for SpatRasters. The default is "terra".

## Value

A named list with aligned spatial data components

## See also

Functions for calculating disturbance:
[`DisturbanceMetrics-class`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md),
[`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/reference/disturbanceMetrics.md),
[`rasterizeLineDensity()`](https://landscitech.github.io/caribouMetrics/reference/rasterizeLineDensity.md),
[`reclassDist()`](https://landscitech.github.io/caribouMetrics/reference/reclassDist.md),
[`results()`](https://landscitech.github.io/caribouMetrics/reference/results.md),
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md)

## Examples

``` r
# create example rasters
lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
                     resolution = 250, crs = "EPSG:5070")
lc[] <- 0
nd <- lc
nd[1:30, 1:30] <- 1
ad <- lc
ad[30:50, 3:50] <- 1
lc[] <- 1
lc[70:100, 70:100] <- 2

# create sf objects
lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
                                                            ncol = 2, byrow = TRUE))),
                              crs = 5070))
esk <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 10000, 10000, 0),
                                                            ncol = 2, byrow = TRUE))),
                              crs = 5070))


projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))

# prepare data
res <- loadSpatialInputs(projectPoly = projPol, refRast = lc,
                         inputsList = list(esker = esk, linFeat = lf, natDist = nd,
                                           anthroDist = ad),
                         convertToRast = c("esker", "linFeat"))
#> cropping esker to extent of projectPoly
#> cropping linFeat to extent of projectPoly
                                           
```
