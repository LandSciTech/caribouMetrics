# Calculate metrics of natural and anthropogenic disturbance

Calculate the predictors described in Table 52 of Environment Canada
(2011) Scientific Assessment to Inform the Identification of Critical
Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal
Population, in Canada:2011 Update. Ottawa, Ontario.The variables
calculated by this function include:

- Fire: % fire

- Anthro: % non-overlapping anthropogenic disturbance.

- Total_dist: Percent total non-overlapping fire and anthropogenic
  disturbance.

- Fire_excl_anthro: % fire not overlapping with anthropogenic
  disturbance.

## Usage

``` r
disturbanceMetrics(
  landCover = NULL,
  linFeat = NULL,
  projectPoly = NULL,
  isPercent = TRUE,
  ...
)
```

## Source

Environment Canada. 2011. Scientific Assessment to Inform the
Identification of Critical Habitat for Woodland Caribou (Rangifer
tarandus caribou), Boreal Population, in Canada:2011 Update. Ottawa,
Ontario.

## Arguments

- landCover:

  filename, SpatRaster or RasterLayer. 0 and NA values are assumed to be
  water and omitted from the tabulated area. Note landCover is also used
  to define the input grid, so must be provided even if all values are
  1.

- linFeat:

  filename, SpatRaster, RasterLayer, sf object or a list of these that
  will be combined. Linear features.

- projectPoly:

  filename or sf object. Polygons defining range boundaries.

- isPercent:

  logical. Should the results be returned as a percentage?

- ...:

  optional arguments:

  - natDist: filename, SpatRaster or RasterLayer. Presence or absence of
    natural disturbance, primarily by fire. Should include 40 years
    cumulative disturbance.

  - anthroDist: filename, SpatRaster or RasterLayer. Anthropogenic
    disturbance including harvest. This can have an effect on any type
    of landcover except water. Should include 40 years cumulative
    disturbance.

  - padProjPoly: logical. Should the area around the `projectPoly` be
    used to avoid edge effects? If FALSE, the default, only data from
    inside the `projectPoly` is used. If TRUE then `projectPoly` is
    buffered and the other variables are clipped to the extent of the
    buffered area. Results are always clipped to the original
    `projectPoly`. It is ideal to set this to TRUE and provide a dataset
    that is larger than the `projectPoly` to avoid edge effects.

  - padFocal: logical. This value is passed to the pad argument in
    [`terra::focal`](https://rspatial.github.io/terra/reference/focal.html),
    if it is FALSE then cells near the edge will return NA, if it is
    TRUE a value will be returned for each cell that assumes cells
    outside the input data are 0 for all resource types. This is not a
    good assumption and should be used with caution.

  - bufferWidth: number. Width of buffer applied to anthropogenic
    disturbance in metres. Default is 500.

  - linBuffMethod: character. The method used to buffer linear features
    if they are supplied as sf lines. The default is "raster" in which
    case they are rasterized and then buffered using a moving window
    method. If "sf" then the lines are buffered with st_buffer and then
    rasterized. Either way points are included in the raster output.

  - saveOutput: character. The filename to save the raster of binary
    disturbances with buffered anthropogenic disturbance. Note this will
    overwrite existing files with the same name. The .grd format is
    recommended because it will preserve layer names when the file is
    reloaded.

  - preppedData: list. A list containing pre-prepared input data sets.
    If not NULL then data checks will be skipped. Names must match
    argument names except that `landCover` should be called `refRast`
    and `projectPoly` should be called `projectPolyOrig` See
    [`loadSpatialInputs()`](https://landscitech.github.io/caribouMetrics/reference/loadSpatialInputs.md).

## Value

A DisturbanceMetrics Object see
[`DisturbanceMetrics-class()`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md)

## Details

Note assume natDist and anthroDist include 40 years of cumulative
disturbance. Note that locations where landCover is NA or 0 are omitted
from the tabulated area. Missing layers are omitted from the output, not
interpreted as 0 disturbance. To update an existing DisturbanceMetrics
object with new data see
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md).

## See also

[`DisturbanceMetrics-class()`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md)
for information on the object returned and
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md)
for updating an existing DisturbanceMetrics object.

Functions for calculating disturbance:
[`DisturbanceMetrics-class`](https://landscitech.github.io/caribouMetrics/reference/DisturbanceMetrics-class.md),
[`loadSpatialInputs()`](https://landscitech.github.io/caribouMetrics/reference/loadSpatialInputs.md),
[`rasterizeLineDensity()`](https://landscitech.github.io/caribouMetrics/reference/rasterizeLineDensity.md),
[`reclassDist()`](https://landscitech.github.io/caribouMetrics/reference/reclassDist.md),
[`results()`](https://landscitech.github.io/caribouMetrics/reference/results.md),
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md)

## Examples

``` r
# create example rasters
lc <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0,
                  ymax = 10, crs = "EPSG:5070")
nd <- lc
nd[1:3, 1:3] <- 1
ad <- lc
ad[3:5, 3:5] <- 1
lc[] <- 1

# create sf objects
lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10, 10), 
                                                            ncol = 2, byrow = TRUE))),
                              crs = 5070))
projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))

# calculate disturbance
disturbanceMetrics(landCover = lc,
                   linFeat = lf,
                   natDist = nd,
                   anthroDist = ad,
                   projectPoly = projPol,
                   padFocal = TRUE,
                   bufferWidth = 1)
#> buffering anthropogenic disturbance
#> calculating disturbance metrics
#> An object of class "DisturbanceMetrics"
#> Slot "landCover":
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   :     1 
#> max value   :     1 
#> 
#> Slot "natDist":
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   :     1 
#> max value   :     1 
#> 
#> Slot "anthroDist":
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   :     1 
#> max value   :     1 
#> 
#> Slot "linFeat":
#> [[1]]
#> Simple feature collection with 1 feature and 0 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 10 ymax: 10
#> Projected CRS: NAD83 / Conus Albers
#>                         x
#> 1 LINESTRING (0 0, 10 10)
#> 
#> 
#> Slot "projectPoly":
#> Simple feature collection with 1 feature and 0 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 10 ymax: 10
#> Projected CRS: NAD83 / Conus Albers
#>   sf..st_as_sfc.sf..st_bbox.ad..
#> 1 POLYGON ((0 0, 10 0, 10 10,...
#> 
#> Slot "processedData":
#> class       : SpatRaster 
#> size        : 10, 10, 4  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070) 
#> source(s)   : memory
#> names       : Anthro, Fire, Total_dist, Fire_excl_anthro 
#> min values  :      0,    0,          0,                0 
#> max values  :      1,    1,          1,                1 
#> 
#> Slot "disturbanceMetrics":
#>   zone Anthro Fire Total_dist Fire_excl_anthro
#> 1    1     52    9         58                6
#> 
#> Slot "attributes":
#> $bufferWidth
#> [1] 1
#> 
#> $padProjPoly
#> [1] FALSE
#> 
#> $padFocal
#> [1] TRUE
#> 
#> 
```
