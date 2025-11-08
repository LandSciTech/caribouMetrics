# Rasterize line density

Rasterize line density in meters per hectare.

## Usage

``` r
rasterizeLineDensity(x, r, ptDensity = 1)
```

## Arguments

- x:

  an sf object containing lines and/or points

- r:

  a SpatRaster or RasterLayer object to be used as a template for the
  output raster

- ptDensity:

  a number giving the density to assign to points, in units of `res(r)`.
  A value of 1 indicates one straight line crossing of the pixel. A
  value of 2+2\*2^0.5 is horizontal, vertical, and diagonal crossings.
  If NULL, points in x will be ignored.

## Value

A SpatRaster object with values representing the density of lines in
meters per hectare.

## See also

Functions for calculating disturbance:
[`DisturbanceMetrics-class`](https://landscitech.github.io/caribouMetrics/dev/reference/DisturbanceMetrics-class.md),
[`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/dev/reference/disturbanceMetrics.md),
[`loadSpatialInputs()`](https://landscitech.github.io/caribouMetrics/dev/reference/loadSpatialInputs.md),
[`reclassDist()`](https://landscitech.github.io/caribouMetrics/dev/reference/reclassDist.md),
[`results()`](https://landscitech.github.io/caribouMetrics/dev/reference/results.md),
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/dev/reference/updateDisturbance.md)

## Examples

``` r
# create example raster
lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
                     resolution = 250, crs = "EPSG:5070")

#' # create line
lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
                                                            ncol = 2, byrow = TRUE)),
                                   sf::st_linestring(matrix(c(1, 10001, 10001, 1),
                                                            ncol = 2, byrow = TRUE)),
                                   sf::st_linestring(matrix(c(5001, 10001, 5001, 1),
                                                            ncol = 2, byrow = TRUE))),
                                   crs = 5070))

rastLines <- rasterizeLineDensity(lf, lc)

plot(rastLines)

#> Error in plot.xy(xy, type, ...): invalid type passed to graphics function
```
