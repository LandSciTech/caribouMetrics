# Disturbance Metrics class

An object containing data and results for predictors described in Table
52 of Environment Canada (2011) Scientific Assessment to Inform the
Identification of Critical Habitat for Woodland Caribou (Rangifer
tarandus caribou), Boreal Population, in Canada: 2011 Update. Ottawa,
Ontario. Output variables include:

- Fire: Percent fire

- Anthro: Percent non-overlapping buffered anthropogenic disturbance.

- Total_dist: Percent total non-overlapping fire and anthropogenic
  disturbance.

- Fire_excl_anthro: Percent fire not overlapping with anthropogenic
  disturbance.

## Details

Note that NA values are omitted from tabulated area. Missing layers are
omitted from the output, not interpreted as 0 disturbance

## Slots

- `landCover`:

  SpatRaster distinguishing land from water. 0 or NA is water.

- `natDist`:

  SpatRaster of natural disturbance.

- `anthroDist`:

  SpatRaster of anthropogenic disturbances.

- `linFeat`:

  Sf object with linear features including roads, utilities, and rail.

- `projectPoly`:

  Sf object with polygon of project boundary.

- `processedData`:

  SpatRaster with named layers for each input variable used in the RSF

- `disturbanceMetrics`:

  Data frame of disturbance metric values

- `attributes`:

  A list of arguments from the `caribouHabtat` call.

## See also

See
[`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/reference/disturbanceMetrics.md)
for options when creating a DisturbanceMetrics object.

Functions for calculating disturbance:
[`disturbanceMetrics()`](https://landscitech.github.io/caribouMetrics/reference/disturbanceMetrics.md),
[`loadSpatialInputs()`](https://landscitech.github.io/caribouMetrics/reference/loadSpatialInputs.md),
[`rasterizeLineDensity()`](https://landscitech.github.io/caribouMetrics/reference/rasterizeLineDensity.md),
[`reclassDist()`](https://landscitech.github.io/caribouMetrics/reference/reclassDist.md),
[`results()`](https://landscitech.github.io/caribouMetrics/reference/results.md),
[`updateDisturbance()`](https://landscitech.github.io/caribouMetrics/reference/updateDisturbance.md)
