setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("missingOrNULLOrChar", c("missing", "NULL","character"))

# NOTE: Constructors for each class are defined in the R file bearing the name
# of the class (lower case).

#' Caribou Habitat class
#'
#' An object containing data and results for the caribou habitat resource
#' selection function
#'
#' @seealso See [caribouHabitat()] for options when creating a
#'   CaribouHabitat object.
#' @slot landCover SpatRaster with land cover classified into resource types. See
#'   `resTypeCode` for a legend.
#' @slot esker SpatRaster of esker density in m^2 per hectare.
#' @slot natDist SpatRaster of natural disturbance.
#' @slot anthroDist SpatRaster of anthropogenic disturbances.
#' @slot linFeat SpatRaster of linear feature density in m^2 per hectare.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData SpatRaster with named layers for each input variable
#'   used in the RSF
#' @slot habitatUse SpatRaster with named layers for each season
#' @slot attributes A list of arguments from the `caribouHabtat` call. 
#' @name CaribouHabitat-class
#' @rdname CaribouHabitat-class
#' @family habitat
#' @export CaribouHabitat
#' @importClassesFrom terra SpatRaster

CaribouHabitat <- setClass("CaribouHabitat",
                           slots = c(landCover = "SpatRaster", 
                                     esker  = "SpatRaster",
                                     natDist  = "SpatRaster",
                                     anthroDist = "SpatRaster",
                                     linFeat  = "SpatRaster", 
                                     projectPoly = "sf", 
                                     processedData = "SpatRaster",
                                     habitatUse = "SpatRaster",
                                     attributes = "list"),
                           prototype = list(
                             landCover = terra::rast(matrix(NA)), 
                             esker  = terra::rast(matrix(NA)),
                             natDist = terra::rast(matrix(NA)),
                             linFeat  = terra::rast(matrix(NA)),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = terra::rast(matrix(NA)),
                             habitatUse = terra::rast(matrix(NA)),
                             attributes = list()
                           ))


#' Disturbance Metrics class
#'
#' An object containing data and results for predictors described in Table 52 of
#' Environment Canada (2011) Scientific Assessment to Inform the Identification
#' of Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal
#' Population, in Canada: 2011 Update. Ottawa, Ontario. Output variables
#' include:
#'  * Fire: Percent fire
#'  * Anthro: Percent non-overlapping buffered anthropogenic disturbance.
#'  * Total_dist: Percent total non-overlapping fire and anthropogenic disturbance.
#'  * fire_excl_anthro: Percent fire not overlapping with anthropogenic disturbance.
#'
#' Note that NA values are omitted from tabulated area.
#' Missing layers are omitted from the output, not interpreted as 0 disturbance
#'
#' @seealso See [disturbanceMetrics()] for options when creating a
#'   DisturbanceMetrics object.
#'   
#' @slot landCover SpatRaster distinguishing land from water. 0 or NA is water.
#' @slot natDist SpatRaster of natural disturbance.
#' @slot anthroDist SpatRaster of anthropogenic disturbances.
#' @slot linFeat Sf object with linear features including roads, utilities, and rail.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData SpatRaster with named layers for each input variable used in the RSF
#' @slot disturbanceMetrics Data frame of disturbance metric values
#' @slot attributes A list of arguments from the `caribouHabtat` call. 
#' @name DisturbanceMetrics-class
#' @rdname DisturbanceMetrics-class
#' @family disturbance
#' @export DisturbanceMetrics
#' @importClassesFrom terra SpatRaster

DisturbanceMetrics <- setClass("DisturbanceMetrics",
                           slots = c(landCover  = "SpatRaster",
                                     natDist  = "SpatRaster",
                                     anthroDist = "SpatRaster", 
                                     linFeat  = "list", #this is a hack - need to allow raster or sf, but can't setClassUnion because sf class is not exported from sf package. 
                                     projectPoly = "sf", 
                                     processedData = "SpatRaster",
                                     disturbanceMetrics = "data.frame",
                                     attributes = "list"),
                           prototype = list(
                             landCover = terra::rast(matrix(NA)),
                             linFeat  = list(),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = terra::rast(matrix(NA)),
                             disturbanceMetrics = data.frame(),
                             attributes = list()
                           ))
