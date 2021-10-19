setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("missingOrNULLOrChar", c("missing", "NULL","character"))

# NOTE: Constructors for each class are defined in the R file bearing the name
# of the class (lower case).

#' Caribou Habitat class
#'
#' An object containing data and results for the caribou habitat resource
#' selection function
#'
#' @seealso See \code{\link{caribouHabitat}} for options when creating a
#'   CaribouHabitat object.
#' @slot landCover Raster with land cover classified into resource types. See
#'   \code{resTypeCode} for a legend.
#' @slot esker Raster of esker density in m^2 per hectare.
#' @slot natDist Raster of natural disturbance.
#' @slot anthroDist Raster of anthropogenic disturbances.
#' @slot linFeat Raster of linear feature density in m^2 per hectare.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData RasterStack with named layers for each input variable
#'   used in the RSF
#' @slot habitatUse RasterStack with named layers for each season
#' @slot attributes A list of arguments from the \code{caribouHabtat} call. 
#' @name CaribouHabitat-class
#' @rdname CaribouHabitat-class
#' @export CaribouHabitat
#' @importClassesFrom raster RasterLayer

CaribouHabitat <- setClass("CaribouHabitat",
                           slots = c(landCover = "RasterLayer", 
                                     esker  = "RasterLayer",
                                     natDist  = "RasterLayer",
                                     anthroDist = "RasterLayer",
                                     linFeat  = "RasterLayer", 
                                     projectPoly = "sf", 
                                     processedData = "Raster",
                                     habitatUse = "Raster",
                                     attributes = "list"),
                           prototype = list(
                             landCover = raster(matrix(NA)), 
                             esker  = raster(matrix(NA)),
                             natDist = raster(matrix(NA)),
                             linFeat  = raster(matrix(NA)),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = raster::stack(raster(matrix(NA))),
                             habitatUse = raster::stack(raster(matrix(NA))),
                             attributes = list()
                           ))


#' Disturbance Metrics class
#'
#' An object containing data and results for predictors described in Table 52 of
#' Environment Canada (2011) Scientific Assessment to Inform the Identification
#' of Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal
#' Population, in Canada: 2011 Update. Ottawa, Ontario. Output variables
#' include:
#' \itemize{
#'   \item {Fire:} {\% fire}
#'   \item {Anthro:} {\% non-overlapping anthropogenic disturbance.}
#'   \item {Total_dist:} {Percent total non-overlapping fire and anthropogenic disturbance.}
#'   \item {fire_excl_anthro:} {\% fire not overlapping with anthropogenic disturbance.}
#' }
#'
#' Note that NA values are omitted from tabulated area.
#' Missing layers are omitted from the output, not interpreted as 0 disturbance
#'
#' @seealso See \code{\link{disturbanceMetrics}} for options when creating a
#'   DisturbanceMetrics object.
#'   
#' @slot landCover Raster distinguishing land from water. 0 or NA is water.
#' @slot natDist Raster of natural disturbance.
#' @slot anthroDist Raster of anthropogenic disturbances.
#' @slot linFeat Sf object with linear features including roads, utilities, and rail.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData RasterStack with named layers for each input variable used in the RSF
#' @slot disturbanceMetrics Data frame of disturbance metric values
#' @slot attributes A list of arguments from the \code{caribouHabtat} call. 
#' @name DisturbanceMetrics-class
#' @rdname DisturbanceMetrics-class
#' @export DisturbanceMetrics
#' @importClassesFrom raster RasterLayer

DisturbanceMetrics <- setClass("DisturbanceMetrics",
                           slots = c(landCover  = "RasterLayer",
                                     natDist  = "RasterLayer",
                                     anthroDist = "RasterLayer", 
                                     linFeat  = "list", #this is a hack - need to allow raster or sf, but can't setClassUnion because sf class is not exported from sf package. 
                                     projectPoly = "sf", 
                                     processedData = "Raster",
                                     disturbanceMetrics = "data.frame",
                                     attributes = "list"),
                           prototype = list(
                             landCover = raster(matrix(NA)),
                             linFeat  = list(),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = raster::stack(raster(matrix(NA))),
                             disturbanceMetrics = data.frame(),
                             attributes = list()
                           ))
