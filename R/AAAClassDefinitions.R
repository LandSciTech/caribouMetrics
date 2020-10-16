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
#'   CaribouHabitatUse object.
#' @slot plc Raster of provincial land cover.
#' @slot esker Raster of esker density in m^2 per hectare.
#' @slot fri Raster of forest resource inventory.
#' @slot age Raster of age in years from forest resource inventory.
#' @slot natDist Raster of natural disturbance.
#' @slot linFeat Raster of linear feature density in m^2 per hectare including
#'   roads, utilities, and rail.
#' @slot hexgrid Sf object with polygons in flat-topped hexagonal grid.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData Sf dataframe object with all data needed to calculate
#'   habitat use
#' @slot habitatUse Sf dataframe object with habitat use for each season and
#'   explanatory variables
#' @name CaribouHabitat-class
#' @rdname CaribouHabitat-class
#' @export CaribouHabitat
#' @importClassesFrom raster RasterLayer

CaribouHabitat <- setClass("CaribouHabitat",
                           slots = c(plc = "RasterLayer", 
                                     esker  = "RasterLayer",
                                     fri  = "RasterLayer", 
                                     age  = "RasterLayer",
                                     natDist  = "RasterLayer",
                                     anthroDist = "RasterLayer",
                                     harv = "RasterLayer", 
                                     linFeat  = "RasterLayer", 
                                     projectPoly = "sf", 
                                     processedData = "Raster",
                                     habitatUse = "Raster"),
                           prototype = list(
                             plc = raster(matrix(NA)), 
                             esker  = raster(matrix(NA)),
                             fri  = raster(matrix(NA)), 
                             age  = raster(matrix(NA)),
                             natDist = raster(matrix(NA)),
                             linFeat  = raster(matrix(NA)),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = raster::stack(raster(matrix(NA))),
                             habitatUse = raster::stack(raster(matrix(NA)))
                           ))

#' Landscape State Class
#'
#' Validated landscape states.
#'
#' @details
#'
#' Use lsBuffer() to fill buffer slot.
#' @examples
#' #TODO - update examples
#'
#' @slot maps named list of SpatialPolygons*, SpatialLines*, or Raster* objects. Each element of the list represents a different data set (fires, polygonal anthropogenic disturbance, linear disturbances, etc), and each layer represents a time point.
#' @slot bufferMaps named list of SpatialPolygons* or Raster* objects. Each element of the list represents a different data set (polygonal anthropogenic disturbance, linear disturbances, etc), and each layer represents a time point. See lsBuffer() for details.
#' @slot bufferWidths named list of numbers. Width of buffer in metres for each element of bufferMaps.
#' @slot params named list. Other parameter values
#' @slot mapType character. 'raster', 'polygon', or 'mixed'. 'mixed' indicates a combination of Raster* and SpatialPolygons* map layers.
#' @name LandscapeState-class
#' @rdname LandscapeState-class
#' @export LandscapeState
#' @importFrom methods new slot slotNames
LandscapeState <- setClass("LandscapeState", representation(maps="list",bufferWidths="list",params="list",mapType="character"))

#' Boreal Caribou Landscape State Class
#'
#' Boreal Caribou LandscapeState.
#'
#' @name LandscapeStateBCaribou-class
#' @rdname LandscapeStateBCaribou-class
#' @export LandscapeStateBCaribou
#' @importFrom methods new slot slotNames
LandscapeStateBCaribou <- setClass("LandscapeStateBCaribou", contains='LandscapeState',representation (fred='character'))
