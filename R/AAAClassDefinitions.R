setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("missingOrNULLOrChar", c("missing", "NULL","character"))

#NOTE: Constructors for each class are defined in the R file bearing the name of the class (lower case).

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
