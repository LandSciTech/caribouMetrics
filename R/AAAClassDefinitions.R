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
#' @slot landCover Raster of provincial land cover.
#' @slot esker Raster of esker density in m^2 per hectare.
#' @slot updatedLC Raster of forest resource inventory.
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
                           slots = c(landCover = "RasterLayer", 
                                     esker  = "RasterLayer",
                                     updatedLC  = "RasterLayer", 
                                     age  = "RasterLayer",
                                     natDist  = "RasterLayer",
                                     anthroDist = "RasterLayer",
                                     harv = "RasterLayer", 
                                     linFeat  = "RasterLayer", 
                                     projectPoly = "sf", 
                                     processedData = "Raster",
                                     habitatUse = "Raster",
                                     attributes = "list"),
                           prototype = list(
                             landCover = raster(matrix(NA)), 
                             esker  = raster(matrix(NA)),
                             updatedLC  = raster(matrix(NA)), 
                             age  = raster(matrix(NA)),
                             natDist = raster(matrix(NA)),
                             linFeat  = raster(matrix(NA)),
                             projectPoly = st_sf(col1 = NA,
                                                 geometry = st_sfc(st_point(x = c(1,1)))),
                             processedData = raster::stack(raster(matrix(NA))),
                             habitatUse = raster::stack(raster(matrix(NA))),
                             attributes = list()
                           ))

