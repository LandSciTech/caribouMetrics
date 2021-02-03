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
#'   \code{resTypeCode} for a list.
#' @slot esker Raster of esker density in m^2 per hectare.
#' @slot updatedLC Raster used to update landCover based on new data and
#'   disturbance.
#' @slot age Raster of age in years from forest resource inventory.
#' @slot natDist Raster of natural disturbance.
#' @slot anthroDist Raster of permanent anthropogenic disturbances.
#' @slot harv Raster of disturbance from forest harvesting.
#' @slot linFeat Raster of linear feature density in m^2 per hectare including
#'   roads, utilities, and rail.
#' @slot projectPoly Sf object with polygon of project boundary.
#' @slot processedData RasterStack with named layers for each input variable
#'   used in the RSF
#' @slot habitatUse RasterStack named with layers for each season
#' @slot attributes A list of arguments from the \code{caribouHabtat} call. 
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

