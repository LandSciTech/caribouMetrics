#' @include AAAClassDefinitions.R
NULL

#' Extract results
#'
#' Extract results from CaribouHabitat or DisturbanceMetrics object.
#'
#' @param x A CaribouHabitat or DisturbanceMetrics object.
#' @param type string. The name of the slot to return. If x is a CaribouHabitat
#'   object the default is "both" and the habitatUse and processedData will be
#'   returned as a RasterStack.
#' @param ... arguments passed to methods
#'
#' @return By default a RasterStack if x is a CaribouHabitat object and a
#'   data.frame if x is a DisturbanceMetrics object. 
#'   
#' @examples 
#' # create example rasters for habitat
#' lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
#'                      resolution = 250, crs = "EPSG:5070")
#' lc[] <- 0
#' nd <- lc
#' nd[1:30, 1:30] <- 1
#' ad <- lc
#' ad[30:50, 3:50] <- 1
#' lc[] <- 1
#' lc[70:100, 70:100] <- 2
#' 
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' esk <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 10000, 10000, 0),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' 
#' 
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#' 
#' # calculate relative probability of use
#' res <- caribouHabitat(landCover = lc,
#'                linFeat = lf,
#'                esker = esk,
#'                natDist = nd,
#'                anthroDist = ad,
#'                projectPoly = projPol,
#'                caribouRange = "Nipigon",
#'                winArea = 1000 #leave as default NULL except for small examples
#' )
#' 
#' # default gets raster of both habitatUse and processedData slots
#' resBoth <- results(res)
#' 
#' # provide a slot name to get one of them
#' resHabUse <- results(res, type = "habitatUse")
#' 
#' # create example rasters for disturbance metrics
#' lc <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0,
#'                   ymax = 10, crs = "EPSG:5070")
#' nd <- lc
#' nd[1:3, 1:3] <- 1
#' ad <- lc
#' ad[3:5, 3:5] <- 1
#' lc[] <- 1
#' 
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10, 10), 
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#' 
#' # calculate disturbance
#' dm <- disturbanceMetrics(landCover = lc,
#'                          linFeat = lf,
#'                          natDist = nd,
#'                          anthroDist = ad,
#'                          projectPoly = projPol,
#'                          padFocal = TRUE,
#'                          bufferWidth = 1)
#' 
#' # default is disturbance metrics table
#' resDM <- results(dm)
#' 
#' # can get other slots as well
#' resDMrasters <- results(dm, type = "processedData")
#' 
#' @family habitat
#' @family disturbance
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname results
setMethod("results", signature(x = "CaribouHabitat"), function(x, type = "both"){
  if(nrow(x@habitatUse) < 2){
    stop("This object has empty @habitatUse. Calculate ",
         "habitatUse first with updateCaribou(x)")
  }
  
  if(type == "both"){
    result <- c(x@habitatUse, x@processedData)
    
    return(result)
  }
  
  return(slot(x, type))

})

#' @rdname results
setMethod("results", signature(x = "DisturbanceMetrics"), function(x, type = "disturbanceMetrics"){
  slot(x, type)
})