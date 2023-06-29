#' @include AAAClassDefinitions.R
NULL

#' Calculate binary use
#'
#' Convert continuous probability of use in each season to binary high use in
#' any season or each season using thresholds given in a table.
#' 
#'
#' @param x CaribouHabitat object where habitat use has been calculated, or a
#'   data.frame with columns PID, Spring, Summer, Fall, Winter, or a numeric
#'   vector containing ID values, in which case spring, summer, fall and winter
#'   must also be supplied
#' @param tholdTable By default `threshTable` which contains thresholds
#'   determimed by Rempel (2021) using Youden's J and a false
#'   negative cost of 5. Change at own risk.
#' @param bySeason logical. If FALSE (the default) the result is a single value
#'   of 1 when the habitat use was >= threshold in any season and 0 if not. If
#'   TRUE then the result is a 0 or 1 for each season.
#' @param ... passed to methods
#'
#'
#' @return
#'  If `bySeason` is `FALSE` and `x` is a CaribouHabitat object then the result is 1 layer SpatRaster.
#'  If `bySeason` is `TRUE` and `x` is a CaribouHabitat object then the result is SpatRaster with a layer per season.  
#'  If `bySeason` is `FALSE` and `x` is a data.frame then the result is a data.frame.
#'  If `bySeason` is `TRUE` and `x` is a data.frame then the result is a data.frame.
#'  If `bySeason` is `FALSE` and `x` is a vector then the result is a vector.
#'  If `bySeason` is `TRUE` and `x` is a vector then the result is a data.frame.  
#'  
#' @examples 
#' # create example rasters
#' lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
#'                      resolution = 250, crs = sf::st_crs(5070)$wkt)
#' lc[] <- 0
#' nd <- lc
#' nd[1:30, 1:30] <- 1
#' ad <- lc
#' ad[30:50, 3:50] <- 1
#' lc[] <- 1
#' lc[70:100, 70:100] <- 2
#' 
#' 
#' 
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = sf::st_crs(5070)))
#' esk <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 10000, 10000, 0),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = sf::st_crs(5070)))
#' 
#' 
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#' projPol <- sf::st_set_crs(projPol, sf::st_crs(5070))
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
#' plot(calcBinaryUse(res))
#' plot(calcBinaryUse(res, bySeason = TRUE))
#' @family habitat
#' @export
setGeneric("calcBinaryUse", function(x, ...) standardGeneric("calcBinaryUse"))

#' @rdname calcBinaryUse
setMethod(
  "calcBinaryUse", signature(x = "CaribouHabitat"), 
  function(x, tholdTable = threshTable, bySeason = FALSE){
    caribouRange <- x@attributes$caribouRange$coefRange
    
    tTable <- threshTable %>% filter(.data$Range == caribouRange) %>% 
      arrange(.data$Season) %>% select("Season", "Threshold")
    
    if(bySeason){
      binUse <- x@habitatUse >= tTable$Threshold
    } else {
      binUse <- any(x@habitatUse >= tTable$Threshold)
      names(binUse) <- "BinaryUse"
    }
    
    return(binUse)
  })
