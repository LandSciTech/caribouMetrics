#' Reclassify natural disturbance and harvest layers
#'
#' Classify disturbances that have occurred in a user defined time
#' period. The output from this function can be used as either the
#' `natDist` or `anthroDist` inputs for `caribouMetrics()` or
#' `disturbanceMetrics()`.
#'
#' @param distYr sf object or rasterLayer. A simple feature collection or
#'   raster layer covering the focal area and including the year of disturbance
#'   or time since disturbance as well as geometry (multipolygon) of the
#'   disturbance.
#' @param endYr numeric default 0. Four digit year indicating the latest
#'   year to include from `distYr`. If 0 assume `dateField` values indicate time since disturbance.
#' @param numCumYrs numeric. Number of years before `endYr` to
#'   include.
#' @param template rasterLayer. A raster of the focal region, used as a
#'   template to which disturbance information is projected. This layers
#'   dimensions determine the dimensions of the output. It is recommended to use
#'   the `landCover` raster layer used in `caribouMetrics()` or
#'   `disturbanceMetrics()` to ensure equal dimensions
#' @param dateField character. Name of the column in which disturbance
#'   year/time since disturbance is recorded.
#'
#'
#' @return Returns a binary `RasterLayer` with the same dimensions as the
#'   `template` input. Values of 1 represent areas of disturbance
#'
#' @examples
#' library(sf)
#' # create template raster
#' lc <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0,
#'                      ymx = 10, crs = 5070)
#' 
#' # create fire polygons 
#' corners <- matrix(c(0,0,10,10,5, 0,10,10,0,5), ncol = 2)
#' 
#' fireYr <- st_sf(FIRE_YEAR = c(1990, 2000, 2009, 2015),
#'                      geometry = st_sfc(st_polygon(list(corners[c(1,2,5, 1),])),
#'                                        st_polygon(list(corners[c(2,3,5, 2),])),
#'                                        st_polygon(list(corners[c(3,4,5, 3),])),
#'                                        st_polygon(list(corners[c(4,1,5, 4),])))) 
#' fireYr <- st_set_crs(fireYr, 5070)
#' 
#' # three polygons should be considered disturbed (1) but the 2015 polygon should
#' # not (0)
#' cumFirePresence <- reclassDist(fireYr,
#'                                endYr = 2010,
#'                                numCumYrs = 40,
#'                                template = lc,
#'                                dateField = "FIRE_YEAR")
#' 
#' plot(cumFirePresence)
#' 
#' # with time since disturbance
#' fireYr$FIRE_YEAR <- c(50, 15, 10, 20)
#' cumFirePresence2 <- reclassDist(fireYr,
#'                                endYr = 0,
#'                                numCumYrs = 40,
#'                                template = lc,
#'                                dateField = "FIRE_YEAR")
#'                                
#' plot(cumFirePresence2)
#'
#' @family disturbance
#' @export
#' 

reclassDist <- function(distYr, endYr = 0, numCumYrs, template, dateField){
  if(inherits(distYr, "RasterLayer")){
    if(endYr == 0){
      out <-  distYr < numCumYrs
    } else {
      startYr <- endYr - numCumYrs
      
      out <- (distYr - startYr) >= 0
    }
  } else {
    # If no endYr is defined assume years relate to time since disturbance, and 
    # thus any value smaller than numCumYrs should be recorded as a disturbance
    if(endYr == 0){
      out <- dplyr::filter(distYr,
                           dplyr::between(distYr[[dateField]], 
                                          endYr, 
                                          numCumYrs))
      
    }else{ # If an endYr is included assume years represent dates of disturbance
      startYr <- endYr - numCumYrs
      
      out <- dplyr::filter(distYr,
                           dplyr::between(distYr[[dateField]], 
                                          startYr, 
                                          endYr))
    }
    
    if(requireNamespace("fasterize", quietly = TRUE)){
      if(st_geometry_type(out, by_geometry = FALSE) == "GEOMETRY"){
        out <- st_cast(out)
      }
      out <- fasterize::fasterize(out, template, background = 0)
    } else {
      message("To speed up install fasterize package")
      out <- raster::rasterize(out, template)
    }
  }
  out <- raster::reclassify(out, cbind(NA, 0))
  return(out)
}
