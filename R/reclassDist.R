#' Reclassify natural disturbance and harvest layers
#'
#' Classify disturbances/harvests that have occurred in a user defined time
#' period. The output from this function can be used as either the
#' \code{natDist} or \code{anthroDist} inputs for \code{caribouMetrics()} or
#' \code{disturbanceMetrics()}.
#'
#' @param distYr \code{sf object or rasterLayer}. A simple feature collection
#'   (preferred) or raster layer covering the focal area and including the year
#'   of disturbance (preferred method) or time since disturbance as well as
#'   geometry (multipolygon) of the disturbance.
#' @param endYr \code{numeric} default 0. Four digit date indicating the latest
#'   year to include from \code{distYr}. If no value is included function
#'   assumes \code{dateField} values indicate time since disturbance.
#' @param numCumYrs \code{numeric}. Number of years before \code{endYr} to
#'   include.
#' @param template \code{rasterLayer}. A raster of the focal region, used as a
#'   template over which disturbance information is projected. Thus this layers
#'   dimensions determine the dimensions of the output. It is recommended to use
#'   the \code{landcover} raster layer used in \code{caribouMetrics()} or
#'   \code{disturbanceMetrics()} to ensure equal dimensions
#' @param dateField \code{character}. Name of the column in which disturbance
#'   dates/time since disturbance are recorded.
#'
#'
#' @return Returns a binary \code{RasterLayer} with the same dimensions as the
#'   \code{template} input. Values of 1 represent areas of disturbance
#'
#' @examples
#' \dontrun{
#' plcD <- raster::raster("/your_data_directory/plc50.tif") %>%
#'                        reclassPLC()
#' fireYr <- sf::st_read("/your_data_directory/fireAFFES.shp") %>%
#'                       st_transform(crs = st_crs(plcD))
#'
#' reclassDist(fireYr,
#'             endYr = 2010,
#'             numCumYrs = 30,
#'             template = plcD,
#'             dateField = "FIRE_YEAR")
#' }
#'
#' @export
#' 

reclassDist <- function(distYr, endYr = 0, numCumYrs, template, dateField){
  # TODO find quicker more efficient version
  # If distYr is supplied as a raster convert it to an sf object
  if(inherits(distYr, "RasterLayer")){
    tmp <- stars::st_as_stars(distYr %>% raster::as.factor())
    distYr <- sf::st_as_sf(tmp,
                           as_points = FALSE, merge = TRUE)
  }
  
  # If no endYr is defined assume years relate to time since disturbance, and 
  # thus any value smaller than numCumYrs should be recorded as a disturbance
  if(endYr == 0){
    out <- dplyr::filter(distYr,
                         dplyr::between(distYr[[dateField]], 
                                        numCumYrs, 
                                        endYr))
    
  }else{ # If an endYr is included assume years represent dates of disturbance
    startYr <- endYr - numCumYrs
    
    out <- dplyr::filter(distYr,
                         dplyr::between(distYr[[dateField]], 
                                        startYr, 
                                        endYr))
  }
  
  out <- fasterize::fasterize(out, template, background = 0)
  
  return(out)
}