#' Reclassify natural disturbance and harvest layers
#'
#'
#' @param distYr \code{sf object}. A simple feature collection covering the focal area and including the year of disturbance (preferred method) or time since disturbance as well as geometry (multipolygon) of the disturbance.
#' @param endYr \code{numeric} default 0. Four digit date indicating the latest year to include from \code{distYr}. If no value is included function assumes \code{dateField} values indicate time since disturbance.
#' @param numCumYrs \code{numeric}. Number of years before \code{endYr} to include.
#' @param landCover \code{rasterLayer}. A raster of the land cover for the focal region, used as a template over which disturbance information is projected. Thus this layers dimensions determine the dimensions of the output.
#' @param dateField \code{character}. Name of the column in which disturbance dates/time since disturbance are recorded.
#' 
#' @description This function classifies disturbances/harvests that have occurred in a user defined time period. The output from this function can be used as either the \code{natDist} or \code{harv} inputs for \code{caribouMetrics()} or \code{disturbanceMetrics()}.
#' 
#' @return Returns a binary \code{RasterLayer} with the same dimensions as the \code{landCover} input. Values of 1 represent areas of disturbance or harvest.
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
#'             landCover = plcD, 
#'             dateField = "FIRE_YEAR")
#' }
#' 
#' @export
#' 

reclassDist <- function(distYr, endYr = 0, numCumYrs, landCover, dateField){
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
  
  out <- fasterize::fasterize(out, landCover, background = 0)
  
  return(out)
}

###############################################################################

#' 
#' @param harvYr \code{sf object}. A simple feature collection covering the focal area and including the year of harvest (preferred method) or time since harvest as well as geometry (multipolygon) of the harvest.
#' @rdname reclassDist
#' @export
#' 

reclassHarv <- function(harvYr, endYr = 0, numCumYrs, landCover, dateField){
  # If no endYr is defined assume years relate to time since disturbance, and 
  # thus any value smaller than numCumYrs should be recorded as a disturbance
  if(endYr == 0){
    out <- dplyr::filter(harvYr,
                         dplyr::between(harvYr[[dateField]], 
                                        numCumYrs, 
                                        endYr))
    
  }else{ # If an endYr is included assume years represent dates of disturbance
    startYr <- endYr - numCumYrs
    
    out <- dplyr::filter(harvYr,
                         dplyr::between(harvYr[[dateField]], 
                                        startYr, 
                                        endYr))
  }
  
  out <- fasterize::fasterize(out, landCover, background = 0)
  
  return(out)
}
