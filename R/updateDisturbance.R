#' Update an existing DisturbanceMetrics Object
#'
#' Update a DisturbanceMetrics object in order to avoid reprocessing parts of
#' the data that have not changed. New data is supplied as a named list and the
#' object is updated depending on the elements provided in the list.
#'
#' If \code{newData} contains only linFeat or anthroDist then the new data will
#' be combined with the existing data in the DisturbanceMetrics object to
#' determine the buffered anthropogenic disturbance
#'
#' If  \code{newData} contains only natDist then only Fire is updated.
#'
#' @param distMet DisturbanceMetrics object
#' @param newData named list of objects to be used to update distMet.
#'   Potential names are: natDist, anthroDist, and linFeat.
#' @param linBuffMethod character. The method used to buffer linear features if
#'   they are supplied as sf lines. The default is "raster" in which case they
#'   are rasterized using the stars package and buffered using a moving window
#'   method. If "sf" then the lines are buffered with st_buffer and then
#'   rasterized. Either way points are included in the raster output.
#' @param isPercent logical. Should the results be returned as a percentage?
#'
#'@return A DisturbanceMetrics Object see \code{\link{DisturbanceMetrics-class}}
#'
#'@seealso \code{\link{DisturbanceMetrics-class}} for information on the object
#'  returned and \code{\link{disturbanceMetrics}} for making a new
#'  DisturbanceMetrics object. 
#' @export
#'
#' @examples
#' # create example rasters
#' lc <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10, crs = 5070)
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
#' # new linear features
#' lf2 <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10, 10),
#'                                                              ncol = 2, byrow = TRUE)),
#'                                     sf::st_linestring(matrix(c(0, 9, 9, 10),
#'                                                              ncol = 2, byrow = TRUE)),
#'                                     sf::st_linestring(matrix(c(10, 0, 5, 5),
#'                                                              ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' dm2 <- updateDisturbance(dm, newData = list(linFeat = lf2))
  
updateDisturbance <- function(distMet, newData, linBuffMethod = "raster",
                              isPercent = TRUE){
  
  # check names match expected
  if(is.null(names(newData)) ||
     !all(names(newData) %in% c("natDist", "anthroDist", "linFeat" ))){
    stop("newData must be a named list with any of the names ",
         "natDist ", "anthroDist ",  "linFeat ", 
         "incorrect names provided are:  ", 
         names(newData)[which(!names(newData) %in% 
                                c("natDist", "anthroDist",
                                  "linFeat" ))], call. = FALSE)
  }
  

  newData2 <- inputDataDisturbance(distMet@landCover, 
                                   linFeat = newData$linFeat,
                                   projectPoly = distMet@projectPoly, 
                                   natDist = newData$natDist, 
                                   anthroDist = newData$anthroDist,
                                   bufferWidth = distMet@attributes$bufferWidth, 
                                   padProjPoly = distMet@attributes$padProjPoly,
                                   padFocal = distMet@attributes$padFocal)
  
  if(!is.null(newData$anthroDist)||!is.null(newData$linFeat)){
    if(!is.null(newData$anthroDist)){
      anthroDist <- newData2@anthroDist
    } else{
      anthroDist <- distMet@anthroDist
    }
    
    if(!is.null(newData$linFeat)){
      linFeat <- newData2@linFeat[[1]]
    } else {
      linFeat <- distMet@linFeat[[1]]
    }

    anthroDist <- processAnthroDM(anthroDist, linFeat,
                                  distMet@landCover, 
                                  linBuffMethod = linBuffMethod,
                                  inData = distMet)
    
  } else {
    anthroDist <- distMet@processedData$Anthro
  }
  
  if(!is.null(newData$natDist)){
    # check natDist is real if not make dummy
    if(raster::ncell(newData2@natDist) == 1){
      newData2@natDist <- raster::init(distMet@landCover, 
                              fun = function(x){rep(0, x)}, 
                              filename = raster::rasterTmpFile())
    }
  } else {
    natDist <- distMet@natDist
  }
  
  natDist <- reclassify(natDist, cbind(NA, 0))
  
  distMet@processedData <- calcDMSpatial(anthroDist, natDist, distMet@landCover,
                                         distMet)
  
  ####### Range summaries
  distMet@disturbanceMetrics <- calcDM(distMet, isPercent)
  
  return(distMet)
}