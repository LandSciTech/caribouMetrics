#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating disturbance metrics.
#' 
#' @noRd
inputDataDisturbance <- function(landCover, projectPoly, linFeat = NULL,
                                 natDist = NULL, anthroDist = NULL,
                                 bufferWidth = 500, padProjPoly = FALSE,
                                 padFocal = FALSE, preppedData = NULL) {
  if(is.null(preppedData)){
    if(padProjPoly){
      bufWidth <- bufferWidth
    } else {
      bufWidth <- NULL
    } 
    
    preppedData <- loadSpatialInputs(projectPoly = projectPoly, refRast = landCover, 
                                inputsList = dplyr::lst(linFeat, natDist, anthroDist),
                                bufferWidth = bufWidth,
                                reclassOptions = list(natDist = cbind(NA, 0),
                                                      anthroDist = cbind(NA, 0))) 
  }
  
  if(raster::isLonLat(preppedData$refRast)){
    stop("landCover must have a projected CRS", call. = FALSE)
  }
  
  if(is.null(preppedData$natDist)){
    preppedData$natDist <- raster(matrix(NA))
  }
  
  if(is.null(preppedData$anthroDist)){
    preppedData$anthroDist <- raster(matrix(NA))
  }
  
  if(is.null(preppedData$linFeat)){
    preppedData$linFeat <- raster(matrix(NA))
  }
  
  return(new("DisturbanceMetrics", preppedData$refRast, preppedData$natDist, 
             preppedData$anthroDist, 
             preppedData$linFeat, preppedData$projectPolyOrig,  
             processedData = raster(matrix(NA)), 
             disturbanceMetrics = data.frame(),
             attributes = list(bufferWidth = bufferWidth,
                               padProjPoly = padProjPoly, 
                               padFocal = padFocal)))
}
