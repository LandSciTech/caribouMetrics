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
                                 padFocal = FALSE) {
  
  inData <- loadSpatialInputs(projectPoly = projectPoly, refRast = landCover, 
                              inputsList = dplyr::lst(linFeat, natDist, anthroDist),
                              bufferWidth = 500,
                              reclassOptions = list(natDist = cbind(NA, 0),
                                                    anthroDist = cbind(NA, 0)))
  
  if(raster::isLonLat(inData$refRast)){
    stop("landCover must have a projected CRS", call. = FALSE)
  }
  
  if(is.null(natDist)){
    inData$natDist <- raster(matrix(NA))
  }
  
  if(is.null(anthroDist)){
    inData$anthroDist <- raster(matrix(NA))
  }
  
  if(is.null(linFeat)){
    inData$linFeat <- raster(matrix(NA))
  }
  
  return(new("DisturbanceMetrics", inData$refRast, inData$natDist, 
             inData$anthroDist, 
             inData$linFeat, inData$projectPolyOrig,  
             processedData = raster(matrix(NA)), 
             disturbanceMetrics = data.frame(),
             attributes = list(bufferWidth = bufferWidth,
                               padProjPoly = padProjPoly, 
                               padFocal = padFocal)))
}
