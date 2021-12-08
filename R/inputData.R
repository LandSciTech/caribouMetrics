#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating caribou habitat use.
#' 

#' @noRd
inputData <- function(landCover, esker, linFeat, projectPoly, caribouRange,
                      coefTable, natDist = NULL, anthroDist = NULL, 
                      eskerSave = NULL, linFeatSave = NULL, 
                      winArea = NULL, padProjPoly = FALSE,
                      padFocal = FALSE, ptDensity = 1, 
                      tmplt = NULL) {
  
  if(padProjPoly){
    bufWidth <- sqrt(winArea*10000/pi)
  } else {
    bufWidth <- NULL
  } 
  
  inData <- loadSpatialInputs(projectPoly = projectPoly, refRast = landCover, 
                              inputsList = dplyr::lst(linFeat, natDist, 
                                                      anthroDist, esker),
                              bufferWidth = bufWidth,
                              convertToRast = c("esker", "linFeat"),
                              useTemplate = c("esker", "linFeat"),
                              altTemplate = tmplt,
                              reclassOptions = list(natDist = cbind(NA, 0),
                                                    anthroDist = cbind(NA, 0)),
                              ptDensity = ptDensity)
  if(raster::isLonLat(inData$refRast)){
    stop("landCover must have a projected CRS", call. = FALSE)
  }
  
  if(is.character(landCover)){
    inData$refRast <- reclassPLC(inData$refRast)
  }
  
  # require column called Range
  if(nrow(inData$projectPolyOrig) == 1 & !"Range" %in% names(inData$projectPolyOrig)){
    inData$projectPolyOrig <- inData$projectPolyOrig %>% mutate(Range = caribouRange$Range)
  } else {
    if(!"Range" %in% names(inData$projectPolyOrig)){
      stop("projectPoly must have a column Range that corresponds to the ",
           "caribouRange$Range column", call. = FALSE)
    } else {
      if(!all(inData$projectPolyOrig$Range %in% caribouRange$Range)){
        stop("All values of in projectPoly$Range must have matching",
             " values in caribouRange$Range", call. = FALSE)
      }
    }
  }
  
  .checkInputs(caribouRange, winArea, inData$refRast, coefTable)
  
  if(is.null(natDist)){
    inData$natDist <- raster(matrix(NA))
  }
  
  if(is.null(anthroDist)){
    inData$anthroDist <- raster(matrix(NA))
  }
  
  if(!is.null(linFeatSave)){
    raster::writeRaster(inData$linFeat, 
                        linFeatSave, overwrite = TRUE)
    inData$linFeat <- raster(linFeatSave)
  }
  
  if(!is.null(eskerSave)){
    raster::writeRaster(inData$esker, eskerSave, overwrite = TRUE)
    inData$esker <- raster(eskerSave)
  }
  
  if(is.null(tmplt)){
    tmplt <- raster(inData$refRast) %>% raster::`res<-`(c(400, 400))
  } 
  
  return(new("CaribouHabitat", inData$refRast, inData$esker, 
             inData$natDist, inData$anthroDist,
             inData$linFeat, inData$projectPolyOrig,  
             processedData = raster(matrix(NA)), 
             habitatUse = raster(matrix(NA)),
             attributes = list(caribouRange = caribouRange, winArea = winArea,
                               padProjPoly = padProjPoly, padFocal = padFocal, 
                               tmplt = tmplt)))
}