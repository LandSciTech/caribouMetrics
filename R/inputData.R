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
                      tmplt = NULL, preppedData = NULL) {
  
  if(is.null(preppedData)){
    if(padProjPoly){
      bufWidth <- sqrt(winArea*10000/pi)
    } else {
      bufWidth <- NULL
    } 
    
    preppedData <- loadSpatialInputs(projectPoly = projectPoly, refRast = landCover, 
                                     inputsList = dplyr::lst(linFeat, natDist, 
                                                             anthroDist, esker),
                                     bufferWidth = bufWidth,
                                     convertToRastDens = c("esker", "linFeat"),
                                     useTemplate = c("esker", "linFeat"),
                                     altTemplate = tmplt,
                                     # reclassOptions = list(natDist = cbind(NA, 0),
                                     #                       anthroDist = cbind(NA, 0)),
                                     ptDensity = ptDensity)
  }
  
  if(terra::is.lonlat(preppedData$refRast)){
    stop("landCover must have a projected CRS", call. = FALSE)
  }
  
  if(is.character(landCover)){
    preppedData$refRast <- reclassPLC(preppedData$refRast)
  }
  
  # require column called Range
  if(nrow(preppedData$projectPolyOrig) == 1 && !"Range" %in% names(preppedData$projectPolyOrig)){
    preppedData$projectPolyOrig <- preppedData$projectPolyOrig %>% mutate(Range = caribouRange$Range)
  } else {
    if(!"Range" %in% names(preppedData$projectPolyOrig)){
      stop("projectPoly must have a column Range that corresponds to the ",
           "caribouRange$Range column", call. = FALSE)
    } else {
      if(!all(preppedData$projectPolyOrig$Range %in% caribouRange$Range)){
        stop("All values of in projectPoly$Range must have matching",
             " values in caribouRange$Range", call. = FALSE)
      }
    }
  }
  
  .checkInputs(caribouRange, winArea, preppedData$refRast, coefTable)
  
  if(is.null(preppedData$natDist)){
    preppedData$natDist <- terra::rast(matrix(NA))
  }
  
  if(is.null(preppedData$anthroDist)){
    preppedData$anthroDist <- terra::rast(matrix(NA))
  }
  
  if(!is.null(linFeatSave)){
    terra::writeRaster(preppedData$linFeat, 
                        linFeatSave, overwrite = TRUE)
    preppedData$linFeat <- terra::rast(linFeatSave)
  }
  
  if(!is.null(eskerSave)){
    terra::writeRaster(preppedData$esker, eskerSave, overwrite = TRUE)
    preppedData$esker <- terra::rast(eskerSave)
  }
  
  if(is.null(tmplt)){
    tmplt <- terra::rast(preppedData$refRast) %>% terra::`res<-`(c(400, 400))
  } 
  
  return(new("CaribouHabitat", preppedData$refRast, preppedData$esker, 
             preppedData$natDist, preppedData$anthroDist,
             preppedData$linFeat, preppedData$projectPolyOrig,  
             processedData = terra::rast(matrix(NA)), 
             habitatUse = terra::rast(matrix(NA)),
             attributes = list(caribouRange = caribouRange, winArea = winArea,
                               padProjPoly = padProjPoly, padFocal = padFocal, 
                               tmplt = tmplt)))
}