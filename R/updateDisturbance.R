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
    natDist <- reclassify(newData2@natDist, cbind(NA, 0))
  } else {
    natDist <- distMet@natDist
  }
  
  distMet@processedData <- calcDMSpatial(anthroDist, natDist, distMet@landCover,
                                         distMet)
  
  ####### Range summaries
  distMet@disturbanceMetrics <- calcDM(distMet, isPercent)
  
  return(distMet)
}