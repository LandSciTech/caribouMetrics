#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating disturbance metrics.
#' 
#' @param landCover 
#' @param natDist 
#' @param anthroDist
#' @param linFeat 
#' @param projectPoly 
#' @param ... 
#'
#' @export
setGeneric("inputDataDisturbance", function(landCover,linFeat, projectPoly, ...) standardGeneric("inputDataDisturbance"))

#' @rdname inputDataDisturbance
setMethod(
  "inputDataDisturbance", signature(landCover = "RasterLayer"), 
  function(landCover, linFeat, projectPoly, 
           natDist=NULL,anthroDist = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {
    
    charIn <-  sapply(list(landCover,natDist, anthroDist,  
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 
    
    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }
    
    
    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
    # prep projectPoly
    projPolyOut <- prepProjPoly(projectPoly, landCover, winArea = bufferWidth,
                                padProjPoly)
    projectPoly <- projPolyOut[["projectPoly"]]
    projectPolyOrig <- projPolyOut[["projectPolyOrig"]]
    rm(projPolyOut)

    # combine linFeat
    if(inherits(linFeat, "list")){
      linFeat <- combineLinFeat(linFeat)
    }

    if(is(linFeat, "Spatial")){
      linFeat <- sf::st_as_sf(linFeat)
    } 
    
    if(is(linFeat, "Raster")){
      rastLst <- lst(natDist, anthroDist, linFeat)

    } else {
      rastLst <- lst(natDist, anthroDist)
      linFeat <- checkAlign(linFeat, landCover, "linFeat", "landCover")
    }

    # check alignment of all raster layers
    rastLst <- prepRasts(rastLst, landCover, projectPoly)
    
    if(is.null(natDist)){
      rastLst$natDist <- raster(matrix(NA))
    }
    
    if(is.null(anthroDist)){
      rastLst$anthroDist <- raster(matrix(NA))
    }
    
    if(is(linFeat, "Raster")){
      linFeat <- rastLst$linFeat
    }
    
    return(new("DisturbanceMetrics", rastLst$landCover, rastLst$natDist, 
               rastLst$anthroDist, 
               linFeat, projectPolyOrig,  
               processedData = raster(matrix(NA)), 
               disturbanceMetrics = data.frame(),
               attributes = list(bufferWidth = bufferWidth,
                                 padProjPoly = padProjPoly, padFocal = padFocal)))
  })

#' @rdname inputDataDisturbance
setMethod(
  "inputDataDisturbance", signature(landCover = "character"), 
  function(landCover, linFeat, projectPoly,
           natDist=NULL,anthroDist = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(natDist, anthroDist, projectPoly)
      
      linFeat <- combineLinFeat(linFeat)
      
    } else {
      indata <- lst(landCover,natDist, anthroDist, linFeat,
                    projectPoly)
    }
    
    indata <- loadFromFile(indata)
    
    if(is.character(linFeat)){
      linFeat <- indata$linFeat
    }
    
    return(inputDataDisturbance(landCover=indata$landCover, 
                                natDist = indata$natDist, anthroDist = indata$anthroDist, 
                                linFeat = linFeat, 
                                projectPoly = indata$projectPoly, 
                                bufferWidth = bufferWidth, 
                                padProjPoly = padProjPoly, padFocal = padFocal))
    
  })