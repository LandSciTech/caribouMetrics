#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating caribou habitat use.
#' 

#' @noRd
setGeneric("inputData", function(landCover, esker, linFeat, projectPoly, ...) standardGeneric("inputData"))

setMethod(
  "inputData", signature(landCover = "RasterLayer"), 
  function(landCover, esker, linFeat, projectPoly, caribouRange,
           coefTable,
           natDist = NULL, 
           anthroDist = NULL, 
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE,
           padFocal = FALSE, ptDensity = 1, 
           tmplt = NULL) {

    charIn <-  sapply(list(landCover, esker, 
                           natDist, anthroDist,  
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 

    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }

    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
    if(!inherits(caribouRange, "data.frame")){
      caribouRange <- data.frame(Range = caribouRange, 
                                 coefRange = caribouRange, 
                                 stringsAsFactors = FALSE)
    } else {
      if(any(names(caribouRange) != c("Range", "coefRange"))){
        stop("If caribouRange is a data.frame the column names",
             " must be Range and coefRange", call. = FALSE)
      }
    }

    # require column called Range
    if(nrow(projectPoly) == 1 & !"Range" %in% names(projectPoly)){
      projectPoly <- projectPoly %>% mutate(Range = caribouRange$Range)
    } else {
      if(!"Range" %in% names(projectPoly)){
        stop("projectPoly must have a column Range that corresponds to the ",
             "caribouRange$Range column", call. = FALSE)
      } else {
        if(!all(projectPoly$Range %in% caribouRange$Range)){
          stop("All values of in projectPoly$Range must have matching",
               " values in caribouRange$Range", call. = FALSE)
        }
      }
    }
    
    .checkInputs(caribouRange, winArea, landCover, coefTable)
    
    # prep projectPoly
    projPolyOut <- prepProjPoly(projectPoly, landCover, winArea, padProjPoly)
    projectPoly <- projPolyOut[["projectPoly"]]
    projectPolyOrig <- projPolyOut[["projectPolyOrig"]]
    rm(projPolyOut)

    # check alignment of all raster layers esker and linFeat added if needed
    rastLst <- lst(natDist, anthroDist, esker, linFeat)
    
    rastLst <- purrr::keep(rastLst, is, "RasterLayer")
    
    rastLst <- prepRasts(rastLst, landCover, projectPoly, tmplt, 
                          useTmplt = c("esker", "linFeat"))
    
    if(is.null(natDist)){
      rastLst$natDist <- raster(matrix(NA))
    }
        
    if(is.null(anthroDist)){
      rastLst$anthroDist <- raster(matrix(NA))
    }
    
    if(is.null(tmplt)){
      tmplt <- raster(rastLst$landCover) %>% raster::`res<-`(c(400, 400))
    } 
    
    
    # rasterize eskers
    if(inherits(esker, "sf")){
      esker <- checkAlign(esker, projectPoly, "esker", "projectPoly")
      
      esker <- rasterizeLineDensity(esker, tmplt)
      
    } 
    
    # rasterize linFeat
    if(inherits(linFeat, "list")){
      linFeat <- combineLinFeat(linFeat)
    }
    
    if(is(linFeat, "Spatial")){
      linFeat <- sf::st_as_sf(linFeat)
    } 
    
    if(inherits(linFeat, "sf")){
      
      linFeat <- checkAlign(linFeat, projectPoly, "linFeat", "projectPoly")
      
      linFeat <- rasterizeLineDensity(linFeat, tmplt, ptDensity)
    } 
    
    if(is.null(rastLst$linFeat)){
      rastLst <- c(rastLst, linFeat = linFeat)
    }
    
    if(is.null(rastLst$esker)){
      rastLst <- c(rastLst, esker = esker)
    }
    
    if(!is.null(linFeatSave)){
      raster::writeRaster(rastLst[["linFeat"]], 
                          linFeatSave, overwrite = TRUE)
      rastLst$linFeat <- raster(linFeatSave)
    }
    
    if(!is.null(eskerSave)){
      raster::writeRaster(rastLst[["esker"]], eskerSave, overwrite = TRUE)
      rastLst$esker <- raster(eskerSave)
    }
    
    return(new("CaribouHabitat", rastLst$landCover, rastLst$esker, 
               rastLst$natDist, rastLst$anthroDist,
               rastLst$linFeat, projectPolyOrig,  
               processedData = raster(matrix(NA)), 
               habitatUse = raster(matrix(NA)),
               attributes = list(caribouRange = caribouRange, winArea = winArea,
                                 padProjPoly = padProjPoly, padFocal = padFocal, 
                                 tmplt = tmplt)))
})

#' @noRd
setMethod(
  "inputData", signature(landCover = "character"), 
  function(landCover, esker, linFeat,  projectPoly, caribouRange, coefTable,
           natDist = NULL, 
           anthroDist = NULL,
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE, 
           padFocal = FALSE, 
           ptDensity = 1) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(landCover, esker, 
                    natDist, anthroDist,  
                    projectPoly)
      
      linFeat <- combineLinFeat(linFeat)
      
    } else {
      indata <- lst(landCover, esker,
                    natDist, anthroDist, linFeat,
                    projectPoly)
    }
    
    indata <- loadFromFile(indata)
    
    if(is.character(linFeat)){
      linFeat <- indata$linFeat
    }
   
    indata$landCover <- indata$landCover %>% reclassPLC()
  
    return(inputData(landCover = indata$landCover, esker = indata$esker, 
                     natDist = indata$natDist, anthroDist = indata$anthroDist, 
                     linFeat = linFeat, 
                     projectPoly = indata$projectPoly, 
                     caribouRange = caribouRange, eskerSave = eskerSave, 
                     linFeatSave = linFeatSave, winArea = winArea, 
                     padProjPoly = padProjPoly, padFocal = padFocal, 
                     ptDensity = ptDensity, coefTable = coefTable))
    
  })