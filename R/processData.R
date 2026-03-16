
#' @include AAAClassDefinitions.R
NULL

#' Process the data after it is loaded and checked by inputData
#' @noRd
setGeneric("processData", function(inData, newData, ...) standardGeneric("processData"))

# Method for initial processing when no new data is provided
#' @noRd
setMethod(
  "processData", signature(inData = "DisturbanceMetrics", newData = "missing"), 
  function(inData, linBuffMethod = "raster", isPercent) {
    #inData=dm
    
    anthroDist <- inData@anthroDist
    natDist <- inData@natDist
    landCover <- inData@landCover
    
    # check natDist is real if not make dummy
    if(terra::ncell(natDist) == 1){
      natDist <- makeDummyRast(landCover)
    }
    
    natDist <- terra::classify(natDist, cbind(NA, 0))

    anthroDist <- processAnthroDM(anthroDist, inData@linFeat[[1]], landCover, 
                                  linBuffMethod = linBuffMethod,
                                  inData = inData)
    
    
    inData@processedData <- calcDMSpatial(anthroDist, natDist, landCover, 
                                          inData)
    
    ####### Range summaries
    inData@disturbanceMetrics <- calcDM(inData, isPercent)
    
    return(inData)
  })

processAnthroDM <- function(anthroDist, linFeat, landCover, 
                            linBuffMethod = "raster", inData){
  
  # check anthroDist is real if not make dummy
  if(terra::ncell(anthroDist) == 1){
    anthroDist <- makeDummyRast(landCover)
  }
  anthroDist <- terra::classify(anthroDist, cbind(NA, 0))
  
  ############ Buffer anthropogenic disturbance
  #To speed calculations, include points from linFeat in raster.
  message("buffering anthropogenic disturbance")
  
  anthroDist <- anthroDist > 0
  
  
  if(!is(linFeat, "SpatRaster")){
    if(any(c("POINT", "MULTIPOINT") %in% sf::st_geometry_type(linFeat))){
      lfLine <- sf::st_collection_extract(linFeat, "LINESTRING")
      lfPt <- sf::st_collection_extract(linFeat, "POINT")
    } else {
      lfPt <- slice(linFeat, 0)
      lfLine <- linFeat
    }
    if(linBuffMethod == "raster"){
      lfRas <- terra::rasterize(terra::vect(lfLine), y = anthroDist, 
                                touches = TRUE, background = 0) 
      
      anthroDist <- terra::mask(anthroDist, lfRas, inverse = TRUE,
                                maskvalue = 0, updatevalue = 1)
    } 
    if(nrow(lfPt) > 0){
      lfRasPt <- terra::rasterize(terra::vect(lfPt), anthroDist, 
                                field = 1, background = 0) 
      
      anthroDist <- terra::mask(anthroDist, lfRasPt, inverse = TRUE,
                                maskvalue = 0, updatevalue = 1)
    } 
  }else{
    lfRas <- linFeat
    
    anthroDist <- terra::mask(anthroDist, lfRas, inverse = TRUE,
                               maskvalue = 0, updatevalue = 1)
  }
  
  # window radius 
  winRad <- (inData@attributes$bufferWidth/terra::res(anthroDist[[1]])[1]) %>% 
    round(digits = 0)*terra::res(anthroDist[[1]])[1] %>% 
    round()
  
  if(!inData@attributes$padProjPoly){
    anthroDist <- terra::mask(anthroDist, terra::vect(inData@projectPoly))
  }
  
  anthroDist <- movingWindowAvg(rast = anthroDist, radius = winRad,
                                nms = "ANTHRO",
                                naExternal = ifelse(inData@attributes$padFocal,
                                                    "expand", "NA"), 
                                naInternal = "interpolate", 
                                offset = FALSE)
  
  anthroDist <- anthroDist > 0 
  
  if(!is(linFeat, "SpatRaster") && linBuffMethod != "raster"){
    ############## Buffer linear features
    message("buffering linear features")
    
    #Note points were included with raster above.
    lf <- linFeat %>% sf::st_collection_extract("LINESTRING")
    
    # simplify could speed this up but using raster method instead
    #lf <- st_simplify(lf, dTolerance = res(anthroDist)[1])
    
    linBuff <- sf::st_buffer(lf, inData@attributes$bufferWidth)
    
    if(sf::st_geometry_type(linBuff, by_geometry = FALSE) == "GEOMETRY"){
      linBuff <- sf::st_collection_extract(linBuff, "POLYGON")
    }

    linBuff <- terra::rasterize(linBuff, anthroDist, field = 1, background = 0)

    anthroDist <- (linBuff + anthroDist) > 0
  }
  return(anthroDist)
}

calcDMSpatial <- function(anthroDist, natDist, landCover, inData){
  Fire_excl_anthro <- terra::lapp(c(natDist, anthroDist),
                                  fun = function(x, y){
                                    (x - y) > 0 
                                  })
  
  all <- (anthroDist + natDist) > 0
  
  outStack <- c(anthroDist, natDist, all, Fire_excl_anthro)
  
  rm(anthroDist, natDist, all, Fire_excl_anthro)
  
  pp = terra::rasterize(inData@projectPoly, outStack[[1]])
  
  #set NAs from landcover
  toNA <- terra::lapp(c(landCover, pp),
                          fun = function(x, y){
                            ifelse(is.na(x) |x == 0| is.na(y), NA, 1)
                          })
  
  outStack <- terra::mask(outStack, toNA)
  rm(toNA, pp)
  
  names(outStack) = c("Anthro", 
                      "Fire", 
                      "Total_dist",
                      "Fire_excl_anthro")
  return(outStack)
}

calcDM <- function(inData, isPercent){
  message("calculating disturbance metrics")

  rr <- terra::extract(inData@processedData, terra::vect(inData@projectPoly), 
                      fun = "mean", bind = TRUE, na.rm = TRUE) %>% 
    as.data.frame() %>% 
    select(all_of(names(inData@processedData)), everything()) %>% 
    mutate(zone = 1:n(), .before = everything())
  
  if (isPercent == TRUE) {
    rr$Anthro <- rr$Anthro * 100
    rr$Fire <- rr$Fire * 100
    rr$Total_dist <- rr$Total_dist * 100
    rr$Fire_excl_anthro <- rr$Fire_excl_anthro * 100
  }
  return(rr)
}