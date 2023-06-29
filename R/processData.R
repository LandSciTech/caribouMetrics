
#' @include AAAClassDefinitions.R
NULL

#' Process the data after it is loaded and checked by inputData
#' @noRd
setGeneric("processData", function(inData, newData, ...) standardGeneric("processData"))

# Method for initial processing when no new data is provided
#' @noRd
#' @importFrom raster plot
setMethod(
  "processData", signature(inData = "CaribouHabitat", newData = "missing"), 
  function(inData) {
 
    tmplt <- inData@attributes$tmplt
    
    # resample to 16ha and flag cells with >35% disturbance
    expVars <- applyDist(inData@landCover, inData@natDist, inData@anthroDist, 
                         tmplt)

    # create raster attribute table
    inData@landCover <- terra::categories(
      inData@landCover, layer = 1,
      left_join(terra::unique(inData@landCover) %>% setNames("landCover"),
                resTypeCode, by = c(landCover = "code")) %>%
        setNames(c("ID", "ResourceType")))
    
    # Resample linFeat and esker to 16 ha
    if(any(terra::res(inData@linFeat) < terra::res(tmplt)[1])){
      message("resampling linFeat to match landCover resolution")
      
      inData@linFeat <- terra::resample(inData@linFeat, tmplt,
                                         method = "bilinear")
    }
    if(!terra::compareGeom(expVars[[1]], inData@linFeat, stopOnError = FALSE)){
      stop("linFeat could not be aligned with landCover.",
           " Please provide a higher resolution raster or a shapefile",
           call. = FALSE)
    }
    
    if(any(terra::res(inData@esker) < terra::res(tmplt)[1])){
      message("resampling esker to match landCover resolution")
      
      inData@esker <- terra::resample(inData@esker, tmplt,
                                       method = "bilinear")
    }
    if(!terra::compareGeom(expVars[[1]], inData@esker, stopOnError = FALSE)){
      stop("esker could not be aligned with landCover.",
           " Please provide a higher resolution raster or a shapefile",
           call. = FALSE)
    }
    
    # Get window area from table b/c some models used different sizes
    # if(is.null(inData@attributes$winArea)){
    #   inData@attributes$winArea <- coefTable %>% 
    #     filter(.data$Range %in% inData@attributes$caribouRange$coefRange) %>% 
    #     pull(.data$WinArea) %>% 
    #     unique()
    # }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/terra::res(expVars[[1]])[1]) %>% 
      round(digits = 0)*terra::res(expVars[[1]])[1] %>% 
      round()
    
    # calculate moving window average for all explanatory variables
    expVars <- c(expVars, inData@linFeat, inData@esker) 
    
    layernames <- c(resTypeCode %>% arrange(.data$code) %>%
                      pull(.data$ResourceType) %>% 
                      as.character(),
                    "TDENLF", "ESK")
    
    if(!inData@attributes$padProjPoly){
      expVars <- terra::mask(expVars, inData@projectPoly)
    }
    
    message("Applying moving window.")
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                               nms = layernames, 
                               pad = inData@attributes$padFocal, usePfocal = FALSE)
    
    inData@processedData <- c(expVars, 
                              expVars[["MIX"]] + expVars[["DEC"]], 
                              makeDummyRast(expVars[[1]],1)) %>% 
      `names<-`(c(names(expVars), "LGMD", "CONST"))
    
    return(inData)
    
  })

# method to update processed data when new data is supplied avoids unnecessary
# steps if new data is not provided for all variables

#' @noRd
setMethod(
  "processData", 
  signature(inData = "CaribouHabitat", newData = "list"), 
  function(inData, newData) {
    
    if(is.null(names(newData)) ||
       !all(names(newData) %in% c("landCover", "natDist", "anthroDist",
                                  "linFeat" ))){
      stop("newData must be a named list with any of the names ",
           "landCover ", "natDist ", "anthroDist ",  "linFeat ", 
           "incorrect names provided are:  ", 
           names(newData)[which(!names(newData) %in% 
                                  c("landCover", "natDist", "anthroDist",
                                    "linFeat" ))], call. = FALSE)
    }
    
    tmplt <- inData@attributes$tmplt
    
    if(!is.null(newData$linFeat)){
      # rasterize linFeat
      if(inherits(newData$linFeat, "list")){
        newData$linFeat <- combineLinFeat(newData$linFeat)
      }
      
      if(inherits(newData$linFeat, "sf")){
        newData$linFeat <- checkAlign(newData$linFeat, inData@projectPoly,
                                      "linFeat", "projectPoly")
        
        inData@linFeat <- rasterizeLineDensity(newData$linFeat, tmplt)
      }
      
    }
    
    newData <- rapply(newData, f = terra::rast, classes = "RasterLayer", how = "replace")
    
    if(!all(sapply(newData[setdiff(names(newData), "linFeat")], is, "SpatRaster"))){
      stop("All data supplied in the newData list must be SpatRaster or RasterLayer objects", 
           call. = FALSE)
    }
    
    # check alignment of new data with projectPoly except for linFeat
    rastCRSMatch <- lapply(newData[which(names(newData) != "linFeat")],
                           terra::compareGeom,
                           y = inData@landCover, crs = TRUE, res = TRUE, ext = FALSE, 
                           rowcol = FALSE, stopOnError = FALSE) %>%
      unlist() %>% all()
    
    if(!rastCRSMatch){
      stop("all raster data sets must have matching resolution", call. = FALSE)
    }
    
    newData <- purrr::map2(newData, names(newData), 
                           ~cropIf(.x, inData@projectPoly, .y, "projectPoly"))
    
    if(any(names(newData) %in% c("landCover", "natDist", "anthroDist" ))){
      
      # TODO: should you be able to update only Nat or Anthro dist?
      if(is.null(newData$landCover)){
        newData$landCover <- inData@landCover
        warning("existing landCover being updated with new disturbance", 
                call. = FALSE)
      }
      
      if(is.null(newData$natDist)){
        newData$natDist <- inData@natDist
        warning("existing natDist being used to update landCover", 
                call. = FALSE)
      }
      if(is.null(newData$anthroDist)){
        newData$anthroDist <- inData@anthroDist
        warning("existing anthroDist being used to update landCover", 
                call. = FALSE)
      }
      
      expVars <- applyDist(newData$landCover, natDist = newData$natDist, 
                           anthroDist = newData$anthroDist, tmplt = tmplt)
    }
    
    # resample linFeat if provided
    if(inherits(newData$linFeat, "SpatRaster")){
      message("resampling linFeat to match landCover resolution")
      newData$linFeat <- terra::resample(newData$linFeat, tmplt, 
                                          method = "bilinear")
      inData@linFeat <- newData$linFeat
    } 
    
    # calculate moving window average for changed explanatory variables
    if(all(c("linFeat", "landCover") %in% names(newData))){
      expVars <- c(expVars, inData@linFeat)
      
      layernames <- c(resTypeCode %>% arrange(.data$code) %>%
                        pull(.data$ResourceType) %>% as.character(), "TDENLF")
    } else if ("landCover" %in% names(newData)){
      layernames <- resTypeCode %>% arrange(.data$code) %>%
        pull(.data$ResourceType) %>% as.character()
    } else {
      expVars <- inData@linFeat
      
      layernames <- "TDENLF"
    }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/terra::res(expVars[[1]])[1]) %>% 
      round(digits = 0)*terra::res(expVars[[1]])[1] %>% 
      round()
    
    message("Applying moving window.")
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                               nms = layernames, 
                               pad = inData@attributes$padFocal, 
                               usePfocal = FALSE)
    
    if(all(c("MIX","DEC") %in% names(expVars))){
      expVars <- c(expVars, (expVars[["MIX"]] + expVars[["DEC"]])) %>% 
        `names<-`(c(names(expVars), "LGMD"))
    }
    
    notUpdated <- which(!names(inData@processedData) %in% names(expVars)) 
    
    expVars <- c(expVars, inData@processedData[[notUpdated]])
    # order output stack to match 
    inData@processedData <- expVars[[names(inData@processedData)]]
    
    
    return(inData)
  })

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
                                pad = inData@attributes$padFocal, 
                                offset = FALSE, usePfocal = FALSE)
  
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
  fire_excl_anthro <- terra::lapp(c(natDist, anthroDist),
                                  fun = function(x, y){
                                    (x - y) > 0 
                                  })
  
  all <- (anthroDist + natDist) > 0
  
  outStack <- c(anthroDist, natDist, all, fire_excl_anthro)
  
  rm(anthroDist, natDist, all, fire_excl_anthro)
  
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
                      "fire_excl_anthro")
  return(outStack)
}

calcDM <- function(inData, isPercent){
  message("calculating disturbance metrics")

  rr <- terra::extract(inData@processedData, terra::vect(inData@projectPoly), 
                      fun = "mean", bind = TRUE) %>% 
    as.data.frame() %>% 
    select(all_of(names(inData@processedData)), everything()) %>% 
    mutate(zone = 1:n(), .before = everything())
  
  if (isPercent == TRUE) {
    rr$Anthro <- rr$Anthro * 100
    rr$Fire <- rr$Fire * 100
    rr$Total_dist <- rr$Total_dist * 100
    rr$fire_excl_anthro <- rr$fire_excl_anthro * 100
  }
  return(rr)
}