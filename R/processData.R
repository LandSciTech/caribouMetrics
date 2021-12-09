
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
    inData@landCover <- raster::ratify(inData@landCover)
    
    levels(inData@landCover)[[1]] <- left_join(raster::levels(inData@landCover)[[1]], 
                                               resTypeCode, by = c(ID = "code"))
    
    # Resample linFeat and esker to 16 ha
    if(any(raster::res(inData@linFeat) < raster::res(tmplt)[1])){
      message("resampling linFeat to match landCover resolution")
      
      inData@linFeat <- raster::resample(inData@linFeat, tmplt,
                                         method = "bilinear")
    }
    if(!compareRaster(expVars[[1]], inData@linFeat, stopiffalse = FALSE)){
      stop("linFeat could not be aligned with landCover.",
           " Please provide a higher resolution raster or a shapefile",
           call. = FALSE)
    }
    
    if(any(raster::res(inData@esker) < raster::res(tmplt)[1])){
      message("resampling esker to match landCover resolution")
      
      inData@esker <- raster::resample(inData@esker, tmplt,
                                       method = "bilinear")
    }
    if(!compareRaster(expVars[[1]], inData@esker, stopiffalse = FALSE)){
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
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1] %>% 
      round()
    
    # calculate moving window average for all explanatory variables
    expVars <- expVars %>% 
      addLayer(inData@linFeat) %>% 
      addLayer(inData@esker) 
    
    layernames <- c(resTypeCode %>% arrange(.data$code) %>%
                      pull(.data$ResourceType) %>% 
                      as.character(),
                    "TDENLF", "ESK")
    
    if(!inData@attributes$padProjPoly){
      expVars <- raster::mask(expVars, inData@projectPoly)
    }
    
    message("Applying moving window.")
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                               nms = layernames, 
                               pad = inData@attributes$padFocal)
    
    inData@processedData <- expVars %>% 
      raster::addLayer(expVars[["MIX"]] + expVars[["DEC"]]) %>% 
      raster::addLayer(raster::init(expVars[[1]], 
                                    fun = function(x){rep(1, x)})) %>% 
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
    
    if(!all(sapply(newData[setdiff(names(newData), "linFeat")], is, "RasterLayer"))){
      stop("All data supplied in the newData list must be RasterLayer objects", 
           call. = FALSE)
    }
    
    # check alignment of new data with projectPoly except for linFeat
    if(!do.call(raster::compareRaster, 
                c(`names<-`(newData[which(names(newData) != "linFeat")], NULL),
                  list(inData@landCover,inData@landCover,
                       res = TRUE, extent = FALSE, 
                       rowcol = FALSE, stopiffalse = FALSE)))){
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
    if(inherits(newData$linFeat, "Raster")){
      message("resampling linFeat to match landCover resolution")
      newData$linFeat <- raster::resample(newData$linFeat, tmplt, 
                                          method = "bilinear")
      inData@linFeat <- newData$linFeat
    } 
    
    # calculate moving window average for changed explanatory variables
    if(all(c("linFeat", "landCover") %in% names(newData))){
      expVars <- expVars %>% addLayer(inData@linFeat)
      
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
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1] %>% 
      round()
    
    message("Applying moving window.")
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                               nms = layernames, 
                               pad = inData@attributes$padFocal)
    
    if(all(c("MIX","DEC") %in% names(expVars))){
      expVars <- expVars %>% 
        raster::addLayer(expVars[["MIX"]] + expVars[["DEC"]]) %>% 
        `names<-`(c(names(expVars), "LGMD"))
    }
    
    notUpdated <- which(!names(inData@processedData) %in% names(expVars)) 
    
    expVars <- raster::stack(expVars, 
                             inData@processedData[[notUpdated]])
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
    if(raster::ncell(natDist) == 1){
      natDist <- raster::init(landCover, 
                              fun = function(x){rep(0, x)}, 
                              filename = raster::rasterTmpFile())
    }
    
    natDist <- reclassify(natDist, cbind(NA, 0))
    
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
  if(raster::ncell(anthroDist) == 1){
    anthroDist <- raster::init(landCover, 
                               fun = function(x){rep(0, x)}, 
                               filename = raster::rasterTmpFile())
  }
  anthroDist <- reclassify(anthroDist, cbind(NA, 0))
  
  ############ Buffer anthropogenic disturbance
  #To speed calculations, include points from linFeat in raster.
  message("buffering anthropogenic disturbance")
  
  anthroDist <- anthroDist > 0
  
  
  if(!is(linFeat, "RasterLayer")){
    if(linBuffMethod == "raster"){
      # rasterize roads to template
      tmplt <- stars::st_as_stars(sf::st_bbox(anthroDist), nx = raster::ncol(anthroDist),
                                  ny = raster::nrow(anthroDist), values = 0)
      
      lfRas <- stars::st_rasterize(linFeat[attr(linFeat, "sf_geometry")],
                                   template = tmplt,
                                   options = "ALL_TOUCHED=TRUE") %>%
        as("Raster")
      
      anthroDist <- raster::mask(anthroDist, lfRas, inverse = TRUE,
                                 maskvalue = 0, updatevalue = 1)
    } else {
      lfPt <- linFeat %>% dplyr::filter(st_is(linFeat , "POINT"))
      
      if(nrow(lfPt) > 0){
        lfPt <- as(lfPt, "Spatial")
        
        lfRas <- raster::rasterize(lfPt, anthroDist, field = 1, background = 0) 
        
        anthroDist <- raster::mask(anthroDist, lfRas, inverse = TRUE,
                                   maskvalue = 0, updatevalue = 1)
      } 
    }
  }else{
    lfRas <- linFeat
    
    anthroDist <- raster::mask(anthroDist, lfRas, inverse = TRUE,
                               maskvalue = 0, updatevalue = 1)
  }
  
  # window radius 
  winRad <- (inData@attributes$bufferWidth/res(anthroDist[[1]])[1]) %>% 
    round(digits = 0)*res(anthroDist[[1]])[1] %>% 
    round()
  
  if(!inData@attributes$padProjPoly){
    anthroDist <- raster::mask(anthroDist, inData@projectPoly)
  }
  
  anthroDist <- movingWindowAvg(rast = anthroDist, radius = winRad,
                                nms = "ANTHRO", 
                                pad = inData@attributes$padFocal, 
                                offset = FALSE)
  
  anthroDist <- anthroDist > 0 
  
  if(!is(linFeat, "RasterLayer") && linBuffMethod != "raster"){
    ############## Buffer linear features
    message("buffering linear features")
    
    #Note points were included with polygons above.
    lf <- linFeat %>% dplyr::filter(!st_is(linFeat , "POINT"))
    
    # simplify could speed this up but using raster method instead
    #lf <- st_simplify(lf, dTolerance = res(anthroDist)[1])
    
    linBuff <- st_buffer(lf, inData@attributes$bufferWidth)
    
    if(st_geometry_type(linBuff, by_geometry = FALSE) == "GEOMETRY"){
      linBuff <- st_collection_extract(linBuff, "POLYGON")
    }
    
    # faster rasterization
    if(requireNamespace("fasterize", quietly = TRUE)){
      linBuff <- fasterize::fasterize(linBuff, anthroDist)
    } else {
      message("To speed up install fasterize package")
      linBuff <- raster::rasterize(linBuff, anthroDist)
    }
    
    linBuff <- linBuff > 0
    linBuff <- reclassify(linBuff, cbind(NA, 0))
    anthroDist <- (linBuff + anthroDist) > 0
  }
  return(anthroDist)
}

calcDMSpatial <- function(anthroDist, natDist, landCover, inData){
  fire_excl_anthro <- raster::overlay(natDist, 
                                      anthroDist,
                                      fun = function(x, y){
                                        x <- x - y 
                                        x <- x > 0
                                        return(x)
                                      })
  
  all <- (anthroDist + natDist) > 0
  
  outStack <- raster::stack(anthroDist, 
                            natDist, 
                            all, 
                            fire_excl_anthro)
  
  rm(anthroDist, natDist, all, fire_excl_anthro)
  
  if(requireNamespace("fasterize", quietly = TRUE)){
    pp = fasterize::fasterize(inData@projectPoly, outStack[[1]])
  } else {
    message("To speed up install fasterize package")
    pp = raster::rasterize(inData@projectPoly, outStack[[1]])
  }
  
  #set NAs from landcover
  toNA <- raster::overlay(landCover, pp,
                          fun = function(x, y){
                            ifelse(is.na(x) |x == 0| is.na(y), NA, 1)
                          })
  
  outStack <- raster::overlay(outStack, toNA, fun = function(x, y){x * y})
  rm(toNA, pp)
  
  names(outStack) = c("Anthro", 
                      "Fire", 
                      "Total_dist",
                      "fire_excl_anthro")
  return(outStack)
}

calcDM <- function(inData, isPercent){
  message("calculating disturbance metrics")
  
  # 4 times faster and 7 times less memory needed vs using raster::zonal
  rr <- exactextractr::exact_extract(inData@processedData, inData@projectPoly, 
                                     fun = "mean", append_cols = TRUE) %>% 
    select(contains("mean."), everything()) %>% 
    rename_with(~gsub("mean\\.", "", .x)) %>% 
    mutate(zone = 1:n(), .before = everything())
  
  if (isPercent == TRUE) {
    rr$Anthro <- rr$Anthro * 100
    rr$Fire <- rr$Fire * 100
    rr$Total_dist <- rr$Total_dist * 100
    rr$fire_excl_anthro <- rr$fire_excl_anthro * 100
  }
  return(rr)
}