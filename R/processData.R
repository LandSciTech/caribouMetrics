
#' @include AAAClassDefinitions.R
NULL

#' Process the data after it is loaded and checked by inputData

#' @export
setGeneric("processData", function(inData, newData, ...) standardGeneric("processData"))

# Method for initial processing when no new data is provided
#' @rdname processData
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

      inData@esker <- raster::resample(inData@esker, tmplt,
                                       method = "bilinear")
    }
    if(!compareRaster(expVars[[1]], inData@esker, stopiffalse = FALSE)){
      stop("esker could not be aligned with landCover.",
           " Please provide a higher resolution raster or a shapefile",
           call. = FALSE)
    }

    # Get window area from table b/c some models used different sizes
    if(is.null(inData@attributes$winArea)){
      inData@attributes$winArea <- coefTable %>% 
        filter(Range %in% inData@attributes$caribouRange$coefRange) %>% 
        pull(WinArea) %>% 
        unique()
    }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1] %>% 
      round()
    
    # calculate moving window average for all explanatory variables
    expVars <- expVars %>% 
      addLayer(inData@linFeat) %>% 
      addLayer(inData@esker) 
    
    layernames <- c(resTypeCode %>% arrange(code) %>%
                      pull(ResourceType) %>% 
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
      raster::addLayer(raster::setValues(expVars[[1]], 1)) %>% 
      `names<-`(c(names(expVars), "LGMD", "CONST"))

    return(inData)
    
  })

# method to update processed data when new data is supplied avoids unnecessary
# steps if new data is not provided for all variables

#' @rdname processData
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
        #tmplt <- raster(inData@landCover) %>% raster::`res<-`(c(400, 400))
        inData@linFeat <- rasterizeLineDensity(newData$linFeat, tmplt)
      }
      
    }

    if(!all(sapply(newData, is, "RasterLayer"))){
      stop("All data supplied in the newData list must be RasterLayer objects", 
           call. = FALSE)
    }
    
    # check alignment of new data with landCover except for linFeat
    if(!do.call(raster::compareRaster, 
                c(`names<-`(newData[which(names(newData) != "linFeat")], NULL),
                  list(inData@landCover,inData@landCover,
                       res = TRUE, extent = FALSE, 
                       rowcol = FALSE, stopiffalse = FALSE)))){
      stop("all raster data sets must have matching resolution", call. = FALSE)
    }
    
    newData <- purrr::map2(newData, names(newData), 
                ~cropIf(.x, inData@landCover, .y, "landCover"))
    
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
      newData$linFeat <- raster::resample(newData$linFeat, inData@linFeat, 
                                          method = "bilinear")
      inData@linFeat <- newData$linFeat
    } 

    # calculate moving window average for changed explanatory variables
    if(all(c("linFeat", "landCover") %in% names(newData))){
      expVars <- expVars %>% addLayer(inData@linFeat)
      
      layernames <- c(resTypeCode %>% arrange(code) %>%
                        pull(ResourceType) %>% as.character(), "TDENLF")
    } else if ("landCover" %in% names(newData)){
      layernames <- resTypeCode %>% arrange(code) %>%
        pull(ResourceType) %>% as.character()
    } else {
      expVars <- inData@linFeat
      
      layernames <- "TDENLF"
    }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(inData@attributes$winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1]
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,  
                               nms = layernames, 
                               pad = inData@attributes$padFocal)
    
    if(all(c("MIX","DEC") %in% names(expVars))){
      expVars <- expVars %>% 
        raster::addLayer(expVars[["MIX"]] + expVars[["DEC"]]) %>% 
        `names<-`(c(names(expVars), "LGMD"))
    }
    
    notUpdated <- which(!names(inData@processedData) %in% names(expVars)) 
    
    inData@processedData <- raster::stack(expVars, 
                                          inData@processedData[[notUpdated]])
    
    return(inData)
  })

# Method for initial processing when no new data is provided
#' @rdname processData
#' @importFrom raster plot
setMethod(
  "processData", signature(inData = "DisturbanceMetrics", newData = "missing"), 
  function(inData) {
    #inData=dm

    anthroDist <- inData@anthroDist
    natDist <- inData@natDist
    landCover <- inData@landCover
    
    # check anthroDist and natDist are real if not make dummy
    if(length(raster::unique(anthroDist)) == 0){
      anthroDist <- landCover
      anthroDist[] <- 0
    }
    if(length(raster::unique(natDist)) == 0){
      natDist <- landCover
      natDist[] <- 0
    }
    
    anthroDist[is.na(anthroDist)] = 0
    natDist[is.na(natDist)] = 0
    ############
    #Buffer anthropogenic disturbance
    #To speed calculations, include points from linFeat in raster.
    message("buffering anthropogenic disturbance")
    
    expVars <- (anthroDist) > 0
    

    if(is(inData@linFeat[[1]], "RasterLayer")){
      
      expVars[inData@linFeat[[1]] > 0] = 1
    }else{
      lfPt <- inData@linFeat[[1]] %>% dplyr::filter(st_is(. , "POINT"))
      
      if(nrow(lfPt) > 0){

        lfPt <- as(lfPt, "Spatial")
        
        lfR = raster::rasterize(lfPt, expVars, field = "ID")    
        expVars[!is.na(lfR)] = 1
      }  
    }

    # window radius 
    winRad <- (inData@attributes$bufferWidth/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1] %>% 
      round()
    
    layernames <- c("ANTHRO")
    
    if(!inData@attributes$padProjPoly){
      expVars <- raster::mask(expVars, inData@projectPoly)
    }

    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                               nms = layernames, 
                               pad = inData@attributes$padFocal, 
                               offset = FALSE)
    
    expVars <- expVars > 0 

    if(!is(inData@linFeat[[1]], "RasterLayer")){
      ##############
      #Buffer linear features
      message("buffering linear features")
      
      #Note points were included with polygons above.
      lf <- inData@linFeat[[1]] %>% dplyr::filter(!st_is(. , "POINT"))
      
      linBuff <- st_buffer(lf,inData@attributes$bufferWidth)
      
      if(st_geometry_type(linBuff, by_geometry = FALSE) == "GEOMETRY"){
        linBuff <- st_collection_extract(linBuff, "POLYGON")
      }
      
      # faster rasterization
      if(requireNamespace("fasterize", quietly = TRUE)){
        linBuff <- fasterize::fasterize(linBuff, expVars)
      } else {
        message("To speed up install fasterize package")
        linBuff <- raster::rasterize(linBuff, expVars)
      }
      
      linBuff <- linBuff > 0 
      linBuff[is.na(linBuff)] = 0
      anthroBuff <- (linBuff + expVars) > 0
    }else{
      anthroBuff <- expVars > 0
    }
    
    fire_excl_anthro <- raster::overlay(natDist, 
                                        anthroBuff,
                                        fun = function(x, y){
                                          x <- x - y 
                                          x <- x > 0
                                          return(x)
                                        })
    
    all <- (anthroBuff + natDist) > 0
    
    
    outStack <- raster::stack(anthroBuff, 
                              natDist, 
                              all, 
                              fire_excl_anthro)

    if(requireNamespace("fasterize", quietly = TRUE)){
      pp = fasterize::fasterize(inData@projectPoly,outStack[[1]])
    } else {
      message("To speed up install fasterize package")
      pp = raster::rasterize(inData@projectPoly,outStack[[1]])
    }
    
    #set NAs from landcover
    outStack[is.na(landCover) | (landCover == 0)|is.na(pp)] = NA
    
    names(outStack) = c("anthroBuff", 
                        "natDist", 
                        "totalDist",
                        "fire_excl_anthro")
    
    inData@processedData <- outStack
    
    #######
    #Range summaries
    message("calculating disturbance metrics")
        
    rr <- raster::zonal(outStack,pp,fun="mean",na.rm=T)
    rr <- as.data.frame(rr)
    polyInfo = as.data.frame(inData@projectPoly)
    polyInfo$geometry=NULL
    polyInfo$zone = seq(1:nrow(polyInfo))
    rr=merge(rr,polyInfo)
    
    inData@disturbanceMetrics <- rr

    return(inData)
  })
