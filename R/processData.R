
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
    if(inData@attributes$updateLC){
      # Update landCover based on updatedLC, age and disturbance history
      updted <- updateLC(landCover = inData@landCover, updatedLC = inData@updatedLC, age = inData@age, 
                         natDist = inData@natDist, anthroDist = inData@anthroDist,
                         harv = inData@harv)
      
      inData@natDist <- updted$natDist
    } else {
      updted <- list(landCover = inData@landCover)
    }

    
    # resample to 16ha and flag cells with >35% disturbance
    expVars <- applyDist(updted$landCover, inData@natDist, inData@anthroDist, 
                         inData@harv)
    
    rm(updted)
    
    # create raster attribute table
    inData@landCover <- raster::ratify(inData@landCover)
    
    levels(inData@landCover)[[1]] <- left_join(raster::levels(inData@landCover)[[1]], 
                                         resTypeCode, by = c(ID = "code"))

    # Resample linFeat and esker to 16 ha
    if(any(raster::res(inData@linFeat) < 400)){
      message("resampling linFeat to match landCover resolution")
      tmplt <- inData@landCover %>% raster::`res<-`(c(400, 400))
      
      inData@linFeat <- raster::resample(inData@linFeat, tmplt,
                                         method = "bilinear")
    }
    if(!compareRaster(expVars[[1]], inData@linFeat, stopiffalse = FALSE)){
      stop("linFeat could not be aligned with landCover.",
           " Please provide a higher resolution raster or a shapefile",
           call. = FALSE)
    }
    
    if(any(raster::res(inData@esker) < 400)){
      # eskFile <- tempfile()
      # raster::writeRaster(inData@esker, eskFile)
      
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
      inData@attributes$winArea <- coefTableHR %>% 
        filter(Range == inData@attributes$caribouRange) %>% 
        pull(WinArea) %>% 
        max()
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
       !all(names(newData) %in% c("updatedLC", "age", "natDist", "anthroDist", "harv", 
                                  "linFeat" ))){
      stop("newData must be a named list with any of the names ",
           "updatedLC ", "age ", "natDist ", "anthroDist ", "harv ", "linFeat ", 
           "incorrect names provided are:  ", 
           names(newData)[which(!names(newData) %in% 
                                  c("updatedLC", "age", "natDist", "anthroDist",
                                    "harv", "linFeat" ))], call. = FALSE)
    }
    
    if(!all(sapply(newData, is, "RasterLayer"))){
      stop("All data supplied in the newData list must be RasterLayer objects")
    }
    
    # check combinations of data in newData make sense
    if(!is.null(newData$updatedLC) && 
       !any(c("age", "natDist", "harv") %in% names(newData))){
      stop("to use updatedLC data either harv or natDist and age must be provided")
    }
    if(any(c("age", "natDist") %in% names(newData)) && 
       !all(c("age", "natDist") %in% names(newData))){
      stop("both natDist and age must be provided to use either")
    }
    
    # check alignment of new data with landCover
    newData <- purrr::map_at(newData, c("harv", "age", "natDist", "updatedLC"),
                              ~aggregateIf(.x, inData@landCover, names(.x), "landCover"))
    
    newData <- purrr::map2(newData, names(newData), 
                ~checkAlign(.x, inData@landCover, .y, "landCover"))
    
    #raster::compareRaster(newData, inData@landCover)
    
    # Get window area from table b/c some models used different sizes
    if(is.null(inData@attributes$winArea)){
      inData@attributes$winArea <- coefTableHR %>% 
        filter(Range == inData@attributes$caribouRange) %>% 
        pull(WinArea) %>% 
        max()
    }
    
   
    LCUpdated <- any(names(newData) %in%
                             c("updatedLC", "age", "natDist", "anthroDist", "harv"))
    
    if(LCUpdated){
      
      if(is.null(newData$updatedLC)){
        newData$updatedLC <- inData@updatedLC
      }
      if(is.null(newData$age)){
        newData$age <- inData@age
      }
      if(is.null(newData$natDist)){
        newData$natDist <- inData@natDist
      }
      if(is.null(newData$anthroDist)){
        newData$anthroDist <- inData@anthroDist
      }
      if(is.null(newData$harv)){
        newData$harv <- inData@harv
      }
      
      # Update landCover based on updatedLC, age and disturbance history
      updted <- updateLC(landCover = inData@landCover, updatedLC = newData$updatedLC, age = newData$age,
                          natDist = newData$natDist, 
                          anthroDist = newData$anthroDist, harv = newData$harv)
      
      inData@natDist <- updted$natDist
      
      expVars <- applyDist(updted$landCover, natDist = newData$natDist, 
                           anthroDist = newData$anthroDist,
                           harv = newData$harv)
      
      
      rm(updted) 
    } 
    
    # resample linFeat if provided
    if(!is.null(newData$linFeat)){
      newData$linFeat <- raster::resample(newData$linFeat, inData@linFeat, 
                                          method = "bilinear")
      inData@linFeat <- newData$linFeat
    } 

    # calculate moving window average for changed explanatory variables
    if("linFeat" %in% names(newData) && LCUpdated){
      expVars <- expVars %>% addLayer(inData@linFeat)
      
      layernames <- c(resTypeCode %>% arrange(code) %>%
                        pull(ResourceType) %>% as.character(), "TDENLF")
    } else if (LCUpdated){
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
        raster::addLayer(raster::setValues(expVars[[1]], 1)) %>% 
        `names<-`(c(names(expVars), "LGMD", "CONST"))
    }
    
    notUpdated <- which(!names(inData@processedData) %in% names(expVars)) 
    
    inData@processedData <- raster::stack(expVars, 
                                          inData@processedData[[notUpdated]])
    
    return(inData)
  })