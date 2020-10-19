
#' @include AAAClassDefinitions.R
NULL


#' @export
setGeneric("processData", function(inData, newData, ...) standardGeneric("processData"))

# Method for initial processing when no new data is provided
#' @rdname processData
#' @importFrom raster plot
setMethod(
  "processData", signature(inData = "CaribouHabitat", newData = "missing"), 
  function(inData, friLU, caribouRange, winArea = NULL, padProjPoly = FALSE) {
    # TODO:
    # get fri lu from raster legend
    
    expectedRanges <- paste0(unique(coefTableHR$Range), collapse = ", ")
    
    if(stringr::str_detect(expectedRanges, caribouRange, negate = TRUE)){
      stop("caribouRange must match one of: ", expectedRanges)
    }
    
    # checks types and match names
    if(!inherits(friLU, "data.frame")){
      stop("friLU must be a data.frame", call. = FALSE)
    }
    if(!is.numeric(friLU[,1])){
      stop("The first column of friLU must be numeric", 
           call. = FALSE)
    }
    
    if(any(is.na(friLU[,2]))){
      warning("friLU contains NA in the second column. ", 
              "Cells with these IDs will be replaced with values from plc ",
              "in the calculation of probability of use", call. = FALSE)
    }
    
    if(!all(unique(friLU[,2]) %in% c(unique(rfuToResType$RegionalForestUnit), NA))){
      stop("The second column of friLU must match a regional forest unit: ", 
           paste0(c(unique(rfuToResType$RegionalForestUnit), NA), sep = ", "),
           call. = FALSE)
    }
    
    if(!all(raster::unique(inData@fri) %in% friLU[,1])){
      stop("All unique values in fri must be present in friLU", call. = FALSE)
    }
    
    friLU <- friLU %>% set_names("ID", "RFU")
    
    # convert from FRI and PLC to ResType using supplied fu to rfu and internal
    # rfu to restype and plc to restype
    friLU <- left_join(friLU, rfuToResType, 
                       by = c("RFU" = "RegionalForestUnit")) %>% 
      select(-RFU)
    
    # reclassify plc and fri to resource types based on look up tables
    rclPLC <- plcToResType %>% 
      mutate(ResourceTypeCode = as.factor(ResourceType) %>% as.numeric()) %>% 
      select(PLCCode, ResourceType, ResourceTypeCode)
    
    resTypeLevels <- rclPLC %>% 
      select(ResourceType, ResourceTypeCode) %>% 
      distinct() %>% 
      set_names(c("resType", "code")) %>% arrange(code)
    
    if(!"DTN" %in% resTypeLevels$resType){
      resTypeLevels <- resTypeLevels %>% add_row(resType = "DTN", code = 99)
    }
    
    rclPLC <- rclPLC %>% select(-ResourceType) %>% 
      as.matrix(rclPLC, rownames.force = FALSE)
    
    rclFRI <-  friLU %>% select(ID, ResourceType) %>%
      left_join(resTypeLevels, by = c("ResourceType" = "resType")) %>%
      select(ID, code) %>% 
      as.matrix()
    
    inData@plc <- reclassify(inData@plc, rclPLC)
    
    inData@fri <- reclassify(inData@fri, rcl = rclFRI)
    

    
    # Update PLC based on FRI, age and disturbance history
    updted <- updatePLC(plc = inData@plc, fri = inData@fri, age = inData@age, 
                        natDist = inData@natDist, anthroDist = inData@anthroDist,
                        harv = inData@harv, resTypeLevels = resTypeLevels)
    
    expVars <- updted$plc
    inData@natDist <- updted$natDist
    rm(updted)
    
    # create raster attribute table
    inData@plc <- raster::ratify(inData@plc)
    
    levels(inData@plc)[[1]] <- left_join(raster::levels(inData@plc)[[1]], 
                                         resTypeLevels, by = c(ID = "code"))

    # Resample linFeat and esker to 16 ha
    tmplt <- inData@plc %>% raster::`res<-`(c(400, 400))
    inData@linFeat <- raster::resample(inData@linFeat, tmplt,
                                       method = "bilinear")

    inData@esker <- raster::resample(inData@esker, tmplt,
                                       method = "bilinear")
    
    # Get window area from table b/c some models used different sizes
    if(is.null(winArea)){
      winArea <- coefTableHR %>% filter(Range == caribouRange) %>% 
        pull(WinArea) %>% 
        max()
    }
    
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1]
    
    # calculate moving window average for all explanatory variables
    expVars <- expVars %>% 
      addLayer(inData@linFeat) %>% 
      addLayer(inData@esker) 
    
    layernames <- c(resTypeLevels %>% arrange(code) %>%
                      pull(resType) %>% as.character(),
                    "TDENLF", "ESK")
    
    if(!padProjPoly){
      expVars <- raster::mask(expVars, inData@projectPoly)
    }
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,
                                 nms = layernames)
    
    inData@processedData <- expVars %>% 
      raster::addLayer(expVars[["MIX"]] + expVars[["DEC"]]) %>% 
      raster::addLayer(raster::setValues(expVars[[1]], 1)) %>% 
      `names<-`(c(names(expVars), "LGMD", "CONST"))

    return(inData)
    
  })

# method to update processed data when new data is supplied 
#' @rdname processData
setMethod(
  "processData", 
  signature(inData = "CaribouHabitat", newData = "list"), 
  function(inData, newData, friLU, caribouRange, winArea = NULL) {
    
    if(is.null(names(newData)) ||
       !all(names(newData) %in% c("fri", "age", "natDist", "anthroDist", "harv", 
                                  "linFeat" ))){
      stop("newData must be a named list with any of the names ",
           "fri ", "age ", "natDist ", "anthroDist ", "harv ", "linFeat ", 
           "incorrect names provided are:  ", 
           names(newData)[which(!names(newData) %in% 
                                  c("fri", "age", "natDist", "anthroDist",
                                    "harv", "linFeat" ))], call. = FALSE)
    }
    
    if(!all(sapply(newData, is, "RasterLayer"))){
      stop("All data supplied in the newData list must be RasterLayer objects")
    }
    
    # check combinations of data in newData make sense
    if(!is.null(newData$fri) && 
       !any(c("age", "natDist", "harv") %in% names(newData))){
      stop("to use fri data either harv or natDist and age must be provided")
    }
    if(any(c("age", "natDist") %in% names(newData)) && 
       !all(c("age", "natDist") %in% names(newData))){
      stop("both natDist and age must be provided to use either")
    }
    
    # check alignment of new data with plc
    newData <- purrr::map2(newData, names(newData), 
                ~checkAlign(.x, inData@plc, .y, "plc"))
    
    #raster::compareRaster(newData, inData@plc)
    
    # get resTypeLevels from inData@plc RAT
    resTypeLevels <- inData@plc %>% raster::levels() %>% .[[1]] %>% 
      set_names(c("code", "resType")) %>% 
      mutate(resType = as.character(resType))
    
    expectedRanges <- paste0(unique(coefTableHR$Range), collapse = ", ")
    
    if(stringr::str_detect(expectedRanges, caribouRange, negate = TRUE)){
      stop("caribouRange must match one of: ", expectedRanges)
    }
    
    # convert from FRI to ResType using supplied fu to rfu and internal
    # rfu to restype 
    # checks types and match names
    if(!is.numeric(friLU[,1])){
      stop("The first column of friLU must be numeric", call. = FALSE)
    }
    
    if(any(is.na(friLU[,2]))){
      warning("friLU contains NA in the second column. ", 
              "Cells with these IDs will be replaced with values from plc ",
              "in the calculation of probability of use", call. = FALSE)
    }
    
    if(!all(unique(friLU[,2]) %in% c(unique(rfuToResType$RegionalForestUnit), NA))){
      stop("The second column of friLU must match a regional forest unit: ", 
           paste0(c(unique(rfuToResType$RegionalForestUnit), NA), sep = ", "), 
           call. = FALSE)
    }
    
    friLU <- friLU %>% set_names("ID", "RFU")
    
    # Get window area from table b/c some models used different sizes
    if(is.null(winArea)){
      winArea <- coefTableHR %>% filter(Range == caribouRange) %>% 
        pull(WinArea) %>% 
        max()
    }
    
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
    
    friLU <- left_join(friLU, rfuToResType, 
                       by = c("RFU" = "RegionalForestUnit")) %>% 
      select(-RFU)
    
    rclFRI <-  friLU %>% select(ID, ResourceType) %>%
      left_join(resTypeLevels, by = c("ResourceType" = "resType")) %>%
      select(ID, code) %>% 
      as.matrix()
    
    if(!is.null(newData$fri)){
      newData$fri <- reclassify(newData$fri, rcl = rclFRI)
    }
    
    PLCUpdated <- any(names(newData) %in%
                             c("fri", "age", "natDist", "anthroDist", "harv"))
    
    if(PLCUpdated){
      
      if(is.null(newData$fri)){
        newData$fri <- inData@fri
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
      
      # Update PLC based on FRI, age and disturbance history
      updted <- updatePLC(plc = inData@plc, fri = newData$fri, age = newData$age,
                          natDist = newData$natDist, 
                          anthroDist = newData$anthroDist, harv = newData$harv,
                          resTypeLevels = resTypeLevels)
      
      expVars <- updted$plc
      inData@natDist <- updted$natDist
      
      rm(updted) 
    } 
    
    # resample linFeat if provided
    if(!is.null(newData$linFeat)){
      newData$linFeat <- raster::resample(newData$linFeat, inData@linFeat, 
                                          method = "bilinear")
      inData@linFeat <- newData$linFeat
    } 

    # calculate moving window average for changed explanatory variables
    if("linFeat" %in% names(newData) && PLCUpdated){
      expVars <- expVars %>% addLayer(inData@linFeat)
      
      layernames <- c(raster::levels(inData@plc)[[1]] %>% arrange(ID) %>%
                        pull(resType) %>% as.character(), "TDENLF")
    } else if (PLCUpdated){
      layernames <- raster::levels(inData@plc)[[1]] %>% arrange(ID) %>%
                        pull(resType) %>% as.character()
    } else {
      expVars <- inData@linFeat
      
      layernames <- "TDENLF"
    }
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    winRad <- (sqrt(winArea*10000/pi)/res(expVars[[1]])[1]) %>% 
      round(digits = 0)*res(expVars[[1]])[1]
    
    expVars <- movingWindowAvg(rast = expVars, radius = winRad,  
                               nms = layernames)
    
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