#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating caribou habitat use.
#' 
#' @param landCover 
#'
#' @param esker 
#' @param updatedLC 
#' @param age 
#' @param natDist 
#' @param anthroDist
#' @param harv
#' @param linFeat 
#' @param projectPoly 
#' @param ... 
#'
#' @export
setGeneric("inputData", function(landCover, esker, linFeat, projectPoly, ...) standardGeneric("inputData"))

#' @rdname inputData
setMethod(
  "inputData", signature(landCover = "RasterLayer"), 
  function(landCover, esker, linFeat, projectPoly, caribouRange,
           caribouRangesCoefs = caribouRange,
           updatedLC = NULL, age = NULL, natDist = NULL, 
           anthroDist = NULL, harv = NULL,
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE,
           padFocal = FALSE, ptDensity = 1) {

    charIn <-  sapply(list(landCover, esker, updatedLC, age, natDist, anthroDist, harv, 
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 

    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }
    
    if(!is.null(updatedLC) && any(is.null(age), is.null(natDist), is.null(harv))){
      stop("harv, age and natDist must be supplied to use updatedLC", 
           call. = FALSE)
    }

    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
    # Get window area from table b/c some models used different sizes
    if(is.null(winArea)){
      winArea <- coefTableHR %>% filter(Range %in% caribouRange) %>% 
        pull(WinArea) %>% 
        unique()
    }
    
    # Error if caribouRanges have different winAreas
    if(length(winArea) > 1){
      stop("If multiple caribouRanges are supplied they ",
           "must all use the same winArea",
           call. = FALSE)
    }
    
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
    
    .checkInputs(caribouRange, winArea, landCover, updatedLC)
    
    if(st_crs(projectPoly) != st_crs(landCover)){
      projectPoly <- st_transform(projectPoly, crs = st_crs(landCover))
      message("projectPoly being transformed to have crs matching landCover")
    }
    
    projectPolyOrig <- projectPoly
    
    if(padProjPoly){

      # window radius is radius of circle with winArea rounded to even number of
      # raster cells based on resolution
      winRad <- (sqrt(winArea*10000/pi)/res(landCover)[1]) %>% 
        round(digits = 0)*res(landCover)[1]
      
      projectPoly <- projectPoly %>% summarise() %>% st_buffer(winRad*3)
    }
    
    landCover <- checkAlign(landCover, projectPoly, "landCover", "projectPoly")

    # rasterize eskers
    esker <- checkAlign(esker, landCover, "esker", "landCover")
    if(inherits(esker, "sf")){
      tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      esker <- rasterizeLineDensity(esker, tmplt)
    }
    if(!is.null(eskerSave)){
      raster::writeRaster(esker, eskerSave, overwrite = TRUE)
      esker <- raster(eskerSave)
    }

    # rasterize linFeat
    if(inherits(linFeat, "list")){
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
    }

    linFeat <- checkAlign(linFeat, landCover, "linFeat", "landCover")
    
    if(inherits(linFeat, "sf")){
      tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      linFeat <- rasterizeLineDensity(linFeat, tmplt, ptDensity)
    }
    if(!is.null(linFeatSave)){
      raster::writeRaster(linFeat, linFeatSave, overwrite = TRUE)
      linFeat <- raster(linFeatSave)
    }
    
    # check alignment of other layers
    if(!is.null(updatedLC)){
      updatedLC <- aggregateIf(updatedLC, landCover, "updatedLC", "landCover") %>%
        checkAlign(landCover, "updatedLC", "landCover")
      compareRaster(landCover, updatedLC)
    } else {
      updatedLC <- raster(matrix(NA))
    }
    
    if(!is.null(age)){
      age <- aggregateIf(age, landCover, "age", "landCover") %>% 
        checkAlign(landCover, "age", "landCover")
      compareRaster(landCover, age)
    } else {
      age <- raster(matrix(NA))
    }
      
    if(!is.null(natDist)){
      natDist <- aggregateIf(natDist, landCover, "natDist", "landCover") %>% 
        checkAlign(landCover, "natDist", "landCover")
      compareRaster(landCover, natDist)
    } else {
      natDist <- raster(matrix(NA))
    }
        
    if(!is.null(anthroDist)){
      anthroDist <- aggregateIf(anthroDist, landCover, "anthroDist", "landCover") %>% 
        checkAlign(landCover, "anthroDist", "landCover")
      compareRaster(landCover, anthroDist)
    } else {
      anthroDist <- raster(matrix(NA))
    }
    
    if(!is.null(harv)){
      harv <- aggregateIf(harv, landCover, "harv", "landCover") %>%
        checkAlign(landCover, "harv", "landCover")
      compareRaster(landCover, harv)
    } else {
      harv <- raster(matrix(NA))
    }


    return(new("CaribouHabitat", landCover, esker, updatedLC, age, 
               natDist, anthroDist, harv,
               linFeat, projectPolyOrig,  
               processedData = raster(matrix(NA)), 
               habitatUse = raster(matrix(NA)),
               attributes = list(caribouRange = caribouRange, winArea = winArea,
                                 padProjPoly = padProjPoly, padFocal = padFocal, 
                                 updateLC = length(raster::unique(updatedLC)) > 0)))
})

#' @rdname inputData
setMethod(
  "inputData", signature(landCover = "character"), 
  function(landCover, esker, linFeat,  projectPoly, caribouRange, 
           updatedLC = NULL, age = NULL, natDist = NULL, 
           anthroDist = NULL, harv = NULL, 
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE, 
           padFocal = FALSE, friLU = NULL, ptDensity = 1) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(landCover, esker, updatedLC, age, natDist, anthroDist, harv, 
                    projectPoly)
      
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
      
    } else {
      indata <- lst(landCover, esker, updatedLC, age, natDist, anthroDist, harv, linFeat,
                    projectPoly)
    }
    

    # remove NULLs from indata
    indata <- indata[which(!vapply(indata, function(x) is.null(x), 
                                   FUN.VALUE = TRUE))]

    
    charIn <-  indata %>% unlist(recursive = FALSE) %>% is.character()
    
    if(!charIn){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }  
    
    filesExist <- sapply(indata, file.exists)
    if(!all(filesExist)){
      stop("Path(s) for ",
           paste0(names(filesExist)[which(!filesExist)], collapse = ", "), 
           " do(es) not exist")
    }
    
    vect <- names(indata)[which(grepl(".shp$", indata))]
    rast <- names(indata)[which(!grepl(".shp$", indata))]
    
    neverVect <- c("landCover", "updatedLC", "age", "natDist", "anthroDist", "harv")
    neverRast <- c("projectPoly")
    
    if(any(vect %in% neverVect)){
      stop("Extension is .shp but a raster file must be provided for: ",
           paste0(vect[which(vect %in% neverVect)], collapse = ", "))
    }
    
    if(any(rast %in% neverRast)){
      stop("Extension is not .shp but a shapefile must be provided for: ",
           paste0(rast[which(rast %in% neverRast)], collapse = ", "))
    }
    
    
    indata[vect] <- lapply(indata[vect], st_read, quiet = TRUE, agr = "constant")
    indata[rast]<- lapply(indata[rast], raster)
   
    if(is.character(linFeat)){
      linFeat <- indata$linFeat
    }
    if(is.null(friLU)){
      stop("friLU is required when using the character method of caribouHabitat", 
           call. = FALSE)
    }
    
    if(!is.null(indata$updatedLC)){
      if(is.null(friLU)){
        stop("friLU is required when updatedLC is supplied")
      }
      indata$updatedLC <- indata$updatedLC %>% reclassFRI(friLU)
    }
    
    indata$landCover <- indata$landCover %>% reclassPLC()
  
    return(inputData(landCover = indata$landCover, esker = indata$esker, 
                     updatedLC = indata$updatedLC, age = indata$age, 
                     natDist = indata$natDist, anthroDist = indata$anthroDist, 
                     harv = indata$harv, linFeat = linFeat, 
                     projectPoly = indata$projectPoly, 
                     caribouRange = caribouRange, eskerSave = eskerSave, 
                     linFeatSave = linFeatSave, winArea = winArea, 
                     padProjPoly = padProjPoly, padFocal = padFocal, 
                     ptDensity = ptDensity))
    
  })