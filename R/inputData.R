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
    
    rastLst <- list(landCover, updatedLC, age, natDist, 
                    anthroDist, harv)
    
    # remove NULLs from rastLst
    rastLst <- rastLst[which(!vapply(rastLst, function(x) is.null(x), 
                                   FUN.VALUE = TRUE))]
    
    if(!do.call(raster::compareRaster, c(rastLst, list(res = TRUE, extent = FALSE, 
                             rowcol = FALSE, stopiffalse = FALSE)))){
      stop("all raster data sets must have matching resolution", call. = FALSE)
    }
    rm(rastLst)
    
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
    if(nrow(projectPoly) == 1){
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
    
    # Get window area from table b/c some models used different sizes
    if(is.null(winArea)){
      winArea <- coefTableHR %>% filter(Range %in% caribouRange$coefRange) %>% 
        pull(WinArea) %>% 
        unique()
    }
    
    # Error if caribouRanges have different winAreas
    if(length(winArea) > 1){
      stop("If multiple caribouRanges are supplied they ",
           "must all use the same winArea",
           call. = FALSE)
    }
    .checkInputs(caribouRange, winArea, landCover, updatedLC)
    
    if(st_crs(projectPoly) != st_crs(landCover)){
      projectPoly <- st_transform(projectPoly, crs = st_crs(landCover))
      message("projectPoly being transformed to have crs matching landCover")
    }
    
    projectPolyOrig <- projectPoly
    
    # union together multiple range polygons for raster processing
    projectPoly <- projectPoly %>% summarise()
    
    if(padProjPoly){

      # window radius is radius of circle with winArea rounded to even number of
      # raster cells based on resolution
      winRad <- (sqrt(winArea*10000/pi)/res(landCover)[1]) %>% 
        round(digits = 0)*res(landCover)[1]
      
      projectPoly <- projectPoly %>% st_buffer(winRad*3)
    }
    
    landCover <- checkOverlap(landCover, projectPoly, "landCover", "projectPoly") %>%
      cropIf(projectPoly, "landCover", "projectPoly")

    # rasterize eskers
    esker <- checkAlign(esker, landCover, "esker", "landCover")
    if(inherits(esker, "sf")){
      tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      esker <- rasterizeLineDensity(esker, tmplt)
    } else {
      esker <- checkOverlap(esker, landCover, "esker", "landCover") %>% 
        cropIf(landCover, "esker", "landCover")
      
      chk <- any(raster::compareRaster(esker, landCover, res = TRUE, 
                                       extent = FALSE, rowcol = FALSE,
                                       stopiffalse = FALSE), 
                 raster::compareRaster(esker,
                                       raster(landCover) %>%
                                         raster::`res<-`(c(400, 400)),
                                       res = TRUE, extent = FALSE, 
                                       rowcol = FALSE,
                                       stopiffalse = FALSE))
      if(!chk){
        stop("esker is not aligned with landCover")
      }
    }
    
    if(!is.null(eskerSave)){
      raster::writeRaster(esker, eskerSave, overwrite = TRUE)
      esker <- raster(eskerSave)
    }

    # rasterize linFeat
    if(inherits(linFeat, "list")){
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
    }

    if(inherits(linFeat, "sf")){
      linFeat <- checkAlign(linFeat, landCover, "linFeat", "landCover")
      
      tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      linFeat <- rasterizeLineDensity(linFeat, tmplt, ptDensity)
    } else {
      linFeat <- checkOverlap(linFeat, landCover, "linFeat", "landCover") %>% 
        cropIf(landCover, "linFeat", "landCover")
      
      chk <- any(raster::compareRaster(linFeat, landCover, res = TRUE, 
                                       extent = FALSE, rowcol = FALSE,
                                       stopiffalse = FALSE), 
                 raster::compareRaster(linFeat,
                                       raster(landCover) %>%
                                         raster::`res<-`(c(400, 400)),
                                       res = TRUE, extent = FALSE, 
                                       rowcol = FALSE,
                                       stopiffalse = FALSE))
      if(!chk){
        stop("linFeat is not aligned with landCover")
      }
    }
    if(!is.null(linFeatSave)){
      raster::writeRaster(linFeat, linFeatSave, overwrite = TRUE)
      linFeat <- raster(linFeatSave)
    }
    
    # check alignment of other layers
    if(!is.null(updatedLC)){
      updatedLC <- cropIf(updatedLC, landCover, "updatedLC", "landCover") 
      
      compareRaster(landCover, updatedLC)
    } else {
      updatedLC <- raster(matrix(NA))
    }
    
    if(!is.null(age)){
      age <- cropIf(age, landCover, "age", "landCover")

      compareRaster(landCover, age)
    } else {
      age <- raster(matrix(NA))
    }
      
    if(!is.null(natDist)){
      natDist <- cropIf(natDist, landCover, "natDist", "landCover")

      compareRaster(landCover, natDist)
    } else {
      natDist <- raster(matrix(NA))
    }
        
    if(!is.null(anthroDist)){
      anthroDist <- cropIf(anthroDist, landCover, "anthroDist", "landCover")

      compareRaster(landCover, anthroDist)
    } else {
      anthroDist <- raster(matrix(NA))
    }
    
    if(!is.null(harv)){
      harv <- cropIf(harv, landCover, "harv", "landCover")
      
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
      
      expectNames = c("roads","rail","utilities")
      if(is.null(names(linFeat))){
        missingNames = expectNames  
      }else{
        missingNames = setdiff(names(linFeat),expectNames)
        
      }
      if(length(missingNames)>0){
        stop(paste0("Expecting names of linFeat list to include ",paste(expectNames,collapse=",")))
      }    
      
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