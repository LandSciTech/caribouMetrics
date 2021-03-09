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
#' @param harv
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
           natDist=NULL,anthroDist = NULL, harv = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {

    charIn <-  sapply(list(landCover,natDist, anthroDist, harv, 
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 

    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }
    

    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
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
      
      projectPoly <- st_buffer(projectPoly, winRad*3)
    }
    
    landCover <- checkAlign(landCover, projectPoly, "landCover", "projectPoly")
    
    # combine linFeat
    if(inherits(linFeat, "list")){
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
    }
    if(!(is(linFeat, "sf") || is(linFeat, "sfc"))){
      if(is(linFeat, "Spatial")){
        linFeat <- sf::st_as_sf(linFeat)
      } else if(is(linFeat, "Raster")){
        # roads <- rasterToLineSegments(roads)
        linFeat <- raster::rasterToPoints(linFeat, fun = function(x){x > 0}, 
                                        spatial = TRUE) %>% 
          sf::st_as_sf()
      }
    }
    
    
    linFeat <- checkAlign(linFeat, landCover, "linFeat", "landCover")
    
    # check alignment of other layers
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
    

    return(new("DisturbanceMetrics", landCover, natDist, anthroDist, harv,
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
           natDist=NULL,anthroDist = NULL, harv = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(natDist, anthroDist, harv, 
                    projectPoly)
      
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
      
    } else {
      indata <- lst(landCover,natDist, anthroDist, harv, linFeat,
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
    
    neverVect <- c("natDist", "anthroDist", "harv")
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

    return(inputDataDisturbance(landCover=indata$landCover, 
                     natDist = indata$natDist, anthroDist = indata$anthroDist, 
                     harv = indata$harv, linFeat = linFeat, 
                     projectPoly = indata$projectPoly, 
                     bufferWidth = bufferWidth, 
                     padProjPoly = padProjPoly, padFocal = padFocal))
    
  })