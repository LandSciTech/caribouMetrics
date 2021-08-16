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
           natDist=NULL,anthroDist = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {
    
    charIn <-  sapply(list(landCover,natDist, anthroDist,  
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 
    
    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }
    
    
    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
    rastLst <- list(landCover, 
                    natDist, anthroDist)
    
    # remove NULLs from rastLst
    rastLst <- rastLst[which(!vapply(rastLst, function(x) is.null(x), 
                                     FUN.VALUE = TRUE))]
    
    if(!do.call(raster::compareRaster, c(rastLst, list(landCover, res = TRUE, extent = FALSE, 
                                                       rowcol = FALSE, stopiffalse = FALSE)))){
      stop("all raster data sets must have matching resolution", call. = FALSE)
    }
    rm(rastLst)
    
    
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
      
      projectPoly <- st_buffer(projectPoly, winRad*3)
    }
    
    landCover <- checkOverlap(landCover, projectPoly, "landCover", "projectPoly") %>%
      cropIf(projectPoly, "landCover", "projectPoly")
    
    # combine linFeat
    if(inherits(linFeat, "list")){
      
      linFeat <- combineLinFeat(linFeat)
    }
    if(!(is(linFeat, "sf") || is(linFeat, "sfc"))){
      if(is(linFeat, "Spatial")){
        linFeat <- sf::st_as_sf(linFeat)
      } 
      # roads <- rasterToLineSegments(roads)
      #linFeat <- raster::rasterToPoints(linFeat, fun = function(x){x > 0}, 
      #                                spatial = TRUE) %>% 
    }
    
    
    
    
    linFeat <- checkAlign(linFeat, landCover, "linFeat", "landCover")
    
    if(is(linFeat, "Raster")){
      
      tt = try(compareRaster(landCover, linFeat), silent = TRUE)
      if(class(tt)=="try-error"){
        stop("landcover and linFeat rasters do not have the",
             " same extent, number of rows and columns, projection, ",
             "resolution, or origin. Use raster::compareRaster() to ",
             "identify the problem.", call. = FALSE)
      }
    }
    
    # check alignment of other layers
    if(!is.null(natDist)){
      natDist <- checkOverlap(natDist, projectPoly, "natDist", "projectPoly") %>%
        cropIf(projectPoly, "natDist", "projectPoly")
      
      tt = try(compareRaster(landCover, natDist),silent=T)
      if(class(tt)=="try-error"){
        stop("landcover and natDist rasters do not have the same have the same extent, number of rows and columns, projection, resolution, or origin. Use raster::compareRaster() to identify the problem.")
      }
    } else {
      natDist <- raster(matrix(NA))
    }
    
    if(!is.null(anthroDist)){
      anthroDist <- checkOverlap(anthroDist, projectPoly, "anthroDist", "projectPoly") %>%
        cropIf(projectPoly, "anthroDist", "projectPoly")
      #anthroDist <- cropIf(anthroDist, landCover, "anthroDist", "landCover")
      
      tt = try(compareRaster(landCover, anthroDist),silent=T)
      if(class(tt)=="try-error"){
        stop("landcover and anthroDist rasters do not have the same have the same extent, number of rows and columns, projection, resolution, or origin. Use raster::compareRaster() to identify the problem.")
      }
    } else {
      anthroDist <- raster(matrix(NA))
    }
    
    return(new("DisturbanceMetrics", landCover, natDist, anthroDist, 
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
           natDist=NULL,anthroDist = NULL,
           bufferWidth = 500, padProjPoly = FALSE,
           padFocal = FALSE) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(natDist, anthroDist, projectPoly)
      
      linFeat <- combineLinFeat(linFeat)
      
    } else {
      indata <- lst(landCover,natDist, anthroDist, linFeat,
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
    
    neverVect <- c("natDist", "anthroDist")
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
                                linFeat = linFeat, 
                                projectPoly = indata$projectPoly, 
                                bufferWidth = bufferWidth, 
                                padProjPoly = padProjPoly, padFocal = padFocal))
    
  })