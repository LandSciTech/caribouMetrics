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
#' @param natDist 
#' @param anthroDist
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
           coefTable,
           natDist = NULL, 
           anthroDist = NULL, 
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE,
           padFocal = FALSE, ptDensity = 1, 
           tmplt =  raster::`res<-`(raster(landCover), c(400, 400))) {

    charIn <-  sapply(list(landCover, esker, 
                           natDist, anthroDist,  
                           linFeat, projectPoly), 
                      function(x) "character" %in% class(x)) 

    if(any(charIn)){
      stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
    }

    if(raster::isLonLat(landCover)){
      stop("landCover must have a projected CRS", call. = FALSE)
    }
    
    
    #Note linear features can also be rasters here. And natDist/anthroDist can be polygons.
    #also - should we be requiring fully aligned rasters, not just same resolution?
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

    # require column called Range
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
    
    .checkInputs(caribouRange, winArea, landCover, coefTable)
    
    if(st_crs(projectPoly) != st_crs(landCover)){
      projectPoly <- st_transform(projectPoly, crs = st_crs(landCover))
      message("projectPoly being transformed to have crs matching landCover")
    }
    
    projectPolyOrig <- projectPoly
    
    # union together multiple range polygons for raster processing
    projectPoly <- projectPoly %>% summarise()
    
    # pad projPoly to 3 times the window radius, using the larger if multiple
    if(padProjPoly){

      # window radius is radius of circle with winArea rounded to even number of
      # raster cells based on resolution
      winRad <- (sqrt(max(winArea)*10000/pi)/res(landCover)[1]) %>% 
        round(digits = 0)*res(landCover)[1]
      
      projectPoly <- projectPoly %>% st_buffer(winRad*3)
    }
    
    
    landCover <- checkOverlap(landCover, projectPoly, "landCover", "projectPoly") %>%
      cropIf(projectPoly, "landCover", "projectPoly")
    

    # rasterize eskers
    if(inherits(esker, "sf")){
      esker <- checkAlign(esker, projectPoly, "esker", "projectPoly")
      
      #tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      esker <- rasterizeLineDensity(esker, tmplt)
    } else {
      esker <- checkOverlap(esker, projectPoly, "esker", "projectPoly") %>% 
        cropIf(projectPoly, "esker", "projectPoly")
      
      chk <- any(raster::compareRaster(esker, landCover, res = TRUE, 
                                       extent = FALSE, rowcol = FALSE,
                                       stopiffalse = FALSE), 
                 raster::compareRaster(esker, tmplt,
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
      
      linFeat <- combineLinFeat(linFeat)
      
    }

    if(inherits(linFeat, "sf")){

      linFeat <- checkAlign(linFeat, projectPoly, "linFeat", "projectPoly")
      
      #tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
      linFeat <- rasterizeLineDensity(linFeat, tmplt, ptDensity)
    } else {
      linFeat <- checkOverlap(linFeat, projectPoly, "linFeat", "projectPoly") %>% 
        cropIf(projectPoly, "linFeat", "projectPoly")
      
      chk <- any(raster::compareRaster(linFeat, landCover, res = TRUE, 
                                       extent = FALSE, rowcol = FALSE,
                                       stopiffalse = FALSE), 
                 raster::compareRaster(linFeat, tmplt,
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
    if(!is.null(natDist)){

      natDist <- checkOverlap(natDist, projectPoly, "natDist", "projectPoly") %>%
        cropIf(projectPoly, "natDist", "projectPoly")
      #natDist <- cropIf(natDist, landCover, "natDist", "landCover")
      #Note: orginally all cropIf calls used landcover, but that led to compareRaster errors.
      #Switching to projectPoly appears to solve this problem, but may introduce others. Will see.

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
    
    return(new("CaribouHabitat", landCover, esker, natDist, anthroDist,
               linFeat, projectPolyOrig,  
               processedData = raster(matrix(NA)), 
               habitatUse = raster(matrix(NA)),
               attributes = list(caribouRange = caribouRange, winArea = winArea,
                                 padProjPoly = padProjPoly, padFocal = padFocal, 
                                 tmplt = tmplt)))
})

#' @rdname inputData
setMethod(
  "inputData", signature(landCover = "character"), 
  function(landCover, esker, linFeat,  projectPoly, caribouRange, coefTable,
           natDist = NULL, 
           anthroDist = NULL,
           eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE, 
           padFocal = FALSE, 
           ptDensity = 1) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(landCover, esker, 
                    natDist, anthroDist,  
                    projectPoly)
      
      linFeat <- combineLinFeat(linFeat)
      
    } else {
      indata <- lst(landCover, esker,
                    natDist, anthroDist, linFeat,
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
    
    neverVect <- c("landCover", "natDist", "anthroDist")
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
   
    indata$landCover <- indata$landCover %>% reclassPLC()
  
    return(inputData(landCover = indata$landCover, esker = indata$esker, 
                     natDist = indata$natDist, anthroDist = indata$anthroDist, 
                     linFeat = linFeat, 
                     projectPoly = indata$projectPoly, 
                     caribouRange = caribouRange, eskerSave = eskerSave, 
                     linFeatSave = linFeatSave, winArea = winArea, 
                     padProjPoly = padProjPoly, padFocal = padFocal, 
                     ptDensity = ptDensity, coefTable = coefTable))
    
  })