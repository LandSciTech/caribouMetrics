#' @include AAAClassDefinitions.R
NULL
#' @include spatialAlignFns.R
NULL

#' Load input data
#' 
#' Load input data for calculating caribou habitat use.
#' 
#' @param plc 
#'
#' @param esker 
#' @param fri 
#' @param age 
#' @param natDist 
#' @param anthroDist
#' @param harv
#' @param linFeat 
#' @param projectPoly 
#' @param ... 
#'
#' @export
setGeneric("inputData", function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, projectPoly, ...) standardGeneric("inputData"))

#' @rdname inputData
setMethod(
  "inputData", signature(plc = "RasterLayer"), 
  function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, projectPoly, 
           caribouRange, eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE, padFocal = FALSE, friLU = NULL) {

    charIn <-  lapply(list(plc, esker, fri, age, natDist, anthroDist, harv, 
                           linFeat, projectPoly),
                      function(x) "character" %in% class(x))

    if(any(charIn)){
      stop("All data must be supplied as sf and raster objects or character
                 paths not a mixture of each",
           call. = FALSE)
    }

    if(raster::isLonLat(plc)){
      stop("plc must have a projected CRS", call. = FALSE)
    }
    
    .checkInputs(fri, caribouRange, friLU, winArea)
    
    if(st_crs(projectPoly) != st_crs(plc)){
      projectPoly <- st_transform(projectPoly, crs = st_crs(plc))
      message("projectPoly being transformed to have crs matching plc")
    }
    
    projectPolyOrig <- projectPoly
    
    if(padProjPoly){
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
      winRad <- (sqrt(winArea*10000/pi)/res(plc)[1]) %>% 
        round(digits = 0)*res(plc)[1]
      
      projectPoly <- st_buffer(projectPoly, winRad*3)
    }
    
    plc <- checkAlign(plc, projectPoly, "plc", "projectPoly")

    # rasterize eskers
    esker <- checkAlign(esker, plc, "esker", "plc")
    if(inherits(esker, "sf")){
      tmplt <- raster(plc) %>% raster::`res<-`(c(400, 400))
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
    
    linFeat <- checkAlign(linFeat, plc, "linFeat", "plc")
    
    if(inherits(linFeat, "sf")){
      tmplt <- raster(plc) %>% raster::`res<-`(c(400, 400))
      linFeat <- rasterizeLineDensity(linFeat, tmplt)
    }
    if(!is.null(linFeatSave)){
      raster::writeRaster(linFeat, linFeatSave, overwrite = TRUE)
      linFeat <- raster(linFeatSave)
    }
    
    # # Resample linFeat and esker to 16 ha
    # tmplt <- plc %>% raster::`res<-`(c(400, 400))
    # if(any(res(linFeat) < 400)){
    #   linFeat <- raster::resample(inData@linFeat, tmplt,
    #                               method = "bilinear")
    # }
    # 
    # if(any(res(esker) < 400)){
    #   esker <- raster::resample(inData@esker, tmplt,
    #                             method = "bilinear")
    # }
    # 
    # check alignment of other layers
    fri <- aggregateIf(fri, plc, "fri", "plc") %>% checkAlign(plc, "fri", "plc")
    age <- aggregateIf(age, plc, "age", "plc") %>% checkAlign(plc, "age", "plc")
    natDist <- aggregateIf(natDist, plc, "natDist", "plc") %>% 
      checkAlign(plc, "natDist", "plc")
    anthroDist <- aggregateIf(anthroDist, plc, "anthroDist", "plc") %>% 
      checkAlign(plc, "anthroDist", "plc")
    harv <- aggregateIf(harv, plc, "harv", "plc") %>%
      checkAlign(plc, "harv", "plc")

    # check crs, alignment, will give error if crs, extent, or res are different
    # Error in compareRaster(plc, fri, age, natDist, linFeat, esker) :
    #   different extent
    # SE: seems clear enough to me for a user to understand
    compareRaster(plc, fri, age, natDist, anthroDist, harv)

    # # these don't need matching extent
    # compareRaster(plc, linFeat, esker, crs = TRUE, res = TRUE,
    #               extent = FALSE, rowcol = FALSE)
    
    return(new("CaribouHabitat", plc, esker, fri, age, natDist, anthroDist, harv,
               linFeat, projectPolyOrig,  
               processedData = raster(matrix(NA)), 
               habitatUse = raster(matrix(NA)),
               attributes = list(caribouRange = caribouRange, winArea = winArea,
                                 padProjPoly = padProjPoly, padFocal = padFocal)))
})

#' @rdname inputData
setMethod(
  "inputData", signature(plc = "character"), 
  function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat,  projectPoly,
           caribouRange, eskerSave = NULL, linFeatSave = NULL, 
           winArea = NULL, padProjPoly = FALSE, padFocal = FALSE, friLU = NULL) {
    
    if(inherits(linFeat, "list")){
      indata <- lst(plc, esker, fri, age, natDist, anthroDist, harv, 
                    projectPoly)
      
      linFeat <- combineLinFeat(linFeat$roads, linFeat$rail, linFeat$utilities)
      
    } else {
      indata <- lst(plc, esker, fri, age, natDist, anthroDist, harv, linFeat,
                    projectPoly)
    }
    
    charIn <-  lapply(indata, 
                      function(x) "character" %in% class(x))
    
    if(!all(charIn)){
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
    
    neverVect <- c("plc", "fri", "age", "natDist", "anthroDist", "harv")
    neverRast <- c("projectPoly")
    
    if(any(vect %in% neverVect)){
      stop("Extension is .shp but a raster file must be provided for: ",
           paste0(vect[which(vect %in% neverVect)], collapse = ", "))
    }
    
    if(any(rast %in% neverRast)){
      stop("Extension is not .shp but a shapefile must be provided for: ",
           paste0(rast[which(rast %in% neverRast)], collapse = ", "))
    }
    
    
    indata[vect] <- lapply(indata[vect], st_read, quiet = TRUE)
    indata[rast]<- lapply(indata[rast], raster)
   
    if(is.character(linFeat)){
      linFeat <- indata$linFeat
    }
    
    return(inputData(indata$plc, indata$esker, indata$fri, indata$age,
                     indata$natDist, indata$anthroDist, indata$harv, linFeat, 
                     indata$projectPoly, caribouRange = caribouRange, 
                     eskerSave, linFeatSave, winArea = winArea, 
                     padProjPoly = padProjPoly, padFocal = padFocal,
                     friLU = friLU))
    
  })