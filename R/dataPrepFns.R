#' Functions used for data prep for disturbanceMetrics and caribouHabitat
#' 
#' @noRd
prepProjPoly <- function(projectPoly, landCover, bufferWidth, padProjPoly){
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
    # winRad <- (sqrt(max(winArea)*10000/pi)/res(landCover)[1]) %>% 
    #   round(digits = 0)*res(landCover)[1]
    
    winRad <- (bufferWidth/terra::res(landCover)[1]) %>% 
      round(digits = 0)*terra::res(landCover)[1] %>% 
      round()
    
    projectPoly <- projectPoly %>% st_buffer(winRad*3)
  }
  return(lst(projectPoly, projectPolyOrig))
}

prepRasts <- function(rastLst, landCover, projectPoly, tmplt = NULL, 
                      useTmplt = NULL){
  
  landCover <- checkAlign(landCover, projectPoly, "landCover", "projectPoly") 
  
  if(is.null(tmplt)){
    tmplt <- terra::rast(landCover) %>% terra::`res<-`(c(400, 400))
  }
    
  # remove NULLs from rastLst
  rastLst <- rastLst[which(!vapply(rastLst, function(x) is.null(x), 
                                   FUN.VALUE = TRUE))]
  rastCRSMatch <- lapply(rastLst, terra::compareGeom,
                         y = landCover, crs = TRUE, res = FALSE, ext = FALSE, 
                         rowcol = FALSE, stopOnError = FALSE) %>%
    unlist() %>% all()
  
  if(!rastCRSMatch){
    stop("all raster data sets must have matching CRS", call. = FALSE)
  }
  
  # check alignment of each raster against projectPoly and compareRaster with
  # landCover
  rastLst <- purrr::map2(rastLst, names(rastLst),
              ~checkAlign(.x, projectPoly, .y, "projectPoly")) 
  
  # # need to crop tmplt too so that it will match extent
  # tmplt <- cropIf(tmplt, projectPoly, "tmplt", "projectPoly")
  
  # tmplt is usually res 400 400 raster that linFeat and esker are rasterized to
  tmpltUse <- rep_len(list(), length(rastLst)) %>% as.list()
  if(!is.null(useTmplt)){
    tmpltUseInd <- which(names(rastLst) %in% useTmplt)
    for (i in tmpltUseInd) {
      tmpltUse[[i]] <- tmplt
    }
    
  }

  purrr::pwalk(list(rastLst, tmpltUse, names(rastLst)),
               ~checkCompRast(x = ..1, y = landCover, nmx = ..3,
                              nmy = "landCover", y2 = ..2))
  
  rastLst$refRast <- landCover
  
  return(rastLst)
}

loadFromFile <- function(indata){
  # remove NULLs from indata
  indata <- indata[which(!vapply(indata, function(x) is.null(x), 
                                 FUN.VALUE = TRUE))]
  
  if(length(indata) == 0){
    return(list())
  }
  
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
  
  neverVect <- c("refRast", "natDist", "anthroDist")
  neverRast <- c("projectPoly")
  
  if(any(vect %in% neverVect)){
    stop("Extension is .shp but a raster file must be provided for: ",
         paste0(vect[which(vect %in% neverVect)], collapse = ", "))
  }
  
  if(any(rast %in% neverRast)){
    stop("Extension is not .shp but a shapefile must be provided for: ",
         paste0(rast[which(rast %in% neverRast)], collapse = ", "))
  }
  
  
  indata[vect] <- lapply(indata[vect], sf::st_read, quiet = TRUE, agr = "constant")
  indata[rast]<- lapply(indata[rast], terra::rast)
  
  return(indata)
}