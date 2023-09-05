#' Load Spatial Input Data
#'
#' Load spatial input data from file and then align the inputs to the
#' `projectPoly` and `refRast`. Inputs can be file names or spatial objects.
#'
#' @param projectPoly character or sf. A polygon delineating the study area.
#' @param refRast character or raster. A raster which will be used as the
#'   template that all other rasters must align to (but see `alttemplate`).
#' @param inputsList list. A named list of inputs that are either spatial
#'   objects or file names of spatial objects. ".shp" is the only extension
#'   accepted for vector data and all other extensions will be passed to
#'   [terra::rast()]. If an element is a list these are assumed to be linear features
#'   and they are combined.
#' @param convertToRast character. Optional. Names of elements of `inputsList`
#'   that should be converted to raster after loading.
#' @param convertToRastDens character. Optional. Names of elements of
#'   `inputsList` that should be converted to raster line density after loading.
#' @param altTemplate raster. Optional template raster for raster inputs that
#'   can have a different resolution from the `refRast`.
#' @param useTemplate character. Optional. Names of elements of `inputsList`
#'   that use `altTemplate`.
#' @param reclassOptions list. An optional named list containing a function, a
#'   list where the first element is a function that takes the corresponding
#'   `inputsList` element as its first argument and the subsequent elements are
#'   named arguments for that function, or a matrix that will be passed to
#'   [terra::classify()].
#' @param bufferWidth numeric. The width of a moving window that will be applied
#'   to the data. If it is supplied a buffer of 3*`bufferWidth` around the
#'   `projectPoly` is used so that rasters will be cropped to a larger area.
#'   This can be used to avoid edge effects in moving window calculations
#' @param ptDensity number. Only used if an element of `inputsList` is a list
#'   that contains a mixture of rasters and lines and is included in
#'   convertToRast. See [rasterizeLineDensity()].
#' @param rastOut character. The format that rasters should be output with.
#'   "raster" for RasterLayers and "terra" for SpatRasters. The default is "terra".
#'
#' @return A named list with aligned spatial data components
#'
#' @family habitat
#' @export
#'
#' @examples
#' # create example rasters
#' lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
#'                      resolution = 250, crs = "EPSG:5070")
#' lc[] <- 0
#' nd <- lc
#' nd[1:30, 1:30] <- 1
#' ad <- lc
#' ad[30:50, 3:50] <- 1
#' lc[] <- 1
#' lc[70:100, 70:100] <- 2
#'
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' esk <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 10000, 10000, 0),
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#'
#'
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#'
#' # prepare data
#' res <- loadSpatialInputs(projectPoly = projPol, refRast = lc,
#'                          inputsList = list(esker = esk, linFeat = lf, natDist = nd,
#'                                            anthroDist = ad),
#'                          convertToRast = c("esker", "linFeat"))
#'                                            
loadSpatialInputs <- function(projectPoly, refRast, inputsList, convertToRast = NULL,
                              convertToRastDens = NULL,
                              altTemplate = NULL, useTemplate = NULL,
                              reclassOptions = NULL, bufferWidth = NULL,
                              ptDensity = 1, rastOut = "terra"){
  
  allData <- list(inputsList, projectPoly = projectPoly, refRast = refRast) %>% 
    purrr::list_flatten()
  
  # check that names match across lists and arguments
  inNames <- names(allData)
  
  refNames <- c(useTemplate, convertToRast, convertToRastDens, names(reclassOptions))
  
  if(any(!refNames %in% inNames)){
    stop(paste0(refNames[which(!refNames %in% inNames)], collapse = ", "), 
         " does not match any names in inputsList, or refRast or projectPoly", 
         call. = FALSE)
  }
  
  # load the data
  #remove NULLs
  nullNames <- names(purrr::keep(allData, is.null))
  allData <- purrr::compact(allData)
  
  loaded <- loadFromFile(purrr::keep(allData, is.character))
  
  combined <- purrr::map(purrr::keep(allData, ~is.list(.x) & !is.data.frame(.x)),
                         combineLinFeat)
  
  spatObjs <- purrr::discard(allData, ~is.character(.x) | 
                               (is.list(.x) & !is.data.frame(.x)))
  
  allData <- purrr::list_flatten(list(loaded, combined, spatObjs))
  
  if(is(allData$refRast, "Raster")){
    allData$refRast <- terra::rast(allData$refRast)
  }
  
  if(!is(allData$refRast, "SpatRaster")){
    stop("refRast must be a SpatRaster or RasterLayer", call. = FALSE)
  }
  
  # convert and RasterLayer inputs to SpatRaster
  allData <- rapply(allData, f = terra::rast, classes = "RasterLayer", how = "replace")
  
  # Align and crop the data
  projPolyOut <- prepProjPoly(allData$projectPoly, allData$refRast, bufferWidth,
                              !is.null(bufferWidth))
  
  allData$projectPoly <- projPolyOut[["projectPoly"]]
  allData$projectPolyOrig <- projPolyOut[["projectPolyOrig"]]
  
  rm(projPolyOut)
  
  rasters <- purrr::keep(allData, ~is(.x, "SpatRaster"))
  rasters <- rasters[-which(names(rasters) == "refRast")]
  rasters <- prepRasts(rasters,
                       landCover = allData$refRast, 
                       projectPoly = allData$projectPoly,
                       tmplt = altTemplate,
                       useTmplt = useTemplate)
  
  notRasters <- purrr::discard(allData, ~is(.x, "SpatRaster"))
  notRasters <- notRasters[-which(names(notRasters) %in% 
                                    c("projectPoly", "projectPolyOrig"))] 
  notRasters <- purrr::map2(notRasters, names(notRasters), 
                  ~checkAlign(.x, allData$projectPoly, .y, "projectPoly"))
  
  allData <- purrr::list_flatten(list(rasters, notRasters, 
                           allData[which(names(allData) %in% c("projectPoly", "projectPolyOrig"))]))
  
  # Process the data
  if(length(nullNames) > 0){
    convertToRastDens <- dplyr::setdiff(convertToRastDens, nullNames)
    convertToRast <- dplyr::setdiff(convertToRast, nullNames)
  }
  
  if(!is.null(convertToRast)){
    # skip ones that are already rasters
    convertToRast <- dplyr::setdiff(convertToRast,
                                        names(purrr::keep(allData, is, "SpatRaster")))
    
    tmplts <- allData$refRast
    
    allData[convertToRast] <- purrr::map(
      allData[convertToRast],  
      ~ terra::rasterize(.x, allData$refRast, background = 0)
    )
  }
  
  if(!is.null(convertToRastDens)){
    # skip ones that are already rasters
    convertToRastDens <- dplyr::setdiff(convertToRastDens,
                                    names(purrr::keep(allData, is, "SpatRaster")))
    
    if(!is.null(useTemplate)){
      useTemplate <- dplyr::setdiff(useTemplate,
                                    names(purrr::keep(allData, is, "SpatRaster")))
      if(is.null(altTemplate)){
        altTemplate <- terra::rast(allData$refRast) %>% terra::`res<-`(c(400, 400))
      }
    } 
    
    tmplts <- purrr::map(convertToRastDens,
                         ~if(.x %in% useTemplate){altTemplate}else{allData$refRast})
    
    allData[convertToRastDens] <- purrr::map2(allData[convertToRastDens], tmplts, 
                                          rasterizeLineDensity, 
                                          ptDensity = ptDensity)
  }
  
  if(length(nullNames) > 0){
  reclassOptions <- reclassOptions[which(!names(reclassOptions) %in% nullNames)]
  }
  
  if(!is.null(reclassOptions)){
    reclassed <- allData[names(reclassOptions)]
    
    reclassed <- purrr::map2(
      reclassed, reclassOptions, 
      function(x, fn){
        if(is.matrix(fn)){
          return(terra::classify(x, fn))
        } else if(is.function(fn)){
          return(fn(x))
        } else if(is.list(fn)){
          args <- fn[-1]
          fn_args <- names(formals(fn[[1]]))
          args <- purrr::list_flatten(list(args, x))
          names(args) <- c(names(fn[-1]), fn_args[1])
          if("template" %in% fn_args){
            args <- purrr::list_flatten(list(args, template = allData$refRast))
          }
          fn <- fn[[1]]
          return(do.call(fn, args))
        }
      })
    
    allData[names(reclassOptions)] <- reclassed
  }
  
  if(rastOut == "raster"){
    allData <- rapply(allData, function(x){as(x, "Raster")},
                      classes = c("SpatRaster"), how = "replace")
  } else if(rastOut == "terra"){
    allData <- rapply(allData, terra::rast, classes = c("RasterLayer"), how = "replace")
  } else {
    stop("rastOut must be 'terra' or 'raster' not: ", rastOut)
  }


  return(allData)
}
