#' @include AAAClassDefinitions.R
NULL

setMethod(f = "initialize", signature = "DisturbanceMetrics",
          definition = function(.Object, landCover, natDist, 
                                anthroDist, linFeat, projectPoly, isPercent = FALSE, 
                                processedData, disturbanceMetrics, attributes){
            .Object@landCover <- landCover
            .Object@natDist <- natDist
            .Object@anthroDist <- anthroDist
            .Object@linFeat <- list(linFeat)
            .Object@projectPoly <- projectPoly
            .Object@processedData <- processedData
            .Object@disturbanceMetrics <- disturbanceMetrics 
            .Object@attributes <- attributes
            return(.Object)
          })

#'Calculate metrics of natural and anthropogenic disturbance
#'
#'Calculate the predictors described in Table 52 of Environment Canada (2011)
#'Scientific Assessment to Inform the Identification of Critical Habitat for
#'Woodland Caribou (Rangifer tarandus caribou), Boreal Population, in
#'Canada:2011 Update. Ottawa, Ontario.The variables calculated by this function
#'include:
#'* Fire: % fire
#'* Anthro: % non-overlapping anthropogenic disturbance.
#'* Total_dist: Percent total non-overlapping fire and anthropogenic disturbance.
#'* fire_excl_anthro: % fire not overlapping with anthropogenic disturbance.
#'
#'Note assume natDist and anthroDist include 40 years of cumulative disturbance.
#'Note that locations where landCover is NA or 0 are omitted from the tabulated
#'area. Missing layers are omitted from the output, not interpreted as 0
#'disturbance. To update an existing DisturbanceMetrics object with new data see
#'[updateDisturbance()].
#'
#'@param landCover filename, SpatRaster or RasterLayer. 0 and NA values are assumed to be
#'  water and omitted from the tabulated area. Note landCover is also used to
#'  define the input grid, so must be provided even if all values are 1.
#'@param linFeat filename, SpatRaster, RasterLayer, sf object or a list of these that will
#'  be combined. Linear features.
#'@param projectPoly filename or sf object. Polygons defining range boundaries.
#'@param isPercent logical. Should the results be returned as a percentage? 
#'@param ... optional arguments:
#'   * natDist: filename, SpatRaster or RasterLayer. Presence or absence of natural
#'   disturbance, primarily by fire. Should include 40 years cumulative
#'   disturbance.
#'   * anthroDist: filename, SpatRaster or RasterLayer. Anthropogenic disturbance including
#'   harvest. This can have an effect on any type of landcover except water.
#'   Should include 40 years cumulative disturbance.
#'   * padProjPoly: logical. Should the area around the `projectPoly` be
#'   used to avoid edge effects? If FALSE, the default, only data from inside the
#'   `projectPoly` is used. If TRUE then `projectPoly` is buffered and
#'   the other variables are clipped to the extent of the buffered area. Results
#'   are always clipped to the original `projectPoly`. It is ideal to set
#'   this to TRUE and provide a dataset that is larger than the
#'   `projectPoly` to avoid edge effects.
#'   * padFocal: logical. This value is passed to the pad argument in
#'   `terra::focal`, if it is FALSE then cells near the edge will return
#'   NA, if it is TRUE a value will be returned for each cell that assumes cells
#'   outside the input data are 0 for all resource types. This is not a good
#'   assumption and should be used with caution.
#'   * bufferWidth: number. Width of buffer applied to anthropogenic
#'   disturbance in metres. Default is 500.
#'   * linBuffMethod: character. The method used to buffer linear features if
#'   they are supplied as sf lines. The default is "raster" in which case they
#'   are rasterized and then buffered using a moving window
#'   method. If "sf" then the lines are buffered with st_buffer and then
#'   rasterized. Either way points are included in the raster output.
#'   * saveOutput: character. The filename to save the raster of binary
#'   disturbances with buffered anthropogenic disturbance. Note this will
#'   overwrite existing files with the same name. The .grd format is recommended
#'   because it will preserve layer names when the file is reloaded.
#'   * preppedData: list. A list containing pre-prepared input data sets. If
#'   not NULL then data checks will be skipped. Names must match argument names
#'   except that `landCover` should be called `refRast` and
#'   `projectPoly` should be called `projectPolyOrig`
#'   See [loadSpatialInputs()].
#' 
#' 
#'@return A DisturbanceMetrics Object see [DisturbanceMetrics-class()]
#'
#'@seealso [DisturbanceMetrics-class()] for information on the object
#'  returned and [updateDisturbance()] for updating an existing
#'  DisturbanceMetrics object.
#'
#'@source Environment Canada. 2011. Scientific Assessment to Inform the
#'  Identification of Critical Habitat for Woodland Caribou (Rangifer tarandus
#'  caribou), Boreal Population, in Canada:2011 Update. Ottawa, Ontario.
#'  
#' @examples 
#' # create example rasters
#' lc <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0,
#'                   ymax = 10, crs = "EPSG:5070")
#' nd <- lc
#' nd[1:3, 1:3] <- 1
#' ad <- lc
#' ad[3:5, 3:5] <- 1
#' lc[] <- 1
#' 
#' # create sf objects
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10, 10), 
#'                                                             ncol = 2, byrow = TRUE))),
#'                               crs = 5070))
#' projPol <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(ad)))
#' 
#' # calculate disturbance
#' disturbanceMetrics(landCover = lc,
#'                    linFeat = lf,
#'                    natDist = nd,
#'                    anthroDist = ad,
#'                    projectPoly = projPol,
#'                    padFocal = TRUE,
#'                    bufferWidth = 1)
#'
#' @family disturbance
#'@export
disturbanceMetrics <- function(landCover = NULL, linFeat = NULL, 
                               projectPoly = NULL, isPercent = TRUE, ...) {

    dots <- list(...)
    
    # check required args
    if(is.null(dots$preppedData)){
      missReqArgs <- purrr::map_lgl(lst(landCover, linFeat, projectPoly),
                                    is.null)
      if(any(missReqArgs)){
        stop("The required arguments ", paste0(names(missReqArgs)[which(missReqArgs)], collapse = ", "),
             " are missing with no default")
      }
    }
    
    # check all optional arguments are in expected names
    expDotArgs <- c("natDist", "anthroDist",  "bufferWidth", 
                    "padProjPoly", "padFocal", "linBuffMethod", 
                    "saveOutput", "preppedData")
    
    if(!all(names(dots) %in% expDotArgs)){
      stop("Argument ", names(dots)[which(!names(dots) %in% expDotArgs)], 
           " does not match an expected argument. See ?disturbanceMetrics for arguments")
    }
    
    inputDataArgs <- dots[c("landCover","natDist","anthroDist","bufferWidth", 
                            "padProjPoly", "padFocal", "preppedData")]
    
    inputDataArgs <- inputDataArgs[which(lapply(inputDataArgs, length) > 0)]
    
    x <- do.call(inputDataDisturbance, c(lst(landCover,linFeat, projectPoly), 
                              inputDataArgs))
    
    linBuffMethod <- dots[["linBuffMethod"]]
    
    if(length(linBuffMethod) == 0){
      linBuffMethod <- "raster"
    }
    
    x <- processData(x, linBuffMethod = linBuffMethod, isPercent = isPercent)

    if(!is.null(dots$saveOutput)){
      
      byLayer <- grepl("\\.asc$|\\.sdat$|\\.rst$", dots$saveOutput)
      if(!byLayer && !grepl("\\.grd", dots$saveOutput)){
        warning("Saving output to ", dots$saveOutput, 
                ". Layernames will not be preserved.",
                " Use .grd format to preserve names")
      }
 
      terra::writeRaster(x@processedData, filename = dots$saveOutput, 
                          overwrite = TRUE, bylayer = byLayer, 
                          suffix = "names")
    }
    
    return(x)
  }
