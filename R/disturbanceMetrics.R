#' @include AAAClassDefinitions.R
NULL

#' @name DisturbanceMetrics
#' @rdname DisturbanceMetrics-class
setMethod(f = "initialize", signature = "DisturbanceMetrics",
          definition = function(.Object, landCover, natDist, 
                                anthroDist, linFeat, projectPoly, is.percent = FALSE, 
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

#'disturbanceMetrics
#'
#'Calculate the predictors described in Table 52 of Environment Canada (2011)
#'Scientific Assessment to Inform the Identification of Critical Habitat for
#'Woodland Caribou (Rangifer tarandus caribou), Boreal Population, in
#'Canada:2011 Update. Ottawa, Ontario.The variables calculated by this function
#'include: \itemize{ \item {Fire:} {\% fire} \item {Anthro:} {\% non-overlapping
#'anthropogenic disturbance.} \item {Total_dist:} {Percent total non-overlapping
#'fire and anthropogenic disturbance.} \item {fire_excl_anthro:} {\% fire not
#'overlapping with anthropogenic disturbance.} }
#'
#'Note assume natDist and anthroDist include 40 years of cumulative disturbance.
#'Note that locations where landCover is NA or 0 are omitted from the tabulated
#'area. Missing layers are omitted from the output, not interpreted as 0
#'disturbance. To update an existing DisturbanceMetrics object with new data see
#'\code{\link[caribouMetrics]{updateDisturbance}}.
#'
#'@param landCover filename or RasterLayer. 0 and NA values are assumed to be
#'  water and omitted from the tabulated area. Note landCover is also used to
#'  define the input grid, so must be provided even if all values are 1.
#'@param natDist filename or RasterLayer. Presence or absence of natural
#'  disturbance, primarily by fire. Should include 40 years cumulative
#'  disturbance. Optional.
#'@param anthroDist filename or RasterLayer. Anthropogenic disturbance including
#'  harvest. This can have an effect on any type of landcover except water.
#'  Should include 40 years cumulative disturbance. Optional.
#'@param linFeat filename, RasterLayer, sf object or a list of these that will
#'  be combined. Linear features.
#'@param projectPoly filename or sf object. Polygons defining range boundaries.
#'@param padProjPoly logical. Should the area around the \code{projectPoly} be
#'  used to avoid edge effects? If FALSE, the default, only data from inside the
#'  \code{projectPoly} is used. If TRUE then \code{projectPoly} is buffered and
#'  the other variables are clipped to the extent of the buffered area. Results
#'  are always clipped to the original \code{projectPoly}. It is ideal to set
#'  this to TRUE and provide a dataset that is larger than the
#'  \code{projectPoly} to avoid edge effects.
#'@param padFocal logical. This value is passed to the pad argument in
#'  \code{raster::focal}, if it is FALSE then cells near the edge will return
#'  NA, if it is TRUE a value will be returned for each cell that assumes cells
#'  outside the input data are 0 for all resource types. This is not a good
#'  assumption and should be used with caution.
#'@param bufferWidth number. Width of buffer applied to anthropogenic
#'  disturbance in metres. Default is 500.
#'@param linBuffMethod character. The method used to buffer linear features if
#'  they are supplied as sf lines. The default is "raster" in which case they
#'  are rasterized using the stars package and buffered using a moving window
#'  method. If "sf" then the lines are buffered with st_buffer and then
#'  rasterized. Either way points are included in the raster output.
#'
#'@return A DisturbanceMetrics Object see \code{\link{DisturbanceMetrics-class}}
#'
#'@seealso \code{\link{DisturbanceMetrics-class}} for information on the object
#'  returned and \code{\link{updateDisturbance}} for updating and existing
#'  DisturbanceMetrics object.
#'
#'@source Environment Canada. 2011. Scientific Assessment to Inform the
#'  Identification of Critical Habitat for Woodland Caribou (Rangifer tarandus
#'  caribou), Boreal Population, in Canada:2011 Update. Ottawa, Ontario.
#'  
#' @examples 
#' # create example rasters
#' lc <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10, crs = 5070)
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
#'
#'@export
setGeneric("disturbanceMetrics", 
           function(landCover, linFeat, projectPoly, is.percent = TRUE, ...) 
             standardGeneric("disturbanceMetrics"))

setMethod(
  "disturbanceMetrics", 
  signature(landCover = "ANY"), 
  function(landCover, linFeat, projectPoly, is.percent = TRUE, ...) {

    dots <- list(...)
    
    inputDataArgs <- dots[c("landCover","natDist","anthroDist","bufferWidth", "padProjPoly", "padFocal")]
    
    inputDataArgs <- inputDataArgs[which(lapply(inputDataArgs, length) > 0)]
    
    x <- do.call(inputDataDisturbance, c(lst(landCover,linFeat, projectPoly), 
                              inputDataArgs))
    
    linBuffMethod <- dots[["linBuffMethod"]]
    
    if(length(linBuffMethod) == 0){
      linBuffMethod <- "raster"
    }
    
    x <- processData(x, linBuffMethod = linBuffMethod)

    if(!is.null(dots$saveOutput)){
      
      byLayer <- grepl("\\.asc$|\\.sdat$|\\.rst$", dots$saveOutput)
      if(!byLayer && !grepl("\\.grd", dots$saveOutput)){
        warning("Saving output to ", dots$saveOutput, 
                ". Layernames will not be preserved.",
                " Use .grd format to preserve names")
      }
 
      raster::writeRaster(x@habitatUse, filename = dots$saveOutput, 
                          overwrite = TRUE, bylayer = byLayer, 
                          suffix = "names")
    }
    
    if (is.percent == TRUE) {
      x@disturbanceMetrics$Anthro <- x@disturbanceMetrics$Anthro * 100
      x@disturbanceMetrics$Fire <- x@disturbanceMetrics$Fire * 100
      x@disturbanceMetrics$Total_dist <- x@disturbanceMetrics$Total_dist * 100
      x@disturbanceMetrics$fire_excl_anthro <- x@disturbanceMetrics$fire_excl_anthro * 100
    }
    
    return(x)
  })
