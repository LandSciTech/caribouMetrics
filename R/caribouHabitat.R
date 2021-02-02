#' @include AAAClassDefinitions.R
NULL

#' @name CaribouHabitat
#' @rdname CaribouHabitat-class
setMethod(f = "initialize", signature = "CaribouHabitat",
          definition = function(.Object, landCover, esker, updatedLC, age, natDist, 
                                anthroDist, harv, linFeat, projectPoly,  
                                processedData, habitatUse, attributes){
            .Object@landCover <- landCover 
            .Object@esker <- esker
            .Object@updatedLC <- updatedLC
            .Object@age <- age
            .Object@natDist <- natDist
            .Object@anthroDist <- anthroDist
            .Object@harv <- harv
            .Object@linFeat <- linFeat
            .Object@projectPoly <- projectPoly
            .Object@processedData <- processedData
            .Object@habitatUse <- habitatUse 
            .Object@attributes <- attributes
            return(.Object)
          })

#'caribouHabitat
#'
#'Calculate the probability of caribou habitat use in spring, summer, fall and
#'winter for caribou ranges in Northern Ontario, based on Hornseth and Rempel,
#'2016.
#'
#'Caribou habitat use is calculated based on the availability of resources and
#'the presence of disturbances on the landscape. The primary source of resource
#'information is the \code{landCover} but this is can be updated based on more
#'recent \code{updatedLC} data and disturbance information. All data sources can
#'be provided either as filenames or as spatial files. The result is a
#'CaribouHabitat object which has methods defined for plotting and extracting
#'the results. To update an existing CaribouHabitat object with new data see
#'\link[caribouMetrics]{updateCaribou}.
#'
#'
#'@param landCover filename or RasterLayer. Provincial landcover class
#'@param esker filename, RasterLayer or sf object. Eskers. If it is a
#'  RasterLayer then it should be esker density in m^2/ha.
#'@param updatedLC filename or RasterLayer. Land cover data used to update the
#'  landCover raster in areas that were disturbed since the landCover data was
#'  created. If NULL, the default the landCover will not be updated
#'@param age filename or RasterLayer. Tree age in years. Used to inform whether
#'  a cell should be updated after disturbance
#'@param natDist filename or RasterLayer. Presence or absence of natural
#'  disturbance, primarily by fire.
#'@param anthroDist filename or RasterLayer. Anthropogenic disturbance other
#'  than harvest. This can have an effect on any type of landcover except water.
#'@param harv filename or RasterLayer. Harvest history. This can only have an
#'  effect on forest landcover types and will not affect wetlands or water.
#'@param linFeat filename, RasterLayer, sf object or named list with elements
#'  roads, rail, and utilities. Linear features. If it is a RasterLayer then it
#'  should be linear feature density in m^2/ha.
#'@param projectPoly filename or sf object. Polygon defining the project area.
#'@param caribouRange character. The range where caribou were located. See
#'  \code{unique(coefTableHR$Range)} for options.
#'@param eskerSave filename to save rasterized esker data.
#'@param linFeatSave filename to save rasterized linear feature data.
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
#'@param saveOutput character. The filename to save the rasterBrick of habitat
#'  use probabilities to. Note this will overwrite any existing files.
#'
#'@source Rempel, R. and M. Hornseth. 2018. Range-specific seasonal resource
#'  selection probability functions for 13 caribou ranges in Northern Ontario.
#'  Ontario Ministry of Natural Resources and Forestry, Science and Research
#'  Branch, Peterborough, ON. Science and Research Internal File Report IFR-01.
#'  42 p. + appends.
#'
#'  Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
#'  woodland caribou (Rangifer tarandus caribou) across a gradient of
#'  anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'  \url{https://doi.org/10.1139/cjz-2015-0101}
#'
#'@export
setGeneric("caribouHabitat", 
           function(landCover, esker, linFeat, projectPoly, caribouRange, ...) 
             standardGeneric("caribouHabitat"))

setMethod(
  "caribouHabitat", 
  signature(landCover = "ANY"), 
  function(landCover, esker, linFeat, projectPoly, caribouRange, ...) {

    dots <- list(...)
    
    inputDataArgs <- dots[c("updatedLC", "age", "natDist", "anthroDist", 
                            "harv","winArea", "eskerSave", "linFeatSave", 
                            "padProjPoly", "friLU", "padFocal")]
    
    inputDataArgs <- inputDataArgs[which(lapply(inputDataArgs, length) > 0)]
    
    x <- do.call(inputData, c(lst(landCover, esker, linFeat, projectPoly,
                                  caribouRange), 
                              inputDataArgs))
    
    x <- processData(x)
    
    x <- updateCaribou(x)
    
    if(!is.null(dots$saveOutput)){
      
      byLayer <- grepl("\\.asc$|\\.sdat$|\\.rst$", dots$saveOutput)
 
      raster::writeRaster(x@habitatUse, filename = dots$saveOutput, 
                          overwrite = TRUE, bylayer = byLayer, 
                          suffix = "names")
    }
    
    return(x)
  })
