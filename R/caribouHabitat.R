#' @include AAAClassDefinitions.R
NULL

#' @name CaribouHabitat
#' @rdname CaribouHabitat-class
setMethod(f = "initialize", signature = "CaribouHabitat",
          definition = function(.Object, plc, esker, fri, age, natDist, 
                                anthroDist, harv, linFeat, projectPoly,  
                                processedData, habitatUse, attributes){
            .Object@plc <- plc 
            .Object@esker <- esker
            .Object@fri <- fri
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
#'Caluculate the probability of caribou habitat use in spring, summer, fall and
#'winter for caribou ranges in Northern Ontario, based on Hornseth and Rempel,
#'2016.
#'
#'Caribou habitat use is calculated based on the availability of resources and
#'the presence of disturbances on the landscape. The primary source of resource
#'information is the \code{plc} but this is can be updated based on more recent
#'\code{fri} data and disturbance information. All data sources can be provided
#'either as filenames or as spatial files. The result is a CaribouHabitat object
#'which has methods defined for plotting and extracting the results. To update
#'an existing CaribouHabitat object with new data see
#'\link[caribouMetrics]{updateCaribou}.
#'
#'
#'@param plc filename or RasterLayer. Provincial landcover class
#'@param esker filename, RasterLayer or sf object. Eskers. If it is a
#'  RasterLayer then it should be esker density in m^2/ha.
#'@param fri filename or RasterLayer. Forest resource inventory class ids must
#'  correspond to the ids in the \code{friLU}.
#'@param age filename or RasterLayer. Tree age in years
#'@param natDist filename or RasterLayer. Presence or absence or natural
#'  disturbance, primarily by fire
#'@param anthroDist filename or RasterLayer. Anthropogenic disturbance
#'@param harv filename or RasterLayer. Harvest history
#'@param linFeat filename, RasterLayer, sf object or named list with elements
#'  roads, rail, and utilities. Linear features. If it is a RasterLayer then it
#'  should be linear feature density in m^2/ha
#'@param projectPoly filename or sf object. Polygon defining the project area
#'@param caribouRange character. The range where caribou were located. See
#'  \code{unique(coefTableHR$Range)}
#'@param eskerSave filename to save rasterized esker data
#'@param linFeatSave filename to save rasterized linear feature data
#'@param padProjPoly logical. Should the area around the \code{projectPoly} be
#'  used to avoid edge effects? If FALSE, the default, only data from inside the
#'  \code{projectPoly} is used. If TRUE then \code{projectPoly} is buffered and
#'  the other variables are clipped to the extent of the buffered area. Results
#'  are always clipped to the original \code{projectPoly}. It is ideal to set
#'  this to TRUE and provide a dataset that is larger than the
#'  \code{projectPoly} to avoid edge effects.
#'@param friLU data.frame. A look up table to convert local forest units to
#'  regional forest units. It should have two columns, the first must contain
#'  all the unique values in the supplied \code{fri} raster and the second must
#'  contain the names of regional forest units matching those provided in the
#'  table \code{rfuToResType}
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
setGeneric("caribouHabitat", function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, projectPoly, caribouRange, ...) standardGeneric("caribouHabitat"))

setMethod(
  "caribouHabitat", 
  signature(plc = "ANY"), 
  function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, 
           projectPoly, caribouRange, ...) {

    dots <- list(...)
    
    inputDataArgs <- dots[c("winArea", "eskerSave", "linFeatSave", "padProjPoly",
                            "friLU")]
    
    inputDataArgs <- inputDataArgs[which(sapply(inputDataArgs, length) > 0)]
    
    processDataArgs <- dots[c("friLU")]
    processDataArgs <- processDataArgs[which(sapply(processDataArgs, length) > 0)]
    
    x <- do.call(inputData, c(lst(plc, esker, fri, age, natDist, anthroDist, harv,
                                  linFeat, projectPoly, caribouRange), 
                              inputDataArgs))
    
    if(is.null(processDataArgs$friLU)){
      warning("Argument friLU is required to process the data. The loaded", 
              " data is being returned. To calculate caribou habitat use pass ",
              "the returned object to updateCaribou with the required arguments",
              call. = FALSE)
      return(x)
    }
    x <- do.call(processData, c(x, processDataArgs))
    
    x <- updateCaribou(x)
    
    if(!is.null(dots$saveOutput)){
      
      byLayer <- grepl("\\.asc$|\\.sdat$|\\.rst$", dots$saveOutput)
 
      raster::writeRaster(x@habitatUse, filename = dots$saveOutput, 
                          overwrite = TRUE, bylayer = byLayer, 
                          suffix = "names")
    }
    
    return(x)
  })
