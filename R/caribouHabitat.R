#' @include AAAClassDefinitions.R
NULL

#' @name CaribouHabitat
#' @rdname CaribouHabitat-class
setMethod(f = "initialize", signature = "CaribouHabitat",
          definition = function(.Object, plc, esker, fri, age, natDist, 
                                anthroDist, harv, linFeat, projectPoly,  
                                processedData, habitatUse){
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
            return(.Object)
          })

#' caribouHabitat
#' 
#' @param x CaribouHabitat object, or named list with required elements: plc,
#'   esker, fri, age, natDist, linFeat,  projectPoly, friLU, and optional
#'   elements: hexgrid = NULL, offset = c(0,0), eskerSave = NULL, hexgridSave =
#'   NULL, linFeatSave = NULL, winArea = NULL see \code{\link{inputData}} and
#'   \code{\link{processData}} for details on these arguments.
#'
#' @param caribouRange character the range where caribou were located. See
#'   \code{unique(coefTableHR$Range)}
#'

#' @export
setGeneric("caribouHabitat", function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, projectPoly, ...) standardGeneric("caribouHabitat"))

#' @rdname caribouHabitat
setMethod(
  "caribouHabitat", 
  signature(plc = "ANY"), 
  function(plc, esker, fri, age, natDist, anthroDist, harv, linFeat, projectPoly, ...) {
    
    dots <- list(...)
    
    inputDataArgs <- dots[c("winArea", "eskerSave", "linFeatSave", "padProjPoly", 
                         "caribouRange")]
    
    inputDataArgs <- inputDataArgs[which(sapply(inputDataArgs, length) > 0)]
    
    processDataArgs <- dots[c("friLU", "caribouRange", "winArea", "padProjPoly")]
    processDataArgs <- processDataArgs[which(sapply(processDataArgs, length) > 0)]
    
    updateArgs <- dots["caribouRange"]
    updateArgs <- updateArgs[which(sapply(updateArgs, length) > 0)]
    
    x <- do.call(inputData, c(lst(plc, esker, fri, age, natDist, anthroDist, harv,
                                  linFeat, projectPoly), inputDataArgs))
    
    if(is.null(processDataArgs$friLU)){
      warning("Argument friLU is required to process the data. The loaded", 
              " data is being returned. To calculate caribou habitat use pass ",
              "the returned object to updateCaribou with the required arguments",
              call. = FALSE)
      return(x)
    }
    x <- do.call(processData, c(x, processDataArgs))
    
    if(is.null(updateArgs$caribouRange)){
      warning("Argument caribouRange is required to calculate use. The processed", 
              " data is being returned. To calculate caribou habitat use pass ",
              "the returned object to updateCaribou with the required arguments")
      return(x)
    }
    
    return(updateCaribou(x, caribouRange = updateArgs$caribouRange))
  })

# TODO: add methods for print, plot, other useful generics
# TODO: deal with warning messages

# JH: add save method