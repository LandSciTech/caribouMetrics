#' @include AAAClassDefinitions.R
NULL

#' Extract results
#'
#' Extract results from CaribouHabitat or DisturbanceMetrics object.
#'
#' @param x A CaribouHabitat or DisturbanceMetrics object.
#' @param type string. The name of the slot to return. If x is a CaribouHabitat
#'   object the default is "both" and the habitatUse and processedData will be
#'   returned as a RasterStack.
#'
#' @return By default a RasterStack if x is a CaribouHabitat object and a
#'   data.frame if x is a DisturbanceMetrics object. 
#'
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname results
setMethod("results", signature(x = "CaribouHabitat"), function(x, type = "both"){
  if(nrow(x@habitatUse) < 2){
    stop("This object has empty @habitatUse. Calculate ",
         "habitatUse first with updateCaribou(x)")
  }
  
  if(type == "both"){
    result <- raster::stack(x@habitatUse, x@processedData)
    
    return(result)
  }
  
  return(slot(x, type))

})

#' @rdname results
setMethod("results", signature(x = "DisturbanceMetrics"), function(x, type = "disturbanceMetrics"){
  slot(x, type)
})