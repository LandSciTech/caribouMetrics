#' @include AAAClassDefinitions.R
NULL

#' Extract results
#'
#' Extract results from caribouHabitat object.
#'
#' @param x A caribouHabitat object.
#' @param type string. Either "both" for the habitat use results and the
#'   processed data or "habitatUse" for just the habitat use results.
#'
#' @return A RasterStack with explanatory varibales and predictions for each
#'   season.
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
  
  if(type =="habitatUse"){
    return(x@habitatUse)
  }

})