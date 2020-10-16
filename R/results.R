#' @include AAAClassDefinitions.R
NULL

#' Extract results
#'
#' Extract results from caribouHabitat object.
#'
#' @param x A caribouHabitat object.
#'
#' @return A RasterStack with explanatory varibales and predictions for each
#' season.
#'
#' @export
setGeneric("results", function(x) standardGeneric("results"))

#' @rdname results
setMethod("results", signature(x = "CaribouHabitat"), function(x){
  if(nrow(x@habitatUse) < 2){
    stop("This object has empty @habitatUse. Calculate ",
         "habitatUse first with updateCaribou(x)")
  }
  
  result <- raster::stack(x@habitatUse, x@processedData)
  
  return(result)
})