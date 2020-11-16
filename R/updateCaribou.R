#' @include AAAClassDefinitions.R
NULL

#' updateCaribou
#'
#' @param x CaribouHabitat object
#'
#' @param fri,age,natDist,linFeat RasterLayer objects to be used to update x.
#'   fri is required, if the others are missing the layers from x will be used.
#'
#'   

#' @export
setGeneric("updateCaribou", function(x, newData, ...) standardGeneric("updateCaribou"))

# method to get final calculation from processed data in CaribouHabitat object
setMethod(
  "updateCaribou", 
  signature(x = "CaribouHabitat", newData = "missing"), 
  function(x, ...) {
    dots <- list(...)

    # process the data if not already done
    if(nrow(x@processedData) < 2){
      
      if(is.null(dots$friLU)){
        stop("friLU is required to process data")
      }
      
      x <- processData(x, friLU = dots$friLU)
    }
    
    # calculate RSP
    coefTable <- coefTableHR %>% 
      filter(stringr::str_detect(Range, x@attributes$caribouRange))
    
    x@habitatUse <- calcRSP(x@processedData, coefTable)
    
    x@habitatUse <- raster::mask(x@habitatUse, x@projectPoly)
    return(x)
  })

# method to update processed data when new data is supplied 
#' @rdname updateCaribou
setMethod(
  "updateCaribou", 
  signature(x = "CaribouHabitat", newData = "list"), 
  function(x, newData, friLU, resultsOnly = FALSE) {
    
    # process the data if not already done
    if(nrow(x@processedData) < 2){
      stop("x@processedData is empty. Run updateCaribou with no additional
           data to process the initial data before updating")
    }
    
    
    x <- processData(x, newData, friLU)
    
    x <- updateCaribou(x)
    
    if(resultsOnly){
      return(x@habitatUse)
    }
    
    return(x)
  
  })
