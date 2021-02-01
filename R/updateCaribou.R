#' @include AAAClassDefinitions.R
NULL

#' updateCaribou
#'
#' @param CarHab CaribouHabitat object
#'
#' @param updatedLC,age,natDist,linFeat RasterLayer objects to be used to update CarHab.
#'   updatedLC is required, if the others are missing the layers from CarHab will be used.
#'
#'   

#' @export
setGeneric("updateCaribou", function(CarHab, newData, ...) standardGeneric("updateCaribou"))

# method to get final calculation from processed data in CaribouHabitat object
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "missing"), 
  function(CarHab, ...) {
    dots <- list(...)

    # process the data if not already done
    if(nrow(CarHab@processedData) < 2){
      
      if(is.null(dots$friLU)){
        stop("friLU is required to process data")
      }
      
      CarHab <- processData(CarHab, friLU = dots$friLU)
    }
    
    # calculate RSP
    coefTable <- coefTableHR %>% 
      filter(stringr::str_detect(Range, CarHab@attributes$caribouRange))
    
    CarHab@habitatUse <- calcRSP(CarHab@processedData, coefTable)
    
    projRas <- raster::rasterize(CarHab@projectPoly, CarHab@habitatUse[[1]], getCover=TRUE)
    projRas[projRas==0] <- NA
    
    CarHab@habitatUse <- raster::mask(CarHab@habitatUse, projRas )
    CarHab@habitatUse <- raster::crop(CarHab@habitatUse, CarHab@projectPoly, snap = "out")
    CarHab@processedData <- raster::mask(CarHab@processedData, projRas )
    CarHab@processedData <- raster::crop(CarHab@processedData, CarHab@projectPoly, snap = "out")
    
    return(CarHab)
  })

# method to update processed data when new data is supplied 
#' @rdname updateCaribou
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "list"), 
  function(CarHab, newData, friLU, resultsOnly = FALSE) {
    
    # process the data if not already done
    if(nrow(CarHab@processedData) < 2){
      stop("CarHab@processedData is empty. Run updateCaribou with no additional
           data to process the initial data before updating")
    }
    
    
    CarHab <- processData(CarHab, newData, friLU)
    
    CarHab <- updateCaribou(CarHab)
    
    if(resultsOnly){
      return(CarHab@habitatUse)
    }
    
    return(CarHab)
  
  })
