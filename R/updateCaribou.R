#' @include AAAClassDefinitions.R
NULL

#' Update an existing CaribouHabitat Object
#'
#' Update a CaribouHabitat object in order to avoid reprocessing parts of the
#' data that have not changed. New data is supplied as a named list and the
#' object is updated depending on the elements provided in the list.
#'
#' If \code{newData} contains any of updatedLC, natDist, age, or harv and
#' \code{updateType} is "disturbed" then the landCover is updated for disturbed
#' areas and if any of updatedLC, natDist, age, or harv are missing the data
#' stored in the CaribouHabitat object is reused. If \code{newData} contains
#' only linFeat then only linear features will be updated.
#'
#' If \code{updateType} is "entire" then \code{newData} must contain updatedLC
#' and the landCover in the CaribouHabitat object will be replaced and new
#' projections made for the whole landscape. If natDist, harv, or anthroDist are
#' not provided the the data stored in the CaribouHabitat object is reused.
#'
#' @param CarHab CaribouHabitat object
#' @param newData a named list of RasterLayer objects to be used to update
#'   CarHab. Potential names are: updatedLC, age, natDist, harv, anthroDist, and
#'   linFeat.
#' @param updateType character. The default is "disturbed" which means that only
#'   disturbed areas are updated to updatedLC based on natDist, age, and harv
#'   following the process used by Rempel. If \code{updateType} is "entire" then
#'   the current landCover is replaced with updatedLC.
#' @param resultsOnly logical. If FALSE the whole CaribouHabitat object is
#'   returned. If TRUE then only the habitatUse RasterStack is returned.
#'
#' @return If \code{resultsOnly} is FALSE an updated CaribouHabitat object. If
#'   \code{resultsOnly} is TRUE a RasterStack with a layer for each season.
#'
#' @export
setGeneric("updateCaribou", function(CarHab, newData, ...) standardGeneric("updateCaribou"))

# method to get final calculation from processed data in CaribouHabitat object
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "missing"), 
  function(CarHab, coefTable = coefTableHR, doScale = FALSE) {
    # process the data if not already done
    if(nrow(CarHab@processedData) < 2){
      
      CarHab <- processData(CarHab)
    }
    

    # which coefficients to use for which range
    rangeCoefLst <- CarHab@projectPoly %>% 
      left_join(CarHab@attributes$caribouRange, by = "Range") %>% 
      split(.$coefRange)
    
    applyCalcRSP <- function(dat, rangeCoef, doScale, coefTable){
      dat <- raster::mask(dat, rangeCoef)
      
      coefT <- coefTable %>% 
        filter(Range %in% rangeCoef$coefRange)
      
      habitatUse <- calcRSP(dat, coefT, doScale = doScale)
    }
    
    habUseLst <- lapply(rangeCoefLst, applyCalcRSP, dat = CarHab@processedData,
                        doScale = doScale, coefTable = coefTable)
    
    if(length(habUseLst)> 1){
      # do.call doesn't work with names
      names(habUseLst) <- NULL
      
      habUseLst$fun <- function(...){sum(..., na.rm = TRUE)}
      
      CarHab@habitatUse <- do.call(raster::overlay, habUseLst)
      
      names(CarHab@habitatUse) <- names(habUseLst[[1]])
      
      # 0 created by sum when all are NA but 0 is impossible as result of RSF
      # becasue it is logistic so safe to set 0 to NA
      CarHab@habitatUse[CarHab@habitatUse == 0] <- NA
    } else {
      CarHab@habitatUse <- habUseLst[[1]]
    }

    
    
    # Takes a long time not sure it is worth it
    # # This keeps cells that are only partially in the polygon 
    # projRas <- raster::rasterize(CarHab@projectPoly, CarHab@habitatUse[[1]],
    #                              getCover=TRUE)
    # projRas[projRas==0] <- NA
    # 
    # CarHab@habitatUse <- raster::mask(CarHab@habitatUse, projRas )
    # CarHab@habitatUse <- raster::crop(CarHab@habitatUse, CarHab@projectPoly,
    #                                   snap = "out")
    # 
    # CarHab@processedData <- raster::mask(CarHab@processedData, projRas )
    # CarHab@processedData <- raster::crop(CarHab@processedData, CarHab@projectPoly,
    #                                      snap = "out")
    
    return(CarHab)
  })

# method to update processed data when new data is supplied 
#' @rdname updateCaribou
setMethod(
  "updateCaribou", 
  signature(CarHab = "CaribouHabitat", newData = "list"), 
  function(CarHab, newData, updateType = "disturbed", resultsOnly = FALSE, 
           coefTable = coefTableHR, doScale = FALSE) {
    
    if(!updateType %in% c("disturbed", "entire")){
      stop("updateType is not recognized please use 'disturbed' or 'entire'",
           call. = FALSE)
    }
    
    if(nrow(CarHab@processedData) < 2){
      stop("CarHab@processedData is empty. Run updateCaribou with no additional
           data to process the initial data before updating")
    }
    
    
    CarHab <- processData(CarHab, newData, updateType)
    
    CarHab <- updateCaribou(CarHab, coefTable = coefTable, doScale = doScale)
    
    if(resultsOnly){
      return(CarHab@habitatUse)
    }
    
    return(CarHab)
  
  })
