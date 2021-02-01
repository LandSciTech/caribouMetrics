

#' Update provincial land cover
#'
#' Update provincial land cover based on FRI, age and disturbance history
#'
#' @param plc 
#' @param fri 
#' @param age 
#' @param natDist 
#' @param anthroDist 
#' @param harv 
#'
#' @return
#' @export
#'
#' @examples
#' 
setGeneric("updatePLC", function(x, newData, ...) standardGeneric("updatePLC"))

setMethod(
  "updatePLC", 
  signature(x = "missing", newData = "missing"), 
  function(plc, fri, age, natDist, anthroDist, harv){

    # get code for DTN (natural disturbance)
    DTNcode <- resTypeCode %>%
      filter(ResourceType == "DTN") %>%
      pull(code)

    # Find cells in PLC that meet condition for update: had natDist but have
    # age > 35 OR where FRI indicates Harvest (since PLC was created)
    toUpdate <- ((age > 35) & (natDist == 1)) | harv == 1 
    
    # only update if there is FRI data available
    toUpdate <- toUpdate - is.na(fri)
    
    # Update cells in toUpdate in PLC to value in FRI
    plc <- mask(plc, mask = toUpdate, maskvalue = 1, 
                updatevalue = NA)
    plc <- cover(plc, fri)
    
    # remove cells that were updated from natDist
    natDist <- mask(natDist, mask = toUpdate, maskvalue = 1, 
                    updatevalue = 0)
    
    # update PLC to DTN if natDist and age <= 35
    yfNatDist <- (age <= 35 | is.na(age)) & natDist == 1
    
    plc <- mask(plc, yfNatDist, maskvalue = 1, 
                updatevalue = DTNcode)
    
    return(lst(plc, natDist))
  }
)

