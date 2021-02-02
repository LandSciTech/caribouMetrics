

#' Update land cover
#'
#' Update land cover based on updatedLC, age and disturbance history
#'
#' @param landCover 
#' @param updatedLC 
#' @param age 
#' @param natDist 
#' @param harv 
#'
#' @return
#' @export
#'
#' @examples
#' 
setGeneric("updateLC", function(x, newData, ...) standardGeneric("updateLC"))

setMethod(
  "updateLC", 
  signature(x = "missing", newData = "missing"), 
  function(landCover, updatedLC, age, natDist, harv){
    
    # get code for DTN (natural disturbance)
    DTNcode <- resTypeCode %>%
      filter(ResourceType == "DTN") %>%
      pull(code)
    
    # Find cells in landCover that meet condition for update: had natDist but have
    # age > 35 OR where updatedLC indicates Harvest (since landCover was created)
    toUpdate <- ((age > 35) & (natDist == 1)) | harv == 1 
    
    # only update if there is updatedLC data available
    toUpdate <- toUpdate - is.na(updatedLC)
    
    # Update cells in toUpdate in landCover to value in updatedLC
    landCover <- mask(landCover, mask = toUpdate, maskvalue = 1, 
                      updatevalue = NA)
    landCover <- cover(landCover, updatedLC)
    
    # remove cells that were updated from natDist
    natDist <- mask(natDist, mask = toUpdate, maskvalue = 1, 
                    updatevalue = 0)
    
    # update landCover to DTN if natDist and age <= 35
    yfNatDist <- (age <= 35 | is.na(age)) & natDist == 1
    
    landCover <- mask(landCover, yfNatDist, maskvalue = 1, 
                      updatevalue = DTNcode)
    
    return(lst(landCover, natDist))
  }
)

