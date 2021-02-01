

#' Check inputs to caribouHabitat
#'
#' Check the inputs that are needed in processData or updateCaribou so the
#' errors happen before data is processed.
#'
#' @param fri fri data raster
#' @param ... other variables passed to caribouHabitat
#'
#'
#' @examples
.checkInputs <- function(caribouRange, winArea, plc, fri) {
  expectedRanges <- paste0(unique(coefTableHR$Range), collapse = ", ")
  
  if(stringr::str_detect(expectedRanges, caribouRange, negate = TRUE)){
    stop("caribouRange must match one of: ", expectedRanges, call. = FALSE)
  }
  
  if(!is.null(winArea)){
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
  }
  
  if(!all(raster::cellStats(plc, unique) %in% resTypeCode$code)){
    stop("plc raster has values outside the range of: ", 
         paste0(resTypeCode$code, collapse = ", "), 
         ". plc must be reclassified to resource types")
  }
  
  if(!all(raster::cellStats(fri, unique) %in% resTypeCode$code)){
    stop("plc raster has values outside the range of: ", 
         paste0(resTypeCode$code, collapse = ", "), 
         ". plc must be reclassified to resource types")
  }
  
}
  