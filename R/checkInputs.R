

#' Check inputs to caribouHabitat
#'
#' Check the inputs that are needed in processData or updateCaribou so the
#' errors happen before data is processed.
#'
#' @param updatedLC updatedLC data raster
#' @param ... other variables passed to caribouHabitat
#'

.checkInputs <- function(caribouRange, winArea, landCover, updatedLC, coefTable) {
  expectedRanges <- paste0(unique(coefTable$Range), collapse = ", ")
  
  if(any(stringr::str_detect(expectedRanges, caribouRange$Range, negate = TRUE))){
    stop("Caribou Range must match one of: ", expectedRanges, call. = FALSE)
  }
  
  if(!is.null(winArea)){
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
  }
  
  if(!all(raster::unique(landCover) %in% resTypeCode$code)){
    stop("landCover raster has values outside the range of: ", 
         paste0(resTypeCode$code, collapse = ", "), 
         ". landCover must be reclassified to resource types",
         call. = FALSE)
  }
  
  if(!is.null(updatedLC)){
    if(!all(raster::unique(updatedLC) %in% c(resTypeCode$code, NA))){
      stop("updatedLC raster has values ", raster::unique(updatedLC),
           " outside the range of: ", 
           paste0(resTypeCode$code, collapse = ", "), 
           ". updatedLC must be reclassified to resource types", 
           call. = FALSE)
    }
  }

}
  