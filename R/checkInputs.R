

#' Check inputs to caribouHabitat
#'
#' Check the inputs that are needed in processData or updateCaribou so the
#' errors happen before data is processed.
#'
#' @noRd

.checkInputs <- function(caribouRange, winArea, landCover, coefTable) {
  expectedRanges <- unique(coefTable$Range)
  
  if(any(!caribouRange$Range %in% expectedRanges)){
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

}
  