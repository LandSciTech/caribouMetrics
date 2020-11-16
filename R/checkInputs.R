

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
.checkInputs <- function(fri, caribouRange, friLU, winArea) {
  expectedRanges <- paste0(unique(coefTableHR$Range), collapse = ", ")
  
  if(stringr::str_detect(expectedRanges, caribouRange, negate = TRUE)){
    stop("caribouRange must match one of: ", expectedRanges, call. = FALSE)
  }
  
  if(!is.null(winArea)){
    if(!is.numeric(winArea)){
      stop("winArea must be a number (in hectares)", call. = FALSE)
    }
  }
  
  if(!is.null(friLU)){
    # checks types and match names
    if(!inherits(friLU, "data.frame")){
      stop("friLU must be a data.frame", call. = FALSE)
    }
    if(!is.numeric(friLU[,1])){
      stop("The first column of friLU must be numeric", 
           call. = FALSE)
    }
    
    if(any(is.na(friLU[,2]))){
      warning("friLU contains NA in the second column. ", 
              "Cells with these IDs will be replaced with values from plc ",
              "in the calculation of probability of use", call. = FALSE)
    }
    
    if(!all(unique(friLU[,2]) %in% c(unique(rfuToResType$RegionalForestUnit), NA))){
      stop("The second column of friLU must match a regional forest unit: ", 
           paste0(unique(friLU[,2])[which(!unique(friLU[,2]) %in%
                                              c(unique(rfuToResType$RegionalForestUnit), NA))],
                   sep = ", "),
           " does not match",
           call. = FALSE)
    }
  
    if(!all(raster::unique(fri) %in% c(friLU[,1], NA))){
      stop("All unique values in fri must be present in friLU",
           paste0(raster::unique(fri)[which(!raster::unique(fri) %in% c(friLU[,1], NA))]),
           call. = FALSE)
    }
    
  }
}
  