#' Reclassify provincial land cover
#'
#' Reclassify provincial land cover classes into resource types use for
#' \code{caribouHabitat}. The default lookup table \code{plcToResType} is
#' provided with the package and assumes Ontario 2001 land cover classes.
#'
#' @param plc raster. Provincial land cover
#' @param plcLU Lookup table to convert plc classes to resource types
#'
#' @return a RasterLayer with classes matching \code{resTypeCode}
#'
#' @rdname reclassPLC
#' @export
reclassPLC <- function(plc, plcLU = plcToResType){
  if(!is.null(plcLU)){
    # checks types and match names
    if(!inherits(plcLU, "data.frame")){
      stop("plcLU must be a data.frame", call. = FALSE)
    }
    if(!is.numeric(plcLU[,1])){
      stop("The first column of plcLU must be numeric", 
           call. = FALSE)
    }
    
    if(any(is.na(plcLU[,2]))){
      stop("plcLU contains NA in the second column. ", 
           "All land cover classes must have a resource type. ",
           "Use \"other\" for classes that do not fit another resource type ", 
           call. = FALSE)
    }
    
    if(!all(unique(plcLU[,2]) %in% unique(resTypeCode$ResourceType))){
      stop("The second column of plcLU must match a resource type: ", 
           paste0(unique(plcLU[,2])[which(!unique(plcLU[,2]) %in%
                                            unique(resTypeCode$ResourceType))],
                  sep = ", "),
           " does not match",
           call. = FALSE)
    }
    
    uniPLC <- raster::unique(plc)
    if(!all(uniPLC %in% c(plcLU[,1], NA))){
      stop("All unique values in plc must be present in plcLU",
           paste0(uniPLC[which(!uniPLC %in% c(plcLU[,1], NA))]),
           call. = FALSE)
    }
    
  }
  # reclassify plc to resource types based on look up tables
  rclPLC <- plcLU %>% 
    left_join(resTypeCode, by = "ResourceType")%>% 
    select(-ResourceType) %>% 
    as.matrix(rclPLC, rownames.force = FALSE)
  
  plc <- reclassify(plc, rclPLC)
}