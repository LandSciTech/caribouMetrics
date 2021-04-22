#' Reclassify plc and fri
#'
#'
#' @param fri raster. Forest resource inventory. Forest resource inventory class
#'   ids must correspond to the ids in the \code{friLU}.
#' @param friLU data.frame. Lookup table to convert numbers in raster to
#'   regional forest unit codes. data.frame.  It should have two columns, the
#'   first must contain all the unique values in the supplied \code{LC}
#'   raster and the second must contain the names of regional forest units
#'   matching those provided in the table \code{rfuToResType}
#' @param rfuLU data.frame. Lookup table to convert RFU codes to resource types

#' @param plc raster. Provincial land cover
#' @param plcLU Lookup table to convert plc numbers to resource types 
#' @export
#' 
reclassFRI <- function(fri, friLU, rfuLU = rfuToResType){
  # TODO: get fri lu from raster legend?
  
  # check look up matches up
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
    
    if(!all(unique(friLU[,2]) %in% c(unique(rfuLU$RegionalForestUnit), NA))){
      stop("The second column of friLU must match a regional forest unit: ", 
           paste0(unique(friLU[,2])[which(!unique(friLU[,2]) %in%
                                            c(unique(rfuToResType$RegionalForestUnit), NA))],
                  sep = ", "),
           " does not match. Options are: ",paste0(unique(rfuLU$RegionalForestUnit),sep =", "),
           call. = FALSE)
    }
    
    if(!all(raster::unique(fri) %in% c(friLU[,1], NA))){
      stop("All unique values in fri must be present in friLU",
           paste0(raster::unique(fri)[which(!raster::unique(fri) %in% c(friLU[,1], NA))]),
           call. = FALSE)
    }
    
  }
  
  friLU <- friLU %>% set_names("ID", "RFU")
  
  # convert from FRI to ResType using supplied fu to rfu and internal
  # rfu to restype
  rclFRI <- left_join(friLU, rfuToResType, 
                     by = c("RFU" = "RegionalForestUnit")) %>% 
    left_join(resTypeCode, by = "ResourceType") %>% 
    select(-RFU, -ResourceType) %>% 
    as.matrix()
  
  fri <- reclassify(fri, rcl = rclFRI)
}

#' 
#' @rdname reclassFRI
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