#'Reclassify provincial land cover
#'
#'Reclassify provincial land cover classes into resource types used for
#'`caribouHabitat`. Resource types are groups of land cover classes
#'expected to be relevant to boreal caribou (Hornseth and Rempel, 2016).
#'
#'The default lookup table `plcToResType` is provided with the package and
#'assumes Ontario 2001 Provincial Land Cover (PLC) classes. The
#'`fnlcToResType` table for the Ontario Far North Land Cover (FNLC) is
#'included with the package or you can supply a custom lookup table. The custom
#'lookup table must have two columns with names "PLCCode" and "ResourceType"
#'where "PLCCode" is the value in the land cover raster and "ResourceType" is
#'the corresponding letter code for the resource type. The possible resource
#'types are included in `resTypeCode` and descriptions of these codes can
#'be found in table S1.3 of [supplementary material](https://osf.io/x6e2q) for
#'Dyson et. al (2022).
#'
#'@param plc raster. Provincial land cover
#'@param plcLU Lookup table to convert land cover classes to resource types.
#'  Options are: `plcToResType` for Provincial Land Cover,
#'  `fnlcToResType` for the Ontario Far North Land Cover, or a custom
#'  lookup table. See Details.
#'
#'@return a RasterLayer with classes matching `resTypeCode`
#'
#' @examples
#' lc <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10, crs = 5070)
#' lc[] <- 13 # conifer
#' lc[1:3, 1:3] <- 1 # open water
#' lc[3:5, 3:5] <- 11 # deciduous
#' lc[8:10, 8:10] <- 19 # treed bog
#' lc[6:7, 6:7] <- 21 #treed fen
#'
#' resTypes <- reclassPLC(lc)
#'plot(resTypes)
#'@rdname reclassPLC
#'
#'@source Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
#'  woodland caribou (Rangifer tarandus caribou) across a gradient of
#'  anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'  <https://doi.org/10.1139/cjz-2015-0101>
#'
#'  Dyson, M., Endicott, S., Simpkins, C., Turner, J. W., Avery-Gomm, S.,
#'  Johnson, C. A., Leblond, M., Neilson, E. W., Rempel, R., Wiebe, P. A.,
#'  Baltzer, J. L., Stewart, F. E. C., & Hughes, J. (in press). Existing caribou
#'  habitat and demographic models need improvement for Ring of Fire impact
#'  assessment: A roadmap for improving the usefulness, transparency, and
#'  availability of models for conservation.
#'  <https://doi.org/10.1101/2022.06.01.494350>
#'
#' @family habitat
#'@export
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
      # check dataType aligns with minValue to avoid errors
      if(raster::minValue(plc) < 0 && grepl("U", raster::dataType(plc))){
        warning("raster::dataType(plc) is unsigned (", raster::dataType(plc),
             ") but the minimum value is negative (", raster::minValue(plc),
             "). Please resave plc with an appropriate datatype.", call. = FALSE)
      }
      stop("All unique values in plc must be present in plcLU: ",
           paste0(uniPLC[which(!uniPLC %in% c(plcLU[,1], NA))],
                  sep = ", "), " are missing.",
           call. = FALSE)
    }
    
  }
  # reclassify plc to resource types based on look up tables
  rclPLC <- plcLU %>% 
    left_join(resTypeCode, by = "ResourceType")%>% 
    select(-"ResourceType") %>% 
    as.matrix(rclPLC, rownames.force = FALSE)
  
  plc <- reclassify(plc, rclPLC)
}
