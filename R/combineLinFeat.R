
#' Combine linear features
#'
#' Combine roads, rail and utilities in to one linear features object. All
#' linear features will be combined into one vector file which will be used to
#' calculate linear feature density. If the linear feature is provided as a
#' raster it will be converted to points which are interpreted based on
#' ptDensity parameter of `rasterizeLineDensity`.
#'
#' @param linFeats a list of linear feature data sets that will be combined.
#'   Data sets can be in the formats, SpatialLines, sf, raster, or a character
#'   vector with extension .shp or any extension accepted by the raster package
#'
#' @return sf object
#' @noRd

combineLinFeat <- function(linFeats){
  
  crsUse <- processLinFeat(linFeats[[1]]) %>% st_crs()
  
  linFeats <- purrr::map(linFeats, processLinFeat, crsUse = crsUse)
  
  names(linFeats) <- NULL
  
  do.call(rbind, linFeats) %>% mutate(linFID = 1:n()) %>% 
    st_set_agr("constant")
}

processLinFeat <- function(x, crsUse = NULL){
  if(is.character(x)){
    if(grepl(".shp$", x)){
      x <- sf::st_read(x, quiet = TRUE, agr = "constant")
    } else {
      x <- terra::rast(x)
    }
  }
  if(is(x, "RasterLayer")){
    x <- terra::rast(x)
  } 
  if(is(x, "SpatRaster")){
    x <- terra::subst(x, from = 0, to = NA) %>% 
      terra::as.points(na.rm = TRUE) %>% 
      sf::st_as_sf() %>% sf::st_set_agr("constant")
  } 
  if(is(x, "Spatial")){
    x <- sf::st_as_sf(x) %>% sf::st_set_agr("constant")
  }
  
  if(!is(x, "sf")){
    stop("linFeat elements must be one of the classes: sf, Spatial, ",
         "RasterLayer, or SpatRaster, not ",
         class(x), call. = FALSE)
  }
  
  if(is.null(crsUse)){
    return(x)
  }
  x <- x %>% 
    transmute(linFID = 1) %>% 
    st_transform(crsUse) 
}
